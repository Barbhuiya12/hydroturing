#!/usr/bin/env python3
"""HydroTuring adapter for the lumped FLEX/HBV model of chrimerss/HydrologicModels.

The equations are those of `lumped_model/HBVMod.py` at commit cc0aa6f
(the semi-distributed directory's copy of the same file, which clips every
outflow at its store), rewritten so that the step is a parameter rather
than the constant `dt = 1` of the original, and with two deviations that
are stated here so they can be argued with:

* Transpiration in the original is `Ep * Su / (Sumax * Ce)` with no upper
  limit, so for `Ce < 1` a wet soil transpires above the potential rate.
  FLEX's published form is `Ep * min(1, Su / (Ce * Sumax))`; that limit is
  applied here. A physical reference model must not evaporate more than
  the atmosphere asks for, and without the limit the et_plausible
  criterion catches it (README.md has the number).
* The catchment's soil and canopy capacities come from static.json rather
  than from the calibrated parameter set, because a physical model is told
  its catchment; the remaining parameters keep their Wark values.

Stores reported: canopy (Si), soil (Su), groundwater (the slow reservoir
Ss), and the fast reservoir plus the water inside the triangular lag as
`channel`, since both hold runoff that has been generated and not yet
released. No snow module: snow is identically zero and precipitation
below freezing is treated as rain, which conserves water and is wrong
about timing, a limitation of the model rather than the adapter.

Rates with a time in their units are rescaled to the step exactly as any
submitted model's must be: per-day fractions as 1 - (1 - k)^dt, per-day
amounts as amount * dt, the lag in days.

A `stage` is also reported, as a diagnostic rather than a store: it is the
depth the reach's stored water makes in the channel cross-section, so that
`momentum/stage-discharge-monotonic` has a gauge to read. It is derived from
`channel` and is never differenced into any budget.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLUMNS = ["time", "pr", "evspsbl", "mrro", "dis", "gwex", "mrso", "snw", "canopy", "gw", "channel", "stage"]
MODEL = {"name": "flex_lumped", "version": "1.0.0"}

# Default reach geometry, used when the catchment does not hand one over.
# These are the same defaults the reference rating adapters carry, so every
# model judged on this probe is read in the same channel.
DEFAULT_WIDTH_M = 18.0
DEFAULT_SLOPE = 0.0015
DEFAULT_MANNING_N = 0.035
DEFAULT_REACH_LENGTH_M = 4500.0
SECONDS_PER_DAY = 86400.0

TIMESTEP_DAYS = {"PT1D": 1.0, "PT1H": 1.0 / 24.0, "PT15M": 1.0 / 96.0, "PT5M": 1.0 / 288.0, "PT1M": 1.0 / 1440.0}

# lumped_model/A_MC_HBV.py's feasible ranges, and a set inside them near the
# repository's best Wark fit (C_run_model_lumped.py): Imax Ce Sumax beta Pmax Tlag Kf Ks.
PARAMS = {
    "Imax": 2.0,      # mm, replaced by canopy_capacity_mm when the catchment gives one
    "Ce": 0.68,       # soil moisture fraction at which transpiration reaches its potential
    "Sumax": 90.0,    # mm, replaced by soil_capacity_mm when the catchment gives one
    "beta": 1.85,     # partition curvature
    "Pmax": 0.09,     # mm/day, percolation at a full unsaturated store
    "Tlag": 1.1,      # days, triangular lag
    "Kf": 0.1,        # 1/day, fast reservoir
    "Ks": 0.008,      # 1/day, slow reservoir
}


def lag_weights(tlag_steps: float) -> list[float]:
    """Weigfun.py: a triangle of base Tlag, discretised to steps, summing to one."""
    nmax = int(-(-tlag_steps // 1))  # ceil
    if nmax <= 1:
        return [1.0]
    w = [0.0] * nmax
    th = tlag_steps / 2.0
    nh = int(th // 1)
    for i in range(nh):
        w[i] = ((i + 1) - 0.5) / th
    i = nh
    w[i] = (1 + ((i + 1) - 1) / th) * (th - (th // 1)) / 2 + (1 + (tlag_steps - (i + 1)) / th) * ((th // 1) + 1 - th) / 2
    for i in range(nh + 1, int(tlag_steps // 1)):
        w[i] = (tlag_steps - (i + 1) + 0.5) / th
    if tlag_steps > tlag_steps // 1:
        w[int(tlag_steps // 1)] = (tlag_steps - (tlag_steps // 1)) ** 2 / (2 * th)
    total = sum(w)
    return [x / total for x in w]


def per_step(fraction_per_day: float, dt: float) -> float:
    return 1.0 - (1.0 - fraction_per_day) ** dt


def stage_of(mrro_mm_per_day: float, gw_mm: float, static: dict, dt: float) -> float:
    """The level a gauge in the reach would read, in metres.

    A stage is a *length* read off a staff gauge in a cross-section, so it is
    built from the flow the reach is carrying rather than from a catchment
    depth. Manning's normal depth is the honest bridge between the two: the
    flow through the channel sets the depth it runs at, and the depth is what
    a gauge reads.

        Q   = w * h * (1/n) * h^(2/3) * S^(1/2)
        h   = ( Q * n / (w * sqrt(S)) )^(3/5)

    The flow is the reach's own drain plus the baseflow its groundwater store
    is releasing: `mrro` is the water in transit that returns the flood, `gw`
    is the water still on its way down the catchment. The two have different
    time constants, and the gauge sees both, which is why the level on the
    rising limb sits below the level the same flow makes once the store behind
    it has drained.

    Both terms are flow *rates*, in mm/day: a stage is a reading a gauge would
    give at an instant, so it cannot depend on how often the model chooses to
    write a row. An earlier version divided the stores — depths — by `dt` to
    make a rate, which made the same reach read 78x deeper at PT1M than at
    PT1D. `mrro` is already a rate and `gw` is turned into one below.
    """
    area_km2 = float(static.get("area_km2", 0.0))
    width_m = float(static.get("width_m", DEFAULT_WIDTH_M))
    slope = float(static.get("slope", DEFAULT_SLOPE))
    manning_n = float(static.get("manning_n", DEFAULT_MANNING_N))
    if area_km2 <= 0.0 or width_m <= 0.0 or slope <= 0.0 or dt <= 0.0:
        return 0.0
    # The reach's own drain is already a rate; the groundwater store contributes
    # at the rate it releases over the step it was reported for.
    gw_mm_per_day = max(gw_mm, 0.0) / dt
    q_m3s = (max(mrro_mm_per_day, 0.0) + gw_mm_per_day) * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY
    if q_m3s <= 0.0:
        return 0.0
    return (q_m3s * manning_n / (width_m * slope ** 0.5)) ** 0.6


def simulate(forcing: list[dict], static: dict, dt: float) -> list[dict]:
    p = dict(PARAMS)
    if "soil_capacity_mm" in static:
        p["Sumax"] = float(static["soil_capacity_mm"])
    if "canopy_capacity_mm" in static:
        p["Imax"] = float(static["canopy_capacity_mm"])
    kf, ks = per_step(p["Kf"], dt), per_step(p["Ks"], dt)
    pmax = p["Pmax"] * dt
    weights = lag_weights(p["Tlag"] / dt)

    si = 0.0
    su = 0.5 * p["Sumax"]
    sf = 0.0
    ss = 0.0
    generated: list[float] = []  # unrouted runoff per step, for the lag
    rows = []
    for step in forcing:
        pr_rate, pet_rate = step["pr"], step["pet"]
        P = pr_rate * dt
        Ep = pet_rate * dt

        # The human term, first call on the store: a prescribed withdrawal
        # (`abstr`, mm/day net of return flow) is taken from the unsaturated
        # store, and whatever the store cannot supply later from the day's
        # generated runoff. Absent the column nothing changes.
        want = max(step.get("abstr", 0.0), 0.0) * dt
        removed = min(su, want)
        su -= removed

        # Interception store; evaporation from it only on rainless steps.
        if P > 0.0:
            si += P
            pe = max(0.0, si - p["Imax"])
            si -= pe
            ei = 0.0
        else:
            pe = 0.0
            ei = min(Ep, si)
            si -= ei

        # Unsaturated store: a share rho of effective rain goes to fast runoff.
        if pe > 0.0:
            rho = (su / p["Sumax"]) ** p["beta"]
            su += (1.0 - rho) * pe
            quf = rho * pe
        else:
            quf = 0.0

        # Transpiration, limited by soil moisture and by demand.
        ep_left = max(0.0, Ep - ei)
        ea = ep_left * min(1.0, su / (p["Sumax"] * p["Ce"]))
        ea = min(ea, su)
        su -= ea

        # Percolation to the slow reservoir.
        qus = min(pmax * su / p["Sumax"], su)
        su -= qus

        # Fast and slow linear reservoirs.
        sf += quf
        qf = min(kf * sf, sf)
        sf -= qf
        ss += qus
        qs = min(ks * ss, ss)
        ss -= qs

        # The human term, second call: the day's generated runoff.
        rest = want - removed
        take = min(qf, rest)
        qf -= take
        removed += take
        rest -= take
        take = min(qs, rest)
        qs -= take
        removed += take

        generated.append(qf + qs)
        n = len(generated)
        routed = sum(weights[k] * generated[n - 1 - k] for k in range(min(len(weights), n)))

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (ei + ea) / dt,
            "mrro": routed / dt,
            "gwex": -removed / dt,
            "mrso": su,
            "snw": 0.0,
            "canopy": si,
            "gw": ss,
            "channel": sf,  # completed below with the lag's contents
        })
    # Water inside the lag: cumulative generated minus cumulative routed.
    cum_gen = 0.0
    cum_out = 0.0
    for row, g in zip(rows, generated):
        cum_gen += g
        cum_out += row["mrro"] * dt
        row["channel"] += cum_gen - cum_out
        # The discharge the gauge reads is the reach's own outflow plus the
        # baseflow still arriving: the same flow `stage_of` is handed, reported
        # in m3/s so `momentum/stage-discharge-monotonic` can score the rating
        # against discharge rather than against the store.
        row["dis"] = (
            (row["mrro"] + max(row["gw"], 0.0) / dt)
            * 1e-3
            * static["area_km2"]
            * 1e6
            / SECONDS_PER_DAY
        )
        row["stage"] = stage_of(row["mrro"], row["gw"], static, dt)
    return rows


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet", "abstr"):
            if key in row:
                row[key] = float(row[key])
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    rows = simulate(forcing, static, TIMESTEP_DAYS[timestep])

    out = io_dir / request["output"]["table"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows),
                    "notes": {"snw": "identically zero; the model has no snow module",
                              "channel": "fast reservoir plus water inside the triangular lag"}}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
