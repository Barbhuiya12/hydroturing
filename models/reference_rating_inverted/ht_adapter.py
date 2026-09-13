#!/usr/bin/env python3
"""HydroTuring adapter for reference_rating_inverted: the routed bucket whose
gauge is wired to the wrong part of the reach.

The routing is the same two-path reach the honest model uses — a fast
floodplain and a slow channel — and the discharge is identical. The fault is
which of them the gauge reads: it stands in the fast path instead of the slow
one.

That inverts the loop rather than destroying it. The fast path fills and
empties within the flood, so the gauge peaks *with* the wave instead of after
it, and at equal water in transit the rise reads above the fall. Stage is
still a strictly increasing function of the quantity it is drawn from, so the
rating never falls against discharge and the monotonicity check passes
outright: this model fails only the loop, which is exactly the separation the
two criteria exist to make.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

COLUMNS = [
    'time', 'pr', 'evspsbl', 'mrro', 'dis',
    'mrso', 'snw', 'canopy', 'channel', 'stage',
]

MODEL = {"name": "reference_rating_inverted", "version": "1.0.0"}

EVAP_SHAPE = 0.5
SECONDS_PER_DAY = 86400.0

FLOODPLAIN_RESIDENCE_D = 0.35
CHANNEL_RESIDENCE_D = 4.0

TIMESTEP_DAYS = {
    "PT1D": 1.0,
    "PT1H": 1.0 / 24.0,
    "PT15M": 1.0 / 96.0,
    "PT5M": 1.0 / 288.0,
    "PT1M": 1.0 / 1440.0,
}


def _stage(store_mm, static):
    area_km2 = static["area_km2"]
    width_m = static.get("width_m", 18.0)
    reach_m = static.get("reach_length_m", 4500.0)
    stored_m3 = max(store_mm, 0.0) * 1e-3 * area_km2 * 1e6
    return float(stored_m3 / (width_m * reach_m))


def simulate(forcing, static, dt_days=1.0):
    soil_cap = static["soil_capacity_mm"]
    canopy_cap = static["canopy_capacity_mm"]
    ddf = static["degree_day_factor_mm_per_C_day"]
    k_base = static["baseflow_coefficient"]
    t_snow = static["snow_threshold_degC"]
    area_km2 = static["area_km2"]

    k_fast = 1.0 - math.exp(-dt_days / FLOODPLAIN_RESIDENCE_D)
    k_slow = 1.0 - math.exp(-dt_days / CHANNEL_RESIDENCE_D)

    soil = 0.5 * soil_cap
    swe = 0.0
    canopy = 0.0
    fast = 0.0
    slow = 0.0
    rows = []

    for step in forcing:
        pr_rate, tas, pet_rate = step["pr"], step["tas"], step["pet"]
        pr = pr_rate * dt_days
        pet = pet_rate * dt_days

        snowfall = pr if tas < t_snow else 0.0
        rain = 0.0 if tas < t_snow else pr

        swe += snowfall
        melt = min(swe, ddf * max(tas - t_snow, 0.0) * dt_days)
        swe -= melt

        water_in = rain + melt
        intercepted = min(canopy_cap - canopy, water_in)
        canopy += intercepted
        throughfall = water_in - intercepted

        canopy_evap = min(canopy, pet)
        canopy -= canopy_evap
        pet_left = pet - canopy_evap

        soil += throughfall
        surface = max(0.0, soil - soil_cap)
        soil -= surface
        baseflow = k_base * soil * dt_days
        soil -= baseflow
        soil_evap = min(soil, pet_left * min(1.0, soil / (EVAP_SHAPE * soil_cap)))
        soil -= soil_evap

        yield_mm = surface + baseflow
        to_fast = 0.7 * yield_mm
        to_slow = yield_mm - to_fast
        fast += to_fast
        slow += to_slow

        q_fast = fast * k_fast
        q_slow = slow * k_slow
        fast -= q_fast
        slow -= q_slow

        q_total = q_fast + q_slow
        dis_m3s = q_total / dt_days * 1e-3 * area_km2 * 1e6 / SECONDS_PER_DAY
        # The fault: the gauge stands in the fast path.
        stage = _stage(fast + q_fast, static)

        rows.append({
            "time": step["time"],
            "pr": pr_rate,
            "evspsbl": (canopy_evap + soil_evap) / dt_days,
            "mrro": q_total / dt_days,
            "dis": dis_m3s,
            "mrso": soil,
            "snw": swe,
            "canopy": canopy,
            "channel": fast + slow,
            "stage": stage,
        })
    return rows


def read_request(path: Path) -> tuple[dict, Path]:
    request = json.loads(path.read_text())
    return request, path.parent


def read_forcing(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("pr", "tas", "pet"):
            if key in row:
                row[key] = float(row[key])
    return rows


def write_result(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()

    request_path = Path(args.request).resolve()
    request, io_dir = read_request(request_path)
    forcing = read_forcing(io_dir / request["input"]["forcing"])
    static = json.loads((io_dir / request["input"]["static"]).read_text())

    timestep = request.get("timestep", "PT1D")
    if timestep not in TIMESTEP_DAYS:
        raise SystemExit(f"unsupported timestep {timestep!r}")
    rows = simulate(forcing, static, TIMESTEP_DAYS[timestep])

    write_result(io_dir / request["output"]["table"], rows)
    (io_dir / request["output"]["run"]).write_text(
        json.dumps({"status": "ok", "model": MODEL, "n_steps": len(rows)}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
