"""Meltwater may not leave a snowpack that still holds cold content.

A pack below freezing refreezes the water produced at its surface: the latent
heat pays down its cold content, and nothing drains until that deficit is gone
and the pore space has filled. That is ripening, and it is how Snow-17 routes
its outflow (Anderson, 2006): water that reaches a pack still carrying `NEGHS`
goes back into the ice, and only the excess over the deficit and the holding
capacity leaves. The criterion holds a model's reported outflow `snm` to its own
reported cold content `csnow`:

    R = sum over steps that end with cold content C_t > C_eps of snm_t dt
    eps = R / W <= threshold,  W the peak pack water `snw`

What it does not forbid is melt. In a layered pack, surface melt percolates
into colder snow beneath it and refreezes (Marsh and Woo, 1984), which is why the
ROADMAP's earlier wording -- no melt while the pack is below freezing -- would
fail a correct layered model. Refreezing inside the pack is invisible here
because it releases no water; only water that leaves is counted.

The residual rests on the model's own cold content, so the criterion also bounds
what that report can say. It is a consistency check, not an energy balance:
it catches a model whose drainage ignores the cold content it itself carries,
and it cannot catch a model that reports less cold content than it has. No bound
built from the forcing alone could: over this case's winter the pack could
absorb several times more radiant energy than snowfall brings in cold, so the
forcing admits almost any cold content. Only a surface energy balance pins it,
which `melt_energy` carries for the models that report one.

Four checks bound the report, each set by this case rather than by snowpacks in
general:

- cold content and outflow cannot be negative, either of which would take a
  step out of the sum or cancel one inside it;
- the pack must open cold: a pack built under air that never reached freezing
  cannot open the scored block within `min_opening_depression_k` of melting;
- cold content cannot fall faster than energy arrives: net radiation plus
  sensible heat from air warmer than the pack, whose mean temperature follows
  from its own cold content and ice mass. Cold content carried away by ice that
  leaves is credited, since that is mass loss rather than warming;
- a pack under `min_peak_pack_mm` is refused, because melt at the snow-ground
  interface drains however cold the pack is and is fixed per day rather than
  proportional to the pack, so on a thin pack it would read as a violation.
"""

from __future__ import annotations

import numpy as np

from hydroturing.criteria.base import (
    FAIL,
    PASS,
    CriterionIncompatibleError,
    CriterionResult,
    criterion,
    make_window,
    segments,
)
from hydroturing.protocol import RunResult
from hydroturing.spec import ProbeSpec

# Specific heat of ice, J kg-1 K-1, the value `melt_energy` and the reference
# snowpack use. With `snw` in mm, which is kg m-2, the pack's mean temperature
# below freezing is -csnow / (C_ICE * ice).
C_ICE = 2100.0
SECONDS_PER_DAY = 86400.0

_KNOWN = {
    "outflow", "pack", "liquid", "cold_content", "driver", "air_temperature",
    "precipitation", "segment_column", "scored_label", "threshold",
    "min_peak_pack_mm", "cold_content_tolerance_j", "min_opening_depression_k",
    "sensible_heat_coefficient", "dry_tolerance_mm",
}


@criterion("snowpack_ripening")
def snowpack_ripening(run: RunResult, probe: ProbeSpec, params: dict) -> CriterionResult:
    """Water leaving a pack that still holds cold content, as a share of the pack."""
    # A misspelled name would otherwise be dropped in silence, and an author who
    # meant to change a bound would never learn it had not moved.
    unknown = set(params) - _KNOWN
    if unknown:
        raise ValueError(f"snowpack_ripening: unknown parameters {sorted(unknown)}")

    outflow = str(params.get("outflow", "snm"))
    pack = str(params.get("pack", "snw"))
    liquid = str(params.get("liquid", "lwsnl"))
    cold_content = str(params.get("cold_content", "csnow"))
    driver = str(params.get("driver", "rn"))
    air = str(params.get("air_temperature", "tas"))
    precipitation = str(params.get("precipitation", "pr"))
    segment_column = str(params.get("segment_column", "_regime"))
    scored_label = str(params.get("scored_label", "melt"))
    threshold = float(params.get("threshold", 0.05))
    min_peak = float(params.get("min_peak_pack_mm", 300.0))
    c_eps = float(params.get("cold_content_tolerance_j", 1.0))
    min_depression = float(params.get("min_opening_depression_k", 2.0))
    k_sensible = float(params.get("sensible_heat_coefficient", 10.0))
    dry_tolerance_mm = float(params.get("dry_tolerance_mm", 1e-6))

    # A non-finite or negative bound would disable a check rather than tighten
    # it: `.nan` compares false against everything, so it would read as passed.
    bounds = (threshold, min_peak, c_eps, min_depression, k_sensible, dry_tolerance_mm)
    if not all(np.isfinite(b) and b >= 0.0 for b in bounds):
        raise ValueError(
            "snowpack_ripening needs finite, non-negative threshold, "
            "min_peak_pack_mm, cold_content_tolerance_j, min_opening_depression_k, "
            f"sensible_heat_coefficient and dry_tolerance_mm; got {bounds}"
        )

    w = make_window(run, probe)
    for var in (outflow, pack, cold_content):
        if var not in w.table.columns:
            raise ValueError(f"snowpack_ripening needs '{var}' in the model result")
    # The forcing columns are the case's, read as bounds on what the weather
    # makes possible rather than as drivers a model must consume -- which is why
    # the probe does not list them under `requires.forcing`, and why a model
    # that never reads net radiation is still judged.
    for var in (driver, air, precipitation):
        if var not in w.forcing.columns:
            raise ValueError(
                f"snowpack_ripening needs forcing column '{var}'; this probe's "
                "generator does not produce it"
            )

    scored = [b for b in segments(w, segment_column) if b[0] == scored_label]
    if not scored:
        raise ValueError(
            f"snowpack_ripening found no '{scored_label}' block in column "
            f"'{segment_column}'; the generator has to label the stretch the pack "
            "ripens over"
        )

    # Rain on a cold pack is its own question -- much of it refreezes, but some
    # models route it straight through the snow module, and `snm` counts rain
    # passing through. The block has to be dry so that what leaves is meltwater.
    rain = w.forcing[precipitation].to_numpy(dtype=float)
    for label, start, stop in scored:
        fell = float(w.volume(rain[start:stop]).sum())
        if fell > dry_tolerance_mm:
            raise ValueError(
                f"snowpack_ripening scores '{label}', but {fell:.3f} mm of "
                f"'{precipitation}' falls within it; the block has to be dry, or "
                "rain passing through the pack is counted as meltwater leaving it"
            )

    q = w.volume(w.table[outflow].to_numpy(dtype=float))
    snw = w.table[pack].to_numpy(dtype=float)
    lw = (
        w.table[liquid].to_numpy(dtype=float)
        if liquid in w.table.columns
        else np.zeros(len(w.table))
    )
    csnow = w.table[cold_content].to_numpy(dtype=float)
    rn = w.forcing[driver].to_numpy(dtype=float)
    tas = w.forcing[air].to_numpy(dtype=float)
    dt_seconds = w.dt_days * SECONDS_PER_DAY

    if not (np.isfinite(q).all() and np.isfinite(snw).all() and np.isfinite(csnow).all()):
        return _fail("non-finite outflow, pack or cold content in the scored record",
                     outcome="non_finite")

    last = max(stop for _, _, stop in scored)
    peak = float(snw[:last].max()) if last else 0.0
    if peak < min_peak:
        raise CriterionIncompatibleError(
            f"the pack peaks at {peak:.1f} mm, under the {min_peak:g} mm this "
            "criterion needs: melt at the snow-ground interface drains however "
            "cold the pack is and does not scale with it, so on a thinner pack it "
            "would be scored as a violation"
        )

    blocks = []
    for label, start, stop in scored:
        if start > 0:
            c0 = float(csnow[start - 1])
            ice0 = float(snw[start - 1] - lw[start - 1])
        else:
            c0 = float(w.state0[cold_content])
            ice0 = float(w.state0[pack]) - (
                float(w.state0[liquid]) if liquid in w.state0.index else 0.0
            )
        c_open = np.concatenate([[c0], csnow[start:stop - 1]])      # at each step's start
        ice_open = np.concatenate([[ice0], (snw - lw)[start:stop - 1]])
        c_end = csnow[start:stop]
        ice_end = (snw - lw)[start:stop]
        qb = q[start:stop]

        # A negative cold content takes its step out of the sum; a negative
        # outflow cancels drainage inside it. Either is a report the contract
        # does not allow, and either would buy the pass this criterion guards.
        if (np.concatenate([[c0], c_end]) < -c_eps).any():
            return _fail(f"'{cold_content}' is negative in the '{label}' block; cold "
                         "content is an energy deficit and cannot be",
                         outcome="negative_cold_content", block=label)
        if (qb < -1e-9).any():
            return _fail(f"'{outflow}' is negative in the '{label}' block; water "
                         "leaving the snow module cannot be",
                         outcome="negative_outflow", block=label)

        # The pack must open cold. This case holds the air at or below -12 C for
        # the whole accumulation and below -14 C into the block, so a pack
        # reported within a few kelvin of melting on the opening row is a state
        # the forcing did not produce -- the plainest way to empty the sum.
        opening_temperature = -c0 / (C_ICE * ice0) if ice0 > 1e-6 else 0.0
        if ice0 > 1e-6 and opening_temperature > -min_depression:
            return _fail(
                f"the pack opens the '{label}' block at a mean {opening_temperature:.2f} C, "
                f"within {min_depression:g} K of melting, after a winter the air "
                "never brought above freezing: a pack reported ripe before it has "
                "warmed empties the check rather than passing it",
                outcome="opens_ripe", block=label,
                opening_pack_temperature_c=opening_temperature,
            )

        # Cold content cannot fall faster than energy arrives. What reaches the
        # pack in a step is at most the net radiation plus sensible heat from air
        # warmer than the pack; the pack's temperature comes from its own cold
        # content and ice. Ice that leaves takes its cold content with it, and
        # that fall is mass loss, not warming, so it is credited.
        ice_safe = np.maximum(ice_open, 1e-6)
        t_pack = -c_open / (C_ICE * ice_safe)
        supply = (np.maximum(rn[start:stop], 0.0)
                  + k_sensible * np.maximum(tas[start:stop] - t_pack, 0.0)) * dt_seconds
        carried = c_open * np.clip((ice_open - ice_end) / ice_safe, 0.0, 1.0)
        warming = c_open - c_end - carried
        excess = warming - supply
        worst_rate = int(np.argmax(excess))
        if excess[worst_rate] > c_eps:
            return _fail(
                f"the pack loses {warming[worst_rate] / 1e6:.2f} MJ/m2 of cold content in "
                f"one step of the '{label}' block, against the "
                f"{supply[worst_rate] / 1e6:.2f} MJ/m2 that net radiation and sensible "
                "heat could have brought it: it ripens faster than energy arrives",
                outcome="ripens_too_fast", block=label,
                step=start + worst_rate,
                warming_mj_m2=float(warming[worst_rate] / 1e6),
                supply_mj_m2=float(supply[worst_rate] / 1e6),
            )

        # A step counts when the pack still holds cold content at its end: the
        # water left and the pack was cold afterwards, so it cannot have left
        # because the pack ripened. A step that ends ripe may have paid its
        # deficit and drained within itself, which is what Snow-17 does when the
        # water reaching it outweighs NEGHS, and counting it would fail a model
        # for ripening on the day it ripened.
        cold = c_end > c_eps
        leaked = float(qb[cold].sum())
        share = leaked / peak
        blocks.append({
            "label": label,
            "outflow_while_cold_mm": leaked,
            "share_of_peak_pack": share,
            "outflow_total_mm": float(qb.sum()),
            "cold_steps": int(cold.sum()),
            "steps": int(stop - start),
            "opening_pack_temperature_c": opening_temperature,
            "opening_cold_content_mj_m2": c0 / 1e6,
            "worst_warming_share_of_supply": float(
                (warming / np.maximum(supply, 1.0)).max()
            ),
        })

    worst = max(blocks, key=lambda b: b["share_of_peak_pack"])
    ok = worst["share_of_peak_pack"] <= threshold
    return CriterionResult(
        name="snowpack_ripening",
        status=PASS if ok else FAIL,
        value=worst["share_of_peak_pack"],
        threshold=threshold,
        message=(
            f"{worst['outflow_while_cold_mm']:.1f} mm left the pack while it held "
            f"cold content, {worst['share_of_peak_pack']:.2%} of its "
            f"{peak:.0f} mm peak (limit {threshold:.0%})"
            + ("" if ok else ": meltwater drained from a pack that had not ripened")
        ),
        diagnostics={"blocks": blocks, "peak_pack_mm": peak},
    )


def _fail(message: str, **diagnostics) -> CriterionResult:
    return CriterionResult(
        name="snowpack_ripening",
        status=FAIL,
        value=None,
        threshold=None,
        message=message,
        diagnostics=diagnostics,
    )
