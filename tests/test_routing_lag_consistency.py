"""Unit tests for the paired routing-lag criteria.

The acceptance gate exercises real subprocess reference models.  These tests
keep the measurement itself pinned: Snyder geometry, baseline removal, event and
peak centroids, discrete-time tolerance, response guards, and malformed paired
cases each have a small counterexample here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hydroturing import registry
from hydroturing.criteria import get, is_paired
from hydroturing.harness import build_case, compatibility_issues
from hydroturing.protocol import Case, RunResult
from hydroturing.spec import Criterion, ProbeSpec

VARIANTS = ("small", "medium", "large", "xlarge")
AREAS = (30.0, 300.0, 3000.0, 10000.0)
LENGTHS = (9.7965015486, 39.0005751285, 155.2640861436, 319.7409482737)
TARGET_LAGS = (1, 1, 2, 3)
SPINUP_DAYS = 5
PERIOD_DAYS = 90
EVENT = slice(SPINUP_DAYS + 20, SPINUP_DAYS + 21)

COMMON = {
    "event_column": "_event_pr",
    "runoff": "mrro",
    "baseline_days": 30,
    "min_response_fraction": 0.01,
    "min_pre_event_days": 10,
    "min_post_event_days": 30,
}

GEOMETRY_KEYS = {
    "area_km2",
    "main_channel_length_km",
    "centroid_channel_length_km",
}
EVALUATED_MODELS = (
    "cwatm",
    "dhbv2",
    "flex_lumped",
    "flex_topo",
    "google_flood_forecast",
    "lisflood",
    "sacsma_snow17",
    "summa",
    "wflow_sbm",
)


@pytest.fixture(scope="module")
def registered_probe() -> ProbeSpec:
    return registry.find_probe("momentum/routing-lag-consistency")


def _probe() -> ProbeSpec:
    return ProbeSpec(
        id="momentum/routing-lag-consistency",
        title="Routing lag follows catchment geometry",
        law="momentum",
        track="synthetic",
        version=1,
        authors=({"name": "Test Author"},),
        citation="",
        requires_fluxes=("mrro",),
        requires_states=(),
        generator="generate.py",
        n_seeds=3,
        timestep="PT1D",
        period_years=PERIOD_DAYS / 365,
        spinup_days=SPINUP_DAYS,
        max_output_mb=2.0,
        max_runtime_s=120.0,
        variants=VARIANTS,
        criteria=(
            Criterion("lag_time_bounds", {}),
            Criterion("scaling_monotonicity", {}),
        ),
        must_pass=("reference_router",),
        must_fail={"reference_bad_router": "lag_time_bounds"},
        provenance="synthetic",
        path=Path("."),
        period_days=PERIOD_DAYS,
        min_window_days=PERIOD_DAYS,
    )


def _forcing() -> pd.DataFrame:
    n = SPINUP_DAYS + PERIOD_DAYS
    extra = np.zeros(n)
    extra[EVENT] = 50.0
    return pd.DataFrame(
        {
            "time": pd.date_range("2000-01-01", periods=n, freq="D").strftime(
                "%Y-%m-%d"
            ),
            "pr": extra.copy(),
            "tas": np.full(n, 15.0),
            "pet": np.zeros(n),
            "_event_pr": extra,
        }
    )


def _run(variant: str, peak_lag_days: int, *, tie: bool = False) -> RunResult:
    forcing = _forcing()
    baseline = 0.2
    runoff = np.full(len(forcing), baseline)
    # The one-row rain pulse has its centroid in that event row.
    centre = EVENT.start + peak_lag_days
    if tie:
        runoff[centre : centre + 2] += 5.0
    else:
        runoff[centre - 1 : centre + 2] += (2.0, 5.0, 2.0)
    index = VARIANTS.index(variant)
    static = {
        "area_km2": AREAS[index],
        "main_channel_length_km": LENGTHS[index],
        "centroid_channel_length_km": 0.5 * LENGTHS[index],
    }
    case = Case(
        probe_id=f"momentum/routing-lag-consistency@{variant}",
        seed=11,
        forcing=forcing,
        static=static,
        spinup_steps=SPINUP_DAYS,
        timestep="PT1D",
    )
    table = pd.DataFrame({"time": forcing["time"], "mrro": runoff})
    return RunResult(case, table, {"status": "ok"}, 0.0)


def _runs(lags=TARGET_LAGS) -> dict[str, RunResult]:
    return {name: _run(name, lag) for name, lag in zip(VARIANTS, lags)}


def test_criteria_are_registered_as_paired():
    assert is_paired("lag_time_bounds")
    assert is_paired("scaling_monotonicity")


def test_registered_probe_declares_every_input_its_verdict_depends_on(registered_probe):
    assert registered_probe.requires_forcing == ("pr",)
    assert registered_probe.requires_static == (
        "area_km2",
        "main_channel_length_km",
        "centroid_channel_length_km",
    )


def test_registered_variants_share_forcing_and_change_only_geometry(registered_probe):
    cases = {
        variant: build_case(registered_probe, 20260914, variant)
        for variant in registered_probe.variants
    }
    control = cases[registered_probe.variants[0]]

    for case in cases.values():
        pd.testing.assert_frame_equal(case.forcing, control.forcing, check_exact=True)
        assert set(case.static) == set(control.static)
        assert {
            key: value for key, value in case.static.items() if key not in GEOMETRY_KEYS
        } == {
            key: value for key, value in control.static.items() if key not in GEOMETRY_KEYS
        }

    geometries = [
        tuple(cases[variant].static[key] for key in GEOMETRY_KEYS)
        for variant in registered_probe.variants
    ]
    assert len(set(geometries)) == len(registered_probe.variants)
    for key in GEOMETRY_KEYS:
        values = [cases[variant].static[key] for variant in registered_probe.variants]
        assert values == sorted(values)
        assert len(set(values)) == len(registered_probe.variants)

    scored = control.after_spinup(control.forcing)
    assert (scored["_event_pr"] > 0.0).sum() == 1
    np.testing.assert_array_equal(scored["pr"], scored["_event_pr"])


@pytest.mark.parametrize(
    "model_name",
    [
        "reference_snyder_router",
        "reference_instant_router",
        "reference_inverse_router",
    ],
)
def test_routing_references_declare_compatible_inputs(registered_probe, model_name):
    case = build_case(registered_probe, 20260914, "small")
    assert compatibility_issues(
        registry.find_model(model_name), registered_probe, case
    ) == []


def test_evaluated_models_are_not_judged_without_geometry_inputs(registered_probe):
    case = build_case(registered_probe, 20260914, "small")
    for model_name in EVALUATED_MODELS:
        model = registry.find_model(model_name)
        issues = compatibility_issues(model, registered_probe, case)
        message = "; ".join(issues)
        assert "main_channel_length_km" in message, model.name
        assert "centroid_channel_length_km" in message, model.name


def test_snyder_bounds_and_area_scaling_pass_for_physical_lags():
    probe = _probe()
    runs = _runs()
    bounds = get("lag_time_bounds")(
        runs,
        probe,
        {**COMMON, "lower_ratio": 0.5, "upper_ratio": 2.0,
         "discretization_tolerance_days": 0.5},
    )
    scaling = get("scaling_monotonicity")(
        runs,
        probe,
        {
            **COMMON,
            "reversal_tolerance_days": 0.5,
            "min_span_days": 2.0,
        },
    )
    assert bounds.passed, bounds.message
    assert scaling.passed, scaling.message
    expected_days = [
        bounds.diagnostics["variants"][name]["expected_lag_days"]
        for name in VARIANTS
    ]
    np.testing.assert_allclose(
        expected_days,
        [0.63110113, 1.12305225, 2.25004715, 3.33515212],
        rtol=1e-8,
    )
    assert scaling.diagnostics["increments_days"] == {
        "small->medium": 0.0,
        "medium->large": 1.0,
        "large->xlarge": 1.0,
    }
    assert scaling.diagnostics["span_days"] == 2.0


def test_instantaneous_runoff_fails_the_lower_lag_bound():
    result = get("lag_time_bounds")(
        _runs((0, 0, 0, 0)),
        _probe(),
        {**COMMON, "lower_ratio": 0.5, "upper_ratio": 2.0,
         "discretization_tolerance_days": 0.5},
    )
    assert not result.passed
    assert result.value > 0
    assert "outside" in result.message


def test_inverse_area_scaling_fails_monotonicity_without_failing_the_span():
    result = get("scaling_monotonicity")(
        _runs((1, 2, 1, 3)),
        _probe(),
        {
            **COMMON,
            "reversal_tolerance_days": 0.5,
            "min_span_days": 2.0,
        },
    )
    assert not result.passed
    assert result.value == pytest.approx(2.0)
    assert result.diagnostics["minimum_adjacent_increment_days"] == -1.0
    assert "medium to large" in result.message


def test_equal_runoff_maxima_use_their_temporal_centroid():
    runs = _runs()
    runs["small"] = _run("small", 1, tie=True)
    result = get("lag_time_bounds")(
        runs,
        _probe(),
        {**COMMON, "lower_ratio": 0.5, "upper_ratio": 2.0,
         "discretization_tolerance_days": 0.5},
    )
    small = result.diagnostics["variants"]["small"]
    assert small["tied_peak_steps"] == 2
    assert small["observed_lag_days"] == pytest.approx(1.5)


@pytest.mark.parametrize("criterion_name", ["lag_time_bounds", "scaling_monotonicity"])
def test_zero_response_fails_as_a_model_answer(criterion_name):
    runs = _runs()
    for name, run in list(runs.items()):
        table = run.table.copy()
        table["mrro"] = 0.2
        runs[name] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get(criterion_name)(runs, _probe(), dict(COMMON))
    assert not result.passed
    assert "no measurable storm response" in result.message
    assert result.diagnostics["response_fraction"] == pytest.approx(0.0)


def test_direct_criterion_call_defensively_fails_nonfinite_runoff():
    runs = _runs()
    table = runs["medium"].table.copy()
    table.loc[SPINUP_DAYS + 25, "mrro"] = np.nan
    run = runs["medium"]
    runs["medium"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("lag_time_bounds")(runs, _probe(), dict(COMMON))
    assert not result.passed
    assert result.diagnostics == {"variant": "medium", "nonfinite_count": 1}


def test_a_peak_on_the_final_row_has_no_measurable_lag():
    runs = _runs()
    run = runs["xlarge"]
    table = run.table.copy()
    table["mrro"] = 0.2
    table.loc[EVENT.start:, "mrro"] = np.linspace(
        0.21, 5.0, len(table) - EVENT.start
    )
    runs["xlarge"] = RunResult(run.case, table, run.meta, run.wall_seconds)
    result = get("lag_time_bounds")(runs, _probe(), dict(COMMON))
    assert not result.passed
    assert "final row" in result.message


def test_forcing_must_be_exactly_the_same_between_catchments():
    runs = _runs()
    run = runs["large"]
    forcing = run.case.forcing.copy()
    forcing.loc[0, "tas"] += 0.01
    runs["large"] = RunResult(
        Case(
            probe_id=run.case.probe_id,
            seed=run.case.seed,
            forcing=forcing,
            static=run.case.static,
            spinup_steps=run.case.spinup_steps,
            timestep=run.case.timestep,
        ),
        run.table,
        run.meta,
        run.wall_seconds,
    )
    with pytest.raises(ValueError, match="exactly the same forcing"):
        get("lag_time_bounds")(runs, _probe(), dict(COMMON))


def test_event_annotation_must_be_contiguous_and_keep_response_tail():
    runs = _runs()
    for name, run in list(runs.items()):
        forcing = run.case.forcing.copy()
        forcing["pr"] = 0.0
        forcing["_event_pr"] = 0.0
        forcing.loc[len(forcing) - 3 :, ["pr", "_event_pr"]] = 10.0
        case = Case(
            probe_id=run.case.probe_id,
            seed=run.case.seed,
            forcing=forcing,
            static=run.case.static,
            spinup_steps=run.case.spinup_steps,
            timestep=run.case.timestep,
        )
        runs[name] = RunResult(case, run.table, run.meta, run.wall_seconds)
    with pytest.raises(ValueError, match="truncates the routing experiment"):
        get("lag_time_bounds")(runs, _probe(), dict(COMMON))


def test_unannotated_rain_in_the_scored_record_is_a_case_error():
    runs = _runs()
    for name, run in list(runs.items()):
        forcing = run.case.forcing.copy()
        forcing.loc[SPINUP_DAYS + 60, "pr"] = 1.0
        runs[name] = RunResult(
            Case(
                probe_id=run.case.probe_id,
                seed=run.case.seed,
                forcing=forcing,
                static=run.case.static,
                spinup_steps=run.case.spinup_steps,
                timestep=run.case.timestep,
            ),
            run.table,
            run.meta,
            run.wall_seconds,
        )
    with pytest.raises(ValueError, match="must be dry outside"):
        get("lag_time_bounds")(runs, _probe(), dict(COMMON))


def test_daily_quantisation_may_flatten_one_adjacent_lag_step():
    result = get("scaling_monotonicity")(
        _runs((0, 0, 1, 2)),
        _probe(),
        {
            **COMMON,
            "reversal_tolerance_days": 0.5,
            "min_span_days": 2.0,
        },
    )
    assert result.passed, result.message


def test_a_fixed_lag_fails_the_full_ladder_span():
    result = get("scaling_monotonicity")(
        _runs((1, 1, 1, 1)),
        _probe(),
        {
            **COMMON,
            "reversal_tolerance_days": 0.5,
            "min_span_days": 2.0,
        },
    )
    assert not result.passed
    assert "lag span" in result.message


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("main_channel_length_km", 0.0, "positive area"),
        ("centroid_channel_length_km", np.nan, "non-finite static field"),
        ("centroid_channel_length_km", 1000.0, "longer than its main channel"),
        ("area_km2", 20.0, "outside Snyder's published"),
    ],
)
def test_invalid_static_geometry_is_a_case_error(field, value, message):
    runs = _runs()
    run = runs["small"]
    static = dict(run.case.static)
    static[field] = value
    runs["small"] = RunResult(
        Case(
            probe_id=run.case.probe_id,
            seed=run.case.seed,
            forcing=run.case.forcing,
            static=static,
            spinup_steps=run.case.spinup_steps,
            timestep=run.case.timestep,
        ),
        run.table,
        run.meta,
        run.wall_seconds,
    )
    with pytest.raises(ValueError, match=message):
        get("lag_time_bounds")(runs, _probe(), dict(COMMON))
