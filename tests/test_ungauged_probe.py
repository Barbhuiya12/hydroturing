import importlib.util
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from hydroturing.harness import build_case, evaluate_criteria, load_generator
from hydroturing.protocol import Case, RunResult
from hydroturing.registry import find_probe
from hydroturing.seeds import gate_seeds
from hydroturing.spec import REPO_ROOT


def load(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def probe():
    return find_probe("mass/ungauged-basin-closure")


@pytest.mark.parametrize("seed", range(8))
def test_single_case_is_deterministic_and_bounded(probe, seed):
    first, again = build_case(probe, seed), build_case(probe, seed)
    pd.testing.assert_frame_equal(first.forcing, again.forcing)
    assert first.static == again.static
    assert len(first.forcing) == 4380 and first.spinup_steps == 730
    assert first.forcing.pr.between(0, 35).all()
    assert first.forcing.tas.between(1, 31).all()
    assert first.forcing.pet.between(.68, 4.6).all()
    for i in range(730, 4380, 365):
        assert first.forcing.pr.iloc[i:i+365].sum() > 0
    assert 80 <= first.static["soil_capacity_mm"] <= 700
    assert .5 <= first.static["canopy_capacity_mm"] <= 2.5


def test_gate_covers_declared_categories_without_stratifying_random_evaluation(probe):
    seeds = gate_seeds(probe.id, probe.n_seeds)
    assert {s % 8 for s in seeds} == set(range(8))
    assert {s % 5 for s in seeds} == set(range(5))
    assert not probe.variants


@pytest.mark.parametrize("fault,criterion", [
    ("loss", "closure"), ("gain", "closure"), ("capacity", "state_bounds"),
    ("forcing", "forcing_fidelity"), ("et", "et_plausible"),
    ("negative", "non_degenerate"), ("frozen", "non_degenerate")])
def test_causal_faults_are_dormant_inside_domain_and_detected_outside(probe, fault, criterion):
    module = load(REPO_ROOT / "models/_spatial_faults/ht_adapter.py")
    for seed in [0, 1, 6]:
        case = build_case(probe, seed)
        records = case.forcing.to_dict("records")
        rows = module.simulate(records, case.static, fault=fault)
        scores = {r.name: r for r in evaluate_criteria({"control": RunResult(case, pd.DataFrame(rows), {}, 0)}, probe)}
        if seed in [0, 1]:
            assert rows == module.bucket.simulate(records, case.static)
            assert all(r.passed for r in scores.values())
        else:
            assert not scores[criterion].passed
            if fault in ["capacity", "et", "negative"]:
                assert scores["closure"].passed


def synthetic(probe, delayed=False):
    n = 2190
    rain = np.where(np.arange(n) % 20 == 0, 10., 0.) if delayed else 1.+np.arange(n)%3
    et = .2*rain if delayed else .4*rain
    runoff = np.r_[np.zeros(70), .8*rain[:-70]] if delayed else .6*rain
    storage = np.cumsum(rain-et-runoff) if delayed else 100.+10*np.sin(np.arange(n)/100)
    if not delayed:
        runoff -= np.diff(storage, prepend=storage[0])
    times = pd.date_range("2000-01-01", periods=n).strftime("%Y-%m-%d")
    forcing = pd.DataFrame(dict(time=times, pr=rain, tas=15., pet=3.))
    table = pd.DataFrame(dict(time=times, pr=rain, mrro=runoff, evspsbl=et,
                              mrso=0. if delayed else storage, snw=0., canopy=0., gw=storage if delayed else 0.))
    case = Case(probe.id, 0, forcing, dict(soil_capacity_mm=320., canopy_capacity_mm=2.), 730)
    return RunResult(case, table, {}, 0.)


def test_delayed_conservative_response_does_not_require_seven_day_correlation(probe):
    from hydroturing.criteria.degeneracy import non_degenerate
    run = synthetic(probe, delayed=True)
    assert all(r.passed for r in evaluate_criteria({"control": run}, probe))
    assert not non_degenerate(run, probe, {"min_response": .05}).passed


def test_annual_diagnostics_carry_prior_states_and_do_not_change_full_verdict(probe):
    diagnostic = load(REPO_ROOT / "scripts/diagnose_ungauged_closure.py")
    run = synthetic(probe)
    result = diagnostic.diagnose(run, probe)
    assert all(year["relative_residual"] < 1e-12 for year in result["years"])
    table = run.table.copy()
    table.loc[730:1094, "mrro"] += .5
    table.loc[1095:1459, "mrro"] -= .5
    altered = replace(run, table=table)
    assert all(r.passed for r in evaluate_criteria({"control": altered}, probe))
    result = diagnostic.diagnose(altered, probe)
    assert result["diagnostic_only"] is True
    assert result["full_window_closure"]["status"] == "pass"
    assert sum(year["exceeds_full_window_tolerance"] for year in result["years"]) == 2


def test_reference_domain_includes_its_boundaries():
    module = load(REPO_ROOT / "models/_spatial_faults/ht_adapter.py")
    for soil in [200., 450.]:
        for canopy in [1., 2.]:
            assert not module.outside(dict(soil_capacity_mm=soil, canopy_capacity_mm=canopy))
