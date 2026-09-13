# Water balance integrity across declared catchment attributes

Shunan Zhou, Dalian University of Technology, Dalian, China. [ORCID](https://orcid.org/0009-0006-4072-6799), GitHub: Davidanan.

## Question and scope

Hold a model, its learned weights and its non-prescribed parameters fixed. Without recalibration against target-catchment discharge observations, does its reported budget remain closed across the catchment attributes supplied to it?

This is a synthetic daily rainfall–runoff experiment. The **experimental reference domain** is soil capacity 200–450 mm and canopy capacity 1–2 mm. Domain membership is defined by these intervals. It is not inferred from a model's training data. Capacity is the size of a store, not the quantity of water currently in it: capacities stay fixed within each independently initialised simulation.

## Native seed-per-catchment design

One seed generates one catchment and one weather record. All five native criteria see that run, exactly as on `mass/catchment-closure`. There are no paired variants, scoring wrappers, schema changes or default seed-stratification changes. Twelve fixed gate seeds cover all eight attribute categories below. Default evaluations use twelve fresh framework seeds; the categories and climates actually sampled are reported, without a guarantee that every category occurs in every random evaluation.

| `seed % 8` | Category | Soil capacity, mm | Canopy capacity, mm |
|---|---|---:|---:|
| 0 | reference | 320 | 2 |
| 1 | reference_medium | 230–280 | 2 |
| 2 | soil_small | 80–140 | 2 |
| 3 | soil_large | 500–700 | 2 |
| 4 | canopy_small | 320 | 0.5–0.8 |
| 5 | canopy_large | 320 | 2.2–2.5 |
| 6 | joint_small | 80–140 | 0.5–0.8 |
| 7 | joint_large | 500–700 | 2.2–2.5 |

Within a band, capacity is uniform. The experiment samples selected single-attribute changes and two joint corners, not every combination or a measured population frequency. It imposes no monotonicity or equality between catchments' hydrographs. Attribute and weather random streams are separate; the model receives capacities and forcing, not the generator seed or category label.

## Generated weather and duration

PT1D: 730 warm-up days followed by 3,650 evaluated days. Evaluation years are consecutive 365-day blocks, not calendar years. `min_window_days: 3650` prevents a submitted model's event-window setting from shortening this experiment. Two-year warm-up is not a claim of equilibrium for every model; actual storage change is always retained.

| `seed % 5` | Weather regime | Wet intercept / seasonal amplitude | Mean T / amplitude, °C | Mean PET / amplitude, mm/day | Latitude |
|---|---|---|---|---|---|
| 0 | temperate winter-wet | 0.28 / 0.10 | 15 / 9 | 2.2 / 1.3 | 38° |
| 1 | warm humid | 0.40 / 0.03 | 24 / 3 | 3.0 / 0.5 | 20° |
| 2 | summer rainfall | 0.30 / −0.20 | 23 / 5 | 3.0 / 1.0 | 25° |
| 3 | seasonally water-limited | 0.18 / 0.12 | 19 / 8 | 2.6 / 1.3 | 35° |
| 4 | cool maritime rain | 0.32 / 0.06 | 10 / 6 | 1.6 / 0.8 | 48° |

Wet-day probability is the intercept plus seasonal amplitude times cos(phase), plus 0.22 if the preceding day was wet. Wet-day amounts are drawn directly as `35 × Beta(0.9,3.2)` mm, without clipping or seed rejection. T uses uniform noise ±3°C and PET a factor uniform in [0.85,1.15]. Resulting bounds are P∈[0,35] mm/day, T∈[1,31]°C and PET∈[0.68,4.6] mm/day. The supplied conditions are rain-only; sustained ice, severe aridity, stronger storms and managed water systems are not represented by this generator. No model is exempted from any generated case on the basis of its familiar climate.

Area is 250 km². Bucket auxiliary defaults are baseflow coefficient 0.006/day, snow threshold 0°C and melt factor 3.2 mm/(°C day). These are declared attributes, not parameters fitted to an observed hydrograph.

The rainfall bound also stays inside a sufficient stability condition for the unchanged explicit FLEX soil partition: with `x=S/C` and beta=1.85, `1−x^beta ≤ beta(1−x)` for x∈[0,1]. Thus `p≤C/beta` prevents a capacity overshoot before subsequent drainage. At C=80 mm, the sufficient bound is 43.24 mm, above the generated 35 mm. This is a declared generator choice, not a repair to FLEX or a promise about stronger storms.

## Budget and native criteria

For the complete evaluation window W:

```text
R = sum_W[(P + G - ET - Q) * dt] - (S_end - S_before_W)
e = abs(R) / sum_W[P * dt]
```

P comes from supplied forcing; G is the signed external `gwex` if reported. ET=`evspsbl`, Q=`mrro`, dt=1 day. S contains required soil/snow/canopy and every reported groundwater/channel store, using the same catchment-area depth. Internal exchanges are not external inputs. A model must report every store it has, with its physical mapping documented in its adapter.

| Native criterion | Requirement |
|---|---|
| `closure` | Full-window e≤0.05; no pooling across seeds |
| `state_bounds` | Soil/canopy within supplied capacities; snow, groundwater and channel nonnegative where reported; native 1e-6 mm tolerance |
| `et_plausible` | ET nonnegative; full-window ET/PET≤1 with native numerical tolerances |
| `non_degenerate` | Runoff ratio against P+G in [0.02,0.98]; Q/ET CV≥0.1; seven-day rainfall/runoff correlation ≥0.05 |
| `forcing_fidelity` | Native maximum echoed-P deviation divided by mean supplied P≤1e-6 |

State bounds are checked over the native scored window; the budget uses the preceding measured inventory. The generated rainfall denominator is positive; the native criterion rejects a zero denominator. There is no custom threshold margin, clipping of model outputs or repair of the residual.

The native seven-day rainfall/runoff correlation screen is enabled (`min_response: 0.05`), consistent with catchment-closure and the rain-only guidance in the spatial-extrapolation template. It is a scoped anti-degeneracy screen, not a consequence of the mass-conservation identity. A conservative, sufficiently delayed response can fail this screen; that limitation does not justify disabling it only for this probe. Thresholds are fixed before evaluating models.

## Annual diagnostics and interpretation

`scripts/diagnose_ungauged_closure.py` reports ten annual residuals, each using the actual inventory immediately before that year. It regenerates the recorded seed, checks input identity and validates the stored output through the native protocol. It never modifies the standard verdict or archive.

Opposite annual errors may cancel in the full budget. Annual exceedances are therefore visible, but full-window PASS does not establish annual or timestep closure. A regression test explicitly exercises this limitation. Other possible escapes include a conservative attribute-blind small store, mutually compensating fabricated terms and errors below the engineering tolerance. Passing does not establish prediction accuracy, correct capacity sensitivity or the physical provenance of self-reported external exchange.

## Positive and negative evidence

The three positive references are `reference_bucket`, `flex_lumped` and `flex_topo`, with unchanged numerical equations. Their existing capacity mappings are declared in `model.yaml`. FLEX-lumped maps the supplied capacities to `Sumax` and `Imax`; FLEX-topo scales the area-weighted soil capacities and sets interception capacity in each unit.

The probe requires consumption of both `soil_capacity_mm` and `canopy_capacity_mm` through the existing `requires.static` mechanism. Declaring an input is not enough: its native mapping must match the reported storage. Wflow, SUMMA, CWatM and LISFLOOD also have existing mappings, now declared without changing their model equations. Models lacking the required input capability receive the framework's native incompatibility outcome. SAC-SMA has no canopy store and is no longer a must-pass reference for this two-capacity experiment; δHBV is not given a dedicated capacity override. Google also lacks budget outputs. All ten previously evaluated models receive a native archive record on this probe; no pending-evaluation exception is introduced.

Seven causal fixtures use one shared implementation of the existing bucket equations. They are inactive inside the experimental reference domain, and read no case IDs, generator seeds, future data or evaluation boundaries. Their failures are:

| Fixture suffix | Out-of-domain fault | Required criterion | Budget may still close? |
|---|---|---|---|
| loss | Export ET and Q at 0.8 of their values | closure | No |
| gain | Add 0.1P to Q without a debit | closure | No |
| capacity | Add fixed soil and canopy reporting offsets, each equal to its supplied capacity + 1 mm | state_bounds | Yes |
| forcing | Simulate and echo 0.8P | forcing_fidelity | Not against supplied P |
| et | ET=1.1PET with compensating Q | et_plausible | Yes |
| negative | Reverse Q and compensate in ET | non_degenerate | Yes |
| frozen | Zero Q/ET with constant states | non_degenerate | No |

The negative-runoff fixture is caught by the native runoff-ratio check; this probe introduces no separate daily flux-sign criterion. The seven attribute-dependent fixtures are the required negative references. The existing `reference_cheater` and `reference_degenerate` do not consume both required capacities and are not assigned a must-fail verdict on this probe; their other gates are unchanged. The seven new references also pass the ordinary fixed-attribute catchment-closure cross-check, demonstrating why sampling attributes adds evidence beyond that probe. For the capacity fixture, normal bucket dynamics run with the supplied capacities. Only the reported absolute soil and canopy stores receive an offset, constant throughout the run including spinup. At full fault severity each store exceeds its own capacity even in the largest-capacity category, while the identical endpoint offsets cancel from storage differences. This verifies absolute-storage bounds independently of closure. The fault is dormant inside the experimental reference domain. It is not a claim that real models must fail outside that domain, and is never applied to submitted-model results. Regression tests cover all eight attribute categories across all five climates (seeds 0–39), including the large-capacity cases. Fault sizes are declared test levels rather than inferred real-model error rates.

## Reproduction

```bash
ht validate
ht gate --probe mass/ungauged-basin-closure
ht run --model flex_lumped --probe mass/ungauged-basin-closure --gate-seeds --workdir /tmp/ungauged --json /tmp/flex.json
python scripts/diagnose_ungauged_closure.py /tmp/ungauged/flex_lumped__mass__ungauged-basin-closure__1679270760 --seed 1679270760 --output /tmp/annual.json
pytest -q
```

Scientific-model demonstrations use the original adapters, serial execution and the standard per-invocation 60-second limit. The probe acceptance gate's trusted references run without Docker. Model runtime and probe overhead are distinct quantities.

## Sources and parameter semantics

- [HydroTuring spatial-extrapolation template](https://github.com/Flood-Lab/HydroTuring/blob/main/templates/probe.extrapolation-space.template.yaml): seed-per-catchment design, 730-day warm-up and broad capacity sampling.
- [van Oorschot et al. (2024)](https://doi.org/10.5194/hess-28-2313-2024): root-zone storage estimates reaching hundreds of millimetres, supporting magnitude rather than identity with every model's soil storage datum.
- [Zhang et al. (2006)](https://doi.org/10.5194/hess-10-65-2006): canopy interception storage magnitudes.
- [Zhong et al. (2022)](https://doi.org/10.5194/hess-26-5647-2022): leaf/canopy/land-area conventions; this probe uses catchment-area depths throughout.

These sources support physically plausible scales and definitions. They do not establish the frequency of the eight synthetic catchment categories or identify an AI training-data boundary.

## Review validation

The revised gate and all 24 suite gates pass. Bucket and both FLEX references also pass all 40 independent attribute/climate combinations (seeds 0–39). Each of the seven fault fixtures is tested across the same categories/climates and passes the ordinary catchment-closure cross-check. The capacity fixture additionally checks each storage bound independently, so one store cannot hide a missing check on the other.

Native archived evaluations, with the original model equations and the 60-second container limit:

| Model | Native result |
|---|---|
| `reference_bucket` | PASS (OK) |
| `flex_lumped` | PASS (OK) |
| `flex_topo` | PASS (OK) |
| `sacsma_snow17` | N/A (INCOMPATIBLE) |
| `dhbv2` | N/A (INCOMPATIBLE) |
| `google_flood_forecast` | N/A (INCOMPLETE) |
| `wflow_sbm` | PASS (OK) |
| `summa` | FAIL (VIOLATION) |
| `cwatm` | PASS (OK) |
| `lisflood` | PASS (OK) |

SUMMA passes closure on all 12 seeds but fails companion canopy-storage, ET/PET and runoff-ratio checks. These are not a finding of water leakage or an established solver/probe defect. The general adapter consumes capacity inputs, but native interception parameters do not guarantee the same total-storage ceiling; PET-derived atmospheric inputs also do not guarantee the same potential-ET definition. The [complete SUMMA diagnostic report](SUMMA_DIAGNOSTIC.md) supplies per-seed evidence, distinguishes facts from hypotheses, and gives a concrete adapter/contract follow-up plan. The unchanged [native JSON report](SUMMA_GATE_REPORT.json) accompanies it. No result or threshold was altered for this documentation review.
