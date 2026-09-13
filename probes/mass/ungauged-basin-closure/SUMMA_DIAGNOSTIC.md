# SUMMA: interpretation of the ungauged-basin-closure evaluation

## Finding and scope

The packaged SUMMA adapter (4.0.0-f787fa5.4) receives **FAIL (VIOLATION)** under the current contract. Its reported water budget passes on all 12 gate seeds; failures concern companion canopy-storage, ET/PET and runoff-ratio checks. This is **not a finding of water leakage or an established defect in SUMMA's solver or this probe's closure calculation**. The leading concern is the translation between generic benchmark inputs and native process definitions. Consuming an input does not establish semantic equivalence.

The verdict is retained unchanged. This diagnostic review changes no model equations, forcing, thresholds, capability decisions or archived outcomes. It does not certify that every component is defect-free: attribution of the individual failures still requires native process diagnostics. A successful probe acceptance gate and a submitted model's FAIL answer different questions.

## Evidence and provenance

The evidence is the native 12-seed report supplied alongside this document as [SUMMA_GATE_REPORT.json](SUMMA_GATE_REPORT.json), copied byte-for-byte from the review run. It is an existing evaluation, not a new simulation. The code baseline is upstream `123e9415476d6b7ae83e288add3a971af7cb5e43`, with the PR #74 corrections described in the probe README. Evaluation uses 730 daily spinup steps and a 3,650-day scored window.

| Check | Recorded finding | Interpretation |
|---|---|---|
| Closure | 12/12 pass; maximum precipitation-normalised cumulative residual approximately 3.90e-10 | The reported complete water budget closes |
| State bounds | 11/12 fail; worst seed 1946872492, canopy capacity about 0.704674 mm, excess 0.03222855 mm over 17 rows | Exceeds the current absolute-storage contract; not a rounding error at tolerance 1e-6 mm |
| ET plausibility | 7/12 fail; worst seed 1559595569, ET 8789.567857 mm / PET 5838.299826 mm = 1.505501 | The supplied PET and native ET need a definition/mapping audit |
| Non-degeneracy | 4/12 fail; aggregate-selected seed 545691341: runoff ratio 0.004411, below 0.02; seven-day response correlation 0.088651, above 0.05 | This case fails the runoff-ratio screen, not the response-correlation screen |

These aggregate findings come from different seeds. The ET-worst seed has runoff ratio about 0.2481 and passes non-degeneracy. Seed 38739227 has an even smaller runoff ratio, 0.0007850, while ET/PET is 0.996476 and passes ET plausibility. Thus 0.004411 is not the minimum runoff ratio, and low runoff cannot universally be attributed to ET exceeding PET.

## 1. Canopy capacity is consumed, but its meaning is not a guaranteed total-storage ceiling

The [adapter](../../../models/summa/ht_adapter.py), in the canopy-parameter block, reads `canopy_capacity_mm`. For positive capacity C it uses the configured mixed-forest vegetation and writes:

```text
refInterceptCapRain = C / max(monthly LAI + SAI)
refInterceptCapSnow = C / max(monthly LAI + SAI)
canopy output = scalarCanopyLiq + scalarCanopyIce
```

The documented maximum vegetation area index is 5, so C about 0.704674 mm produces each reference parameter about 0.140935 kg m^-2 per unit vegetation area. This is an actual input mapping, not an ignored field. Soil capacity is separately mapped to column depth through saturated water content; the reported worst bounds failure concerns canopy.

The existing [SUMMA adapter README](../../../models/summa/README.md), particularly “Rain on a freezing canopy” and its later warm-day diagnostics, records source tracing in which liquid above interception capacity drains at a finite rate rather than being immediately clipped. A threshold that initiates drainage is not necessarily an inviolable storage maximum. Continuing inflow plus finite drainage is the leading specific candidate for the small excess here, supported by earlier warm-day diagnostics (drainage coefficient 0.005 s^-1), **but not established for this seed**.

Two additional distinctions matter: separate liquid and snow interception parameters do not by themselves constrain the sum of reported liquid and ice; the snow-interception limit need not bound ice formed by intercepted rain freezing. The generator supplies positive air temperatures, but that alone cannot establish the canopy surface phase. Historical freezing diagnostics must not be presented as the demonstrated cause of this warm-air experiment.

This review independently checks the adapter mapping and criterion code. Native drainage/ice mechanisms above are attributed to the repository's existing source-trace notes; the pinned Fortran source was not independently retrieved in this review. The saved gate report lacks the phase-resolved time series needed to identify the mechanism at the offending rows.

## 2. PET-to-weather reconstruction does not guarantee the same potential-ET definition

The [generator](generate.py) supplies precipitation, air temperature, prescribed PET and dates. PET is a seasonal curve with random variation; the experiment does not supply a complete radiation, humidity and wind record.

In `mock_atmosphere`, the adapter converts PET to a reference net-radiation target through Priestley–Taylor (coefficient 1.26, reference temperature 20 degrees C). It constructs shortwave/longwave inputs with reference albedo 0.23 and supplies relative humidity 70%, wind 2 m/s and elevation-dependent standard pressure. SUMMA then solves its own surface energy balance with its vegetation, surface temperature, albedo and exchange/resistance processes.

The chain is **prescribed PET → reference radiation → reconstructed atmosphere → native surface ET**. It is not an inverse conversion: actual net radiation need not equal the reference target, and the prescribed PET is not imposed as a native ET ceiling. The strongest supported explanation is a reference-surface/atmospheric-definition mismatch in this generic translation. Current evidence cannot apportion the 50.55% excess among radiation, sensible heat, assumed humidity/wind and vegetation processes, or establish that its magnitude is physically reasonable.

[FAO-56, chapter 5](https://www.fao.org/4/X0490E/x0490e0a.htm) distinguishes reference-surface evapotranspiration from vegetation-specific evapotranspiration, whose ratio can exceed one. This supports checking definitions; it does **not** identify HydroTuring PET as FAO reference ET or justify this measured ratio. Earlier sensitivity experiments in the SUMMA README likewise motivate a mapping audit but do not diagnose this seed.

## 3. What this reveals about the common contract

These checks are inherited from the existing [catchment-closure probe](../catchment-closure/probe.yaml) and common [bounds criteria](../../../src/hydroturing/criteria/bounds.py), rather than introduced as SUMMA-specific rules.

- **Input capability is necessary but insufficient.** `uses_static` / `needs_static` establish declared consumption; they cannot verify that a native parameter has the same physical meaning as the requested quantity.
- **A scalar capacity needs a common definition.** A total-storage hard maximum, a phase-specific interception capacity and a drainage threshold are not interchangeable merely because they share units.
- **A PET upper bound needs a consistent surface and atmosphere.** Completing sparse forcing with assumed weather can produce a different demand from the supplied reference index. The resulting mismatch can recur in other probes using this adapter and contract.
- **A runoff floor is a scoped anti-degeneracy assumption, not conservation itself.** Low runoff can coexist with a closed budget and ET below PET. The floor remains enabled here; its applicability to different native surface/soil physics deserves review independently of whether a particular model passes.
- **A contract FAIL is not a solver diagnosis.** The current result correctly records the checks that failed. Its scientific interpretation must also disclose mapping confounds; neither automatic solver blame nor automatic exoneration follows from the label.

No probe-local tolerance relaxation, new framework mechanism or unilateral change to N/A is proposed in this PR. The unresolved issue is common-contract interpretation and adapter validation, not evidence that ungauged water-budget closure is an invalid question.

## 4. Concrete follow-up and decision criteria

1. **Canopy diagnostic, adapter maintainer:** reproduce seed 1946872492 with the same native configuration. Preserve liquid storage, ice storage, canopy temperature, phase-specific capacities, interception and drainage at all offending times. Check units, vegetation-area scaling and the mapping of interval-end output. This separates a translation error from finite drainage or phase semantics.
2. **ET diagnostic, adapter maintainer:** reproduce seed 1559595569 with original forcing. Preserve target and actual net radiation, latent and sensible heat, canopy/surface temperature and ET components. Check both water and energy balances before attributing the discrepancy to one forcing assumption.
3. **Low-runoff diagnostic:** compare seeds 545691341 and 38739227 using precipitation, native ET, runoff, all storage changes and declared exchanges. Test whether low runoff follows the native partitioning under these inputs, rather than assuming the ET-upper-bound failure explains it.
4. **Contract decision, framework and adapter maintainers:** explicitly settle whether canopy capacity is a hard maximum of total liquid plus ice, and whether PET is a reference demand or potential ET for the modeled surface. Confirm how the existing applicability mechanism should handle a model that consumes a parameter without implementing that meaning. Do not relabel results before this decision.
5. **Correction only after diagnosis:** if a unit, area or parameter-write error is found, fix the general adapter and rerun affected closure, state and PET/energy probes. If the native process legitimately uses a different definition, resolve the common contract or documented applicability with maintainers. Do not force compliance by clipping storage/ET, increasing drainage, scaling radiation, changing humidity/wind to fit this seed, or raising this probe's ET limit to 1.6.

For this PR, the requested review decision is whether the unchanged native FAIL, accompanied by this disclosed mapping limitation and follow-up plan, is an acceptable archived benchmark outcome. Completing the native attribution work is a separate general-adapter investigation; this report does not claim it has been completed. The probe's acceptance evidence remains its compatible physical references, targeted negative controls and regression checks.

## Per-seed audit

The following values are extracted directly from the accompanying native report. Bounds values are the criterion's maximum violation; a zero denotes no violation. Non-degeneracy values below are its reported runoff ratio, not response correlation.

| Seed | Closure residual | Bounds status / excess (mm) | ET/PET | Non-degeneracy status / runoff ratio |
|---|---:|---|---:|---|
| 1679270760 | 2.25309e-10 | fail / 0.012763194 | 1.072288 | pass / 0.15376538 |
| 38739227 | 1.71903e-10 | fail / 0.01070357 | 0.996476 | fail / 0.00078502583 |
| 545691341 | 1.73177e-10 | fail / 0.0068755413 | 1.318471 | fail / 0.0044113652 |
| 1052643455 | 3.20955e-10 | fail / 0.012219339 | 1.175137 | pass / 0.064581236 |
| 1559595569 | 2.7681e-10 | fail / 0.01316234 | 1.505501 | pass / 0.24807346 |
| 2066547683 | 3.56669e-10 | fail / 0.013773495 | 0.720278 | fail / 0.0044039468 |
| 426016150 | 1.12979e-10 | fail / 0.02481174 | 1.006134 | pass / 0.19090867 |
| 932968264 | 3.26741e-10 | fail / 0.026163769 | 1.445418 | pass / 0.23454851 |
| 1439920378 | 2.25441e-10 | pass / 0 | 0.585813 | pass / 0.079108235 |
| 1946872492 | 2.30531e-10 | fail / 0.032228553 | 0.984039 | fail / 0.00099074529 |
| 306340959 | 3.90078e-10 | fail / 0.017485047 | 1.444601 | pass / 0.21261855 |
| 813293073 | 2.41738e-10 | fail / 0.00060855482 | 0.639326 | pass / 0.034301617 |
