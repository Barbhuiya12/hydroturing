# momentum/stage-discharge-monotonic

A staff gauge is the oldest instrument in hydrology and it imposes two
constraints that no water budget does. This probe asks a model that reports
the stage of its channel to satisfy both.

The first is monotonicity. Stage is single-valued in discharge: read the reach
at the same flow twice and the gauge reads the same level, and a higher flow
reads a higher level. A model that computes its stage from something other
than the store it claims to describe — a running maximum, a smoothed flow, a
depth from a different variable entirely — decouples the gauge from the flow
and fails this.

The second is the loop. A flood wave steepens as it arrives, so at the same
discharge the rising limb sits *lower* than the falling limb. The sign of that
loop is physical rather than conventional: a rating read backwards, high while
the flow is still climbing and low once it is leaving, inverts it. A
monotonicity check alone cannot see that, which is why there are two criteria.
The limbs are paired on the channel store — the water in transit — which is
the quantity the hysteresis is about; a reach that holds nothing in transit is
read on its discharge instead.

Four years of temperate weather with five multi-day storms, in a single reach
whose width, slope and roughness are handed to every model through the static
file. `stage` is a diagnostic the probe asks the model to report in metres. It
is deliberately not one of the storages the suite differences its budget over:
a stage is a reading, not a volume, and summing it into `reported_states`
would corrupt every closure test in the suite.

| Criterion | Asserts |
| --- | --- |
| `rating_monotonic` | stage does not fall against the running maximum as discharge rises |
| `rating_loop` | at matching discharge, the rising limb sits below the falling limb |
| `non_degenerate` | discharge and stage vary with the weather |

Must-fail: `reference_rating_drift`, whose stage is a running maximum of its
own discharge; `reference_rating_inverted`, which reads the loop backwards;
and `reference_flat_stage`, whose gauge reports a constant.
