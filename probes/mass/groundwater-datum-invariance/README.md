# mass/groundwater-datum-invariance

A common translation of the numerical vertical datum must not alter the
river–aquifer exchange trajectory. `control`, `datum_up` and `datum_down`
share the same stage differences, storage properties and conductance; only the
absolute stage, initial head and elevations are shifted.

The paired `datum_flux_invariance` criterion compares the non-cancelling,
absolute exchange difference against the control gross exchange, with a small
absolute floor. It also requires at least 0.1 mm of activity. Groundwater
balance and signed component checks run independently for every variant.
