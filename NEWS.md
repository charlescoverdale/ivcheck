# ivcheck 0.1.0

* Initial release.
* Adds `iv_kitagawa()`: Kitagawa (2015) variance-weighted Kolmogorov-Smirnov test of instrument validity under binary treatment and discrete instrument.
* Adds `iv_mw()`: Mourifie and Wan (2017) conditional-inequality reformulation; supports covariates.
* Adds `iv_testjfe()`: Frandsen, Lefgren, and Leslie (2023) joint test for exclusion and monotonicity in judge-fixed-effects designs.
* Adds `iv_check()`: one-shot wrapper that runs all applicable tests on a fitted instrumental variable model.
* Adds `iv_power()`: Monte Carlo power simulator.
* S3 dispatch for `fixest` and `ivreg` model objects.
* Optional `modelsummary` glue registered on load.
