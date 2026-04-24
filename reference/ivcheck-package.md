# ivcheck: Tests for Instrumental Variable Validity

Implements tests for the identifying assumptions of instrumental
variable models, the local exclusion restriction and monotonicity
conditions required for local average treatment effect identification.
Covers Kitagawa (2015)
[doi:10.3982/ECTA11974](https://doi.org/10.3982/ECTA11974) , Mourifie
and Wan (2017)
[doi:10.1162/REST_a_00622](https://doi.org/10.1162/REST_a_00622) , and
Frandsen, Lefgren, and Leslie (2023)
[doi:10.1257/aer.20201860](https://doi.org/10.1257/aer.20201860) .
Includes a one-shot wrapper that runs all applicable tests on a fitted
instrumental variable model. Dispatches on 'fixest' and 'ivreg' model
objects.

## See also

Useful links:

- <https://github.com/charlescoverdale/ivcheck>

- Report bugs at <https://github.com/charlescoverdale/ivcheck/issues>

## Author

**Maintainer**: Charles Coverdale <charlesfcoverdale@gmail.com>
