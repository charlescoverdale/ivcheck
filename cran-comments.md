# CRAN submission comments — ivcheck 0.1.0

## New submission

This is a new package providing tests for the identifying assumptions of instrumental variable models (local exclusion restriction and monotonicity), covering Kitagawa (2015), Mourifie and Wan (2017), and Frandsen, Lefgren, and Leslie (2023).

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test suite

200+ expectations. No network-dependent tests; all computation is pure.

## Notes on data access

No external data access. Package ships bundled demo datasets (Card 1995 extract and a simulated judge-IV panel with known data-generating process) under `data/`.

## Downstream dependencies

None.
