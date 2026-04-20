# CRAN submission comments - ivcheck 0.1.0

## New submission

This is a new package providing tests for the identifying assumptions of
instrumental variable models (the local exclusion restriction and
monotonicity conditions required for LATE identification). The package
implements Kitagawa (2015, Econometrica), Mourifie and Wan (2017, Review
of Economics and Statistics), and Frandsen, Lefgren, and Leslie (2023,
American Economic Review), plus a one-shot wrapper that runs every
applicable test on a fitted IV model.

## R CMD check results

Local `R CMD check --as-cran`: 0 errors, 0 warnings, 1 NOTE.

The single NOTE ("unable to verify current time") is an environmental
issue on the author's machine and does not reproduce on CRAN
infrastructure.

## Test suite

110+ testthat expectations covering structure, invariants, known-value
cases, edge cases, published-number reproduction (Card 1995), and
end-to-end S3 dispatch against `fixest` and `ivreg` IV models. All tests
run offline; no network-dependent tests.

## Notes for the reviewer

- All DOIs in DESCRIPTION and `@references` have been verified against
  the CrossRef API (`https://api.crossref.org/works/<doi>` returns 200
  for every DOI used).
- Any doi.org links in README.md that `urlchecker::url_check()` flags
  as 403 are a known publisher-side bot-block false positive (AEA and
  MIT Press routinely 403 automated HEAD requests). The CrossRef
  verification confirms the DOIs resolve.
- Package uses `parallel::mclapply` for the multiplier bootstrap with
  an explicit cap at 2 cores when `_R_CHECK_LIMIT_CORES_` is set, per
  CRAN parallel-computing policy.
- No data is downloaded at runtime. A bundled `card1995` dataset (3,003
  rows, cleaned extract of Card 1995 National Longitudinal Survey of
  Young Men data) is included under `data/` via `usethis::use_data()`
  with xz compression. Provenance documented in `data-raw/`.
- S3 dispatch for `fixest` and `ivreg` uses `.onLoad` to conditionally
  register methods; no hard dependency on either package. Both are
  Suggests.

## Downstream dependencies

None. This is a new package.
