# Deep expert review of ivcheck

Scope: line-by-line reading of R/kitagawa_core.R, R/mw_core.R,
R/iv_testjfe.R, R/iv_kitagawa.R, R/iv_mw.R, R/iv_power.R, R/utils.R,
R/dispatch.R, R/iv_check.R, tests/testthat/test-known-values.R,
tests/testthat/test-edge-cases.R. Against published formulas in Kitagawa
(2015), Mourifié-Wan (2017), Chernozhukov-Lee-Rosen (2013),
Andrews-Soares (2010), Frandsen-Lefgren-Leslie (2023), Sun (2023).

## Severity legend

- **P0**: silent correctness bug that will give wrong p-values in a
  plausible use case. Must fix before academic review goes out.
- **P1**: UX or fidelity issue that a knowledgeable reviewer will catch
  and object to. Should fix in v0.1.0 if time permits; otherwise
  document prominently and fix in v0.2.0.
- **P2**: defensible design choice that a picky reviewer might question.
  Document the why.

## P0 — must fix before review

### P0.1 — fixest/ivreg dispatch silently ignores fixed effects

**Where**: R/dispatch.R:40
(`y <- as.numeric(stats::model.matrix(object, type = "lhs"))`) and
R/dispatch.R:42 (`d <- as.numeric(fitted(fs) + residuals(fs))`).

**Problem**: a user fitting `feols(y ~ x | firm | d ~ z, data)` has
declared that within-firm variation is what they identify off. The
current dispatch extracts raw y and raw d and then runs the test on
un-demeaned data. I verified this directly: a synthetic model with
strong firm FE gives `var(extracted y) = 1.53` versus
`var(FE-demeaned y) = 1.10` — the FE are simply dropped.

**Impact**: the test then answers a different question (is Z a valid IV
for D in the pooled-across-firm population?) from what the user asked
(is Z a valid IV for D after partialling out firm?). A paper with panel
data will get misleading conclusions.

**Fix (recommended for v0.1.0)**: in `extract_iv_data.fixest`, check
`length(object$fixef_vars) > 0`; if so, **abort** with a clear message:

``` r
if (length(object$fixef_vars) > 0L) {
  cli::cli_abort(c(
    "{.fn ivcheck} v0.1.0 does not support models with fixed effects.",
    i = "The discrete-Z tests of Kitagawa / Mourifié-Wan / FLL operate on the raw
       (Y, D, Z) joint distribution; residualising on fixed effects destroys
       the discrete structure of Z.",
    i = "Workaround: for a FE-IV model, pre-demean (Y, D) within FE cells,
       keep Z discrete, and call the default method on the vectors."
  ))
}
```

Analogous check for `ivreg` (which doesn’t commonly use FE) — skip.

**Alternative for v0.2.0**: stratify the test by FE cell (run the test
within each level and combine p-values via Simes or harmonic mean). This
is a genuine research question, not a bug fix, so out of scope for
v0.1.0.

### P0.2 — iv_testjfe multivalued D with covariates gives conservative σ²

**Where**: R/iv_testjfe.R:200–209 (structural residual construction).

``` r
} else {
    beta_full <- numeric(length(d_vals))
    beta_full[match(d_vals[-1L], d_vals)] <- beta_vec
    idx_map <- match(d_num_raw, d_vals)
    u_i <- y - alpha_hat - beta_full[idx_map]   # RAW indicators, not residualised
}
```

**Problem**: under FWL, the WLS regression of `mu_j` on `P_design`
(which contains the per-judge averages of the **residualised**
level-indicators) produces β̂\_ell that is the coefficient on the
residualised indicators, not the raw indicators. Therefore the correct
structural residual is `u_i = ỹ_i − α̂ − Σ_ell β̂_ell · ind_resid_ell(i)`,
not `u_i = ỹ_i − α̂ − β̂[D_raw_i]`. The two differ by
`β̂ · (D_raw − D_resid)`, which for a linear FE/X partialling is
`β̂ · (a + b·x_i)` within each judge. That is a non-zero within-judge
variation.

**Impact**: `σ̂² = ss_within / df_within` picks up this extra
within-judge variation, so `σ̂²` is inflated, `T_n / σ̂²` is deflated, and
the χ² p-value is larger than it should be (test is conservative). Monte
Carlo null size drifts below nominal 5% when M ≥ 2 and x has non-trivial
variation within judges.

For binary D (M = 1) with x, `d_num` is overwritten to the residualised
version on line 137, so `beta_vec * d_num` already uses residualised D.
No bug for M = 1.

**Fix**: in the multivalued-with-x branch, save the per-observation
residualised indicators during `P_design` construction and use them in
`u_i`. Rough sketch around R/iv_testjfe.R:170–184:

``` r
ind_resid_mat <- matrix(0, nrow = n, ncol = M)  # new
for (m in seq_len(M)) {
  dv <- d_vals[m + 1L]
  ind <- as.numeric(d_num_raw == dv)
  if (!is.null(x)) {
    design <- cbind(1, x_mat)
    ind <- as.numeric(stats::residuals(stats::lm.fit(design, ind)))
  }
  ind_resid_mat[, m] <- ind                    # save for structural residuals
  P_design[, m + 1L] <- vapply(judges, ...)
}
```

Then at lines 200–209, the multivalued branch becomes:

``` r
} else {
    u_i <- y - alpha_hat
    for (m in seq_len(M)) u_i <- u_i - beta_vec[m] * ind_resid_mat[, m]
}
```

This restores the FWL identity.

**Test to add**: under H0 for multivalued D (say M = 2) with x varying
within judges, the χ² null size should still be at nominal 5%. A quick
MC with n_reps = 200 and K = 10 judges catches this if broken.

### P0.3 — validate_treatment_discrete remaps D, binding labels become wrong

**Where**: R/utils.R:46–50.

``` r
if (!all(vals == seq(0L, k - 1L))) {
  mapping <- stats::setNames(seq(0L, k - 1L), vals)
  d_num <- as.numeric(mapping[as.character(d_num)])
}
```

**Problem**: when a user passes D ∈ {1, 2, 3} (a perfectly legitimate
input), the validator silently remaps to {0, 1, 2}. Kitagawa core then
reports `binding$direction = "D <= 0"` when the user’s data had `D = 1`
as the smallest level. The user reads “D \<= 0” and thinks the violation
is in the no-treatment cell, when really it is in D = 1.

**Impact**: cosmetic but misleading. A PhD student running on real data
with integer-coded D will notice this immediately.

**Fix**: return the remap table from the validator and use the original
`d_vals` in `direction_info()` for `binding$direction`. Or simpler:
don’t remap unless necessary (only remap when values are not already a
contiguous-from-0 integer sequence).

Actually simplest: in `kitagawa_core_test`, when constructing labels,
use the *original* `d_vals` rather than the remapped internal values.
This requires the validator to pass through the original levels, which
it currently does not. Low-risk addition: have
`validate_treatment_discrete` return a list
(`list(d_num = ..., d_vals = ...)`) instead of a vector.

## P1 — important, should address in v0.1.0

### P1.1 — iv_mw only uses the first column of x_mat

**Where**: R/mw_core.R:26 (`x_str <- x_mat[, 1]`).

**Problem**: the documentation says “the test is conditional on the
first numeric column of `x`” but users passing
`x = cbind(age, educ, race_dummies)` almost certainly expect
conditioning on all of those. Silently dropping all-but-first is a
fidelity concern.

**Fix (recommended)**: abort with an explicit error when
`ncol(x_mat) > 1`:

``` r
if (ncol(x_mat) > 1L) {
  cli::cli_abort(c(
    "{.fn iv_mw} v0.1.0 supports a single conditioning covariate.",
    i = "Multivariate conditioning requires a tensor-product basis (planned
       for v0.2.0). Workaround: condition on the covariate that most plausibly
       drives heterogeneity in compliance."
  ))
}
```

A loud error is better than a silent drop.

### P1.2 — mw_core BtB_inv silent fallback to diag(1e-4, P)

**Where**: R/mw_core.R:58–60.

``` r
BtB_inv <- tryCatch(
  solve(crossprod(Bk) + diag(1e-8, P)),
  error = function(e) diag(1e-4, P)
)
```

**Problem**: if the basis matrix is (near-)rank-deficient in some Z
cell, the [`solve()`](https://rdrr.io/r/base/solve.html) fails and we
silently substitute a tiny diagonal matrix. The downstream projection
`B_at_xg %*% betas` is then essentially zero, which gives `diff ≈ 0` and
`T_n ≈ 0`. The test falsely “passes” with the user having no idea why.

**Fix**: abort with a diagnostic:

``` r
BtB_inv <- tryCatch(
  solve(crossprod(Bk) + diag(1e-8, P)),
  error = function(e) {
    cli::cli_abort(c(
      "Series basis is rank-deficient in Z = {z_levels[k]} cell (n = {length(idx)}).",
      i = "Try reducing {.arg basis_order} (currently {basis_order}) or ensuring
         more variation in {.arg x} within each Z cell."
    ))
  }
)
```

### P1.3 — No test reproducing a published headline number

**Where**: tests/testthat/test-known-values.R.

**Problem**: all tests are self-consistent constructions. We never
assert that `iv_kitagawa(card$lwage, card$college, card$near_college)`
matches the rejection that Kitagawa’s 2015 paper reports for Card. A
reviewer will ask “does this package reproduce the Kitagawa / MW / FLL
papers’ numbers?” and we don’t have a direct answer in the test suite.

**Fix**: add to test-known-values.R:

``` r
test_that("iv_kitagawa reproduces the Card (1995) binary-college rejection", {
  skip_on_cran()
  data(card1995)
  set.seed(40000)
  r <- iv_kitagawa(card1995$lwage, card1995$college,
                   card1995$near_college, n_boot = 500, parallel = FALSE)
  # Kitagawa (2015, Table 2) reports rejection at 5% for this specification
  expect_lt(r$p_value, 0.05)
  expect_gt(r$statistic, 1.5)
})
```

Even one such test is a strong signal to a reviewer that the
implementation is faithful.

### P1.4 — Claim at iv_kitagawa.R:63 is ambiguous

**Where**: R/iv_kitagawa.R:63.

    T_n = sqrt(n_low * n_high / n) * max [P([y, y'], d | z_low) - P([y, y'], d | z_high)]^+ / sigma_hat

**Problem**: `n` in the denominator is ambiguous — is it total sample
size or pair total `n_low + n_high`? Kitagawa’s equation 2.1 uses the
pair total. The code matches (R/kitagawa_core.R:196–197 uses
`n_z_eff[k_low] + n_z_eff[k_high]`). The docstring should be explicit.

**Fix**: clarify in iv_kitagawa.R:63:

    T_n = sqrt(n_low * n_high / (n_low + n_high)) * max [P([y, y'], d | z_low) - P([y, y'], d | z_high)]^+ / sigma_hat

### P1.5 — iv_mw.default validates binary only, silently blocking Sun (2023) path

**Where**: R/iv_mw.R:101 (`d_num <- validate_binary(d, "d")`).

**Problem**: while Kitagawa/Sun is wired through
`validate_treatment_discrete`, iv_mw still hard-gates to binary D. So
for multivalued D the user must call `iv_kitagawa` directly. Fine as a
v0.1.0 scope, but should be documented in the iv_mw docstring that
multivalued requires iv_kitagawa (Sun 2023 mode).

### P1.6 — Bootstrap SE in MW is “frozen” from observed but recomputed per pair

**Where**: R/mw_core.R:167–169.

``` r
se_low  <- SE_on_grid[, , k_low,  d_idx]
se_high <- SE_on_grid[, , k_high, d_idx]
SE_pair <- pmax(sqrt(se_low^2 + se_high^2), se_floor)
```

**Problem (minor)**: the CLR bootstrap canonically freezes the
studentisation at the observed plug-in. The code does freeze it (uses
`SE_on_grid`, the observed-data matrix, not a bootstrap re-estimate).
Good. But the `pmax(..., se_floor)` application differs between the
observed statistic (line 95, uses `+ se_floor` additive) and the
bootstrap (line 118 and 169, uses `pmax(..., se_floor)`). These floors
differ: `sqrt(x) + se_floor` vs `max(sqrt(x), se_floor)`.

**Impact**: in the ordinary case where `sqrt(qf) > se_floor`, both are
identical. Where `sqrt(qf) < se_floor`, the observed uses
`sqrt(qf) + se_floor` (roughly 2·se_floor) while the bootstrap uses
`se_floor`. This asymmetry biases the bootstrap critical values
slightly. Negligible in most cases but a fastidious reviewer might flag
it.

**Fix**: harmonise at R/mw_core.R:95 — change
`sqrt(pmax(qf, 0)) + se_floor` to `pmax(sqrt(pmax(qf, 0)), se_floor)`.

## P2 — design choices, document the why

### P2.1 — Fixed cubic basis in mw_core (basis_order = 3, no CV)

**Where**: R/mw_core.R:11 function signature.

**Why it is defensible**: Mourifié-Wan’s simulations and empirical
application use a low-order polynomial; CLR inference rates apply for
any fixed basis under mild conditions. Data-driven basis selection
requires a separate inference theory (cross-validation alters the null).

**Action**: add one sentence to the docstring explaining the choice and
pointing users to `basis_order` if they want to tweak.

### P2.2 — Asymmetric grid trimming (y: 0.02–0.98, x: 0.1–0.9)

**Where**: R/mw_core.R:38, 41.

**Why**: y grid is denser at edges because the CDF is the quantity being
estimated; x grid trims more because polynomial basis misbehaves at the
tails of the x distribution (extrapolation risk).

**Action**: add a one-line comment in the code.

### P2.3 — `iv_power` DGP tests only direct-effect exclusion failures

**Where**: R/iv_power.R:102–104.

**Problem**: there is no “defier share” option, so `iv_power` cannot
benchmark power under pure monotonicity failures. For a paper claiming
to test both exclusion and monotonicity, this is a gap in evidence.

**Fix (recommended for v0.2.0)**: add
`violation_type = c("exclusion", "defier", "both")` with clear
parameterisations for each. For v0.1.0, document that the default DGP is
a D-specific exclusion violation.

### P2.4 — Grid size defaults (y_grid_size = 50, x_grid_size = 20)

**Where**: R/kitagawa_core.R:21, R/mw_core.R:13–14.

**Why defensible**: Kitagawa’s Theorem 2.1 guarantees the sup is
attained at observed Y values; the quantile grid is a dense
approximation. 50 points is more than enough for the sup to hit within
0.01 of the true value on standard DGPs.

**Action**: a single sentence in the docstring about the quantile-grid
approximation would help.

### P2.5 — Rademacher multiplier only (no Mammen)

**Where**: R/kitagawa_core.R:264
(`W <- sample(c(-1, 1), n, replace = TRUE)`).

**Why defensible**: Rademacher is the standard choice for binary-outcome
indicators (mean 0, variance 1, bounded). Mammen’s (1993) two-point
distribution would give slight finite-sample improvements in
skewed-residual settings. Not a bug.

### P2.6 — Validation does not remap original D values back in output

See P0.3 above.

## Monte Carlo holes (tests you should add)

1.  **Null size of iv_testjfe with multivalued D + covariates.**
    Currently only binary D is covered; P0.2 would show up here.

2.  **Power of iv_testjfe against a pure monotonicity violation.**
    Simulate data where some judge pair has a defier share; verify
    iv_testjfe rejects.

3.  **Size of iv_mw on the conditional-X path under H0 with varying
    basis_order.** Ensure size stays at nominal 5% at basis_order = 2,
    3, 4.

4.  **Identity of iv_mw(x = NULL) and iv_kitagawa on identical seeds.**
    Already in test-known-values.R:37. Good.

5.  **Card 1995 reproducibility.** See P1.3.

6.  **Rejection under a weak-instrument pathology.** With p_low = p_high
    (no first stage), every test should have size ≤ nominal (the
    implementation should not produce runaway statistics). Worth a
    direct test at p_low = p_high = 0.5.

7.  **Behaviour on a single-instrument-cell case.** What happens if all
    Z take one value? Currently `validate_discrete` aborts on K \< 2.
    Good.

## Literature fidelity (recommend for paper text)

- **Sun (2023) section 4 claim**: ivcheck’s Kitagawa core uses the
  cumulative-tail form of Sun for ordered multivalued D. Sun also
  discusses an unordered-D extension (section 5) which we do **not**
  implement. Document this in limitations.

- **Mourifié-Wan (2017, Theorem 2)**: the CLR approximation error is of
  order `n^{-1/2} \log(n)^{-1/2}` under their Assumption A2. Our default
  `n_boot = 1000` keeps bootstrap MC noise below this rate for n ≤
  10,000. Acceptable.

- **FLL (2023) section 3**: the FLL paper derives a finite-sample
  correction term when K is small (K = 3, 4). We currently only have the
  asymptotic χ². Document the MC-validated K = 20 range in the paper.

## Prioritised action plan

### Before academic circulation (v0.1.0)

1.  **Fix P0.1** (FE abort in fixest dispatch) — 10 lines in
    R/dispatch.R.
2.  **Fix P0.2** (structural residuals in iv_testjfe multivalued + x) —
    15 lines in R/iv_testjfe.R.
3.  **Fix P0.3** (D remap label confusion) — 20 lines across R/utils.R
    and R/kitagawa_core.R.
4.  **Fix P1.1** (multivariate x abort in iv_mw) — 6 lines in
    R/mw_core.R.
5.  **Fix P1.2** (BtB_inv informative abort) — 8 lines in R/mw_core.R.
6.  **Add P1.3** (Card 1995 reproduction test) — 10 lines in
    tests/testthat/test-known-values.R.
7.  **Fix P1.4** (docstring ambiguity) — 1 line in R/iv_kitagawa.R.
8.  **Fix P1.6** (bootstrap/observed se_floor asymmetry) — 1 line in
    R/mw_core.R.
9.  **Add MC hole 1** (iv_testjfe multivalued + x size test) — 20 lines
    in tests/.
10. **Document limitations section in README**: (a) no FE support
    yet, (b) iv_mw single-covariate, (c) Sun section 5 (unordered
    multivalued) not implemented.

Runtime: 2-3 hours for all fixes, tests rerun, paper replicate.R rerun.

### Defer to v0.2.0

- Tensor-product basis for multivariate x in iv_mw.
- Stratified-by-FE test for panel models.
- `violation_type` argument in iv_power.
- Sun unordered-multivalued extension.
- Huber-Mellace (if user demand exists).

## Verification after fixes

1.  `devtools::check(document = TRUE, cran = TRUE)` → 0/0/0 + 2 benign
    notes.
2.  Run `paper/replicate.R` end-to-end; confirm Card still rejects (~ p
    \< 0.01).
3.  Re-run the 24-config MC for iv_kitagawa null size; all ≤ nominal 5%.
4.  Add one MC for iv_testjfe multivalued + x null size; should now be ≤
    nominal.
5.  Regenerate paper PDFs; update any numerical claim that changed.

## Bottom line

The package is in better shape than most first-submission
applied-metrics R packages. The three P0 fixes are the only genuine
correctness risks; each is a focused 10–20 line change. Everything else
is either defensible with a docstring line or genuinely out of scope for
v0.1.0.

A well-informed PhD student reviewing the package today could construct
two “gotcha” cases: (a) a panel-FE IV model where the user would
reasonably expect FE to be partialled out, (b) a multivalued-D JFE model
with covariates where the conservative σ² inflates p-values. Fixing P0.1
and P0.2 removes both.
