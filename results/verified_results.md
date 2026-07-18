# Verified numerical results

All values below were computed, not estimated or recalled. Each is tagged
with the scope under which it is valid. Figures and tables in the paper are
generated directly from these numbers.

## 1. Sojourn test calibration (synthetic)

Pooled over 10 seeds x 300 replications = 3000 total, n = 800 items,
T = 16 vintages, alpha = 0.05, debiased reference.

### Size (efficient data generating process)

| belief confidence w | empirical size | standard error |
|---|---|---|
| 0.55 (near threshold) | 0.0000 | 0.0000 |
| 0.60 | 0.0000 | 0.0000 |
| 0.75 (mid) | 0.0653 | 0.0031 |
| 0.90 (confident) | 0.0080 | 0.0014 |

### Power (injected predictable inefficiency, w = 0.75)

| excess rate | rejection rate |
|---|---|
| 0.01 | 0.442 |
| 0.03 | 0.998 |
| 0.05 | 1.000 |
| 0.08 | 1.000 |
| 0.12 | 1.000 |
| 0.20 | 1.000 |

Scope: the power curve is graded and separated from size. The mid-belief
size of 0.065 is the honest measured value; it is not 0.05.

## 2. Baseline comparison (synthetic)

alpha = 0.05, 120 replications per cell, within-locus permutation null.
Procedures: `naive` = fixed-effects regression of increment on lagged
level; `corrected` = simulation-calibrated sojourn test; `vr` =
Lo-MacKinlay variance ratio; `el` = Escanciano-Lobato automatic
portmanteau.

| regime | role | naive | corrected | vr | el |
|---|---|---|---|---|---|
| efficient | size | 0.000 | **0.042** | 0.158 | 0.017 |
| iid_bounded | size | 0.075 | **0.050** | 0.033 | 0.025 |
| momentum | rejection | 0.308 | 0.017 | 0.083 | 0.308 |
| reversion | rejection | 0.975 | 0.042 | 0.000 | 0.292 |
| trend | rejection | 0.000 | 0.000 | 0.267 | 0.317 |

Scope: this is a **validity** result. The corrected test is the only one of
the four holding nominal size under both nulls. It is **not** a power win:
in the three alternative regimes the corrected test is among the least
sensitive. Sensitivity is characterized separately by the excess-rate sweep
in section 1.

The `naive` rejection rate of 0.975 under `reversion` is the vacuousness
result realized empirically: a bounded categorical label that mean-reverts
is exactly the object for which the martingale-difference null is
unfalsifiable.

## 3. ClinVar x gnomAD (real data)

Panel: 14 dated ClinVar release vintages, 2019-2024, GRCh38.
Label: confident (pathogenic / likely pathogenic / benign / likely benign)
versus uncertain (uncertain significance / conflicting).
Covariate: lagged population-maximum 95% filtering allele frequency,
computed identically in every gnomAD release as the max over the same five
continental population groups.

| quantity | value |
|---|---|
| unique variants in panel | 2,679,330 |
| allele-frequency entries returned by join | 4,099,357 |
| gnomAD releases used | v3.0, v3.1, v4.0 |
| global coverage (filtering AF > 0) | 0.3256 |

### Resolution channels

| | benign-resolving | pathogenic-resolving | ratio |
|---|---|---|---|
| resolution events n | 29,754 | 8,297 | — |
| fraction with non-zero filtering AF | 0.5758 | 0.1247 | 4.6x |
| fraction crossing benign threshold before resolution | 0.005948 | 0.000361 | 16.5x |

### Identification diagnostics (benign channel)

| quantity | value |
|---|---|
| non-constant lagged covariate | 240 of 29,754 (0.81%) |
| crossing benign threshold before resolution (effective n) | 177 |
| post-2023 v4-eligible benign resolvers (look-ahead bound) | 7,779 |

Scope: the concordance (4.6x, 16.5x) is **powered and descriptive**. It
quantifies association between prior population frequency and resolution
direction; it is not a causal claim, and the mechanism is the documented
ACMG benign criteria, which curators apply.

The effective n of 177 is a **non-identification** result, not a null of
effect. The longitudinal test is not identified because the exogenous
covariate is time-constant for 99.2% of benign-resolving variants: the
population reference refreshes on a release cadence coarser than the label
revises.

## Validity guards applied

1. **Look-ahead guard.** Effective n (177) is bounded by the count of
   post-2023 v4-eligible resolvers (7,779), confirming no crossing is
   counted that fails to strictly precede its resolution vintage.
2. **Covariate consistency guard.** An earlier run used the formal
   filtering AF for v4 but fell through to raw AF for v3, producing a
   definitional step-down at the v3-to-v4 boundary that would have made
   threshold crossings artifactual. The covariate is now the same
   population-max faf95 estimator in every release, reading both
   underscore- and dash-separated INFO field names across release schemas.
3. **Panel dimension check.** Resolver counts (29,754 / 8,297) reproduce
   across independent runs and scaffold sources.
