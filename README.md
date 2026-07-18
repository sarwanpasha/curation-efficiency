# Efficiency of Categorical Curation is Partially Identified

Theory, methods, and analysis code for a belief-martingale treatment of
efficiency in vintaged categorical gold standards, with an application to
variant reclassification in ClinVar.

## Overview

Curated gold-standard labels — variant pathogenicity calls, gene
annotations, the versioned resources underpinning genomic benchmarking —
are revised at every release. Asking whether that revision is *efficient*
(each release the curator's best estimate, so changes are unpredictable)
turns out to be harder than it looks.

**The naive formalization is vacuous.** The only `{0,1}`-valued martingale
is a constant path. A categorical label that ever changes therefore cannot
satisfy a martingale-difference efficiency null, so any test of that null on
moving labels is unfalsifiable. This is not a technicality: in the baseline
comparison here, a fixed-effects predictability test rejects efficiency in
97.5% of samples drawn from a bounded mean-reverting process, and a
variance-ratio test rejects a genuinely efficient process at more than three
times its nominal level.

**The fix is to relocate efficiency to a latent belief.** A curator's
posterior about a fixed truth is a `[0,1]`-valued martingale by
construction; the observed categorical label is a *threshold report* of it.
Efficiency is a property of the unobserved belief, and the label is a
quantization. This yields two exact, nuisance-free constraints on
observables:

- **Flip-rate inequality.** `P(flip | p) <= 2 * min(p, 1-p)`. A confident
  belief cannot flip its label often.
- **Sojourn absorption invariant.** `P(run is terminal | entry belief p) =
  |2p - 1|`, *independent of the belief's volatility*. Raw run lengths are
  contaminated by how fast the belief moves; the absorption probability is
  not.

**Efficiency is only partially identified.** For a range of observed flip
behaviours there exists an efficient belief generating them, so no
label-path statistic separates efficiency from inefficiency. Identification
requires an exogenous covariate. This makes the limits on any test's power a
theorem about the problem rather than a defect of a procedure.

## Results at a glance

| | result | scope |
|---|---|---|
| Size control | only procedure of four holding nominal size under both nulls (0.042, 0.050) | validity, not power |
| Sensitivity | graded power 0.442 → 1.000 across a narrow effect-size range | separate excess-rate sweep |
| ClinVar concordance | benign resolutions cross the benign frequency threshold at **16.5x** the pathogenic rate; **4.6x** more likely to be observed in gnomAD | powered, descriptive |
| Longitudinal test | effective n = 177; covariate time-constant for 99.2% of resolvers | **non-identification**, not a null of effect |

Full numbers, with scope tags and the validity guards applied, are in
[`results/verified_results.md`](results/verified_results.md).

## Repository layout

```
src/
  sojourn_calib.py      simulation-calibrated sojourn test + size/power grids
  aggregate_calib.py    pools array-task JSON output into size/power tables
  prescreen.py          admissibility diagnostics (A/B gates, covariate
                        classification, panel typing)
  battery.py            test battery: naive FE, variance ratio, portmanteau,
                        and the corrected test when an admissible covariate
                        is confirmed
scripts/
  run_calibration.sbatch  SLURM array driver for the calibration grid
results/
  verified_results.md   all computed numbers with scope tags
paper/
  main.tex, references.bib, main.pdf
```

## Quick start

```bash
pip install -r requirements.txt

# single-seed calibration run (fast smoke test)
python src/sojourn_calib.py 0 50

# full grid on a cluster: 10 array tasks x 300 replications
sbatch scripts/run_calibration.sbatch
python src/aggregate_calib.py          # pools results/calib_result_*.json
```

`sojourn_calib.py` takes `<seed> <replications>` and emits JSON with size
across belief-confidence levels and a fine near-zero power grid. It is pure
NumPy and CPU-only; no GPU is used.

## Method notes

Two design choices carry most of the statistical weight.

**The belief proxy must not come from the label path.** The test estimates
`p_hat(W) = P(terminal label = 1 | W)` from *terminal outcomes* and the
exogenous covariate, never from the flip behaviour under test. A
persistence-based proxy would be circular: a stable label would imply a
confident belief, which would imply flips should be rare, which is what was
observed in the first place.

**The reference must carry the same estimation noise as the statistic.**
The efficient reference law is simulated across a range of belief
volatilities and taken conservatively, so the volatility nuisance cannot
re-enter. Critically, each simulated replicate estimates *its own* belief
proxy exactly as the observed statistic does. Without this debiasing the
mid-belief regime is over-sized, because the reference is otherwise
dominated by proxy-estimation noise absent from a fixed-proxy null.

## Reproducing the ClinVar analysis

The ETL and join code is cluster- and site-specific and is not included
here. To reproduce:

1. **Assemble vintages.** Fetch dated ClinVar release VCFs and build a long
   panel keyed by `(chrom, pos, ref, alt)` on GRCh38, one row per
   (variant, vintage).
2. **Label.** Map clinical significance to confident (pathogenic / likely
   pathogenic / benign / likely benign) versus uncertain (uncertain
   significance / conflicting). Check uncertainty and conflict markers
   *before* confident markers: `Conflicting_interpretations_of_pathogenicity`
   contains the substring `pathogenic` and is otherwise misclassified.
3. **Join the covariate.** Query versioned gnomAD sites files by remote
   indexed access (no bulk download required). For each vintage use the
   latest release *strictly preceding* the vintage date.
4. **Keep the covariate definition constant across releases.** Compute the
   population-maximum 95% filtering allele frequency as the max over the
   same continental population groups in every release. Field naming drifts
   between schemas (underscore- versus dash-separated), and mixing a formal
   filtering AF in one release with a raw AF in another creates a
   definitional step-down that manufactures spurious threshold crossings.
5. **Prescreen before testing.** Run `src/prescreen.py` to obtain
   admissibility gates, panel type, and covariate classification. Confirm a
   covariate as exogenous only on constructional grounds — that it is
   computable without reading the label. Review-status and submitter-count
   fields fail this test: they describe the classification process itself.
6. **Run the battery.** `src/battery.py` runs the naive, variance-ratio, and
   portmanteau tests unconditionally, and the corrected test only when a
   confirmed admissible covariate exists.

## Interpreting results honestly

Three distinctions are load-bearing and easy to blur.

**Size validity is not power.** The baseline table establishes that standard
predictability tests do not control size on bounded categorical labels, so
their rejections cannot be read as evidence of inefficiency. It does not
establish that the corrected test is more sensitive; in those regimes it is
among the least. Sensitivity is a separate measurement.

**Concordance is not causal discovery.** The 16.5-fold frequency separation
shows that resolution direction is strongly associated with prior population
frequency. The mechanism is the ACMG benign criteria, which curators
consult. The finding quantifies a documented criterion operating at scale.

**Non-identification is not a null of effect.** An effective n of 177 does
not show that frequency fails to drive resolution. It shows that the
longitudinal causal effect is not identifiable, because the covariate is
time-constant for almost every variant. Denser population references would
raise the effective sample and bring the test within reach.

## Requirements

Python 3.9+, NumPy, pandas, SciPy. See `requirements.txt`. A SLURM cluster
is optional and used only to parallelize the calibration grid.

## License

MIT. See [`LICENSE`](LICENSE).
