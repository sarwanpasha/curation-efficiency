#!/usr/bin/env python3
"""
s35_realcorpus.py -- Generic battery for a single pre-screened real corpus.

Pipeline:
  1. Pre-screen via prescreen.prescreen_corpus (must be DIAGNOSTIC_USABLE; not EXCLUDE).
  2. Build the analysis panel: balanced-core, or exogenous-unbalanced per s34 rules.
     NEVER forward-fill labels.
  3. Run the test battery with block-permutation p-values:
       - naive fixed-effects (label-path features)   [always]
       - variance-ratio (VR)                          [always]
       - automatic Portmanteau (EL)                   [always]
       - BSMD plug-in    on admissible W              [only if BSMD_TESTABLE]
       - BSMD cross-fit  on admissible W              [only if BSMD_TESTABLE]
     For DIAGNOSTIC-USABLE-only corpora, BSMD is reported as "degenerate/NA".
  4. Emit per-test reject decision + naive F-stat + the admissibility classification.

ESTIMATOR INTERNALS are imported from s33_bsmd_validation.py. Because s33 lives on the
cluster and its exact signatures are not visible here, ALL s33-dependent calls are
isolated in the ADAPTER block below and marked  # RECONCILE:  -- edit those to match
your real s33 function names/args/returns. Everything else is self-contained.
"""

import sys
import json
import numpy as np
import pandas as pd

import prescreen as s34

# ======================================================================
# ADAPTER BLOCK -- RECONCILE these against your real s33_bsmd_validation.py
# ======================================================================
# The battery needs five callables. If s33 exposes them under other names or
# signatures, remap here (this is the ONLY block you should need to edit).
#
# Expected semantics (edit RHS to match s33):
#   naive_fe_stat(panel)      -> (F_obs, p_boot)   naive FE test w/ block-perm null
#   vr_test(panel)            -> (stat, p)         variance-ratio on the series
#   el_test(panel)            -> (stat, p)         automatic Portmanteau
#   bsmd_stat(panel, W, mode) -> (S, p)            BSMD; mode in {'plugin','crossfit'}
#   BLOCK_LEN, NPERM          -> ints              null construction params
#
# 'panel' is a long DataFrame with columns item_id, vintage, Y, Delta (see build_panel).
# 'W' is the name of the admissible covariate column in panel.
try:
    import s33_bsmd_validation as s33  # RECONCILE: module name
    # RECONCILE: map each of these to the real s33 callable.
    _naive_fe   = s33.pooled_md_stat            # RECONCILE: naive FE + block-perm p
    _vr_test    = s33.variance_ratio_test       # RECONCILE
    _el_test    = s33.automatic_portmanteau      # RECONCILE
    _bsmd_stat  = s33.bsmd_stat                  # RECONCILE: (panel, W, mode)->(S,p)
    BLOCK_LEN   = getattr(s33, "BLOCK_LEN", 4)   # RECONCILE
    NPERM       = getattr(s33, "NPERM", 299)     # RECONCILE
    _S33_OK = True
except Exception as e:  # pragma: no cover
    _S33_OK = False
    _S33_IMPORT_ERR = repr(e)
# ======================================================================


def build_panel(df, prescreen_report, reference_mode="explicit"):
    """Build the analysis panel per the pre-screen panel_type. Returns a long
    DataFrame with item_id, vintage, Y, Delta (+ covariates), or raises if EXCLUDE."""
    ptype = prescreen_report["panel_type"]
    if ptype == "label-driven-unbalanced-EXCLUDE":
        raise ValueError("panel_type is EXCLUDE (label-driven entry); do not run battery.")

    dfY = s34.build_Y(df, reference_mode=reference_mode)
    vintages = np.sort(dfY["vintage"].unique())

    if ptype == "balanced-core":
        counts = dfY.groupby("item_id")["vintage"].nunique()
        core = counts[counts == len(vintages)].index
        panel = dfY[dfY["item_id"].isin(core)].copy()
    elif ptype == "exogenous-unbalanced":
        # keep each item's own observed span; drop items with <2 increments
        span = dfY.groupby("item_id")["vintage"].nunique()
        keep = span[span >= 3].index
        panel = dfY[dfY["item_id"].isin(keep)].copy()
    else:
        raise ValueError(f"unknown panel_type {ptype!r}")

    panel = panel.sort_values(["item_id", "vintage"]).reset_index(drop=True)
    return panel


def run_battery(df, corpus_name, covariates, reference_mode="explicit",
                confirmed_exogenous=None):
    """Full pre-screen + battery for one corpus. Returns a result dict."""
    rep = s34.prescreen_corpus(
        df, covariates=covariates, reference_mode=reference_mode,
        confirmed_exogenous=confirmed_exogenous, corpus_name=corpus_name)

    result = {
        "corpus": corpus_name,
        "prescreen": {
            "DIAGNOSTIC_USABLE": rep["DIAGNOSTIC_USABLE"],
            "BSMD_TESTABLE": rep["BSMD_TESTABLE"],
            "panel_type": rep["panel_type"],
            "admissible_covariates": rep["admissible_covariates"],
            "predictor_classes": {c: v["class"] for c, v in rep["predictor_classes"].items()},
        },
        "tests": {},
        "status": "ok",
    }

    if not rep["DIAGNOSTIC_USABLE"] or rep["panel_type"] == "label-driven-unbalanced-EXCLUDE":
        result["status"] = "skipped: not diagnostic-usable or EXCLUDE panel"
        return result

    if not _S33_OK:
        result["status"] = f"cannot run tests: s33 import failed ({_S33_IMPORT_ERR}). " \
                           f"Reconcile the ADAPTER block."
        return result

    panel = build_panel(df, rep, reference_mode=reference_mode)
    result["n_items_panel"] = int(panel["item_id"].nunique())
    result["n_vintages"] = int(panel["vintage"].nunique())

    # --- always-run tests (naive, VR, EL) ---
    try:
        F_obs, p_naive = _naive_fe(panel)                 # RECONCILE signature
        result["tests"]["naive_FE"] = {"F": float(F_obs), "p": float(p_naive),
                                       "reject_0.05": bool(p_naive < 0.05)}
    except Exception as e:
        result["tests"]["naive_FE"] = {"error": repr(e)}
    for name, fn in [("VR", _vr_test), ("EL", _el_test)]:
        try:
            stat, p = fn(panel)                            # RECONCILE signature
            result["tests"][name] = {"stat": float(stat), "p": float(p),
                                     "reject_0.05": bool(p < 0.05)}
        except Exception as e:
            result["tests"][name] = {"error": repr(e)}

    # --- BSMD: only if a CONFIRMED admissible covariate exists ---
    if rep["BSMD_TESTABLE"] and rep["admissible_covariates"]:
        W = rep["admissible_covariates"][0]  # primary admissible covariate
        result["bsmd_covariate"] = W
        for mode in ("plugin", "crossfit"):
            try:
                S, p = _bsmd_stat(panel, W, mode)          # RECONCILE signature
                result["tests"][f"BSMD_{mode}"] = {"S": float(S), "p": float(p),
                                                   "reject_0.05": bool(p < 0.05)}
            except Exception as e:
                result["tests"][f"BSMD_{mode}"] = {"error": repr(e)}
    else:
        result["tests"]["BSMD_plugin"] = {"status": "degenerate/NA (no admissible W)"}
        result["tests"]["BSMD_crossfit"] = {"status": "degenerate/NA (no admissible W)"}

    return result


def main():
    """CLI: python s35_realcorpus.py corpus.parquet corpus_name [--ref-mode terminal]
              [--covars a,b,c] [--confirm a,b]
       Writes <corpus_name>_battery.json."""
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("corpus_name")
    ap.add_argument("--ref-mode", default="explicit", choices=["explicit", "terminal"])
    ap.add_argument("--covars", default=None, help="comma-separated covariate columns")
    ap.add_argument("--confirm", default=None,
                    help="comma-separated covariates the ANALYST confirms are exogenous "
                         "(computable without the label). REQUIRED for BSMD to run.")
    args = ap.parse_args()

    if args.path.endswith(".parquet"):
        df = pd.read_parquet(args.path)
    else:
        df = pd.read_csv(args.path)

    covars = args.covars.split(",") if args.covars else None
    confirm = args.confirm.split(",") if args.confirm else None

    res = run_battery(df, args.corpus_name, covariates=covars,
                      reference_mode=args.ref_mode, confirmed_exogenous=confirm)
    out = f"{args.corpus_name}_battery.json"
    with open(out, "w") as fh:
        json.dump(res, fh, indent=2)
    print(json.dumps(res, indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
