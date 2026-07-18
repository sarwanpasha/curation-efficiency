#!/usr/bin/env python3
"""
s34_prescreen.py -- Admissibility pre-screen for real vintaged corpora.

Implements the A1-A5 (diagnostic-usable) and B1-B4 (BSMD-testable) checklist,
panel-type classification, and per-predictor admissibility classification.

Self-contained (numpy/pandas only). Importable: the battery script s35 calls
prescreen_corpus(df, ...) directly.

INPUT: long-format DataFrame with columns:
    item_id   : stable entity id, persistent across vintages (locus/variant/article)
    vintage   : ordered release index (int or sortable); the time axis
    label     : the revisable label at (item, vintage). May be non-binary; see `ref`.
    ref       : (optional) the fixed reference value per item. If provided, Y = 1[label==ref].
                If absent, pass reference_mode to derive Y (e.g. 'terminal' uses each
                item's last-vintage label as its fixed reference).
    <covariates...> : any number of candidate predictor columns to classify.

OUTPUT: a dict report with DIAGNOSTIC_USABLE, BSMD_TESTABLE, panel_type, the A/B
gate details with numbers, and per-predictor classification.

Design notes:
- Y_{i,t} = 1[ label_{i,t} agrees with fixed reference ].  Reference is FIXED across
  the window (A2); if derived as 'terminal', it is each item's final-vintage label.
- We NEVER forward-fill labels to balance (would fabricate Delta=0). Unbalanced
  handling is by balanced-core or exogenous-unbalanced only.
- C1 (exogeneity): a covariate is boundary-coupled if it is (near-)deterministically
  a function of the label path. We test this empirically: regress/predict the
  covariate from {lagged label, run-length, cumulative label} within item; if the
  covariate is (almost) perfectly determined by these, it FAILS C1. This is a
  heuristic screen, flagged as such -- final C1 judgement is the analyst's, informed
  by how the covariate was CONSTRUCTED (the decisive test is "computable without Y?").
- C2 (temporal variation): covariate varies within item across vintages.
"""

import numpy as np
import pandas as pd

RESERVED = {"item_id", "vintage", "label", "ref", "Y", "Delta"}


# ----------------------------------------------------------------------
# Panel construction
# ----------------------------------------------------------------------
def build_Y(df, reference_mode="explicit"):
    """Add Y = 1[label agrees with fixed reference] and the increment Delta.
    reference_mode:
      'explicit' : use df['ref'] (must be present); fixed per item.
      'terminal' : reference = each item's label at its max vintage (fixed per item).
    Returns df sorted by (item_id, vintage) with columns Y, Delta added.
    Delta is defined only for t>=1 within each item (first vintage has NaN Delta).
    """
    df = df.copy()
    df = df.sort_values(["item_id", "vintage"]).reset_index(drop=True)
    if reference_mode == "explicit":
        if "ref" not in df.columns:
            raise ValueError("reference_mode='explicit' requires a 'ref' column.")
        ref = df["ref"]
    elif reference_mode == "terminal":
        term = df.groupby("item_id").apply(
            lambda g: g.loc[g["vintage"].idxmax(), "label"]
        )
        ref = df["item_id"].map(term)
    else:
        raise ValueError(f"unknown reference_mode {reference_mode!r}")
    df["Y"] = (df["label"].values == np.asarray(ref)).astype(int)
    df["Delta"] = df.groupby("item_id")["Y"].diff()  # NaN at each item's first vintage
    return df


# ----------------------------------------------------------------------
# A1-A5 : diagnostic-usable gates
# ----------------------------------------------------------------------
def check_diagnostic(df, min_vintages=3, min_match_frac=0.5, min_move_frac=1e-6):
    """Return (passed: bool, details: dict). df must already have Y, Delta."""
    d = {}

    # A1: reconstructable panel + persistent item identity (match rate t vs t-1)
    vintages = np.sort(df["vintage"].unique())
    d["n_items"] = int(df["item_id"].nunique())
    d["n_vintages"] = int(len(vintages))
    # match fraction: of items present at t, fraction also present at t-1
    present = df.groupby("vintage")["item_id"].apply(set)
    match_fracs = []
    for a, b in zip(vintages[1:], vintages[:-1]):
        cur, prev = present.get(a, set()), present.get(b, set())
        if len(cur) > 0:
            match_fracs.append(len(cur & prev) / len(cur))
    d["min_adjacent_match_frac"] = float(np.min(match_fracs)) if match_fracs else 0.0
    A1 = (d["n_vintages"] >= min_vintages) and (d["min_adjacent_match_frac"] >= min_match_frac)
    d["A1_panel_identity"] = bool(A1)

    # A2: fixed reference -- structurally enforced by build_Y (ref fixed per item).
    # We flag whether ref was constant within item (it must be).
    if "ref" in df.columns:
        ref_const = df.groupby("item_id")["ref"].nunique().max() <= 1
    else:
        ref_const = True  # terminal-derived ref is fixed by construction
    d["A2_fixed_reference"] = bool(ref_const)

    # A3: binary Y
    yvals = set(np.unique(df["Y"].dropna()))
    A3 = yvals.issubset({0, 1})
    d["A3_binary_label"] = bool(A3)

    # A4: >=3 vintages (>=2 increments)
    A4 = d["n_vintages"] >= min_vintages
    d["A4_min_vintages"] = bool(A4)

    # A5: nonzero movement
    moved = df["Delta"].fillna(0).abs()
    d["move_frac"] = float((moved > 0).mean())
    A5 = d["move_frac"] > min_move_frac
    d["A5_nonzero_movement"] = bool(A5)

    passed = all([A1, ref_const, A3, A4, A5])
    return passed, d


# ----------------------------------------------------------------------
# Panel-type classification (balanced-core / exogenous-unbalanced / EXCLUDE)
# ----------------------------------------------------------------------
def classify_panel(df, balanced_core_min_frac=0.30, entry_label_bias_tol=0.15):
    """Return (panel_type, details). Detects label-driven entry (selection on outcome)."""
    d = {}
    vintages = np.sort(df["vintage"].unique())
    n_items = df["item_id"].nunique()

    # balanced core: items present in ALL vintages
    counts = df.groupby("item_id")["vintage"].nunique()
    core = counts[counts == len(vintages)].index
    d["balanced_core_size"] = int(len(core))
    d["balanced_core_frac"] = float(len(core) / n_items) if n_items else 0.0

    # entry-exogeneity check: do items ENTERING at vintage t have a different label
    # rate than items CONTINUING? Large gap => entry correlated with label => EXCLUDE.
    entry_bias = []
    present = df.groupby("vintage")["item_id"].apply(set)
    ylook = df.set_index(["item_id", "vintage"])["Y"]
    for a, b in zip(vintages[1:], vintages[:-1]):
        entrants = present.get(a, set()) - present.get(b, set())
        continuing = present.get(a, set()) & present.get(b, set())
        if len(entrants) >= 5 and len(continuing) >= 5:
            ye = np.mean([ylook.get((i, a), np.nan) for i in entrants])
            yc = np.mean([ylook.get((i, a), np.nan) for i in continuing])
            if np.isfinite(ye) and np.isfinite(yc):
                entry_bias.append(abs(ye - yc))
    d["max_entry_label_bias"] = float(np.max(entry_bias)) if entry_bias else 0.0

    if d["balanced_core_frac"] >= balanced_core_min_frac:
        panel_type = "balanced-core"
    elif d["max_entry_label_bias"] <= entry_label_bias_tol:
        panel_type = "exogenous-unbalanced"
    else:
        panel_type = "label-driven-unbalanced-EXCLUDE"
    d["panel_type"] = panel_type
    return panel_type, d


# ----------------------------------------------------------------------
# B1-B4 : per-predictor admissibility classification
# ----------------------------------------------------------------------
def _within_item_var(df, col):
    """Mean within-item variance of a covariate (0 => locus-constant => fails C2)."""
    g = df.groupby("item_id")[col]
    # fraction of items where the covariate varies at all
    varies = g.nunique(dropna=True) > 1
    return float(varies.mean())


def _label_path_determination(df, col):
    """Heuristic C1 screen: how well is the covariate determined by label-path features
    within item? Returns R^2-like determination in [0,1]; high => the covariate is a
    function of the label path => FAILS C1.

    IMPORTANT: this screen includes the CONTEMPORANEOUS label Y_{i,t} and its identity
    among the predictors, because a covariate that equals (or is a function of) the
    current label is the most dangerous boundary-coupled case and must be caught. We
    also short-circuit to R2=1.0 if the covariate is (near-)exactly equal to Y or to a
    simple label-path feature.

    This remains a HEURISTIC. The DECISIVE C1 test is constructional: 'could this
    covariate be computed without ever reading the label column?' The tool therefore
    reports the screen but REQUIRES analyst confirmation (see confirmed_exogenous).
    """
    d = df.sort_values(["item_id", "vintage"]).copy()
    d["Ycur"] = d["Y"]                                   # contemporaneous label
    d["lagY"] = d.groupby("item_id")["Y"].shift(1)
    d["cumY"] = d.groupby("item_id")["Y"].cumsum().shift(1)
    def runlen(y):
        out = np.ones(len(y))
        for k in range(1, len(y)):
            out[k] = out[k-1] + 1 if y.iloc[k] == y.iloc[k-1] else 1
        return pd.Series(out, index=y.index)
    d["run"] = d.groupby("item_id")["Y"].apply(runlen).reset_index(level=0, drop=True)

    # exact-identity short circuit: covariate equals Y (or 1-Y) almost everywhere.
    # NaN-safe, aligned on non-missing covariate rows only.
    mask = d[col].notna() & d["Ycur"].notna()
    if mask.sum() > 0:
        cv = d.loc[mask, col].to_numpy(float)
        yv = d.loc[mask, "Ycur"].to_numpy(float)
        eq_Y = np.mean(np.isclose(cv, yv))
        eq_negY = np.mean(np.isclose(cv, 1.0 - yv))
        if max(eq_Y, eq_negY) > 0.98:
            return 1.0  # it IS the label -> boundary-coupled

    feats = ["Ycur", "lagY", "cumY", "run"]
    sub = d[feats + [col]].dropna()
    if len(sub) < 10 or sub[col].nunique() <= 1:
        return np.nan
    X = np.column_stack([np.ones(len(sub))] + [sub[f].to_numpy(float) for f in feats])
    y = sub[col].to_numpy(float)
    try:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        yhat = X @ beta
        ss_res = np.sum((y - yhat) ** 2)
        ss_tot = np.sum((y - y.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        return float(max(0.0, min(1.0, r2)))
    except Exception:
        return np.nan


def classify_predictors(df, covariates, c1_r2_threshold=0.95, c2_var_frac_min=0.05,
                        confirmed_exogenous=None):
    """Classify each covariate: boundary-coupled / degenerate / admissible.
    c1_r2_threshold: label-path determination above this => FAILS C1 (boundary-coupled).
    c2_var_frac_min: fraction of items with within-item variation below this => FAILS C2.
    confirmed_exogenous: set/list of covariate names the ANALYST has confirmed are
        constructionally exogenous (computable without the label). A covariate can only
        be 'admissible' if it is BOTH heuristically C1-clean AND analyst-confirmed. If
        confirmed_exogenous is None, no covariate is confirmed, and the best any can
        reach is 'admissible-UNCONFIRMED' (which s35 must NOT feed to BSMD until
        confirmed). This makes the constructional C1 test mandatory, not optional.
    """
    confirmed = set(confirmed_exogenous) if confirmed_exogenous is not None else set()
    out = {}
    for col in covariates:
        if col in RESERVED:
            continue
        c2_frac = _within_item_var(df, col)
        c2_pass = c2_frac >= c2_var_frac_min
        r2 = _label_path_determination(df, col)
        c1_heuristic_pass = (not np.isfinite(r2)) or (r2 < c1_r2_threshold)
        c1_confirmed = col in confirmed

        if not c1_heuristic_pass:
            klass = "boundary-coupled"                 # heuristic caught label-coupling
        elif not c2_pass:
            klass = "degenerate"                       # C1 ok but locus-constant
        elif not c1_confirmed:
            klass = "admissible-UNCONFIRMED"           # passes screens, NOT analyst-confirmed
        else:
            klass = "admissible"                       # C1 (heuristic+confirmed) and C2
        out[col] = {
            "class": klass,
            "C1_heuristic_pass": bool(c1_heuristic_pass),
            "C1_labelpath_R2": (None if not np.isfinite(r2) else round(r2, 4)),
            "C1_analyst_confirmed": bool(c1_confirmed),
            "C2_temporal_variation_pass": bool(c2_pass),
            "C2_within_item_var_frac": round(c2_frac, 4),
            "note": "C1 R2 is heuristic. 'admissible' requires BOTH heuristic pass AND "
                    "analyst confirmation that the covariate is computable WITHOUT the "
                    "label. Pass confirmed_exogenous=[...] once verified by construction.",
        }
    return out


# ----------------------------------------------------------------------
# Top-level pre-screen
# ----------------------------------------------------------------------
def prescreen_corpus(df, covariates=None, reference_mode="explicit",
                     min_vintages=3, min_match_frac=0.5,
                     balanced_core_min_frac=0.30, entry_label_bias_tol=0.15,
                     c1_r2_threshold=0.95, c2_var_frac_min=0.05,
                     confirmed_exogenous=None, corpus_name="corpus"):
    """Full pre-screen. Returns a report dict.
    confirmed_exogenous: list of covariate names the analyst has verified are
        constructionally exogenous (computable without the label). REQUIRED for a
        covariate to count as 'admissible' and for the corpus to be BSMD_TESTABLE.
    """
    if covariates is None:
        covariates = [c for c in df.columns if c not in RESERVED]

    dfY = build_Y(df, reference_mode=reference_mode)

    diag_pass, diag = check_diagnostic(
        dfY, min_vintages=min_vintages, min_match_frac=min_match_frac)
    panel_type, panel = classify_panel(
        dfY, balanced_core_min_frac=balanced_core_min_frac,
        entry_label_bias_tol=entry_label_bias_tol)
    preds = classify_predictors(
        dfY, covariates, c1_r2_threshold=c1_r2_threshold,
        c2_var_frac_min=c2_var_frac_min, confirmed_exogenous=confirmed_exogenous)

    admissible_cols = [c for c, v in preds.items() if v["class"] == "admissible"]
    unconfirmed_cols = [c for c, v in preds.items() if v["class"] == "admissible-UNCONFIRMED"]

    # Advisory warnings: the heuristic is best-effort, NOT a guarantee. The binding
    # C1 test is the analyst's constructional confirmation. Warn loudly if any
    # analyst-confirmed covariate shows non-trivial label-path determination, since
    # that is a sign the confirmation may be mistaken.
    warnings = []
    for c, v in preds.items():
        r2 = v["C1_labelpath_R2"]
        if v["C1_analyst_confirmed"] and (r2 is not None) and r2 > 0.5:
            warnings.append(
                f"CONFIRMED covariate {c!r} has label-path R2={r2}: analyst confirmed "
                f"it as exogenous, but it is substantially predicted by the label path. "
                f"RE-CHECK the constructional C1 test ('computable without the label?').")
    if warnings:
        print("!!! C1 ADVISORY WARNINGS (heuristic is not a guarantee) !!!")
        for w in warnings:
            print("  " + w)
    bsmd_testable = bool(diag_pass
                         and panel_type != "label-driven-unbalanced-EXCLUDE"
                         and len(admissible_cols) > 0)

    report = {
        "corpus": corpus_name,
        "DIAGNOSTIC_USABLE": bool(diag_pass and panel_type != "label-driven-unbalanced-EXCLUDE"),
        "BSMD_TESTABLE": bsmd_testable,
        "panel_type": panel_type,
        "admissible_covariates": admissible_cols,
        "admissible_unconfirmed_covariates": unconfirmed_cols,
        "A_gates": diag,
        "panel_details": panel,
        "predictor_classes": preds,
        "reference_mode": reference_mode,
    }
    return report


def print_report(rep):
    print("=" * 70)
    print(f"CORPUS: {rep['corpus']}")
    print(f"  DIAGNOSTIC_USABLE = {rep['DIAGNOSTIC_USABLE']}")
    print(f"  BSMD_TESTABLE     = {rep['BSMD_TESTABLE']}")
    print(f"  panel_type        = {rep['panel_type']}")
    print(f"  admissible W      = {rep['admissible_covariates'] or 'NONE (BSMD degenerate)'}")
    print("  A1-A5 gates:")
    for k, v in rep["A_gates"].items():
        print(f"    {k}: {v}")
    print("  predictor classes:")
    for c, v in rep["predictor_classes"].items():
        print(f"    {c:>24s} -> {v['class']:22s} "
              f"C1_heur={v['C1_heuristic_pass']}(R2={v['C1_labelpath_R2']}) "
              f"C1_conf={v['C1_analyst_confirmed']} "
              f"C2={v['C2_temporal_variation_pass']}(varfrac={v['C2_within_item_var_frac']})")
    print("=" * 70)


if __name__ == "__main__":
    # self-test on a tiny synthetic corpus
    rng = np.random.default_rng(0)
    rows = []
    for i in range(200):
        ref = rng.integers(0, 2)
        y = rng.integers(0, 2)
        w_exo = rng.uniform(-1, 1)  # will be made time-varying below
        for t in range(6):
            y = y if rng.random() < 0.8 else 1 - y
            rows.append(dict(item_id=i, vintage=t, label=y, ref=ref,
                             W_exo=rng.uniform(-1, 1),      # time-varying exogenous (admissible)
                             W_static=w_exo,                 # locus-constant (degenerate)
                             W_lag=y))                        # label-path (boundary-coupled)
    df = pd.DataFrame(rows)
    rep = prescreen_corpus(df, covariates=["W_exo", "W_static", "W_lag"],
                           reference_mode="explicit", corpus_name="selftest")
    print_report(rep)
