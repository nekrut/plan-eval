#!/usr/bin/env python3
"""
Regenerate the 5 paper figures from results.csv + per-run score.json.

usage:
  python3 scripts/make_figures.py --all
  python3 scripts/make_figures.py --fig 1
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

BENCH = Path(__file__).resolve().parent.parent
RESULTS_CSV = BENCH / "results.csv"
FIGS = BENCH / "figures"
FIGS.mkdir(exist_ok=True)

# plan_version (in results.csv) -> (display plan, hardware)
PLAN_MAP = {
    "current":    ("v2",    "jetson"),
    "v1":         ("v1",    "jetson"),
    "5080_v1":    ("v1",    "5080"),
    "5080_v2":    ("v2",    "5080"),
    "5080_v1p25": ("v1.25", "5080"),
    "5080_v1p5":  ("v1.5",  "5080"),
    "5080_v0p5":  ("v0.5",  "5080"),
    "5080_v1g":   ("v1g",   "5080"),
    # Jetson recipe-variant stages (added 2026-06-09 for gemma4:12b — the full
    # Figure-1 row run on the Jetson, mirroring the 5080 sweep in matrix_5080.py).
    "jetson_v1p25": ("v1.25", "jetson"),
    "jetson_v1p5":  ("v1.5",  "jetson"),
    "jetson_v1g":   ("v1g",   "jetson"),
    "jetson_v0p5":  ("v0.5",  "jetson"),
    "jetson_b":     ("v2",    "jetson"),   # Track-B (no-plan) runs collapse to the "B" column
}

# Plan-detail order, lean → detailed; "B" prepended for the no-plan column.
PLAN_COLS = ["B", "v0.5", "v1", "v1g", "v1.25", "v1.5", "v2"]

ANTHROPIC = {"claude-opus-4-7", "claude-sonnet-4-6", "claude-haiku-4-5"}

# For pretty model labels: convert internal `qwen3.6_35b-a3b` -> `qwen3.6:35b-a3b`.
def pretty(model: str) -> str:
    if model in ANTHROPIC:
        return model.replace("claude-", "Claude ").replace("-", " ").title().replace("4 7", "4.7").replace("4 6", "4.6").replace("4 5", "4.5")
    # Restore the colon that aggregate stripped
    parts = model.split("_", 1)
    return parts[0] + ":" + parts[1] if len(parts) == 2 else model


def load_data() -> pd.DataFrame:
    df = pd.read_csv(RESULTS_CSV)
    df["plan"]     = df["plan_version"].map(lambda v: PLAN_MAP.get(v, (v, "?"))[0])
    df["hardware"] = df["plan_version"].map(lambda v: PLAN_MAP.get(v, (v, "?"))[1])
    # Heatmap column: Track B (no plan) collapses into a single "B" column except
    # for v0.5 which is its own designed Track-B variant.
    def col(r):
        if r["plan"] == "v0.5":
            return "v0.5"
        if r["track"] == "B":
            return "B"
        return r["plan"]
    df["plan_col"] = df.apply(col, axis=1)
    return df


def model_order(models):
    """Anthropic top, then by descending mean M3 across all v2 cells."""
    return sorted(models, key=lambda m: (0 if m in ANTHROPIC else 1, m))


def load_em_baseline_M3() -> pd.DataFrame:
    """Pull v2 + v2_defensive baseline (none-injection) cells from the
    error_matrix_*.jsonl files, return as a DataFrame with the same columns
    used by the headline heatmap (model, plan_col, hardware, M3).
    These are the cells where the PATH-shim was not active, so M3 reflects
    pure model-on-recipe performance — directly comparable to results.csv."""
    rows = []
    for jsonl in BENCH.glob("error_matrix_*.jsonl"):
        name = jsonl.stem
        if "_m4" in name:    hw = "m4"
        elif "_a5000" in name: hw = "a5000"
        else:                  hw = "jetson"
        for line in jsonl.read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            cell = r.get("cell", "")
            parts = cell.split("/")
            if len(parts) < 4:
                continue
            model_raw, plan, pattern_at_target, _seed = parts
            if "@" not in pattern_at_target:
                continue
            pattern, _ = pattern_at_target.split("@", 1)
            if pattern != "none":
                continue
            if plan not in ("v2", "v2_defensive"):
                continue
            m3 = r.get("M3")
            if m3 is None:
                continue
            rows.append({
                "model":    model_raw.replace(":", "_"),
                "plan_col": plan,
                "hardware": hw,
                "track":    "A",
                "M3":       m3,
            })
    return pd.DataFrame(rows)


# ────────────────────────────────────────────────────────────────────────────
# Figure 1 — headline matrix heatmap
# ────────────────────────────────────────────────────────────────────────────
def fig1_heatmap(df: pd.DataFrame, out: Path):
    pivot = df.pivot_table(index="model", columns="plan_col", values="M3", aggfunc="mean")
    cnt   = df.pivot_table(index="model", columns="plan_col", values="M3", aggfunc="count")
    cols  = [c for c in PLAN_COLS if c in pivot.columns]
    pivot = pivot[cols]
    cnt   = cnt[cols]
    rows  = model_order(pivot.index)
    pivot = pivot.loc[rows]
    cnt   = cnt.loc[rows]

    fig, ax = plt.subplots(figsize=(8.5, 9.5))
    cmap = mpl.cm.RdYlGn
    cmap.set_bad(color="#dddddd")
    im = ax.imshow(np.ma.masked_invalid(pivot.values), aspect="auto",
                   cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([pretty(m) for m in rows], fontsize=8)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            n = cnt.values[i, j]
            if pd.isna(v):
                continue
            ax.text(j, i, f"{v:.2f}\n(n={int(n)})",
                    ha="center", va="center", fontsize=6.5,
                    color="black" if 0.25 < v < 0.85 else "white")
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("mean M3 (Jaccard)")
    ax.set_xlabel("Plan variant (lean → detailed; B = no plan)")
    ax.set_title("Figure 1. Mean M3 by model × plan variant\n"
                 "(both hardware platforms pooled; grey = untested cell)")
    plt.tight_layout()
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()


def fig1_per_platform(df: pd.DataFrame, out_dir: Path):
    """One heatmap per hardware platform, using results.csv data plus
    error_matrix baseline cells. All four panels share the same y-axis
    (union of all models tested on any platform, in the same order) so
    the per-platform coverage gaps line up row-for-row."""
    PLAN_COLS_FULL = PLAN_COLS + ["v2_defensive"]
    PLATFORM_LABEL = {
        "jetson": "NVIDIA Jetson AGX Orin",
        "5080":   "RTX 5080 desktop",
        "m4":     "MacBook Pro M4 Pro",
        "a5000":  "2x RTX A5000 workstation",
    }

    em = load_em_baseline_M3()
    full = pd.concat([df[["model","plan_col","hardware","M3"]], em[["model","plan_col","hardware","M3"]]],
                     ignore_index=True)

    # Shared y-axis: union of all models tested on ANY platform, sorted
    # Anthropic-first, then alphabetical.
    all_models = sorted(full["model"].unique().tolist(),
                        key=lambda m: (0 if m in ANTHROPIC else 1, m))

    for hw in ["jetson", "5080", "m4", "a5000"]:
        sub = full[full["hardware"] == hw]
        if sub.empty:
            print(f"[fig1_per_platform] {hw}: no data — skipping")
            continue
        pivot = sub.pivot_table(index="model", columns="plan_col", values="M3", aggfunc="mean")
        cnt   = sub.pivot_table(index="model", columns="plan_col", values="M3", aggfunc="count")
        # Reindex to the union model list — missing rows become NaN (grey)
        pivot = pivot.reindex(index=all_models)
        cnt   = cnt.reindex(index=all_models)
        cols  = [c for c in PLAN_COLS_FULL if c in pivot.columns]
        pivot = pivot[cols]
        cnt   = cnt[cols]

        n_present = pivot.notna().any(axis=1).sum()
        fig, ax = plt.subplots(figsize=(max(5, 1.2*len(cols)+2), max(3, 0.45*len(all_models)+1.5)))
        cmap = mpl.cm.RdYlGn
        cmap.set_bad(color="#dddddd")
        im = ax.imshow(np.ma.masked_invalid(pivot.values), aspect="auto",
                       cmap=cmap, vmin=0, vmax=1)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols)
        ax.set_yticks(range(len(all_models)))
        ax.set_yticklabels([pretty(m) for m in all_models], fontsize=8)
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v = pivot.values[i, j]
                n = cnt.values[i, j]
                if pd.isna(v):
                    continue
                ax.text(j, i, f"{v:.2f}\n(n={int(n)})",
                        ha="center", va="center", fontsize=6.5,
                        color="black" if 0.25 < v < 0.85 else "white")
        cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label("mean M3 (Jaccard)")
        ax.set_xlabel("Plan variant (lean → detailed; B = no plan)")
        ax.set_title(f"{PLATFORM_LABEL[hw]}\n"
                     f"({n_present} of {len(all_models)} models tested × {len(cols)} plan variants; grey = untested)")
        plt.tight_layout()
        out_path = out_dir / f"fig1_headline_heatmap_{hw}.png"
        plt.savefig(out_path, dpi=130, bbox_inches="tight")
        plt.close()
        print(f"[fig1_per_platform] wrote {out_path}")

    # Combined 2x2 figure with identical axes per panel — all 26 models
    # × all 8 plan variants on each panel; grey = untested.
    panels = [
        ("jetson", "A"), ("5080",  "B"),
        ("m4",     "C"), ("a5000", "D"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(14, 13), sharex=True, sharey=True)
    cmap = mpl.cm.RdYlGn
    cmap.set_bad(color="#dddddd")

    last_im = None
    for ax, (hw, letter) in zip(axes.flat, panels):
        sub = full[full["hardware"] == hw]
        pivot = sub.pivot_table(index="model", columns="plan_col", values="M3", aggfunc="mean")
        cnt   = sub.pivot_table(index="model", columns="plan_col", values="M3", aggfunc="count")
        # Reindex BOTH axes so every panel is identical shape and order
        pivot = pivot.reindex(index=all_models, columns=PLAN_COLS_FULL)
        cnt   = cnt.reindex(index=all_models, columns=PLAN_COLS_FULL)
        n_present = pivot.notna().any(axis=1).sum()
        last_im = ax.imshow(np.ma.masked_invalid(pivot.values), aspect="auto",
                            cmap=cmap, vmin=0, vmax=1)
        ax.set_xticks(range(len(PLAN_COLS_FULL)))
        ax.set_xticklabels(PLAN_COLS_FULL, rotation=30, ha="right", fontsize=10)
        ax.set_yticks(range(len(all_models)))
        ax.set_yticklabels([pretty(m) for m in all_models], fontsize=9)
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v = pivot.values[i, j]
                if pd.isna(v):
                    continue
                ax.text(j, i, f"{v:.2f}",
                        ha="center", va="center", fontsize=7.5,
                        color="black" if 0.25 < v < 0.85 else "white")
        ax.set_title(f"({letter}) {PLATFORM_LABEL[hw]} — "
                     f"{n_present}/{len(all_models)} models tested",
                     fontsize=11)

    cb = fig.colorbar(last_im, ax=axes.ravel().tolist(),
                      fraction=0.025, pad=0.02)
    cb.set_label("mean M3 (Jaccard)", fontsize=10)
    fig.suptitle("Headline matrix per hardware platform "
                 "(grey = untested combination)", fontsize=13, y=1.00)
    out_combined = out_dir / "fig1_headline_heatmap_per_platform.png"
    plt.savefig(out_combined, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"[fig1_per_platform] wrote {out_combined}")


# ────────────────────────────────────────────────────────────────────────────
# Figure 2 — the v1 cliff and its single-line repair
# ────────────────────────────────────────────────────────────────────────────
def fig2_intermediates(df: pd.DataFrame, out: Path):
    # 5080 only; subset of representative models
    subset = [
        "qwen3.6_27b",       # dense ≥27B that already passes v1
        "qwen3.5_27b",       # dense ≥27B that fails v1, recovered by v1.25
        "qwen3.6_35b-a3b",   # MoE, fails v1, recovered by v1.25
        "qwen3-coder_30b",   # coder
        "gemma4_26b",        # generalist
        "qwen3_14b",         # small dense
        "qwen3.5_9b",        # small dense
        "gemma4_e4b",        # tiny
    ]
    plans = ["v1", "v1.25", "v1.5", "v2"]
    sub = df[(df["hardware"] == "5080") & (df["track"] == "A") & (df["model"].isin(subset)) & (df["plan"].isin(plans))]
    g = sub.groupby(["model", "plan"])["M3"].agg(["mean", "std", "count"]).reset_index()

    fig, ax = plt.subplots(figsize=(10, 5))
    width = 0.2
    x = np.arange(len(subset))
    colors = {"v1": "#d73027", "v1.25": "#fdae61", "v1.5": "#abd9e9", "v2": "#1a9850"}
    for i, p in enumerate(plans):
        means = []
        stds  = []
        for m in subset:
            row = g[(g["model"] == m) & (g["plan"] == p)]
            means.append(row["mean"].iloc[0] if len(row) else np.nan)
            stds.append(row["std"].iloc[0] if len(row) and not pd.isna(row["std"].iloc[0]) else 0.0)
        ax.bar(x + (i - 1.5) * width, means, width, yerr=stds, capsize=2,
               label=p, color=colors[p], edgecolor="black", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([pretty(m) for m in subset], rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("M3 (mean ± std, n=3)")
    ax.set_ylim(0, 1.05)
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle=":")
    ax.legend(title="Plan variant", loc="lower left")
    ax.set_title("Figure 2. The v1 cliff and its single-line repair (RTX 5080, Track A)\n"
                 "v1.25 = v1 + literal lofreq line; v1.5 = v2 with prose stripped")
    plt.tight_layout()
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()


# ────────────────────────────────────────────────────────────────────────────
# Figure 3 — plan value (Track A vs Track B)
# ────────────────────────────────────────────────────────────────────────────
def fig3_track_b(df: pd.DataFrame, out: Path):
    # Three points per model: B (mean across any plan_version where track=B and plan != v0.5),
    # v1 Track A, v2 Track A. Use 5080 where present, else Jetson.
    rows = []
    for m, g in df.groupby("model"):
        b_vals  = g[(g["track"] == "B") & (g["plan"] != "v0.5")]["M3"].values
        a1_vals = g[(g["track"] == "A") & (g["plan"] == "v1")]["M3"].values
        a2_vals = g[(g["track"] == "A") & (g["plan"] == "v2")]["M3"].values
        if len(b_vals) and len(a2_vals):
            rows.append({"model": m, "B": np.mean(b_vals),
                         "v1": np.mean(a1_vals) if len(a1_vals) else np.nan,
                         "v2": np.mean(a2_vals)})
    sub = pd.DataFrame(rows)

    def family(m):
        if m in ANTHROPIC: return "Anthropic frontier"
        if m.endswith("-a3b") or m.startswith("gpt-oss") or m.startswith("glm-4.7-flash") or m.startswith("deepseek-coder"):
            return "MoE local"
        return "Dense local"

    sub["family"] = sub["model"].map(family)
    palette = {"Anthropic frontier": "#1f77b4", "Dense local": "#2ca02c", "MoE local": "#d62728"}

    fig, ax = plt.subplots(figsize=(8, 5.5))
    xs = [0, 1, 2]
    # Add small jitter so converging lines at v2=1.0 are individually visible.
    rng = np.random.default_rng(0)
    for _, r in sub.iterrows():
        ys = [r["B"], r["v1"], r["v2"]]
        jit = rng.uniform(-0.012, 0.012, size=3)
        ax.plot(xs, np.array(ys) + jit, "-o", color=palette[r["family"]], alpha=0.55,
                linewidth=1.0, markersize=4, label=r["family"])
    # de-dupe legend
    handles, labels = ax.get_legend_handles_labels()
    by = dict(zip(labels, handles))
    ax.legend(by.values(), by.keys(), loc="center right")
    ax.set_xticks(xs)
    ax.set_xticklabels(["Track B (no plan)", "Track A v1 (lean)", "Track A v2 (detailed)"])
    ax.set_ylabel("mean M3 (Jaccard)")
    ax.set_ylim(-0.05, 1.1)
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle=":")
    ax.set_title("Figure 3. Plan value: Track B vs Track A v1 vs Track A v2\n"
                 "(per-model means; pooled across hardware where both tested)")
    plt.tight_layout()
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()


# ────────────────────────────────────────────────────────────────────────────
# Figure 4 — plan-source robustness (v1g)
# ────────────────────────────────────────────────────────────────────────────
def collect_per_sample(plan_version: str):
    """Return DataFrame of per-sample Jaccard rows for one plan_version on 5080.
    Failed runs (no per_sample) contribute 4 rows of M3_sample=run-level m3_jaccard
    (typically 0.0) so models that fail this plan still appear in box plots.
    """
    runs_dir = BENCH / f"runs_{plan_version}"
    rows = []
    SAMPLES = ["M117-bl", "M117-ch", "M117C1-bl", "M117C1-ch"]
    for d in runs_dir.glob("*_track-A_seed-*"):
        sj = d / "score.json"
        if not sj.exists():
            continue
        s = json.loads(sj.read_text())
        per = s.get("m3_detail", {}).get("per_sample", {})
        meta = json.loads((d / "meta.json").read_text())
        model = meta.get("model", "").replace("/", "_").replace(":", "_")
        if per:
            for sample, v in per.items():
                rows.append({"model": model, "sample": sample, "M3_sample": v})
        else:
            run_m3 = float(s.get("m3_jaccard") or 0.0)
            for sample in SAMPLES:
                rows.append({"model": model, "sample": sample, "M3_sample": run_m3})
    return pd.DataFrame(rows)


def fig4_v1g(df: pd.DataFrame, out: Path):
    """Grouped bar: per-model M3 mean ± std on v1.25 (hand-authored) vs v1g (IUC mechanical).
    Same plan otherwise; the only delta is who wrote the lofreq snippet."""
    sub = df[(df["hardware"] == "5080") & (df["track"] == "A") & (df["plan"].isin(["v1.25", "v1g"]))]
    g = sub.groupby(["model", "plan"])["M3"].agg(["mean", "std", "count"]).reset_index()
    # Models present in both
    in_both = set(g[g["plan"] == "v1.25"]["model"]) & set(g[g["plan"] == "v1g"]["model"])
    rows = model_order([m for m in in_both])

    fig, ax = plt.subplots(figsize=(11, 5))
    width = 0.38
    x = np.arange(len(rows))
    colors = {"v1.25": "#1a9850", "v1g": "#d73027"}
    for i, p in enumerate(["v1.25", "v1g"]):
        means = []; stds = []
        for m in rows:
            row = g[(g["model"] == m) & (g["plan"] == p)]
            means.append(row["mean"].iloc[0] if len(row) else np.nan)
            stds.append(row["std"].iloc[0] if len(row) and not pd.isna(row["std"].iloc[0]) else 0.0)
        ax.bar(x + (i - 0.5) * width, means, width, yerr=stds, capsize=2,
               color=colors[p], edgecolor="black", linewidth=0.4,
               label=("v1.25 (hand-authored)" if p == "v1.25" else "v1g (Galaxy IUC mechanical)"))
    ax.set_xticks(x)
    ax.set_xticklabels([pretty(m) for m in rows], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("M3 (mean ± std, n=3)")
    ax.set_ylim(0, 1.05)
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle=":")
    ax.legend(loc="lower left")
    ax.set_title("Figure 4. Plan-source robustness: hand-authored vs mechanically-extracted lofreq (RTX 5080, Track A)\n"
                 "Same plan, only the lofreq snippet source differs. Strong models repair noisy IUC extraction; weak models copy it literally and fail.")
    plt.tight_layout()
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()


# ────────────────────────────────────────────────────────────────────────────
# Figure 5 — hardware comparison
# ────────────────────────────────────────────────────────────────────────────
def fig5_hardware(df: pd.DataFrame, out: Path):
    # Panel A: gen_seconds vs M3 on Jetson and 5080 for v2 Track A overlapping models
    v2A = df[(df["track"] == "A") & (df["plan"] == "v2") & (df["model"] != "qwen3.6_35b-a3b") & (df["model"] != "qwen3.6_35b")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    rng = np.random.default_rng(1)
    for hw, marker, color in [("jetson", "o", "#377eb8"), ("5080", "s", "#e41a1c")]:
        sub = v2A[v2A["hardware"] == hw].groupby("model").agg({"gen_seconds": "mean", "M3": "mean"}).reset_index()
        # Jitter the M3 by ±0.01 so models that converge at 1.0 don't fully overlap
        jit = rng.uniform(-0.012, 0.012, size=len(sub))
        ax.scatter(sub["gen_seconds"], sub["M3"] + jit, s=60, marker=marker, color=color,
                   alpha=0.75, edgecolor="black", linewidth=0.4, label=hw)
        # Only label genuine outliers (M3 < 0.95). Top cluster is too crowded.
        for _, r in sub.iterrows():
            if r["M3"] < 0.95:
                ax.annotate(pretty(r["model"]), (r["gen_seconds"], r["M3"]),
                            fontsize=7, alpha=0.9,
                            xytext=(4, 4), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_xlabel("mean generation time (s, log scale)")
    ax.set_ylabel("mean M3")
    ax.set_ylim(-0.05, 1.1)
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle=":")
    ax.legend(title="Hardware", loc="lower right")
    ax.set_title("Panel A — gen_seconds vs M3 on v2 Track A")

    # Panel B: Jetson 30W vs MAXN.  We can't separate from results.csv directly,
    # but the *retried* models are identifiable by tag.  Mark them on the chart.
    MAXN_RETRIED = {
        "glm-4.7-flash", "gpt-oss_20b", "granite4",
        "llama3.3_70b-instruct-q3_K_M", "nemotron-3-nano", "olmo-3.1_32b",
    }
    jet_v2 = df[(df["hardware"] == "jetson") & (df["track"] == "A") & (df["plan"] == "v2")]
    g = jet_v2.groupby("model")["M3"].mean().reset_index()
    g["power"] = g["model"].map(lambda m: "MAXN retry" if m in MAXN_RETRIED else "30 W (original)")
    g = g.sort_values(["power", "M3"], ascending=[False, False])
    ax = axes[1]
    bars = ax.bar(range(len(g)), g["M3"],
                  color=["#984ea3" if p == "MAXN retry" else "#4daf4a" for p in g["power"]],
                  edgecolor="black", linewidth=0.4)
    ax.set_xticks(range(len(g)))
    ax.set_xticklabels([pretty(m) for m in g["model"]], rotation=40, ha="right", fontsize=7.5)
    ax.set_ylim(0, 1.1)
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle=":")
    ax.set_ylabel("mean M3 (Track A v2)")
    handles = [plt.Rectangle((0,0),1,1, color="#4daf4a"),
               plt.Rectangle((0,0),1,1, color="#984ea3")]
    ax.legend(handles, ["30 W (original sweep)", "MAXN (retry of failures)"], loc="lower right")
    ax.set_title("Panel B — Jetson v2 Track A: 30 W (originals) vs MAXN (retried failures)")

    fig.suptitle("Figure 5. Hardware is not the bottleneck", y=1.02)
    plt.tight_layout()
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()


# ────────────────────────────────────────────────────────────────────────────
# Figure 6 — error-handling robustness
# ────────────────────────────────────────────────────────────────────────────
HANDLE_ORDER = ["crash", "propagate", "partial", "recover"]
HANDLE_COLOR = {
    "crash":     "#d73027",   # red — script died with no detection
    "propagate": "#fdae61",   # orange — set -e fired, no structured detection
    "partial":   "#abd9e9",   # light blue — defensive, some/all samples skipped
    "recover":   "#1a9850",   # green — full recovery
}


PLAN_DISPLAY = {
    "PLAN":            "v2",
    "PLAN_v2_defensive": "v2_defensive",
}


def collect_inject_runs(runs_dir: Path):
    """Return DataFrame of all injection runs with handle/recover/diagnose.
    Reads per-run dirs under runs_dir, then augments with rows from any
    error_matrix_*.jsonl files at the repo root that aren't already covered
    (lets us pick up cross-platform results contributed via PR — e.g. an
    M4 run that only landed the JSONL log, not the per-run dirs)."""
    rows = []
    seen_keys = set()  # (run_id, hardware) — Anthropic run_ids repeat across hardware platforms

    # Primary source: per-run dirs (have full meta.json + score.json)
    for d in runs_dir.glob("*_track-A_seed-*"):
        sj = d / "score.json"
        if not sj.exists():
            continue
        s = json.loads(sj.read_text())
        eh = s.get("error_handling") or {}
        meta = json.loads((d / "meta.json").read_text())
        if not eh.get("inject_pattern") or eh["inject_pattern"] == "none":
            continue
        plan_name = meta.get("plan_name") or ""
        rows.append({
            "model":   meta["model"].replace(":", "_"),
            "pattern": eh["inject_pattern"],
            "target":  eh.get("inject_target", ""),
            "plan":    PLAN_DISPLAY.get(plan_name, plan_name),
            "seed":    meta.get("seed"),
            "hardware": "jetson",
            "handle":  eh.get("m_handle"),
            "recover": eh.get("m_recover"),
            "diagnose": eh.get("m_diagnose"),
            "M3":      s.get("m3_jaccard"),
        })
        seen_keys.add((meta.get("run_id"), "jetson"))

    # Augment from JSONL summary logs (covers cross-platform PRs)
    for jsonl in BENCH.glob("error_matrix_*.jsonl"):
        # Hardware tag from filename — match the platform anywhere in the
        # stem so suffixes like _m4_r2 or _a5000_addendum are recognised.
        name = jsonl.stem
        if "_m4" in name:
            hw = "m4"
        elif "_a5000" in name:
            hw = "a5000"
        else:
            hw = "jetson"
        for line in jsonl.read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if "error" in r:
                continue
            if (r.get("run_id"), hw) in seen_keys:
                continue
            cell = r.get("cell", "")
            parts = cell.split("/")
            if len(parts) < 4:
                continue
            model, plan, pattern_at_target, seed_str = parts[0], parts[1], parts[2], parts[3]
            if "@" not in pattern_at_target:
                continue
            pattern, target = pattern_at_target.split("@", 1)
            if pattern == "none":
                continue
            try:
                seed = int(seed_str.replace("seed-", ""))
            except ValueError:
                seed = None
            if r.get("m_handle") is None:
                continue
            rows.append({
                "model":   model.replace(":", "_"),
                "pattern": pattern,
                "target":  target,
                "plan":    plan,
                "seed":    seed,
                "hardware": hw,
                "handle":  r.get("m_handle"),
                "recover": r.get("m_recover"),
                "diagnose": r.get("m_diagnose"),
                "M3":      r.get("M3"),
            })
            if r.get("run_id"):
                seen_keys.add((r["run_id"], hw))
    return pd.DataFrame(rows)


def _modal_handle(series):
    """Most common handle category in a small (n=3) seed group."""
    if len(series) == 0:
        return None
    return series.value_counts().idxmax()


def fig6_error_handling(out: Path, runs_dir: Path = BENCH / "runs_inject"):
    """Heatmap: rows=model, cols=(pattern, target), color=modal handle category.
    One panel per recipe variant (v2 vs v2_defensive)."""
    df = collect_inject_runs(runs_dir)
    if df.empty:
        print(f"[fig6] no runs in {runs_dir} yet — skipping")
        return
    df = df[df["plan"].isin(["v2", "v2_defensive"])]
    if df.empty:
        print("[fig6] no runs have plan-tagging yet")
        return

    plans = ["v2", "v2_defensive"]
    patterns = sorted(df["pattern"].unique())
    targets  = sorted(df["target"].unique())
    cells    = [(p, t) for p in patterns for t in targets if not df[(df["pattern"]==p) & (df["target"]==t)].empty]
    models   = model_order(df["model"].unique())

    fig, axes = plt.subplots(1, len(plans), figsize=(max(8, 2*len(cells)), 1+0.45*len(models)),
                             sharey=True)
    if len(plans) == 1:
        axes = [axes]
    for ax, plan in zip(axes, plans):
        sub = df[df["plan"] == plan]
        grid = np.full((len(models), len(cells)), -1, dtype=int)
        for i, m in enumerate(models):
            for j, (p, t) in enumerate(cells):
                ms = sub[(sub["model"]==m) & (sub["pattern"]==p) & (sub["target"]==t)]
                if len(ms) == 0:
                    continue
                modal = _modal_handle(ms["handle"])
                if modal in HANDLE_ORDER:
                    grid[i, j] = HANDLE_ORDER.index(modal)
        # Plot as a categorical image
        cmap = mpl.colors.ListedColormap([HANDLE_COLOR[h] for h in HANDLE_ORDER])
        masked = np.ma.masked_less(grid, 0)
        im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=0, vmax=len(HANDLE_ORDER)-1)
        ax.set_xticks(range(len(cells)))
        ax.set_xticklabels([f"{p}\n@{t}" for p, t in cells], rotation=30, ha="right", fontsize=7)
        ax.set_yticks(range(len(models)))
        ax.set_yticklabels([pretty(m) for m in models], fontsize=8)
        ax.set_title(plan)

    handles = [plt.Rectangle((0,0),1,1, color=HANDLE_COLOR[h]) for h in HANDLE_ORDER]
    fig.legend(handles, HANDLE_ORDER, loc="lower center", ncol=len(HANDLE_ORDER), fontsize=8,
               bbox_to_anchor=(0.5, -0.02), title="modal m_handle across 3 seeds")
    fig.suptitle("Figure 6. Error-handling robustness across 7 injection patterns × {v2, v2_defensive}", y=1.02)
    plt.tight_layout()
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()


# ────────────────────────────────────────────────────────────────────────────
# Manuscript figures (rewritten narrative): 5080 implementer gradient and
# qwen3.6:27b cross-platform error injection.
# ────────────────────────────────────────────────────────────────────────────

# The six "latest open-weight" implementers tested on RTX 5080 with the
# full plan gradient (Track B → v0.5 → v1 → v1g → v1.25 → v1.5 → v2).
LINEUP_5080 = [
    "qwen3.6_27b",
    "qwen3.6_35b-a3b",
    "gemma4_26b",
    "gemma4_e4b",
    "glm-4.7-flash",
    "gpt-oss_20b",
]


def fig_5080_implementer_gradient(df: pd.DataFrame, out: Path):
    """Headline figure: 6 open-weight implementers × 7 plan variants on RTX 5080.
    Cells are mean M3 across 3 seeds; grey = untested."""
    sub = df[(df["hardware"] == "5080") & (df["model"].isin(LINEUP_5080))]
    pivot = sub.pivot_table(index="model", columns="plan_col", values="M3", aggfunc="mean")
    cnt   = sub.pivot_table(index="model", columns="plan_col", values="M3", aggfunc="count")
    cols  = [c for c in PLAN_COLS if c in pivot.columns]
    pivot = pivot.reindex(index=LINEUP_5080, columns=cols)
    cnt   = cnt.reindex(index=LINEUP_5080, columns=cols)

    fig, ax = plt.subplots(figsize=(1.4*len(cols)+2, 0.55*len(LINEUP_5080)+2))
    cmap = mpl.cm.RdYlGn
    cmap.set_bad(color="#dddddd")
    im = ax.imshow(np.ma.masked_invalid(pivot.values), aspect="auto",
                   cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, fontsize=11)
    ax.set_yticks(range(len(LINEUP_5080)))
    ax.set_yticklabels([pretty(m) for m in LINEUP_5080], fontsize=10)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            n = cnt.values[i, j]
            if pd.isna(v):
                continue
            ax.text(j, i, f"{v:.2f}\n(n={int(n)})",
                    ha="center", va="center", fontsize=8.5,
                    color="black" if 0.25 < v < 0.85 else "white")
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("mean M3 (Jaccard)", fontsize=10)
    ax.set_xlabel("Plan variant (lean → detailed; B = no plan)", fontsize=11)
    ax.set_title("Six latest open-weight implementers × seven Opus-authored plans on RTX 5080",
                 fontsize=11)
    plt.tight_layout()
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"[fig_5080_implementer_gradient] wrote {out}")


def fig_qwen3p6_27b_error_per_platform(out: Path,
                                       runs_dir: Path = BENCH / "runs_inject"):
    """qwen3.6:27b vs claude-opus-4-7 across Jetson, M4, A5000 on v2 and
    v2_defensive plans. Two heatmap rows per panel (model x plan), one
    panel per hardware platform."""
    df_inject = collect_inject_runs(runs_dir)
    if df_inject.empty:
        print(f"[fig_qwen3p6_27b_error_per_platform] no runs in {runs_dir} — skipping")
        return
    df_inject = df_inject[df_inject["plan"].isin(["v2", "v2_defensive"])]
    df_inject = df_inject[df_inject["model"].isin(["qwen3.6_27b", "claude-opus-4-7"])]

    # M4 r2 supersedes M4 r1 (the r1 batch had a wrapper bug)
    df_inject = df_inject.copy()
    # The collect_inject_runs function tags hardware as "m4" for both r1 and r2.
    # We prefer the r2 numbers since they were taken after the wrapper fix.
    # Use unique (model, plan, pattern, target, seed, hardware) — keep last.
    df_inject = df_inject.drop_duplicates(
        subset=["model","plan","pattern","target","seed","hardware"], keep="last")

    platforms = ["jetson", "m4", "a5000"]
    plat_label = {"jetson": "NVIDIA Jetson AGX Orin",
                  "m4":     "MacBook Pro M4 Pro",
                  "a5000":  "2× RTX A5000 workstation"}
    models = ["claude-opus-4-7", "qwen3.6_27b"]
    plans  = ["v2", "v2_defensive"]
    patterns = sorted(df_inject["pattern"].unique())
    targets  = sorted(df_inject["target"].unique())
    cells    = [(p, t) for p in patterns for t in targets
                if not df_inject[(df_inject["pattern"]==p) & (df_inject["target"]==t)].empty]

    fig, axes = plt.subplots(len(platforms), len(plans),
                             figsize=(max(7, 0.6*len(cells)+2), 1.3*len(platforms)*len(models)+1.5),
                             sharex=True, sharey=True)
    cmap = mpl.colors.ListedColormap([HANDLE_COLOR[h] for h in HANDLE_ORDER])

    for row, hw in enumerate(platforms):
        for col, plan in enumerate(plans):
            ax = axes[row, col]
            sub = df_inject[(df_inject["hardware"] == hw) & (df_inject["plan"] == plan)]
            grid = np.full((len(models), len(cells)), -1, dtype=int)
            for i, m in enumerate(models):
                for j, (p, t) in enumerate(cells):
                    ms = sub[(sub["model"]==m) & (sub["pattern"]==p) & (sub["target"]==t)]
                    if len(ms) == 0:
                        continue
                    modal = _modal_handle(ms["handle"])
                    if modal in HANDLE_ORDER:
                        grid[i, j] = HANDLE_ORDER.index(modal)
            masked = np.ma.masked_less(grid, 0)
            ax.imshow(masked, aspect="auto", cmap=cmap, vmin=0, vmax=len(HANDLE_ORDER)-1)
            ax.set_yticks(range(len(models)))
            ax.set_yticklabels([pretty(m) for m in models], fontsize=8)
            if row == 0:
                ax.set_title(plan, fontsize=10)
            if col == 0:
                ax.set_ylabel(plat_label[hw], fontsize=9)
            if row == len(platforms) - 1:
                ax.set_xticks(range(len(cells)))
                ax.set_xticklabels([f"{p}\n@{t}" for p, t in cells],
                                   rotation=30, ha="right", fontsize=7)

    handles = [plt.Rectangle((0,0),1,1, color=HANDLE_COLOR[h]) for h in HANDLE_ORDER]
    fig.legend(handles, HANDLE_ORDER, loc="lower center", ncol=len(HANDLE_ORDER),
               fontsize=8, bbox_to_anchor=(0.5, -0.04),
               title="modal m_handle across 3 seeds")
    fig.suptitle("qwen3.6:27b matches claude-opus-4-7 across platforms on the error matrix",
                 fontsize=11, y=1.005)
    plt.tight_layout()
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"[fig_qwen3p6_27b_error_per_platform] wrote {out}")


# ────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--fig", type=int, choices=[1, 2, 3, 4, 5, 6])
    ap.add_argument("--fig1-per-platform", action="store_true",
                    help="generate one heatmap per hardware platform (Jetson, 5080, M4, A5000)")
    ap.add_argument("--manuscript", action="store_true",
                    help="generate the new manuscript figures (5080 implementer gradient + qwen3.6:27b cross-platform error injection)")
    args = ap.parse_args()

    if not args.all and args.fig is None and not args.fig1_per_platform and not args.manuscript:
        ap.error("specify --all, --fig N, --fig1-per-platform, or --manuscript")

    df = load_data()
    funcs = {
        1: lambda: fig1_heatmap(df,        FIGS / "fig1_headline_heatmap.png"),
        2: lambda: fig2_intermediates(df,  FIGS / "fig2_v1_cliff_repair.png"),
        3: lambda: fig3_track_b(df,        FIGS / "fig3_plan_value.png"),
        4: lambda: fig4_v1g(df,            FIGS / "fig4_v1g_robustness.png"),
        5: lambda: fig5_hardware(df,       FIGS / "fig5_hardware.png"),
        6: lambda: fig6_error_handling(    FIGS / "fig6_error_handling.png"),
    }
    if args.all:
        for i in (1, 2, 3, 4, 5, 6):
            print(f"[fig{i}] generating…")
            funcs[i]()
        fig1_per_platform(df, FIGS)
        print("done")
    elif args.fig1_per_platform:
        fig1_per_platform(df, FIGS)
    elif args.manuscript:
        fig_5080_implementer_gradient(df, FIGS / "ms_fig2_5080_gradient.png")
        fig_qwen3p6_27b_error_per_platform(FIGS / "ms_fig4_qwen3p6_27b_error.png")
    else:
        funcs[args.fig]()
        print(f"wrote figures/fig{args.fig}_*.png")


if __name__ == "__main__":
    main()
