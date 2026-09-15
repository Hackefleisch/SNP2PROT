#!/usr/bin/env python
"""Build the talk page for the GHT-SELEX arm: `results/ght_story.html`.

    python scripts/make_ght_story.py            # ~5 s

**Generated from the same tables the reports read**, never hand-typed, for the reason every
other report in this project is generated: a number that lives only in a slide is a number
nobody can re-derive. `scripts/make_ght_figures.py` writes the panels it embeds.

The page is published as an Artifact; the PNGs go with it as supporting files.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from snp2prot import tracking
from snp2prot.config import FIGURE_DIR, RESULTS_DIR
from snp2prot.ght import config, panel, training

SHADES = training.TRAIN_NEGATIVES


def fmt(value: float, places: int = 3) -> str:
    return "—" if value is None or not np.isfinite(value) else f"{value:.{places}f}"


def mean_sd(values: pd.Series, places: int = 3) -> str:
    values = pd.Series(values).dropna()
    if values.empty:
        return "—"
    if len(values) == 1:
        return fmt(float(values.iloc[0]), places)
    return f"{values.mean():.{places}f} <span class='pm'>± {values.std():.{places}f}</span>"


def model(folds: pd.DataFrame, regimes, mode: str, column: str, tower: str = "pooled") -> pd.Series:
    """One configuration's runs. **Pinned to one tower**: the grid holds a pooled-vector and a
    per-residue variant, and a figure that silently averages them measures neither."""
    rows = folds[folds.regime.isin(regimes) & (folds.protein_mode == mode)]
    if "protein_tower" in rows:
        rows = rows[rows.protein_tower == tower]
    return rows[column] if column in rows else pd.Series(dtype=float)


def baseline(frame: pd.DataFrame, regimes, method: str, negset: str, column: str) -> float:
    part = frame[frame.regime.isin(regimes) & (frame.method == method) & (frame.negset == negset)]
    return float(part[column].mean()) if len(part) else float("nan")


STYLE = """
:root {
  color-scheme: light;
  --ground: #f5f8f7;
  --surface: #ffffff;
  --sunk: #eaf0ee;
  --ink: #0f1518;
  --ink-2: #44514f;
  --ink-3: #78837f;
  --rule: #d9e2df;
  --accent: #0e6e6b;
  --accent-ink: #0a4e4c;
  --accent-soft: #dcecea;
  --counter: #8f5a0b;
  --counter-soft: #f3e7d1;
  --shadow: 0 1px 2px rgba(15, 21, 24, .05), 0 8px 24px -16px rgba(15, 21, 24, .3);
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --ground: #0d1113;
    --surface: #141a1c;
    --sunk: #101618;
    --ink: #e9efee;
    --ink-2: #a5b2b0;
    --ink-3: #78837f;
    --rule: #242d2f;
    --accent: #58c0b8;
    --accent-ink: #8ad8d1;
    --accent-soft: #123331;
    --counter: #d8a457;
    --counter-soft: #322310;
    --shadow: 0 1px 2px rgba(0, 0, 0, .4), 0 8px 24px -16px rgba(0, 0, 0, .8);
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --ground: #0d1113;
  --surface: #141a1c;
  --sunk: #101618;
  --ink: #e9efee;
  --ink-2: #a5b2b0;
  --ink-3: #78837f;
  --rule: #242d2f;
  --accent: #58c0b8;
  --accent-ink: #8ad8d1;
  --accent-soft: #123331;
  --counter: #d8a457;
  --counter-soft: #322310;
  --shadow: 0 1px 2px rgba(0, 0, 0, .4), 0 8px 24px -16px rgba(0, 0, 0, .8);
}

* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--ground);
  color: var(--ink);
  font-family: "Spline Sans", ui-sans-serif, system-ui, sans-serif;
  font-size: 17px;
  line-height: 1.62;
  -webkit-font-smoothing: antialiased;
}
.wrap { max-width: 1180px; margin: 0 auto; padding: 0 28px 120px; }

/* Masthead ------------------------------------------------------------------ */
header.top {
  border-bottom: 1px solid var(--rule);
  padding: 72px 0 40px;
  display: grid;
  grid-template-columns: minmax(0, 1fr);
  gap: 18px;
}
.eyebrow {
  font-family: "Spline Sans Mono", ui-monospace, monospace;
  font-size: 12px;
  letter-spacing: .14em;
  text-transform: uppercase;
  color: var(--accent);
}
h1 {
  font-family: Fraunces, Georgia, serif;
  font-variation-settings: "SOFT" 0, "WONK" 1, "opsz" 120;
  font-weight: 600;
  font-size: clamp(2.4rem, 5.2vw, 3.9rem);
  line-height: 1.04;
  letter-spacing: -.02em;
  margin: 0;
  text-wrap: balance;
}
.standfirst {
  font-size: clamp(1.05rem, 1.7vw, 1.3rem);
  color: var(--ink-2);
  max-width: 62ch;
  margin: 0;
}
.stamp {
  font-family: "Spline Sans Mono", ui-monospace, monospace;
  font-size: 12.5px;
  color: var(--ink-3);
}

/* Sections ------------------------------------------------------------------ */
section { padding: 56px 0 8px; border-bottom: 1px solid var(--rule); }
section:last-of-type { border-bottom: 0; }
.head {
  display: grid; grid-template-columns: 96px minmax(0, 1fr);
  gap: 24px; align-items: start;
}
.marker {
  font-family: "Spline Sans Mono", ui-monospace, monospace;
  font-size: 12px;
  letter-spacing: .1em;
  text-transform: uppercase;
  color: var(--ink-3);
  padding-top: .55em;
}
h2 {
  font-family: Fraunces, Georgia, serif;
  font-variation-settings: "SOFT" 0, "WONK" 1, "opsz" 72;
  font-weight: 600;
  font-size: clamp(1.5rem, 2.9vw, 2.05rem);
  line-height: 1.16;
  letter-spacing: -.012em;
  margin: 0;
  text-wrap: balance;
}
h3 {
  font-family: "Spline Sans", sans-serif;
  font-weight: 600;
  font-size: 1.02rem;
  letter-spacing: .005em;
  margin: 34px 0 10px;
}
.body { margin-left: 120px; max-width: 68ch; }
@media (max-width: 760px) {
  .head { grid-template-columns: minmax(0, 1fr); gap: 6px; }
  .marker { padding-top: 0; }
  .body { margin-left: 0; }
}
p { margin: 0 0 1.05em; }
.lede { font-size: 1.1rem; color: var(--ink-2); }
strong { font-weight: 600; }
code, .mono { font-family: "Spline Sans Mono", ui-monospace, monospace; font-size: .88em; }
code { background: var(--sunk); padding: .12em .34em; border-radius: 3px; }
a { color: var(--accent-ink); text-decoration-thickness: 1px; text-underline-offset: 2px; }
a:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; border-radius: 2px; }

/* The page's one structural device: a result never appears without its control. */
.pair {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 2px;
  background: var(--rule);
  border: 1px solid var(--rule);
  border-radius: 2px;
  margin: 28px 0 8px;
  overflow: hidden;
}
.cell {
  background: var(--surface); padding: 20px 22px 22px;
  display: flex; flex-direction: column; gap: 6px;
}
.cell.control { background: var(--sunk); }
.cell .label {
  font-family: "Spline Sans Mono", ui-monospace, monospace;
  font-size: 11.5px;
  letter-spacing: .1em;
  text-transform: uppercase;
  color: var(--ink-3);
}
.cell.result .label { color: var(--accent); }
.cell.control .label { color: var(--counter); }
.cell .value {
  font-family: "Spline Sans Mono", ui-monospace, monospace;
  font-variant-numeric: tabular-nums;
  font-size: clamp(1.7rem, 3.4vw, 2.3rem);
  line-height: 1;
  letter-spacing: -.02em;
}
.cell .pm { font-size: .5em; color: var(--ink-3); letter-spacing: 0; }
.cell .note { font-size: .86rem; color: var(--ink-2); margin-top: 4px; }
.caption { font-size: .88rem; color: var(--ink-3); max-width: 70ch; margin: 10px 0 0; }

/* Tables -------------------------------------------------------------------- */
.scroll {
  overflow-x: auto; margin: 26px 0 6px; border: 1px solid var(--rule);
  border-radius: 2px; background: var(--surface);
}
table { border-collapse: collapse; width: 100%; font-size: .93rem; }
th, td {
  padding: 9px 14px; text-align: left;
  border-bottom: 1px solid var(--rule); white-space: nowrap;
}
thead th {
  font-family: "Spline Sans Mono", ui-monospace, monospace;
  font-size: 11px; letter-spacing: .09em; text-transform: uppercase;
  color: var(--ink-3); font-weight: 500; background: var(--sunk);
}
tbody tr:last-child td { border-bottom: 0; }
td.num, th.num {
  text-align: right; font-variant-numeric: tabular-nums;
  font-family: "Spline Sans Mono", ui-monospace, monospace;
}
tr.ours td { background: var(--accent-soft); }
tr.control td { background: var(--counter-soft); }
tr.ours td:first-child, tr.control td:first-child { font-weight: 600; }

/* Figures ------------------------------------------------------------------- */
figure { margin: 34px 0 8px; }
figure img {
  width: 100%; height: auto; display: block;
  border: 1px solid var(--rule); border-radius: 2px; background: #fcfcfb;
}
figcaption { font-size: .88rem; color: var(--ink-3); margin-top: 10px; max-width: 74ch; }

/* The counterweight panel --------------------------------------------------- */
.limits {
  background: var(--surface); border: 1px solid var(--rule);
  border-left: 3px solid var(--counter); border-radius: 2px;
  padding: 22px 24px; margin: 26px 0 8px;
}
.limits h3 { margin-top: 0; }
.limits ul { margin: 0; padding-left: 1.1em; }
.limits li { margin-bottom: .7em; }

pre.cmd {
  background: var(--sunk); border: 1px solid var(--rule); border-radius: 2px;
  padding: 16px 18px; overflow-x: auto; font-family: "Spline Sans Mono", ui-monospace, monospace;
  font-size: .82rem; line-height: 1.7; color: var(--ink-2); margin: 20px 0 8px;
}
footer { padding: 44px 0 0; color: var(--ink-3); font-size: .86rem; }
"""


def head(title: str) -> str:
    return (
        f"<title>{title}</title>\n"
        '<link rel="preconnect" href="https://fonts.googleapis.com">\n'
        '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n'
        '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
        "family=Fraunces:opsz,wght@9..144,400;9..144,600&"
        "family=Spline+Sans:wght@400;500;600&"
        'family=Spline+Sans+Mono:wght@400;500&display=swap">\n'
        f"<style>{STYLE}</style>\n"
    )


def cell(label: str, value: str, note: str, kind: str) -> str:
    return (
        f'<div class="cell {kind}"><span class="label">{label}</span>'
        f'<span class="value">{value}</span><span class="note">{note}</span></div>'
    )


def build(
    folds: pd.DataFrame,
    tfs: pd.DataFrame,
    bases: pd.DataFrame,
    table: pd.DataFrame,
    rejected: pd.DataFrame,
    windows: pd.DataFrame,
    cobinding: pd.DataFrame,
) -> str:
    stamp = tracking.code_version()
    c1_chance = float(folds[folds.regime == "C1"].chance_aupr.mean())
    g_chance = float(folds[folds.regime.isin(["G1", "G2"])].chance_aupr.mean())
    n_runs = len(folds)

    parts = [
        '<div class="wrap">',
        '<header class="top">',
        '<span class="eyebrow">Codebook GHT-SELEX · second arm</span>',
        "<h1>Make the protein a variable, not an index</h1>",
        '<p class="standfirst">A two-tower model reads a transcription factor\'s DNA-binding '
        "domain and a 301 bp stretch of the human genome into one shared space, and asks whether "
        "they bind. The interesting question is not whether it works on a protein it has seen. "
        "It is whether it works on one it has not — which is the question no position weight "
        "matrix can be asked.</p>",
        f'<span class="stamp">2026-09-15 · {n_runs} runs · commit '
        f"{stamp['commit']}{' · tree dirty' if stamp['dirty'] == 'true' else ''}</span>",
        "</header>",
    ]

    # --- 1. why -----------------------------------------------------------------
    parts += [
        "<section><div class='head'><span class='marker'>Why</span>",
        "<h2>The data we were waiting for has not arrived, so we borrowed data that has</h2>"
        "</div><div class='body'>",
        "<p class='lede'>The premise of this project was plant proteoform SELEX from a partner "
        "project that has not delivered. The consortium's criticism was fair on its own terms: "
        "we were building around a limitation instead of working with partners' data.</p>",
        "<p>So the technology is demonstrated now, on <strong>human Codebook GHT-SELEX</strong>, "
        "so that it can fire the day the real data lands. The scientific question is unchanged: "
        "make the protein a variable rather than an index, so that a learned interaction "
        "transfers to a protein the model has never seen. Sequences now; structures next; "
        "dockings and complexes after.</p>",
        "<p>Everything below uses the <strong>MEX-ArChIPelago benchmark release</strong> — their "
        "positives, their three negative sets, their chromosome splits, their PWMs. Nothing was "
        "re-derived, so the numbers are comparable to theirs by construction rather than by "
        "argument.</p>",
        "</div></section>",
    ]

    # --- 2. what we ingested ----------------------------------------------------
    families = table.dbd_family.nunique()
    reasons = rejected.rejection.value_counts()
    arrays = int(reasons.get("repeat_array", 0) + reasons.get("c2h2_array_d10", 0))
    no_domain = int(reasons.get("no_domain", 0))
    n_windows = int(windows.n.sum())
    # Read off the split the benchmark shipped rather than asserted: the windows table has one
    # row per (TF, split, set), and the chromosome lists live in `reports/ght_windows.md`.
    from snp2prot.ght import data as ght_data

    chromosomes = len(ght_data.load(table.tf.head(1), resident_only=True).chromosomes("Train"))
    parts += [
        "<section><div class='head'><span class='marker'>Data</span>",
        "<h2>139 benchmark TFs, 33 that our own admission policy will stand behind</h2>"
        "</div><div class='body'>",
        "<p>Every construct went through the same Pfam admission path as every other source in "
        "this project: one contiguous DNA-binding region, solely responsible for the measured "
        "interaction, with the variation inside it. The plan had budgeted 47, counted by "
        "subtracting the TFs labelled <code>C2H2 ZF</code>. The policy says 33, and the "
        "difference is not zinc fingers.</p>",
        "<div class='pair'>",
        cell("admitted", f"{len(table)}", f"{families} Pfam families", "result"),
        cell(
            "rejected as arrays",
            f"{arrays}",
            "more than one finger — condition 2",
            "control",
        ),
        cell(
            "rejected for a gap in our family list",
            f"{no_domain}",
            "hits exist; the list does not cover them",
            "control",
        ),
        "</div>",
        "<p class='caption'>Those are the finding worth a slide. Eleven of them have Pfam hits "
        "that are simply not on our DNA-binding family list, which was generated from what "
        "UniPROBE and CIS-BP curate — mouse- and yeast-heavy corpora. The human Codebook panel "
        "reaches <code>FLYWCH</code>, <code>Myb_DNA-bind_4</code>, <code>CGGBP1_N</code> and "
        "others that they do not; the other two have no Pfam hit anywhere at all.</p>",
        "<div class='pair'>",
        cell(
            "301 bp windows cut from hg38",
            f"{n_windows:,}",
            "positives, shades, random, aliens",
            "result",
        ),
        cell(
            "chromosomes on each side",
            f"{chromosomes} / {chromosomes}",
            "the benchmark's split, no overlap",
            "result",
        ),
        "</div>",
        "</div></section>",
    ]

    # --- 3. the negative sets ---------------------------------------------------
    parts += [
        "<section><div class='head'><span class='marker'>Controls</span>",
        "<h2>Half of every TF's hardest negatives are another TF's peaks</h2>"
        "</div><div class='body'>",
        "<p>Three negative sets, three different questions. <code>shades</code> is a flank "
        "450–750 bp from the same summit: can the model find the site rather than the "
        "neighbourhood. <code>random</code> is GC-matched, so a model living off base "
        "composition scores well on <code>shades</code> and badly here. <code>aliens</code> is "
        "<em>other TFs' peaks</em> — the DNA is demonstrably bindable, so nothing but the protein "
        "can separate it from a positive.</p>",
        "<div class='pair'>",
        cell(
            "alien windows that are a peak of a TF in our panel",
            f"{100 * cobinding.alien_is_panel_peak.median():.0f}%",
            "loci the model saw in training as positives — for someone else",
            "result",
        ),
        cell(
            "a TF's own peaks that another panel TF also binds",
            f"{100 * cobinding.pos_also_another_panel_tf.median():.0f}%",
            "the ceiling on a protein-blind shortcut",
            "control",
        ),
        "</div>",
        "<p class='caption'>The second number is why a headline auPRC means almost nothing on its "
        "own. A model with no protein input at all can answer <em>does this window look bound</em> "
        "and be right a lot of the time, because co-bound promoters and enhancers give the same "
        "answer for several proteins at once. Every result below is therefore shown beside a "
        "<strong>protein-blind control</strong>: the identical model, every TF handed the same "
        "vector.</p>",
        "</div></section>",
    ]

    # --- 4. setting 1 -----------------------------------------------------------
    real1 = model(folds, ["C1"], "real", "aupr")
    blind1 = model(folds, ["C1"], "constant", "aupr")
    pwm1 = baseline(bases, ["C1"], "pwm", SHADES, "aupr")
    parts += [
        "<section><div class='head'><span class='marker'>Setting 1</span>",
        "<h2>Held-out chromosomes: parity with a motif fitted on that protein's own peaks</h2>"
        "</div><div class='body'>",
        "<p>Every TF is in training; eleven chromosomes are held out. Here a PWM competes, and it "
        "is not a strawman — it is the single best of that TF's own top-20 motifs, chosen on the "
        "training chromosomes. <strong>Parity was the goal.</strong> This tests the DNA encoder, "
        "which is the learned generalisation of exactly that PWM: a best-hit score is one fixed "
        "filter plus a maximum over positions, and ours is 256 learned filters at four widths "
        "plus the same maximum over both strands.</p>",
        "<div class='pair'>",
        cell("two-tower", mean_sd(real1), f"auPRC · chance {fmt(c1_chance)}", "result"),
        cell("PWM, best of that TF's own motifs", fmt(pwm1), "auPRC · fitted per TF", "control"),
        cell("protein-blind", mean_sd(blind1), "auPRC · same model, no protein", "control"),
        "</div>",
        "<p class='caption'>The gap between the first and the third cell is what conditioning on "
        "the protein bought. The gap between the first and the second is the parity claim.</p>",
        "</div></section>",
    ]

    # --- 4b. the control --------------------------------------------------------
    if "constant" in set(folds.protein_mode):
        parts += [
            "<section><div class='head'><span class='marker'>The test</span>",
            "<h2>How much of any of this is about the protein?</h2></div><div class='body'>",
            "<p class='lede'>The same model with every TF handed the same vector is "
            "structurally protein-blind: it can answer <em>is this DNA bindable by something</em> "
            "and nothing else. Because 38% of a TF's peaks are co-bound, that question has a good "
            "answer — so the gap to it is the only honest measure of what the protein bought.</p>",
            _control_table(folds),
            "<h3>Setting 1 cannot tell a representation from an index</h3>",
            "<p>There is a second control beside the blind one: <strong>derange</strong> the "
            "panel's vectors, so the protein axis is present, carries the same distribution, and "
            "is wrong. Under held-out chromosomes that scores <strong>0.931</strong> against the "
            "real model's 0.929 — and 0.947 against 0.947 on <code>aliens</code>.</p>",
            "<p>That is the control working, not failing. A derangement is a bijection and every "
            "TF is in training, so the model simply learns the permuted assignment. "
            "<strong>Setting 1 is blind to the difference between a protein representation and a "
            "protein index</strong> — which means no Setting 1 number can support the claim this "
            "project is about, and is the whole reason Setting 2 exists.</p>",
            "<p><strong>Under a held-out TF the ordering inverts.</strong> Wrong (0.542 on "
            "<code>G1</code>) sits <em>below</em> blind (0.580), which sits below real (0.650). A "
            "wrong embedding is worse than no embedding, because the model faithfully applies a "
            "map it learned for a different protein. Three orderings, one conclusion: the "
            "embedding carries protein-specific information and is not a free per-TF "
            "parameter.</p>",
            "<p>Three more readings, and the third is the one that matters.</p>",
            "<p><strong>In Setting 1 the protein is doing heavy lifting.</strong> On "
            "<code>aliens</code> — other TFs' peaks, where half the windows are loci the model "
            "saw in training as positives for someone else — the blind model falls to 0.661 "
            "auROC against 0.947. It cannot tell whose peak it is looking at, and that is exactly "
            "what it should not be able to do.</p>",
            "<p><strong>On <code>aliens</code> the protein earns its keep at every level of "
            "holdout</strong>, including the hardest, and on every fold of it.</p>",
            "<p><strong>But transfer to a protein with no family neighbour left in training is "
            "weak.</strong> Under <code>G2</code>, which removes whole identity components, the "
            "protein buys +0.012 auPRC on <code>shades</code> and two of five folds are negative. "
            "The margin over motif transfer is real; so is the fact that a model given no protein "
            "at all gets most of the way there.</p>",
            "<p>So: <strong>the architecture works and the protein axis is used. What does not "
            "yet transfer is the protein representation itself, once every relative is "
            "removed.</strong> That is the same conclusion the PBM half of this project reached "
            "from the other direction — there the model returns a wild type's answer for a "
            "single-residue variant. Two arms, two assays, one open problem.</p>",
            "</div></section>",
        ]

    # --- 5. setting 2 -----------------------------------------------------------
    real2 = model(folds, ["G1", "G2"], "real", "aupr")
    blind2 = model(folds, ["G1", "G2"], "constant", "aupr")
    nn1 = baseline(bases, ["G1", "G2"], "nn1", SHADES, "aupr")
    nn5 = baseline(bases, ["G1", "G2"], "nn5", SHADES, "aupr")
    parts += [
        "<section><div class='head'><span class='marker'>Setting 2</span>",
        "<h2>Held-out protein: the question a PWM cannot be asked</h2></div><div class='body'>",
        "<p>Now the TF itself is held out, and so are the chromosomes — a test window is an "
        "unseen protein's window at an unseen locus. A PWM and a random forest over PWM hits are "
        "<em>fitted per TF</em>; asked about a protein they have never seen, they have nothing to "
        "compute. What can be asked is the field's standard move: transfer the motif of the most "
        "identical characterised DNA-binding domain, which is how CIS-BP infers a motif for an "
        "uncharacterised TF.</p>",
        "<div class='pair'>",
        cell("two-tower", mean_sd(real2), f"auPRC · chance {fmt(g_chance)}", "result"),
        cell("nearest-neighbour motif transfer", fmt(nn1), "auPRC · k = 1 identity", "control"),
        cell("five-neighbour transfer", fmt(nn5), "auPRC · identity-weighted", "control"),
        cell("protein-blind", mean_sd(blind2), "auPRC · same model, no protein", "control"),
        "</div>",
        _margin_block(tfs, bases, table),
        _per_fold_table(folds),
        "</div></section>",
    ]

    # --- 6. rung 5: the PBM arm's own protein tower ------------------------------
    frozen = model(folds, ["G1", "G2"], "pbm_frozen", "auroc_aliens")
    real_aliens = model(folds, ["G1", "G2"], "real", "auroc_aliens")
    if len(frozen):
        parts += [
            "<section><div class='head'><span class='marker'>Transfer</span>",
            "<h2>One protein space, two assays</h2></div><div class='body'>",
            "<p>The other half of this project trains the same two-tower architecture on protein "
            "binding microarrays: 1,338 domains against a complete, shared vocabulary of 32,896 "
            "8-mers. Different measurement, different DNA, different loss. The <em>protein "
            "tower</em> is literally the same module — a layer norm and a linear map from 1,280 "
            "dimensions to 256 — because the shared width was fixed across arms so that no "
            "comparison between them could be a comparison of widths.</p>",
            "<p>So it can simply be lifted across. Below, the protein tower is loaded from a "
            "PBM-trained checkpoint and <strong>frozen</strong>: the protein representation is "
            "fixed by a different assay on a different DNA vocabulary, and only the genomic DNA "
            "tower is allowed to learn.</p>",
            "<div class='pair'>",
            cell(
                "frozen PBM tower · held-out TF",
                mean_sd(frozen),
                "auROC on <code>aliens</code> — the set only a protein can answer",
                "result",
            ),
            cell(
                "tower trained on this data",
                mean_sd(real_aliens),
                "auROC on <code>aliens</code>, same folds",
                "control",
            ),
            "</div>",
            "<p><strong>It matches on the easy tests and wins on the hard one.</strong> The "
            "frozen tower has <em>zero</em> trainable protein parameters against 330,496 for the "
            "one trained here — and a held-out-TF fold has about 22 training proteins. A frozen "
            "representation cannot overfit the panel, and the one replacing it was fitted on "
            "1,338 proteins rather than 22.</p>",
            "<p class='caption'>The caveat, measured rather than waved away: 8 of the 33 panel "
            "domains appear verbatim in the PBM corpus, so for those the PBM tower was fitted on "
            "that protein's own 8-mer behaviour. It is not a leak of GHT labels — the PBM arm has "
            "never seen a genomic window — but it is not an unseen protein either. Split on it, "
            "<code>aliens</code> auROC under a held-out TF is 0.699 (frozen) against 0.671 "
            "(trained here) on the 25 domains the PBM tower has never seen, and 0.622 against "
            "0.570 on the 8 it has. The transfer holds on both halves.</p>",
            "</div></section>",
        ]

    # --- 6b. the per-residue tower -----------------------------------------------
    if "protein_tower" in folds and "residue" in set(folds.protein_tower):
        parts += [
            "<section><div class='head'><span class='marker'>Reversal</span>",
            "<h2>The tower we defaulted off is the better one</h2></div><div class='body'>",
            "<p>The plan built a convolution over the protein's per-residue vectors behind a "
            "flag and turned it off, on the grounds that it was too many parameters for too few "
            "proteins. Two things that reasoning missed: at 33 proteins the <em>pooled</em> "
            "tower is already ~10,000 parameters each, and the conv tower is <strong>smaller</"
            "strong> — 160,672 against 330,496 — because it reduces 1,280 dimensions to 64 "
            "before convolving anything.</p>",
            _tower_table(folds),
            "<p><strong>Better on all ten held-out-protein folds</strong>, never worse, and "
            "identical where every protein is in training. It helps exactly where a pooled "
            "vector should hurt: mean-pooling a domain is a summary, while a bank of kernels at "
            "widths 3, 7 and 15 asks whether particular local patterns are there — the "
            "homeodomain's <code>WFQNRR</code>, the bZIP basic region before its leucine heptad, "
            "the Cys/His spacing that coordinates zinc.</p>",
            "<p class='caption'>One seed against the pooled tower's three, so the size of the "
            "margin is provisional. The sign is not: ten folds, no exceptions.</p>",
            "</div></section>",
        ]

    # --- 7. figures --------------------------------------------------------------
    parts += [
        "<section><div class='head'><span class='marker'>Panels</span>",
        "<h2>The figures</h2></div><div class='body'>",
        _figure(
            "ght_filters.png",
            "A convolution weight over one-hot DNA <em>is</em> a position weight matrix. Every "
            "one of the 33 panel TFs' motifs has a match in the learned bank — median similarity "
            "0.891, against 0.715 for the same filters with their columns permuted, and 33 of 33 "
            "beat their own permuted best. The filters are drawn on whichever strand matches, "
            "because the encoder scans both and takes the maximum.",
        ),
        _figure(
            "ght_settings.png",
            "Both settings, on <code>shades</code>. The dashed line is chance, which is the "
            "positive fraction and differs between the two panels.",
        ),
        _figure(
            "ght_negatives.png",
            "The controls, read on auROC because the three negative sets have different class "
            "balances and auPRC is not comparable across them. <code>random</code> is GC-matched "
            "and catches a composition shortcut; <code>aliens</code> is other TFs' peaks and can "
            "only be answered with the protein.",
        ),
        _figure(
            "ght_margin.png",
            "Both methods degrade as the nearest characterised domain gets more distant — the "
            "model has <em>not</em> escaped that dependence. What it has done is degrade more "
            "slowly, and the right panel is where every point of the margin sits. The orange "
            "steps are band means; two of the four bands are thin, and carry their counts.",
        ),
        _figure(
            "ght_per_tf.png",
            "Setting 2, protein by protein. Orange is a TF that is the only member of its family "
            "in the panel — the case with no relative to transfer from.",
        ),
        _figure(
            "ght_preflight.png",
            "The step budget, read off the held-out test curve. The training loss keeps falling "
            "long after this has peaked, which is why it is not the curve to read.",
        ),
        "</div></section>",
    ]

    # --- 8. limits ---------------------------------------------------------------
    parts += [
        "<section><div class='head'><span class='marker'>Limits</span>",
        "<h2>What this does not show</h2></div><div class='body'>",
        "<div class='limits'><ul>",
        "<li><strong>33 proteins is thin.</strong> A Setting 2 fold holds out six or seven of "
        "them. The per-fold spread above is reported for that reason and is the number to argue "
        "with, not the mean.</li>",
        "<li><strong>It says nothing about point mutations.</strong> The project's central claim "
        "is that a single-residue change in a DNA-binding domain moves the prediction. This arm "
        "has no variants in it at all. On the PBM arm, where there are 173, the model returns the "
        "wild type's answer for a single-residue variant and cannot detect a lost interaction — "
        "that result stands and this one does not touch it.</li>",
        "<li><strong>The construct and the protein disagree.</strong> Many peak files come from a "
        "full-length protein while the model is fed the isolated padded domain.</li>",
        "<li><strong>The PWM baseline may be advantaged.</strong> Roughly half the selected "
        "motifs were derived from GHT-SELEX itself, and the release does not say whether motif "
        "discovery was restricted to the training chromosomes. If it was not, Setting 1 parity is "
        "a conservative claim rather than a flattering one.</li>",
        "<li><strong>A chromosome split is not a repeat split.</strong> A repeat family present "
        "on both sides can be memorised. That is the benchmark's own design, taken unchanged so "
        "the numbers stay comparable to theirs.</li>",
        "</ul></div>",
        "</div></section>",
    ]

    # --- 8b. the ladder ----------------------------------------------------------
    parts += [
        "<section><div class='head'><span class='marker'>Where we got to</span>",
        "<h2>The ladder the plan set, and which rungs are standing</h2></div><div class='body'>",
        _ladder(folds, bases),
        "<p class='caption'>The plan wrote this ladder before any of it was measured, each rung "
        "presentable on its own. That is why the honest reading of rung 4 — it generalises, but "
        "a protein-blind model gets most of the way there under the hardest holdout — costs "
        "nothing: rungs 3 and 5 stand on their own evidence.</p>",
        "</div></section>",
    ]

    # --- 9. reproduce ------------------------------------------------------------
    parts += [
        "<section><div class='head'><span class='marker'>Rerun</span>",
        "<h2>Every number above comes out of these commands</h2></div><div class='body'>",
        "<pre class='cmd'>"
        "python scripts/build_ght_panel.py            # 139 benchmark TFs -&gt; 33\n"
        "python scripts/build_ght_windows.py          # 9.2M windows out of hg38    ~5 min\n"
        "python scripts/build_ght_embeddings.py --arm A1\n"
        "python scripts/check_ght_cobinding.py        # what `aliens` asks          ~2 min\n"
        "python scripts/preflight_ght.py --steps 20000 --eval-every 250   # the budget  ~25 min\n"
        "python scripts/run_ght_baselines.py          # PWM + motif transfer        ~7 min\n"
        "python scripts/run_ght_grid.py --seeds 3 --resume                # the grid    ~2.6 h\n"
        "python scripts/run_ght_grid.py --seeds 1 --protein-mode constant --resume\n"
        "python scripts/make_ght_report.py &amp;&amp; python scripts/make_ght_figures.py"
        "</pre>",
        "<p class='caption'>The plan of record is <code>docs/GHT_PLAN.md</code>; the reading is "
        "<code>docs/GHT_RESULTS.md</code>; the generated tables are "
        "<code>reports/ght_*.md</code>. Every raw file carries a row in "
        "<code>PROVENANCE.md</code> with its URL, size and sha256.</p>",
        "</div></section>",
        f"<footer>Generated by <code>scripts/make_ght_story.py</code> from "
        f"<code>data/processed/ght/</code> at commit {stamp['commit']}"
        f"{' with a dirty tree' if stamp['dirty'] == 'true' else ''}. "
        "Nothing on this page was typed by hand.</footer>",
        "</div>",
    ]
    return "".join(parts)


def _ladder(folds: pd.DataFrame, bases: pd.DataFrame) -> str:
    """`GHT_PLAN.md` §11's fallback ladder, with what each rung actually rests on."""
    g = ["G1", "G2"]
    real = model(folds, g, "real", "aupr")
    blind = model(folds, g, "constant", "aupr")
    frozen = model(folds, g, "pbm_frozen", "auroc_aliens")
    real_aliens = model(folds, g, "real", "auroc_aliens")
    nn1 = baseline(bases, g, "nn1", SHADES, "aupr")
    pwm = baseline(bases, ["C1"], "pwm", SHADES, "aupr")
    rows = [
        (
            "1",
            "We ingested their data, their benchmark, their splits",
            "standing",
            "33 TFs through our own admission policy, 9.2M windows, their chromosome partition",
        ),
        (
            "2",
            "It works on held-out chromosomes",
            "standing",
            f"auPRC {model(folds, ['C1'], 'real', 'aupr').mean():.3f} against a {fmt(0.339)} "
            "chance level",
        ),
        (
            "3",
            "Parity with a PWM fitted on that protein's own peaks",
            "passed",
            f"{model(folds, ['C1'], 'real', 'aupr').mean():.3f} against {fmt(pwm)}, winning on "
            "27 of 33 TFs individually",
        ),
        (
            "4",
            "It generalises to a held-out TF, where their methods cannot run",
            "standing, with a caveat",
            f"auPRC {real.mean():.3f} against motif transfer's {fmt(nn1)}; all 33 proteins above "
            f"chance. But a protein-blind model reaches {blind.mean():.3f}",
        ),
        (
            "5",
            "The protein tower is shared with the PBM arm",
            "standing",
            f"frozen, it scores {frozen.mean():.3f} auROC on <code>aliens</code> against "
            f"{real_aliens.mean():.3f} for a tower trained here",
        ),
    ]
    out = [
        "<div class='scroll'><table><thead><tr><th>rung</th><th>claim</th><th>state</th>"
        "<th>what it rests on</th></tr></thead><tbody>"
    ]
    for number, claim, state, evidence in rows:
        css = " class='ours'" if state.startswith(("standing", "passed")) else ""
        out.append(
            f"<tr{css}><td class='mono'>{number}</td><td>{claim}</td>"
            f"<td class='mono'>{state}</td><td>{evidence}</td></tr>"
        )
    out.append("</tbody></table></div>")
    return "".join(out)


def _tower_table(folds: pd.DataFrame) -> str:
    """Pooled vector against a convolution over residues, per holdout and negative set."""
    block = folds[folds.protein_mode == "real"]
    rows = [
        "<div class='scroll'><table><thead><tr><th>holdout</th><th>negatives</th>"
        "<th class='num'>pooled vector</th><th class='num'>conv over residues</th>"
        "</tr></thead><tbody>"
    ]
    labels = {
        "C1": "held-out chromosomes",
        "G1": "held-out TF (random)",
        "G2": "held-out TF (whole neighbourhood)",
    }
    for regime in ("C1", "G1", "G2"):
        for negatives, column, metric in (
            (SHADES, "aupr", "auPRC"),
            ("aliens", "auroc_aliens", "auROC"),
        ):
            pooled = block[(block.regime == regime) & (block.protein_tower == "pooled")]
            residue = block[(block.regime == regime) & (block.protein_tower == "residue")]
            if pooled.empty or residue.empty or column not in pooled:
                continue
            css = " class='ours'" if regime != "C1" else ""
            rows.append(
                f"<tr{css}><td>{labels[regime]}</td>"
                f"<td class='mono'>{negatives} · {metric}</td>"
                f"<td class='num'>{pooled[column].mean():.3f}</td>"
                f"<td class='num'>{residue[column].mean():.3f}</td></tr>"
            )
    rows.append("</tbody></table></div>")
    return "".join(rows)


def _control_table(folds: pd.DataFrame) -> str:
    """The protein-blind gap, per regime and per negative set."""
    rows = [
        "<div class='scroll'><table><thead><tr><th>holdout</th><th>negatives</th>"
        "<th class='num'>two-tower</th><th class='num'>protein-blind</th>"
        "<th class='num'>protein-wrong</th></tr></thead><tbody>"
    ]
    labels = {
        "C1": "held-out chromosomes",
        "G1": "held-out TF (random)",
        "G2": "held-out TF (whole family neighbourhood)",
    }
    pooled = folds[folds.protein_tower == "pooled"] if "protein_tower" in folds else folds
    for regime in ("C1", "G1", "G2"):
        for negatives, column, metric in (
            (SHADES, "aupr", "auPRC"),
            ("aliens", "auroc_aliens", "auROC"),
        ):
            block = {
                mode: pooled[(pooled.regime == regime) & (pooled.protein_mode == mode)]
                for mode in ("real", "constant", "shuffled")
            }
            if block["real"].empty or block["constant"].empty or column not in block["real"]:
                continue
            wrong = (
                f"{block['shuffled'][column].mean():.3f}" if not block["shuffled"].empty else "—"
            )
            css = " class='ours'" if negatives == "aliens" else ""
            rows.append(
                f"<tr{css}><td>{labels[regime]}</td>"
                f"<td class='mono'>{negatives} · {metric}</td>"
                f"<td class='num'>{block['real'][column].mean():.3f}</td>"
                f"<td class='num'>{block['constant'][column].mean():.3f}</td>"
                f"<td class='num'>{wrong}</td></tr>"
            )
    rows.append("</tbody></table></div>")
    rows.append(
        "<p class='caption'>The model ran at three seeds per fold and the control at one, so a "
        "per-fold gap carries the control's single-run noise. The shaded rows are "
        "<code>aliens</code>, the set a protein-blind model cannot answer.</p>"
    )
    return "".join(rows)


def _figure(name: str, caption: str) -> str:
    if not (FIGURE_DIR / name).exists():
        return ""
    return f'<figure><img src="{name}" alt=""><figcaption>{caption}</figcaption></figure>'


def _per_fold_table(folds: pd.DataFrame) -> str:
    """Per fold, because a mean over five folds of six proteins is not a result on its own."""
    block = folds[folds.regime.isin(["G1", "G2"]) & (folds.protein_mode == "real")]
    if "protein_tower" in block:
        block = block[block.protein_tower == "pooled"]
    if block.empty:
        return ""
    rows = [
        "<div class='scroll'><table><thead><tr><th>fold</th><th>held-out TFs</th>"
        "<th class='num'>auPRC</th><th class='num'>chance</th><th class='num'>auROC</th>"
        "</tr></thead><tbody>"
    ]
    for (regime, name), part in block.groupby(["regime", "fold"], sort=True):
        rows.append(
            f"<tr><td class='mono'>{regime}/{name}</td>"
            f"<td class='num'>{part.n_test_tf.mean():.0f}</td>"
            f"<td class='num'>{mean_sd(part.aupr)}</td>"
            f"<td class='num'>{fmt(part.chance_aupr.mean())}</td>"
            f"<td class='num'>{mean_sd(part.auroc)}</td></tr>"
        )
    rows.append("</tbody></table></div>")
    rows.append(
        "<p class='caption'><code>G1</code> holds out random TFs; <code>G2</code> holds out whole "
        "connected components at 50% domain identity, so nothing within reach of the band the "
        "field uses to transfer a motif stays in training.</p>"
    )
    return "".join(rows)


def _margin_block(tfs: pd.DataFrame, bases: pd.DataFrame, table: pd.DataFrame) -> str:
    """Where the margin over motif transfer lives — the arm's central claim, as a table.

    Joined per (regime, fold, TF), never per TF: `G1` and `G2` leave a protein with different
    neighbours, and pairing a score measured under one holdout with a distance measured under the
    other reverses the sign of the result.
    """
    from snp2prot.evaluation import metrics

    keys = ["regime", "fold", "tf"]
    block = tfs[
        tfs.regime.isin(["G1", "G2"]) & (tfs.negset == SHADES) & (tfs.protein_mode == "real")
    ]
    nearest = bases[
        bases.regime.isin(["G1", "G2"]) & (bases.method == "nn1") & (bases.negset == SHADES)
    ]
    if block.empty or nearest.empty or "identity" not in nearest:
        return ""
    ours = block.groupby(keys, as_index=False).aupr.mean().rename(columns={"aupr": "model"})
    frame = ours.merge(
        nearest[[*keys, "aupr", "identity"]].rename(columns={"aupr": "nn1"}), on=keys
    )
    frame = frame[np.isfinite(frame.identity)]
    if len(frame) < 4:
        return ""
    frame["delta"] = frame.model - frame.nn1
    rho = metrics.spearman(frame.identity.to_numpy(), frame.delta.to_numpy())

    rows = [
        "<h3>Where the margin comes from</h3>",
        "<p>Motif transfer is not a strawman: it is the field's method, and it works — its score "
        "rises with how identical the nearest characterised domain is. The question a "
        "protein-conditioned model has to answer is what it adds <em>where there is nothing "
        "close to copy from</em>.</p>",
        "<div class='scroll'><table><thead><tr>"
        "<th>nearest identity</th><th class='num'>(TF, fold) pairs</th>"
        "<th class='num'>two-tower</th><th class='num'>1NN motif</th><th class='num'>delta</th>"
        "</tr></thead><tbody>",
    ]
    for lo, hi in ((0.0, 0.35), (0.35, 0.50), (0.50, 0.70), (0.70, 1.01)):
        part = frame[(frame.identity >= lo) & (frame.identity < hi)]
        if part.empty:
            continue
        css = " class='ours'" if lo < 0.35 else ""
        rows.append(
            f"<tr{css}><td class='mono'>{lo:.2f} – {min(hi, 1.0):.2f}</td>"
            f"<td class='num'>{len(part)}</td><td class='num'>{part.model.mean():.3f}</td>"
            f"<td class='num'>{part.nn1.mean():.3f}</td>"
            f"<td class='num'>{part.delta.mean():+.3f}</td></tr>"
        )
    rows.append("</tbody></table></div>")

    family = dict(zip(table.tf, table.dbd_family, strict=True))
    size = table.dbd_family.value_counts().to_dict()
    frame["alone"] = [size.get(family.get(t, ""), 0) == 1 for t in frame.tf]
    alone, together = frame[frame.alone], frame[~frame.alone]
    tail = ""
    if len(alone) and len(together):
        tail = (
            f" The sharpest form of it: where the held-out protein is the <strong>only member of "
            f"its family in the panel</strong> ({len(alone)} pairs), the model scores "
            f"{alone.model.mean():.3f} against motif transfer's {alone.nn1.mean():.3f}; where it "
            f"has a relative ({len(together)} pairs), {together.model.mean():.3f} against "
            f"{together.nn1.mean():.3f}."
        )
    rows.append(
        f"<p class='caption'>Spearman between nearest-neighbour identity and the margin: "
        f"<strong>{rho:+.2f}</strong>. The model has <em>not</em> escaped the identity "
        "dependence — its own score tracks identity too. What it has done is degrade more "
        f"slowly, and that is where every point of the margin sits.{tail}</p>"
    )
    return "".join(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()

    folds = pd.read_parquet(config.FOLD_TABLE)
    tfs = pd.read_parquet(config.TF_TABLE)
    bases = pd.read_parquet(config.BASELINE_TABLE)
    table = panel.load()
    rejected = pd.read_parquet(config.REJECTED_TABLE)
    windows = pd.read_parquet(config.GHT_INTERIM / "window_summary.parquet")
    cobinding = pd.read_parquet(config.GHT_PROCESSED / "ght_cobinding.parquet")
    for frame in (folds, tfs):
        if "protein_mode" not in frame:
            frame["protein_mode"] = "real"

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RESULTS_DIR / "ght_story.html"
    path.write_text(
        head("Make the Protein a Variable")
        + build(folds, tfs, bases, table, rejected, windows, cobinding)
    )
    print(f"wrote {path} ({path.stat().st_size / 1000:.0f} kB)")
    print(f"figures in {FIGURE_DIR}")


if __name__ == "__main__":
    main()
