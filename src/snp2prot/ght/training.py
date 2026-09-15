"""Training the genomic two-tower: one arm, one fold, end to end.

`GHT_PLAN.md` §8. What this implements from it:

**A batch is balanced across proteins, not across windows.** `batch_proteins` TFs are drawn and
`windows_per_protein` of each TF's training windows are drawn under them, so every protein in the
batch contributes the same number of rows whether it has 714 peaks or 26,000. The PBM arm
shuffles the protein axis and keeps the DNA axis complete; here the DNA axis is per-TF and
disjoint, so it is the protein axis that has to be equalised by hand or KDM2A alone supplies 12%
of every step.

**Windows are drawn from the TF's own pool at its own base rate.** Not balanced positive to
negative: `shades` is already about 1:2, the calibration head is initialised at the training
pool's measured rate, and re-balancing would hand the head a prior that is not the data's.

**Positives and `shades` live on the device; `random` and `aliens` are streamed.** 356 MB against
2.4 GB, and the second pair is touched once per run.

**The step budget is measured, not inherited** (`GHT_PLAN.md` §8.1). `scripts/preflight_ght.py`
runs one fold well past where it looks done with `patience: 0` and logs validation *and* test at
every evaluation; the budget is read off the **test** curve. The PBM arm got this backwards once
and the lesson is in `configs/experiment.yaml`: its training loss reached ~0.005 by step 6,000
while held-out AUPR kept climbing for another 9,000 steps.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from snp2prot import tracking
from snp2prot.evaluation import metrics
from snp2prot.ght import config as ght_config
from snp2prot.ght.data import WindowCorpus
from snp2prot.ght.model import GenomicTwoTower, genomic_loss
from snp2prot.ght.splits import GHTFold

#: The negative set every run trains on and selects on (`GHT_PLAN.md` D4). The others are scored
#: at the end as controls and never fitted.
TRAIN_NEGATIVES = "shades"
#: Negative sets scored at the final evaluation, in report order.
EVALUATION_NEGATIVES = ("shades", "random", "aliens")


@dataclass
class GHTResult:
    """What one fold's run produced. Two models, as in the PBM arm, for the same reason."""

    steps_run: int
    best_step: int
    best_validation: float
    history: list[dict] = field(default_factory=list)
    parameter_counts: dict[str, int] = field(default_factory=dict)
    temperature: float = float("nan")
    seconds: float = 0.0
    seconds_per_step: float = float("nan")
    seconds_per_eval: float = float("nan")
    peak_memory_mb: float = float("nan")
    best_state: dict = field(default_factory=dict)
    final_state: dict = field(default_factory=dict)


def device_for(requested: str | None = None) -> torch.device:
    if requested:
        return torch.device(requested)
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


class GHTTrainer:
    """Holds the resident tensors and runs one fold."""

    def __init__(
        self,
        corpus: WindowCorpus,
        protein_vectors: np.ndarray,
        config: dict,
        device: torch.device | None = None,
        seed: int = 0,
        protein_mask: np.ndarray | None = None,
    ):
        self.corpus = corpus
        self.config = config
        self.device = device or device_for()
        torch.manual_seed(seed)

        if len(protein_vectors) != len(corpus.tfs):
            raise ValueError(
                f"{len(protein_vectors)} protein vectors for {len(corpus.tfs)} TFs — the two are "
                "indexed by the same positions for the whole run and must line up"
            )
        self.proteins = torch.from_numpy(np.ascontiguousarray(protein_vectors)).to(self.device)
        # Present only for the per-residue tower, where `proteins` is (n_tf, L, 1280) instead of
        # (n_tf, 1280) and a maximum over a pad position would invent a motif.
        self.protein_mask = (
            torch.from_numpy(np.ascontiguousarray(protein_mask)).to(self.device)
            if protein_mask is not None
            else None
        )
        self.labels = torch.from_numpy(corpus.label.astype(np.int64))
        self.tf_index = torch.from_numpy(corpus.tf_index.astype(np.int64))
        # Host copy of every window; the resident subset is mirrored on the device below.
        self.tokens_cpu = torch.from_numpy(corpus.tokens)
        resident = np.isin(corpus.negset, np.asarray(["", TRAIN_NEGATIVES]))
        self.slot = torch.full((corpus.n_windows,), -1, dtype=torch.int64)
        self.slot[torch.from_numpy(np.flatnonzero(resident))] = torch.arange(int(resident.sum()))
        mirrored = torch.from_numpy(np.flatnonzero(resident))
        self.tokens_gpu = self.tokens_cpu[mirrored].to(self.device)
        self.labels_gpu = self.labels.to(self.device)
        self.tf_index_gpu = self.tf_index.to(self.device)
        self.slot_gpu = self.slot.to(self.device)

        model_cfg = config["model"]
        self.model = GenomicTwoTower(
            protein_features=int(self.proteins.shape[-1]),
            width=int(model_cfg["width"]),
            dna_widths=tuple(int(w) for w in model_cfg["dna"]["widths"]),
            dna_channels=int(model_cfg["dna"]["channels_per_width"]),
            dna_dropout=float(model_cfg["dna"]["dropout"]),
            protein_hidden=int(model_cfg["protein"]["hidden"]),
            protein_dropout=float(model_cfg["protein"]["dropout"]),
            protein_tower=str(model_cfg["protein"].get("tower", "pooled")),
            residue_reduced=int(model_cfg["protein"].get("reduced", 64)),
            residue_channels=int(model_cfg["protein"].get("channels_per_width", 32)),
            temperature=float(model_cfg["temperature"]),
            learn_temperature=bool(model_cfg["learn_temperature"]),
            max_logit_scale=float(model_cfg["max_logit_scale"]),
        ).to(self.device)
        self.contrastive_weight = float(model_cfg.get("contrastive_weight", 0.0))
        self.load_protein_tower(
            model_cfg["protein"].get("init_from"),
            bool(model_cfg["protein"].get("freeze", False)),
        )

    def load_protein_tower(self, checkpoint: str | Path | None, freeze: bool = False) -> None:
        """Start the protein tower from a **PBM** checkpoint's tower, optionally frozen.

        Rung 5 of `GHT_PLAN.md` §11's fallback ladder — *protein tower shared with the PBM arm* —
        and it is a one-line load rather than an adapter because the two towers are literally the
        same module: `ProteinTower(1280, 256)`, four tensors, same names, same shapes. That they
        are is not a coincidence; `ML_PLAN.md` §4.2 fixed `D = 256` across arms so a comparison
        between them could never be a comparison of widths, and the same decision is what makes
        the weights portable.

        **Frozen is the interesting setting.** A tower that is allowed to move can rediscover
        whatever GHT needs and the PBM initialisation becomes a warm start; frozen, the protein
        representation is *fixed by a different assay on a different DNA vocabulary* and only the
        DNA tower may adapt. If that still works, one protein space serves both — which is the
        architectural claim the talk is making.

        **The caveat is measured, not waved away:** 8 of the 33 panel domains appear verbatim in
        the PBM corpus, so for those the PBM tower was fitted on that protein's own 8-mer
        behaviour. It is not a GHT label leak — the PBM arm never saw a genomic window — but it
        is not an unseen protein either, and the report splits the result on it.
        """
        if not checkpoint:
            return
        state = torch.load(Path(checkpoint), map_location="cpu", weights_only=False)
        weights = state.get("state_dict", state)
        tower = {k[len("protein.") :]: v for k, v in weights.items() if k.startswith("protein.")}
        if not tower:
            raise ValueError(f"{checkpoint} carries no `protein.` weights")
        self.model.protein.load_state_dict(tower, strict=True)
        self.model.protein.to(self.device)
        if freeze:
            for parameter in self.model.protein.parameters():
                parameter.requires_grad_(False)

    #: Parameters weight decay must not touch — the same list the PBM arm keeps, minus the null
    #: anchor it does not have. `calibration_bias` starts at the logit of the base rate and
    #: decaying it toward 0 is decaying it toward `p = 0.5`.
    NOT_WEIGHTS = ("logit_scale", "calibration_scale", "calibration_bias")

    def _parameter_groups(self, weight_decay: float) -> list[dict]:
        special = [p for n, p in self.model.named_parameters() if n in self.NOT_WEIGHTS]
        rest = [p for n, p in self.model.named_parameters() if n not in self.NOT_WEIGHTS]
        return [
            {"params": rest, "weight_decay": weight_decay},
            {"params": special, "weight_decay": 0.0},
        ]

    def tokens_for(self, rows: torch.Tensor) -> torch.Tensor:
        """`(n, 301)` int64 tokens on the device, from the GPU mirror where possible.

        A `random` or `aliens` row is not mirrored, so it is gathered on the host and copied. The
        branch is on the rows actually asked for rather than on a flag, because the final
        evaluation asks for both kinds and a caller should not have to know which.
        """
        slots = self.slot_gpu[rows]
        if bool((slots >= 0).all()):
            return self.tokens_gpu[slots].long()
        return self.tokens_cpu[rows.cpu()].to(self.device).long()

    def sample(self, pools: list[torch.Tensor], generator: torch.Generator) -> torch.Tensor:
        """One step's rows: `batch_proteins` TFs, `windows_per_protein` windows under each."""
        cfg = self.config["training"]
        per_protein = int(cfg["windows_per_protein"])
        n_proteins = min(int(cfg["batch_proteins"]), len(pools))
        order = torch.randperm(len(pools), generator=generator, device=self.device)[:n_proteins]
        picked = []
        for p in order.tolist():
            pool = pools[p]
            draw = torch.randint(len(pool), (per_protein,), generator=generator, device=self.device)
            picked.append(pool[draw])
        return torch.cat(picked)

    def _grouped(self, rows: torch.Tensor):
        """`(protein_vectors, mask, group)` for one batch — one entry per DISTINCT protein.

        The protein axis of a batch has at most 33 values and the window axis has thousands, so
        gathering the protein per window is wasted work for the pooled tower and impossible for
        the per-residue one (`GenomicTwoTower.cosine_grouped`).
        """
        of = self.tf_index_gpu[rows]
        unique, group = torch.unique(of, return_inverse=True)
        mask = None if self.protein_mask is None else self.protein_mask[unique]
        return self.proteins[unique], mask, group

    def loss_for(self, rows: torch.Tensor) -> tuple[torch.Tensor, dict[str, float]]:
        tokens = self.tokens_for(rows)
        proteins, mask, group = self._grouped(rows)
        cosine = self.model.cosine_grouped(proteins, group, tokens, mask)
        logits = self.model.calibrate(cosine)
        return genomic_loss(
            logits,
            self.labels_gpu[rows],
            cosine,
            self.tf_index_gpu[rows],
            self.model.logit_scale,
            self.contrastive_weight,
        )

    @torch.no_grad()
    def predict(self, rows: np.ndarray, chunk: int = 8192) -> np.ndarray:
        """`(n,)` calibrated logits for the given window rows, in the order given."""
        self.model.eval()
        index = torch.from_numpy(np.asarray(rows, dtype=np.int64)).to(self.device)
        out = np.empty(len(rows), dtype=np.float32)
        for start in range(0, len(rows), chunk):
            block = index[start : start + chunk]
            tokens = self.tokens_for(block)
            proteins, mask, group = self._grouped(block)
            cosine = self.model.cosine_grouped(proteins, group, tokens, mask)
            out[start : start + chunk] = self.model.calibrate(cosine).float().cpu().numpy()
        self.model.train()
        return out

    def score(self, rows: np.ndarray) -> tuple[dict, pd.DataFrame]:
        """Per-TF auPRC and auROC over `rows`, macro-averaged, each against its chance level."""
        predicted = self.predict(rows)
        tf_of = self.corpus.tf_index[rows]
        labels = self.corpus.label[rows]
        entries = []
        for t in np.unique(tf_of):
            mask = tf_of == t
            truth = (labels[mask] == 1).astype(np.int64)
            if truth.sum() == 0 or truth.sum() == len(truth):
                continue
            entries.append(
                {
                    "tf": str(self.corpus.tfs[t]),
                    "n": int(mask.sum()),
                    "n_pos": int(truth.sum()),
                    "positive_rate": float(truth.mean()),
                    "aupr": metrics.average_precision(truth, predicted[mask]),
                    "auroc": metrics.auroc(truth, predicted[mask]),
                }
            )
        frame = pd.DataFrame(entries)
        # The SAME keys either way. A negative set that is not on disk — the pre-flight loads a
        # resident-only corpus — must produce a row of nans rather than a shorter dict, or the
        # report reads a column that exists for one run and not for the next.
        empty = {
            k: float("nan")
            for k in (
                "aupr",
                "aupr_median",
                "aupr_min",
                "aupr_max",
                "auroc",
                "auroc_median",
                "chance_aupr",
                "aupr_over_chance",
            )
        }
        if frame.empty:
            return empty | {"n_tf": 0.0}, frame
        return (
            {
                "aupr": float(frame.aupr.mean()),
                "aupr_median": float(frame.aupr.median()),
                "aupr_min": float(frame.aupr.min()),
                "aupr_max": float(frame.aupr.max()),
                "auroc": float(frame.auroc.mean()),
                "auroc_median": float(frame.auroc.median()),
                "chance_aupr": float(frame.positive_rate.mean()),
                "aupr_over_chance": float((frame.aupr / frame.positive_rate).mean()),
                "n_tf": float(len(frame)),
            },
            frame,
        )

    def train(
        self,
        train_rows: np.ndarray,
        validation_rows: np.ndarray,
        seed: int = 0,
        on_eval=None,
    ) -> GHTResult:
        """Run one fold to its step budget. Selection by validation, termination by budget."""
        cfg = self.config["training"]
        steps = int(cfg["steps"])
        eval_every = int(cfg["eval_every"])
        patience = int(cfg["patience"])
        warmup = int(cfg["warmup_steps"])

        rate = float((self.corpus.label[train_rows] == 1).mean())
        self.model.set_calibration_prior(min(max(rate, 1e-6), 1 - 1e-6))

        optimiser = torch.optim.AdamW(
            self._parameter_groups(float(cfg["weight_decay"])),
            lr=float(cfg["learning_rate"]),
        )
        generator = torch.Generator(device=self.device).manual_seed(seed)
        pools = self._pools(train_rows)

        result = GHTResult(0, 0, -float("inf"))
        result.parameter_counts = self.model.parameter_counts()
        best_state: dict | None = None
        since_best = 0
        started = time.time()
        step_seconds = 0.0
        eval_seconds = 0.0
        if self.device.type == "cuda":
            torch.cuda.reset_peak_memory_stats(self.device)

        self.model.train()
        for step in range(1, steps + 1):
            for group in optimiser.param_groups:
                group["lr"] = float(cfg["learning_rate"]) * min(1.0, step / max(warmup, 1))
            tick = time.time()
            rows = self.sample(pools, generator)
            loss, parts = self.loss_for(rows)
            optimiser.zero_grad(set_to_none=True)
            loss.backward()
            optimiser.step()
            self.model.clamp_temperature()
            step_seconds += time.time() - tick
            result.steps_run = step

            if step % eval_every == 0 or step == steps:
                tick = time.time()
                summary, _ = self.score(validation_rows) if len(validation_rows) else ({}, None)
                eval_seconds += time.time() - tick
                score = float(summary.get("aupr", float("nan")))
                entry = {
                    "step": step,
                    "loss": float(loss.detach()),
                    "loss_bce": parts["bce"],
                    "loss_infonce": parts["infonce"],
                    "validation_aupr": score,
                    "validation_auroc": float(summary.get("auroc", float("nan"))),
                    "temperature": self.model.temperature,
                    "temperature_clamped": float(self.model.temperature_is_clamped),
                    "calibration_scale": float(self.model.calibration_scale.detach().exp()),
                    "calibration_bias": float(self.model.calibration_bias.detach()),
                }
                result.history.append(entry)
                if on_eval is not None:
                    on_eval(entry, self)

                improved = not math.isnan(score) and score > result.best_validation
                if not len(validation_rows):
                    improved = True
                if improved:
                    result.best_validation = score
                    result.best_step = step
                    best_state = {k: v.detach().clone() for k, v in self.model.state_dict().items()}
                    since_best = 0
                else:
                    since_best += 1
                    if patience and since_best >= patience:
                        break

        result.final_state = self._cpu_state()
        if best_state is not None:
            self.model.load_state_dict(best_state)
        else:
            result.best_step = result.steps_run
        result.best_state = self._cpu_state()
        result.temperature = self.model.temperature
        result.seconds = time.time() - started
        result.seconds_per_step = step_seconds / max(result.steps_run, 1)
        result.seconds_per_eval = eval_seconds / max(len(result.history), 1)
        if self.device.type == "cuda":
            result.peak_memory_mb = torch.cuda.max_memory_allocated(self.device) / 1e6
        return result

    def _pools(self, train_rows: np.ndarray) -> list[torch.Tensor]:
        """One tensor of row positions per training TF — the thing a step draws from."""
        tf_of = self.corpus.tf_index[train_rows]
        pools = []
        for t in np.unique(tf_of):
            rows = np.asarray(train_rows)[tf_of == t]
            pools.append(torch.from_numpy(rows.astype(np.int64)).to(self.device))
        return pools

    def _cpu_state(self) -> dict:
        return {k: v.detach().cpu().clone() for k, v in self.model.state_dict().items()}

    def load_state(self, state: dict) -> None:
        self.model.load_state_dict({k: v.to(self.device) for k, v in state.items()})

    def save(self, path: Path, result: GHTResult, meta: dict) -> Path:
        """Both trained models plus everything needed to rebuild either around them."""
        path.parent.mkdir(parents=True, exist_ok=True)
        torch.save(
            {
                "state_dict": result.best_state or self._cpu_state(),
                "final_state_dict": result.final_state or self._cpu_state(),
                "best_step": result.best_step,
                "final_step": result.steps_run,
                "model_config": self.config["model"],
                "protein_features": int(self.proteins.shape[-1]),
                "parameter_counts": self.model.parameter_counts(),
                **meta,
            },
            path,
        )
        return path


def fold_rows(corpus: WindowCorpus, fold: GHTFold, negset: str = TRAIN_NEGATIVES) -> dict:
    """The four row sets a fold defines, as positions into the corpus.

    `train` and `validation` are disjoint by construction — `C1` carves chromosomes and the TF
    regimes carve TFs, so neither can overlap — and both are checked here rather than trusted,
    because a silent overlap would make the selection signal a training metric.
    """
    train = corpus.rows(fold.fitting_tfs, fold.fitting_chromosomes, negset)
    if fold.validation_chromosomes:
        validation = corpus.rows(fold.fitting_tfs, fold.validation_chromosomes, negset)
    elif fold.validation_tfs:
        validation = corpus.rows(fold.validation_tfs, fold.fitting_chromosomes, negset)
    else:
        validation = np.array([], dtype=np.int64)
    test = corpus.rows(fold.test_tfs, fold.test_chromosomes, negset)
    if np.intersect1d(train, validation).size:
        raise ValueError(f"{fold.label}: training and validation rows overlap")
    if np.intersect1d(train, test).size:
        raise ValueError(f"{fold.label}: training and test rows overlap")
    return {"train": train, "validation": validation, "test": test}


def run_fold(
    fold: GHTFold,
    arm: str,
    corpus: WindowCorpus,
    protein_vectors: np.ndarray,
    protein_model: str,
    config: dict,
    device: torch.device | None = None,
    track: bool = True,
    model_seed: int | None = None,
    on_eval=None,
    protein_mask: np.ndarray | None = None,
    protein_mode: str = "real",
) -> tuple[dict, pd.DataFrame]:
    """Train and score one (arm, fold) pair, and return its summary and per-TF rows.

    The run records **both holdout axes** in its digest (`GHT_PLAN.md` §8.2) and the code stamp
    beside it, because a digest pins the data and not the procedure.
    """
    rows = fold_rows(corpus, fold)
    split_cfg = config["splits"]
    seed = int(split_cfg["seed"])
    model_seed = int(config["model"]["seed"]) if model_seed is None else int(model_seed)

    digests = {
        "test": fold.digest(),
        "train": f"{len(fold.fitting_tfs)}tf/{len(fold.fitting_chromosomes)}chr",
        "validation": (
            f"{len(fold.validation_tfs)}tf"
            if fold.validation_tfs
            else f"{len(fold.validation_chromosomes)}chr"
        ),
    }
    code = tracking.code_version()
    params = {
        "arm": arm,
        "code.commit": code["commit"],
        "code.dirty": code["dirty"],
        "protein_model": protein_model,
        "protein_mode": protein_mode,
        "regime": fold.regime,
        "fold": fold.name,
        "seed": seed,
        "model_seed": model_seed,
        "n_train": len(rows["train"]),
        "n_validation": len(rows["validation"]),
        "n_test": len(rows["test"]),
        "n_train_tf": len(fold.fitting_tfs),
        "n_test_tf": len(fold.test_tfs),
        **{f"model.{k}": v for k, v in config["model"].items()},
        **{f"training.{k}": v for k, v in config["training"].items()},
    }

    trainer = GHTTrainer(
        corpus,
        protein_vectors,
        config,
        device=device,
        seed=model_seed,
        protein_mask=protein_mask,
    )
    tower = str(config["model"]["protein"].get("tower", "pooled"))
    with tracking.start_run(
        f"ght/{arm}/{fold.label}/seed{model_seed}/{protein_mode}/{tower}",
        experiment_name="snp2prot-ght",
        params=params,
        digests=digests,
        enabled=track,
    ) as run:
        result = trainer.train(
            rows["train"],
            rows["validation"],
            seed=model_seed,
            on_eval=lambda entry, tr: (
                run.log_metrics(
                    {k: v for k, v in entry.items() if k != "step"}, step=entry["step"]
                ),
                on_eval(entry, tr, fold) if on_eval is not None else None,
            ),
        )

        # The model at the budget, scored beside the validation-selected one so the selection
        # loss is a reported number per fold rather than an assumption (`TrainingResult`).
        trainer.load_state(result.final_state)
        final_summary, _ = trainer.score(rows["test"])
        trainer.load_state(result.best_state)

        summary: dict = {}
        per_tf = []
        for negatives in EVALUATION_NEGATIVES:
            scored_rows = corpus.rows(fold.test_tfs, fold.test_chromosomes, negatives)
            block, frame = trainer.score(scored_rows)
            if frame is not None and not frame.empty:
                per_tf.append(frame.assign(negset=negatives))
            suffix = "" if negatives == TRAIN_NEGATIVES else f"_{negatives}"
            summary |= {f"{k}{suffix}": v for k, v in block.items()}
        summary["aupr_final"] = final_summary.get("aupr", float("nan"))
        summary["selection_gain"] = summary.get("aupr", float("nan")) - summary["aupr_final"]

        run.log_metrics({f"test.{k}": v for k, v in summary.items() if isinstance(v, float)})
        run.log_metrics({f"params.{k}": v for k, v in result.parameter_counts.items()})

        checkpoint = trainer.save(
            ght_config.checkpoint_file(
                arm,
                fold.regime,
                fold.name,
                model_seed,
                tag=ght_config.run_tag(protein_mode, tower),
            ),
            result,
            {
                "arm": arm,
                "protein_model": protein_model,
                "protein_mode": protein_mode,
                "regime": fold.regime,
                "fold": fold.name,
                "held_out": fold.held_out,
                "test_tfs": list(fold.test_tfs),
                "test_chromosomes": list(fold.test_chromosomes),
                "digests": digests,
                "code": code,
                "model_seed": model_seed,
                "best_step": result.best_step,
                "best_validation": result.best_validation,
            },
        )
        run.log_artifact(checkpoint)

    summary |= {
        "arm": arm,
        # The two knobs that decide what a row MEANS, on every row: a table assembled from a
        # model run and its protein-blind control must never be readable as one experiment.
        "protein_mode": protein_mode,
        "protein_tower": tower,
        "regime": fold.regime,
        "fold": fold.name,
        "seed": seed,
        "model_seed": model_seed,
        "n_train": float(len(rows["train"])),
        "n_validation": float(len(rows["validation"])),
        "n_test": float(len(rows["test"])),
        "n_train_tf": float(len(fold.fitting_tfs)),
        "n_test_tf": float(len(fold.test_tfs)),
        "steps_run": float(result.steps_run),
        "best_step": float(result.best_step),
        "validation_aupr": result.best_validation,
        "temperature": result.temperature,
        "temperature_clamped": float(trainer.model.temperature_is_clamped),
        "seconds": result.seconds,
        "seconds_per_step": result.seconds_per_step,
        "seconds_per_eval": result.seconds_per_eval,
        "peak_memory_mb": result.peak_memory_mb,
        "digest": digests["test"],
        "code_commit": code["commit"],
        "code_dirty": code["dirty"],
        "checkpoint": str(checkpoint),
        **{f"n_params_{k}": float(v) for k, v in result.parameter_counts.items()},
    }
    frame = pd.concat(per_tf, ignore_index=True) if per_tf else pd.DataFrame()
    if not frame.empty:
        frame = frame.assign(
            arm=arm,
            protein_mode=protein_mode,
            protein_tower=tower,
            regime=fold.regime,
            fold=fold.name,
            model_seed=model_seed,
        )
    return summary, frame
