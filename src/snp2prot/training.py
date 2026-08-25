"""The training loop: one arm, one fold, end to end.

`docs/TRAINING.md` is the design document. What this module implements from it:

**Everything is resident** (`ML_PLAN.md` §3.1). The label matrix is 42 MiB of `int8`, the
protein embeddings 6.5 MiB, the 8-mer tokens under a megabyte — all of it goes to the device
once and stays there. There is no dataloader, no shuffling of rows and no I/O in the loop.

**A step draws `B` domains and scores each against the complete 8-mer axis.** Shuffling is on
the protein axis only; the 8-mer axis is never subsampled, so the InfoNCE denominator is exact
rather than estimated. That also makes the training objective structurally identical to the
evaluation metric — §5.3 ranks one protein's 32,896 8-mers and macro-averages over proteins, and
the loss is a smooth surrogate for exactly that.

**The DNA table is computed once per step** and shared across the batch, which is what keeps a
step under a GFLOP. Its activations, not the logits, are what fill the GPU: the whole vocabulary
passes through the convolution twice per step, once per strand.

**Validation is carved from the training domains**, never from the test set, and is used only to
decide when to stop. `docs/TRAINING.md` §7a records the one case it cannot model faithfully.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import torch

from snp2prot import splits, tracking
from snp2prot.data.matrix import KmerMatrix
from snp2prot.embeddings import DomainEmbeddings
from snp2prot.evaluation import metrics
from snp2prot.models import TwoTower, multi_positive_infonce, tokenise


@dataclass
class TrainingResult:
    """What one fold's run produced."""

    steps_run: int
    best_step: int
    best_validation: float
    history: list[dict] = field(default_factory=list)
    parameter_counts: dict[str, int] = field(default_factory=dict)
    temperature: float = float("nan")
    seconds: float = 0.0


def device_for(requested: str | None = None) -> torch.device:
    if requested:
        return torch.device(requested)
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


class Trainer:
    """Holds the resident tensors and runs one fold.

    Constructed once per run rather than once per step, because the point of the design is that
    nothing moves between host and device inside the loop.
    """

    def __init__(
        self,
        matrix: KmerMatrix,
        embeddings: DomainEmbeddings,
        config: dict,
        device: torch.device | None = None,
        seed: int = 0,
    ):
        self.config = config
        self.device = device or device_for()
        torch.manual_seed(seed)

        self.labels = torch.from_numpy(matrix.label).to(self.device)
        self.proteins = torch.from_numpy(embeddings.vectors).to(self.device)
        self.tokens = torch.from_numpy(tokenise(matrix.kmers)).to(self.device)
        self.escore = matrix.escore

        model_cfg = config["model"]
        self.model = TwoTower(
            protein_features=embeddings.width,
            width=int(model_cfg["width"]),
            dna_channels=int(model_cfg["dna"]["channels"]),
            dna_layers=int(model_cfg["dna"]["layers"]),
            protein_hidden=int(model_cfg["protein"]["hidden"]),
            protein_dropout=float(model_cfg["protein"]["dropout"]),
            temperature=float(model_cfg["temperature"]),
            learn_temperature=bool(model_cfg["learn_temperature"]),
            max_logit_scale=float(model_cfg["max_logit_scale"]),
        ).to(self.device)

    def loss_for(self, rows: torch.Tensor, dna_table: torch.Tensor) -> tuple:
        logits = self.model.score(self.proteins[rows], dna_table)
        return multi_positive_infonce(logits, self.labels[rows])

    @torch.no_grad()
    def predict(self, rows: np.ndarray, chunk: int = 256) -> np.ndarray:
        """`(n, K)` predicted scores over the real 8-mers, the same shape the baseline returns.

        The null anchor's column is dropped: it is a training device and a threshold, not an
        8-mer, so it has no place in a ranking of the vocabulary.
        """
        self.model.eval()
        table = self.model.dna_table(self.tokens)
        out = np.empty((len(rows), self.labels.shape[1]), dtype=np.float32)
        index = torch.from_numpy(np.asarray(rows, dtype=np.int64)).to(self.device)
        for start in range(0, len(rows), chunk):
            block = index[start : start + chunk]
            logits = self.model.score(self.proteins[block], table)[:, :-1]
            out[start : start + chunk] = logits.float().cpu().numpy()
        self.model.train()
        return out

    @torch.no_grad()
    def macro_aupr(self, rows: np.ndarray) -> float:
        """Mean per-protein AUPR over `rows`, skipping domains with no positive 8-mer.

        The selection criterion is the reported metric rather than the validation loss, so that
        early stopping optimises the thing the run is judged on.
        """
        if not len(rows):
            return float("nan")
        predicted = self.predict(rows)
        labels = self.labels[torch.from_numpy(rows).to(self.device)].cpu().numpy()
        scores = []
        for i in range(len(rows)):
            scored = labels[i] != -1
            binary = (labels[i][scored] == 1).astype(np.int64)
            if binary.sum():
                scores.append(metrics.average_precision(binary, predicted[i][scored]))
        return float(np.mean(scores)) if scores else float("nan")

    def train(
        self,
        train_rows: np.ndarray,
        validation_rows: np.ndarray,
        seed: int = 0,
        on_eval=None,
    ) -> TrainingResult:
        """Run one fold to its step budget or its patience, whichever comes first."""
        cfg = self.config["training"]
        steps = int(cfg["steps"])
        batch = min(int(cfg["batch_domains"]), len(train_rows))
        eval_every = int(cfg["eval_every"])
        patience = int(cfg["patience"])
        warmup = int(cfg["warmup_steps"])

        optimiser = torch.optim.AdamW(
            self.model.parameters(),
            lr=float(cfg["learning_rate"]),
            weight_decay=float(cfg["weight_decay"]),
        )
        pool = torch.from_numpy(np.asarray(train_rows, dtype=np.int64)).to(self.device)
        generator = torch.Generator(device=self.device).manual_seed(seed)

        result = TrainingResult(0, 0, -float("inf"))
        result.parameter_counts = self.model.parameter_counts()
        best_state = None
        since_best = 0
        started = time.time()

        self.model.train()
        for step in range(1, steps + 1):
            for group in optimiser.param_groups:
                group["lr"] = float(cfg["learning_rate"]) * min(1.0, step / max(warmup, 1))

            order = torch.randperm(len(pool), generator=generator, device=self.device)
            picked = pool[order[:batch]]
            # Recomputed every step: the DNA tower's weights move, so the table does too.
            table = self.model.dna_table(self.tokens)
            loss, _ = self.loss_for(picked, table)

            optimiser.zero_grad(set_to_none=True)
            loss.backward()
            optimiser.step()
            result.steps_run = step

            if step % eval_every == 0 or step == steps:
                score = self.macro_aupr(validation_rows)
                entry = {
                    "step": step,
                    "loss": float(loss.detach()),
                    "validation_aupr": score,
                    "temperature": self.model.temperature,
                    "temperature_clamped": float(self.model.temperature_is_clamped),
                }
                result.history.append(entry)
                if on_eval is not None:
                    on_eval(entry)

                improved = math.isnan(result.best_validation) or (
                    not math.isnan(score) and score > result.best_validation
                )
                if not len(validation_rows):
                    improved = True  # no validation: the last step is the chosen one
                if improved:
                    result.best_validation = score
                    result.best_step = step
                    best_state = {k: v.detach().clone() for k, v in self.model.state_dict().items()}
                    since_best = 0
                else:
                    since_best += 1
                    if patience and since_best >= patience:
                        break

        if best_state is not None:
            self.model.load_state_dict(best_state)
        result.temperature = self.model.temperature
        result.seconds = time.time() - started
        return result


def run_fold(
    fold,
    arm: str,
    domains,
    matrix: KmerMatrix,
    embeddings: DomainEmbeddings,
    distances,
    config: dict,
    trainable: np.ndarray,
    excluded: np.ndarray,
    device: torch.device | None = None,
    track: bool = True,
) -> tuple[dict, pd.DataFrame]:
    """Train and score one (arm, fold) pair, and return its summary and per-domain rows.

    Two exclusions apply to the **training pool** and to nothing else:

    - `trainable` drops the `no_evidence` records, because a silent record with no control
      cannot be told from a failed assay (`T21`);
    - `excluded` drops the C1 evaluation set, which is never trained on by construction (`D6`).

    Both are masks over the corpus, applied to the fold's training side only — a fold's *test*
    set is whatever the regime says it is, and is not quietly narrowed here.
    """
    split_cfg = config["splits"]
    validation_cfg = split_cfg["validation"]
    seed = int(split_cfg["seed"])

    pool = fold.train[trainable[fold.train] & ~excluded[fold.train]]
    inner = splits.Fold(fold.regime, fold.name, fold.test, pool, fold.held_out)
    train_rows, validation_rows = splits.validation_split(
        inner,
        domains,
        distances,
        float(validation_cfg["fraction"]),
        seed,
        str(validation_cfg["grouping"]),
    )

    digests = {
        "test": fold.digest(domains),
        "train": splits.digest_of(domains, train_rows),
        "validation": (
            splits.digest_of(domains, validation_rows) if len(validation_rows) else "none"
        ),
    }
    params = {
        "arm": arm,
        "model_name": embeddings.model,
        "regime": fold.regime,
        "fold": fold.name,
        "seed": seed,
        "n_train": len(train_rows),
        "n_validation": len(validation_rows),
        "n_test": len(fold.test),
        **{f"model.{k}": v for k, v in config["model"].items()},
        **{f"training.{k}": v for k, v in config["training"].items()},
        "splits.s2_min_identity": split_cfg["s2_min_identity"],
        "splits.validation": validation_cfg,
    }

    trainer = Trainer(matrix, embeddings, config, device=device, seed=seed)
    with tracking.start_run(
        f"{arm}/{fold.label}", params=params, digests=digests, enabled=track
    ) as run:
        result = trainer.train(
            train_rows,
            validation_rows,
            seed=seed,
            on_eval=lambda entry: run.log_metrics(
                {k: v for k, v in entry.items() if k != "step"}, step=entry["step"]
            ),
        )

        predicted = trainer.predict(fold.test)
        precision_at = tuple(int(k) for k in config["metrics"]["precision_at"])
        target = float(config["metrics"]["recall_at_precision"])
        scored = [
            metrics.score_domain(
                matrix.label[row],
                matrix.escore[row],
                predicted[i],
                domain=str(matrix.domains[row]),
                precision_at=precision_at,
                precision_target=target,
            )
            for i, row in enumerate(fold.test)
        ]
        summary = metrics.macro_average(scored)
        run.log_metrics({f"test.{k}": v for k, v in summary.items()})
        run.log_metrics({f"params.{k}": v for k, v in result.parameter_counts.items()})

    summary |= {
        "arm": arm,
        "regime": fold.regime,
        "fold": fold.name,
        "n_train": float(len(train_rows)),
        "n_validation": float(len(validation_rows)),
        "steps_run": float(result.steps_run),
        "best_step": float(result.best_step),
        "validation_aupr": result.best_validation,
        "temperature": result.temperature,
        "temperature_clamped": float(trainer.model.temperature_is_clamped),
        "seconds": result.seconds,
        "digest": digests["test"],
        **{f"n_params_{k}": float(v) for k, v in result.parameter_counts.items()},
    }

    per_domain = pd.DataFrame(
        {
            "arm": arm,
            "regime": fold.regime,
            "fold": fold.name,
            "domain": [s.domain for s in scored],
            "family": domains.dbd_family.to_numpy()[fold.test],
            "gene": domains.gene.to_numpy()[fold.test],
            "is_variant": domains.is_variant.to_numpy()[fold.test],
            "verdict": domains.verdict.to_numpy()[fold.test],
            "n_pos": [s.n_pos for s in scored],
            "aupr": [s.aupr for s in scored],
            "auroc": [s.auroc for s in scored],
            "spearman": [s.spearman for s in scored],
        }
    )
    return summary, per_domain
