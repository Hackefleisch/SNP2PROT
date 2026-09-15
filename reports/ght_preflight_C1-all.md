# GHT pre-flight: the step and wall-clock budgets

`scripts/preflight_ght.py --fold C1/all --steps 20000`, arm `A1`, seed 20260915. `docs/GHT_PLAN.md` §8.1.

## (a) The step budget

Read off the **held-out test curve**, not the training loss. Test is logged at every evaluation as a measurement and is never used to select anything — the grid selects on validation like every other run in this project.

| step | train loss | validation auPRC | **test auPRC** | test auROC |
|---:|---:|---:|---:|---:|
| 250 | 0.5307 | 0.6872 | **0.6894** | 0.8131 |
| 1,000 | 0.2097 | 0.9092 | **0.9110** | 0.9545 |
| 1,750 | 0.2305 | 0.9216 | **0.9223** | 0.9601 |
| 2,500 | 0.1517 | 0.9249 | **0.9260** | 0.9627 |
| 3,250 | 0.1772 | 0.9291 | **0.9289** | 0.9638 |
| 4,000 | 0.1530 | 0.9287 | **0.9289** | 0.9633 |
| 4,750 | 0.1505 | 0.9275 | **0.9273** | 0.9627 |
| 5,500 | 0.1110 | 0.9290 | **0.9288** | 0.9628 |
| 6,250 | 0.1324 | 0.9297 | **0.9281** | 0.9625 |
| 7,000 | 0.1328 | 0.9287 | **0.9274** | 0.9620 |
| 7,750 | 0.1360 | 0.9284 | **0.9269** | 0.9616 |
| 8,500 | 0.1058 | 0.9285 | **0.9270** | 0.9617 |
| 9,250 | 0.1124 | 0.9280 | **0.9275** | 0.9618 |
| 10,000 | 0.1122 | 0.9279 | **0.9268** | 0.9618 |
| 10,750 | 0.1015 | 0.9278 | **0.9260** | 0.9611 |
| 11,500 | 0.1101 | 0.9265 | **0.9256** | 0.9611 |
| 12,250 | 0.0992 | 0.9275 | **0.9272** | 0.9616 |
| 13,000 | 0.0891 | 0.9269 | **0.9253** | 0.9609 |
| 13,750 | 0.0987 | 0.9263 | **0.9246** | 0.9603 |
| 14,500 | 0.0913 | 0.9261 | **0.9246** | 0.9603 |
| 15,250 | 0.0744 | 0.9266 | **0.9254** | 0.9605 |
| 16,000 | 0.1158 | 0.9259 | **0.9258** | 0.9607 |
| 16,750 | 0.0635 | 0.9254 | **0.9250** | 0.9602 |
| 17,500 | 0.0830 | 0.9252 | **0.9251** | 0.9602 |
| 18,250 | 0.1216 | 0.9272 | **0.9267** | 0.9609 |
| 19,000 | 0.0955 | 0.9247 | **0.9246** | 0.9600 |
| 19,750 | 0.0871 | 0.9254 | **0.9241** | 0.9599 |

Best test auPRC **0.9296 at step 4,250**; training loss there 0.1507.

Mean test auPRC over the last fifth of the run: **0.9252**, against **0.9255** over the fifth before it.

> The curve has flattened: the last fifth is not above the fifth before it by more than 0.002. A budget at **4,250** steps captures the peak, and the value adopted in `configs/experiment.yaml` is rounded to a round number at or above it.

## (b) The wall-clock budget

| quantity | value |
|---|---:|
| ms per step, steady state | 40.0 |
| seconds per evaluation pass | 1.69 |
| one fold end to end (20,000 steps, eval every 250) | 25.4 min |
| peak GPU memory | 4237 MB |
| training windows | 507,863 |
| validation windows | 125,323 |
| test windows (`shades`) | 550,377 |
| **one fold at the adopted 5,000-step budget** | **4.7 min** |

The fold cost above is the pre-flight's own, at four times the budget it set. The row that matters is the last one: at the adopted budget the full grid is **11 folds x 3 seeds = 2.6 h**, plus the final `random` and `aliens` passes, which this pre-flight does not pay because it loads a resident-only corpus.

## Parameter counts

| tower | parameters | per training protein |
|---|---:|---:|
| protein_tower | 330,496 | 10,015 |
| dna_tower | 80,384 | 2,436 |
| calibration | 2 | 0 |
| **total** | **410,883** | 12,451 |

`docs/TRAINING.md` §5 calls the PBM arm's **245 per protein** the central engineering constraint. At 33 training proteins even the bare linear tower is two orders of magnitude past that, so the panel size — not the tower architecture — is what sets the capacity risk here (`GHT_PLAN.md` §12).

Produced from `d256b62` (tree dirty).
