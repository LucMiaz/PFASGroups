# Experiment A -- Top-3 Feature Sets for Predictable Endpoints

> Endpoints with max ROC-AUC >= 0.65 (best feature set, both models avg).
> Ranked by mean ROC-AUC (both models averaged). All metrics shown as mean +/- SD.

## AR_agonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | TxP_PFAS+g51g52_total | 0.7610 +/- 0.0530 | 0.5067 +/- 0.0987 | 0.4935 +/- 0.0616 | 0.7545 +/- 0.0383 |
| 2nd | PFG_EGR+ring+mol | 0.7590 +/- 0.0380 | 0.3045 +/- 0.0567 | 0.2716 +/- 0.0503 | 0.6848 +/- 0.0359 |
| 3rd | PFG_EGR+spacer+ring+mol | 0.7584 +/- 0.0371 | 0.3041 +/- 0.0553 | 0.2705 +/- 0.0481 | 0.6840 +/- 0.0344 |

## AR_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR_4X+mol | 0.6747 +/- 0.0350 | 0.4968 +/- 0.0479 | 0.2493 +/- 0.0614 | 0.6265 +/- 0.0298 |
| 2nd | PFG_binary+mol | 0.6741 +/- 0.0394 | 0.4881 +/- 0.0539 | 0.2537 +/- 0.0657 | 0.6274 +/- 0.0330 |
| 3rd | PFG_total_component+mol | 0.6723 +/- 0.0411 | 0.4840 +/- 0.0541 | 0.2526 +/- 0.0693 | 0.6266 +/- 0.0345 |

## Aromatase_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_min_dist_to_barycenter+mol | 0.6620 +/- 0.0276 | 0.3337 +/- 0.0389 | 0.1912 +/- 0.0614 | 0.6041 +/- 0.0305 |
| 2nd | PFG_total_component+mol | 0.6616 +/- 0.0336 | 0.3408 +/- 0.0428 | 0.1952 +/- 0.0632 | 0.6049 +/- 0.0332 |
| 3rd | PFG_EGR+spacer+ring+mol | 0.6579 +/- 0.0323 | 0.3254 +/- 0.0393 | 0.1843 +/- 0.0606 | 0.6012 +/- 0.0307 |

## CYP2C19_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR_4X+mol | 0.6922 +/- 0.0331 | 0.5399 +/- 0.0408 | 0.3101 +/- 0.0633 | 0.6577 +/- 0.0325 |
| 2nd | PFG_binary+mol | 0.6906 +/- 0.0282 | 0.5234 +/- 0.0361 | 0.2858 +/- 0.0599 | 0.6459 +/- 0.0309 |
| 3rd | PFG_total_component+mol | 0.6900 +/- 0.0281 | 0.5235 +/- 0.0371 | 0.2813 +/- 0.0602 | 0.6433 +/- 0.0311 |

## CYP2C9_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_total_component+mol | 0.7519 +/- 0.0357 | 0.7378 +/- 0.0399 | 0.3818 +/- 0.0692 | 0.6903 +/- 0.0346 |
| 2nd | PFG_binary+mol | 0.7480 +/- 0.0362 | 0.7331 +/- 0.0444 | 0.3714 +/- 0.0645 | 0.6851 +/- 0.0320 |
| 3rd | PFG_EGR_4X+mol | 0.7452 +/- 0.0270 | 0.7215 +/- 0.0262 | 0.3870 +/- 0.0587 | 0.6930 +/- 0.0293 |

## CYP2D6_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR_4X+mol | 0.7193 +/- 0.0214 | 0.7243 +/- 0.0269 | 0.3216 +/- 0.0454 | 0.6583 +/- 0.0227 |
| 2nd | PFG_binary+mol | 0.7088 +/- 0.0262 | 0.7132 +/- 0.0287 | 0.2745 +/- 0.0455 | 0.6351 +/- 0.0226 |
| 3rd | PFG_total_component+mol | 0.7074 +/- 0.0263 | 0.7137 +/- 0.0292 | 0.2707 +/- 0.0493 | 0.6328 +/- 0.0242 |

## CYP3A4_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR+spacer+mol | 0.7447 +/- 0.0245 | 0.6184 +/- 0.0358 | 0.3347 +/- 0.0642 | 0.6667 +/- 0.0326 |
| 2nd | PFG_total_component+mol | 0.7443 +/- 0.0242 | 0.6255 +/- 0.0358 | 0.3244 +/- 0.0618 | 0.6616 +/- 0.0307 |
| 3rd | PFG_min_dist_to_barycenter+mol | 0.7441 +/- 0.0243 | 0.6174 +/- 0.0346 | 0.3327 +/- 0.0609 | 0.6658 +/- 0.0309 |

## Caspase3_HEPG2

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR+branch+mol | 0.7906 +/- 0.0491 | 0.3152 +/- 0.0801 | 0.3103 +/- 0.1051 | 0.6992 +/- 0.0695 |
| 2nd | PFG_EGR+mol | 0.7897 +/- 0.0505 | 0.3166 +/- 0.0824 | 0.3067 +/- 0.1049 | 0.6967 +/- 0.0691 |
| 3rd | PFG_EGR+spacer+ring+mol | 0.7875 +/- 0.0514 | 0.3169 +/- 0.0834 | 0.3045 +/- 0.1094 | 0.6950 +/- 0.0707 |

## DT40_genotoxicity

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_min_dist_to_barycenter+mol | 0.7195 +/- 0.0274 | 0.7268 +/- 0.0332 | 0.3146 +/- 0.0527 | 0.6560 +/- 0.0262 |
| 2nd | PFG_total_component+mol | 0.7188 +/- 0.0263 | 0.7260 +/- 0.0327 | 0.3155 +/- 0.0508 | 0.6567 +/- 0.0253 |
| 3rd | PFG_binary+mol | 0.7184 +/- 0.0253 | 0.7247 +/- 0.0309 | 0.3068 +/- 0.0508 | 0.6524 +/- 0.0255 |

## MMP_ratio

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_total_component | 0.6827 +/- 0.0350 | 0.4762 +/- 0.0372 | 0.2862 +/- 0.0488 | 0.6526 +/- 0.0291 |
| 2nd | TxP_PFAS+g51g52_total | 0.6803 +/- 0.0379 | 0.4844 +/- 0.0463 | 0.3078 +/- 0.0440 | 0.6580 +/- 0.0262 |
| 3rd | TxP_PFAS | 0.6790 +/- 0.0341 | 0.4677 +/- 0.0389 | 0.3022 +/- 0.0562 | 0.6560 +/- 0.0305 |

## TR_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR_4X+mol | 0.7222 +/- 0.0278 | 0.6171 +/- 0.0457 | 0.3663 +/- 0.0634 | 0.6827 +/- 0.0302 |
| 2nd | PFG_binary+mol | 0.7220 +/- 0.0343 | 0.6363 +/- 0.0439 | 0.3614 +/- 0.0666 | 0.6789 +/- 0.0304 |
| 3rd | PFG_total_component+mol | 0.7215 +/- 0.0343 | 0.6384 +/- 0.0475 | 0.3703 +/- 0.0735 | 0.6838 +/- 0.0342 |
