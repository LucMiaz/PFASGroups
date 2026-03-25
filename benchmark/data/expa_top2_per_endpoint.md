# Experiment A -- Top-2 Feature Sets for Predictable Endpoints

> Endpoints with max ROC-AUC >= 0.72 (best feature set, both models avg).
> Ranked by mean ROC-AUC (both models averaged). All metrics shown as mean +/- SD.

## AR_agonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR+ring+mol | 0.7590 +/- 0.0380 | 0.3045 +/- 0.0567 | 0.2716 +/- 0.0503 | 0.6848 +/- 0.0359 |
| 2nd | PFG_EGR+spacer+ring+mol | 0.7584 +/- 0.0371 | 0.3041 +/- 0.0553 | 0.2705 +/- 0.0481 | 0.6840 +/- 0.0344 |
| 3 | PFG_min_dist_to_barycenter+mol | 0.7582 +/- 0.0364 | 0.2999 +/- 0.0593 | 0.2687 +/- 0.0476 | 0.6842 +/- 0.0350 |

## CYP2C9_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_total_component+mol | 0.7519 +/- 0.0357 | 0.7378 +/- 0.0399 | 0.3818 +/- 0.0692 | 0.6903 +/- 0.0346 |
| 2nd | PFG_binary+mol | 0.7480 +/- 0.0362 | 0.7331 +/- 0.0444 | 0.3714 +/- 0.0645 | 0.6851 +/- 0.0320 |
| 3 | PFG_EGR_4X+mol | 0.7452 +/- 0.0270 | 0.7215 +/- 0.0262 | 0.3870 +/- 0.0587 | 0.6930 +/- 0.0293 |

## CYP3A4_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR+spacer+mol | 0.7447 +/- 0.0245 | 0.6184 +/- 0.0358 | 0.3347 +/- 0.0642 | 0.6667 +/- 0.0326 |
| 2nd | PFG_total_component+mol | 0.7443 +/- 0.0242 | 0.6255 +/- 0.0358 | 0.3244 +/- 0.0618 | 0.6616 +/- 0.0307 |
| 3 | PFG_min_dist_to_barycenter+mol | 0.7441 +/- 0.0243 | 0.6174 +/- 0.0346 | 0.3327 +/- 0.0609 | 0.6658 +/- 0.0309 |

## Caspase3_HEPG2

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR+branch+mol | 0.7906 +/- 0.0491 | 0.3152 +/- 0.0801 | 0.3103 +/- 0.1051 | 0.6992 +/- 0.0695 |
| 2nd | PFG_EGR+mol | 0.7897 +/- 0.0505 | 0.3166 +/- 0.0824 | 0.3067 +/- 0.1049 | 0.6967 +/- 0.0691 |
| 3 | PFG_EGR+spacer+ring+mol | 0.7875 +/- 0.0514 | 0.3169 +/- 0.0834 | 0.3045 +/- 0.1094 | 0.6950 +/- 0.0707 |

## TR_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR_4X+mol | 0.7222 +/- 0.0278 | 0.6171 +/- 0.0457 | 0.3663 +/- 0.0634 | 0.6827 +/- 0.0302 |
| 2nd | PFG_binary+mol | 0.7220 +/- 0.0343 | 0.6363 +/- 0.0439 | 0.3614 +/- 0.0666 | 0.6789 +/- 0.0304 |
| 3 | PFG_total_component+mol | 0.7215 +/- 0.0343 | 0.6384 +/- 0.0475 | 0.3703 +/- 0.0735 | 0.6838 +/- 0.0342 |
