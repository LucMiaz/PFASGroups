# Experiment A -- Top-3 Feature Sets per Endpoint (summary)

> Ranked by mean ROC-AUC (both models averaged). All metrics shown as mean +/- SD.

## AR_agonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR+spacer+mol | 0.7773 +/- 0.0421 | 0.3374 +/- 0.0631 | 0.2951 +/- 0.0697 | 0.7010 +/- 0.0510 |
| 2nd | PFG_EGR+branch+mol | 0.7772 +/- 0.0430 | 0.3363 +/- 0.0631 | 0.2995 +/- 0.0685 | 0.7038 +/- 0.0504 |
| 3rd | PFG_EGR+mol | 0.7772 +/- 0.0414 | 0.3371 +/- 0.0621 | 0.2957 +/- 0.0693 | 0.7012 +/- 0.0506 |

## CYP2C9_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR_4X+mol | 0.7345 +/- 0.0249 | 0.7163 +/- 0.0273 | 0.3782 +/- 0.0544 | 0.6886 +/- 0.0273 |
| 2nd | PFG_binary+mol | 0.7306 +/- 0.0304 | 0.7196 +/- 0.0305 | 0.3440 +/- 0.0599 | 0.6713 +/- 0.0298 |
| 3rd | PFG_EGR+spacer+mol | 0.7271 +/- 0.0297 | 0.7177 +/- 0.0314 | 0.3423 +/- 0.0612 | 0.6704 +/- 0.0307 |

## CYP3A4_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_binary+mol | 0.7423 +/- 0.0255 | 0.6259 +/- 0.0407 | 0.3244 +/- 0.0600 | 0.6625 +/- 0.0308 |
| 2nd | PFG_EGR+ring+mol | 0.7381 +/- 0.0286 | 0.6190 +/- 0.0375 | 0.3194 +/- 0.0647 | 0.6599 +/- 0.0333 |
| 3rd | PFG_EGR+spacer+mol | 0.7372 +/- 0.0297 | 0.6195 +/- 0.0374 | 0.3234 +/- 0.0607 | 0.6619 +/- 0.0315 |

## Caspase3_HEPG2

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR+branch+mol | 0.7965 +/- 0.0492 | 0.3324 +/- 0.0856 | 0.3171 +/- 0.0960 | 0.7046 +/- 0.0605 |
| 2nd | PFG_EGR+mol | 0.7964 +/- 0.0502 | 0.3277 +/- 0.0807 | 0.3172 +/- 0.0945 | 0.7040 +/- 0.0600 |
| 3rd | PFG_EGR+ring+mol | 0.7947 +/- 0.0499 | 0.3286 +/- 0.0843 | 0.3226 +/- 0.1044 | 0.7079 +/- 0.0634 |
