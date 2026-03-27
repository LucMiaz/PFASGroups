# Experiment A -- Top-3 Feature Sets per Endpoint (summary)

> Ranked by mean ROC-AUC (both models averaged). All metrics shown as mean +/- SD.

## AR_agonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR+spacer+mol | 0.7791 +/- 0.0393 | 0.3318 +/- 0.0647 | 0.2996 +/- 0.0663 | 0.7042 +/- 0.0480 |
| 2nd | PFG_EGR+ring+mol | 0.7789 +/- 0.0388 | 0.3326 +/- 0.0649 | 0.2993 +/- 0.0657 | 0.7041 +/- 0.0476 |
| 3rd | PFG_EGR+branch+mol | 0.7787 +/- 0.0395 | 0.3359 +/- 0.0640 | 0.2975 +/- 0.0663 | 0.7027 +/- 0.0479 |

## CYP2C9_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR_4X+mol | 0.7345 +/- 0.0249 | 0.7163 +/- 0.0273 | 0.3782 +/- 0.0544 | 0.6886 +/- 0.0273 |
| 2nd | PFG_binary+mol | 0.7342 +/- 0.0307 | 0.7246 +/- 0.0303 | 0.3569 +/- 0.0594 | 0.6777 +/- 0.0296 |
| 3rd | PFG_EGR+branch+mol | 0.7308 +/- 0.0317 | 0.7214 +/- 0.0310 | 0.3571 +/- 0.0608 | 0.6779 +/- 0.0305 |

## CYP3A4_antagonist

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_binary+mol | 0.7433 +/- 0.0260 | 0.6286 +/- 0.0414 | 0.3199 +/- 0.0552 | 0.6603 +/- 0.0285 |
| 2nd | PFG_EGR+spacer+mol | 0.7377 +/- 0.0293 | 0.6172 +/- 0.0378 | 0.3230 +/- 0.0607 | 0.6617 +/- 0.0309 |
| 3rd | PFG_EGR+branch+mol | 0.7373 +/- 0.0286 | 0.6191 +/- 0.0372 | 0.3217 +/- 0.0635 | 0.6613 +/- 0.0325 |

## Caspase3_HEPG2

| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |
|:---:|---|---:|---:|---:|---:|
| **1st** | PFG_EGR+branch+mol | 0.7962 +/- 0.0478 | 0.3289 +/- 0.0834 | 0.3213 +/- 0.1016 | 0.7059 +/- 0.0625 |
| 2nd | PFG_EGR+spacer+mol | 0.7953 +/- 0.0505 | 0.3277 +/- 0.0824 | 0.3213 +/- 0.1008 | 0.7054 +/- 0.0616 |
| 3rd | PFG_EGR+mol | 0.7940 +/- 0.0500 | 0.3271 +/- 0.0842 | 0.3223 +/- 0.1009 | 0.7062 +/- 0.0622 |
