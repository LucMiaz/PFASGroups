
# PFASGroups Classification Analysis Report

## Overview
This report analyzes the performance of the PFASGroups classification system based on specificity test results.

## Summary Statistics
- Total molecules analyzed: 1832
- Molecules with detections: 1832 (100.0%)
- Molecules with expected groups: 1832
- Specificity rate: 46.3%
- Detection accuracy: 97.3%

## Performance Metrics
- **Specific correct**: 849 molecules correctly classified with high specificity
- **Expected detected**: 1782 molecules where expected groups were detected
- **Coverage**: 100.0% of molecules had at least one functional group detected

## Top Detected Functional Groups
- **Perfluoroalkane-per** (ID: 48): 1689 detections (92.2%)
- **Group_31** (ID: 31): 308 detections (16.8%)
- **Group_22** (ID: 22): 249 detections (13.6%)
- **Group_51** (ID: 51): 249 detections (13.6%)
- **alkane-per** (ID: 23): 190 detections (10.4%)
- **Group_7** (ID: 7): 158 detections (8.6%)
- **Group_36** (ID: 36): 158 detections (8.6%)
- **Group_33** (ID: 33): 156 detections (8.5%)
- **Group_47** (ID: 47): 155 detections (8.5%)
- **halide** (ID: 21): 154 detections (8.4%)

## Top Expected Functional Groups
- **Perfluoroalkane-per** (ID: 48): 775 expected (42.3%)
- **Group_31** (ID: 31): 278 expected (15.2%)
- **Group_36** (ID: 36): 158 expected (8.6%)
- **Group_47** (ID: 47): 155 expected (8.5%)
- **Group_51** (ID: 51): 125 expected (6.8%)
- **Group_29** (ID: 29): 124 expected (6.8%)
- **Group_7** (ID: 7): 124 expected (6.8%)
- **Group_2** (ID: 2): 124 expected (6.8%)
- **Group_33** (ID: 33): 94 expected (5.1%)
- **Perfluoroalkene-per** (ID: 49): 93 expected (5.1%)

## Detection vs Expected Comparison
- **Group_2** (ID: 2): 0/124 detected (0.0% accuracy)
- **Group_7** (ID: 7): 158/124 detected (127.4% accuracy)
- **Group_29** (ID: 29): 0/124 detected (0.0% accuracy)
- **Group_31** (ID: 31): 308/278 detected (110.8% accuracy)
- **Group_33** (ID: 33): 156/94 detected (166.0% accuracy)
- **Group_36** (ID: 36): 158/158 detected (100.0% accuracy)
- **Group_47** (ID: 47): 155/155 detected (100.0% accuracy)
- **Perfluoroalkane-per** (ID: 48): 1689/775 detected (217.9% accuracy)
- **Perfluoroalkene-per** (ID: 49): 0/93 detected (0.0% accuracy)
- **Group_51** (ID: 51): 249/125 detected (199.2% accuracy)

## Molecules with Specificity Issues
- 983 molecules (53.7%) had specificity issues

### Issues by Functional Group Type:
- iodide-per: 31 molecules
- iodide-poly: 31 molecules
- azole-poly: 31 molecules
- azole-per: 31 molecules
- amine-per: 31 molecules
- sulfonamide-per: 31 molecules
- azine-per: 31 molecules
- azine-poly: 31 molecules
- amide-per: 30 molecules
- ether-per: 30 molecules

## Data Quality
- Valid SMILES: 1832 / 1832 (100.0%)
- Molecules with errors: 0 (0.0%)

## Notes
- Specificity issues indicate that the classification detected more groups than expected
- High specificity is desired to avoid false positive functional group assignments
- Detection accuracy shows how well the system finds expected functional groups
