# Mismatch Analysis Report

## Executive Summary

**Total Mismatches: 112 instances (52 unique molecules)**

**Type: ALL FALSE POSITIVES** (HTML detects extra groups that Python doesn't)

**NO FALSE NEGATIVES** (HTML never misses groups that Python detects)

## Root Cause

The **HTML version is NOT filtering groups by pathway** (Perfluoroalkyl/Polyfluoroalkyl), while the **Python test script IS filtering** them.

### In test_html_vs_python.py (line 223):
```python
group_ids = [int(g[0].id) for g in detected_groups if len(g[3])==0 or any([k['SMARTS'] in ['Perfluoroalkyl','Polyfluoroalkyl'] for k in g[3]])]
```

This filters to include only:
1. Groups with no chains (`len(g[3])==0`) OR
2. Groups where at least one chain has SMARTS 'Perfluoroalkyl' or 'Polyfluoroalkyl'

### The HTML version does NOT have this filtering logic!

## Top 5 False Positive Groups

| Group ID | Group Name | Count | Requires Pathway? |
|----------|-----------|-------|-------------------|
| 19 | Hydrofluoroolefins (HFOs) | 56 | NO (smartsPath: null) |
| 5 | Polyfluoroalkylether carboxylic acid | 40 | YES (smartsPath: "poly") |
| 4 | Perfluoroalkylether carboxylic acids | 16 | YES (smartsPath: "per") |
| 21 | Semi-fluorinated alkanes (SFAs) | 8 | YES (smartsPath: "poly") |
| 3 | Perfluoroalkyl dicarboxylic acids | 8 | YES (smartsPath: "per") |

## Detailed Analysis

### Issue 1: Group 19 (Hydrofluoroolefins) - 56 FALSE POSITIVES

**Definition:**
- `smartsPath`: `null` (no pathway requirement)
- `smarts1`: `[#6X3$([#6]=[#6])]` (sp2 carbon in C=C)
- `smarts2`: `[#6X4$([#6]F)]` (sp3 carbon with F)
- `constraints`: `{"gte":{"F":1}, "lte":{}, "eq":{}, "only":["C","F","H"]}`

**Example Mismatch:**
- SMILES: `FC(=C(F)C(C(F)(F)F)(C(F)(F)C(F)(F)F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)F)C(F)(F)F`
- Formula: `C12H2F24`
- Python: detects groups 18, 41, 48, 49 (NO group 19)
- HTML: detects groups 18, **19**, 41, 48, 49 (EXTRA group 19)

**Why the mismatch?**
Group 19 has `smartsPath: null`, meaning it doesn't require a specific pathway type. HOWEVER, the Python test script is filtering it OUT because the chains don't have 'Perfluoroalkyl' or 'Polyfluoroalkyl' SMARTS.

**The HTML should match Python behavior!**

### Issue 2: Groups 4 & 5 (Ether carboxylic acids) - 56 FALSE POSITIVES COMBINED

**Group 4 (Perfluoroalkylether carboxylic acids):**
- `smartsPath`: `"per"` (REQUIRES Perfluoroalkyl pathway)
- Constraint: `{"rel":{"C":{"atoms":["F"],"add":0.5,"div":2}}}`
- Requires Perfluoroalkyl chain between smarts1 and smarts2

**Group 5 (Polyfluoroalkylether carboxylic acid):**
- `smartsPath`: `"poly"` (REQUIRES Polyfluoroalkyl pathway)
- Constraint: `{"rel":{"C":{"atoms":["F","H","Cl","Br","I"],"add":0,"div":2}}}`
- Requires Polyfluoroalkyl chain between smarts1 and smarts2

**Example Mismatch:**
- SMILES: `O=C(O)C(F)(F)C(OC(F)(F)F)(C(F)(F)C(F)(F)F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)F`
- Formula: `C11H3F21O3`
- Python: detects groups 1, 2, 31, 33, 48 (NO groups 4 or 5)
- HTML: detects groups 1, 2, **4**, **5**, 31, 33, 48 (EXTRA groups 4 and 5)

**Why the mismatch?**
These groups require specific pathway types (`"per"` or `"poly"`), but the HTML is detecting them even when the pathway between the functional groups doesn't meet the requirements OR when the Python filtering logic excludes them.

### Issue 3: Group 3 (Perfluoroalkyl dicarboxylic acids) - 8 FALSE POSITIVES

**Definition:**
- `smartsPath`: `"per"` (REQUIRES Perfluoroalkyl pathway)
- `smarts1`: `[#6$([#6](=O)([OH1]))]` (carboxylic acid)
- `smarts2`: `[#6$([#6](=O)([OH1]))]` (carboxylic acid)
- Constraint: `{"rel":{"C":{"atoms":["F"],"add":2,"div":2}}}`

**Example Mismatch:**
- SMILES: `O=C(O)C(F)(F)C(C(F)(F)C(=O)O)(C(F)(F)C(F)(F)F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)F`
- Formula: `C12H2F20O4`
- Python: detects groups 33, 48 (NO group 3)
- HTML: detects groups **3**, 33, 48 (EXTRA group 3)

**Why the mismatch?**
The pathway between the two carboxylic acids likely doesn't meet the Perfluoroalkyl requirements, OR the Python filtering excludes it.

### Issue 4: Group 21 (Semi-fluorinated alkanes) - 8 FALSE POSITIVES

**Definition:**
- `smartsPath`: `"poly"` (REQUIRES Polyfluoroalkyl pathway)
- `smarts1`: `[#6X4$([#6](F)(F))]` (C with 2 F atoms)
- `smarts2`: `[#6X4$([#6][H,I,Br,Cl])]` (C with H/I/Br/Cl)
- Constraint: `{"gte":{"F":1}, "rel":{"C":{"atoms":["Cl","I","Br","F","H"],"add":-1,"div":2}}}`

**Example Mismatch:**
- SMILES: `OCCCC(F)(F)C(F)(C(F)(F)F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F`
- Formula: `C10H7F15O`
- Python: detects groups 15, 48 (NO group 21)
- HTML: detects groups 15, **21**, 48 (EXTRA group 21)

**Why the mismatch?**
Group 21 requires a Polyfluoroalkyl pathway, but HTML detects it without proper filtering.

## Solution

The HTML version needs to implement the SAME filtering logic as the Python test script:

**In the HTML detection logic, after detecting groups, filter them to include only:**
1. Groups where no chains are required (`smartsPath === null` AND no smarts2), OR
2. Groups where at least one detected chain has SMARTS type 'Perfluoroalkyl' or 'Polyfluoroalkyl'

### Pseudo-code for the fix:

```javascript
// After detecting all groups
const filteredGroups = detectedGroups.filter(group => {
    // If group has no chains requirement, include it
    if (group.chains.length === 0) {
        return true;
    }
    
    // If group has chains, check if any chain is Perfluoroalkyl or Polyfluoroalkyl
    return group.chains.some(chain => 
        chain.SMARTS === 'Perfluoroalkyl' || chain.SMARTS === 'Polyfluoroalkyl'
    );
});
```

## Expected Impact

Fixing this pathway filtering issue should resolve **ALL 112 mismatches** because:
- Group 19: 56 instances (will be filtered out when chains aren't Per/Poly)
- Group 5: 40 instances (requires "poly" pathway)
- Group 4: 16 instances (requires "per" pathway)
- Group 21: 8 instances (requires "poly" pathway)
- Group 3: 8 instances (requires "per" pathway)

**Total: 128 false positive instances (some molecules have multiple false positives)**

## Verification Steps

After implementing the fix:
1. Re-export HTML results as `html_results_fixed.csv`
2. Run: `python test_html_vs_python.py --compare html_results_fixed.csv`
3. Expected result: **0 mismatches**
