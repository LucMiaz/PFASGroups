# PFAS Groups Taxonomy

## Overview
This document provides a comprehensive taxonomy of PFAS groups defined in `PFAS_groups_smarts.json`, organized by OECD definitions (IDs 1-28) and generic functional groups (IDs 29-51).

---

## OECD PFAS GROUPS (IDs 1-28)

### 1. PERFLUOROALKYL ACIDS

#### 1.1 Carboxylic Acids
**Group 1: Perfluoroalkyl carboxylic acids (PFCAs)**
- **ID**: 1
- **SMARTS1**: `[#6$([#6](=O)([OH1]))]` (carboxylic acid)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F) ≥ C*0.5/2` (perfluorinated)
- **Dependencies**: 
  - Parent: carboxylic acid (33)
  - Children: Perfluoroalkylether carboxylic acids (4)

**Group 2: Polyfluoroalkyl carboxylic acid (PolyFCAs)**
- **ID**: 2
- **SMARTS1**: `[#6$([#6](=O)([OH1]))]` (carboxylic acid)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F,H,Cl,Br,I) ≥ C*0/2` (polyfluorinated)
- **Dependencies**:
  - Parent: carboxylic acid (33)
  - Children: Polyfluoroalkylether carboxylic acid (5)

**Group 3: Perfluoroalkyl dicarboxylic acids (PFdiCAs)**
- **ID**: 3
- **SMARTS1**: `[#6$([#6](=O)([OH1]))]` (first carboxylic acid)
- **SMARTS2**: `[#6$([#6](=O)([OH1]))]` (second carboxylic acid)
- **Path**: `per` (perfluorinated path between SMARTS)
- **Constraints**: `C_rel: (F) ≥ (C+2)*0.5/2`
- **Dependencies**:
  - Parent: carboxylic acid (33)

**Group 4: Perfluoroalkylether carboxylic acids (PFECAs)**
- **ID**: 4
- **SMARTS1**: `[#6$([#6](=O)([OH1]))]` (carboxylic acid)
- **SMARTS2**: `[#6$([#6]O[#6])]` (ether)
- **Path**: `per` (perfluorinated path)
- **Constraints**: `C_rel: (F) ≥ C*0.5/2`
- **Dependencies**:
  - Parents: Perfluoroalkyl carboxylic acids (1), ether (31)

**Group 5: Polyfluoroalkylether carboxylic acid (PolyFECAs)**
- **ID**: 5
- **SMARTS1**: `[#6$([#6](=O)([OH1]))]` (carboxylic acid)
- **SMARTS2**: `[#6$([#6]O[#6])]` (ether)
- **Path**: `poly` (polyfluorinated path)
- **Constraints**: `C_rel: (F,H,Cl,Br,I) ≥ C*0/2`
- **Dependencies**:
  - Parents: Polyfluoroalkyl carboxylic acid (2), ether (31)

#### 1.2 Sulfonic Acids
**Group 6: Perfluoroalkyl sulfonic acids (PFSAs)**
- **ID**: 6
- **SMARTS1**: `[#6$([#6][#16](=O)(=O)[OH1])]` (sulfonic acid)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F) ≥ (C-0.5)*0.5/2`
- **Dependencies**:
  - Parent: sulfonic acid (36)

**Group 7: Polyfluoroalkyl sulfonic acid (PolyFSAs)**
- **ID**: 7
- **SMARTS1**: `[#6$([#6][#16](=O)(=O)[OH1])]` (sulfonic acid)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F,H,Cl,Br,I) ≥ (C-1)*0.5/2`
- **Dependencies**:
  - Parent: sulfonic acid (36)
  - Children: Polyfluoroalkylether sulfonic acid (10)

**Group 8: Perfluoroalkyl disulfonic acids (PFdiSAs)**
- **ID**: 8
- **SMARTS1**: `[#6$([#6][#16](=O)(=O)[OH1])]` (first sulfonic acid)
- **SMARTS2**: `[#6$([#6][#16](=O)(=O)[OH1])]` (second sulfonic acid)
- **Path**: `per` (perfluorinated path)
- **Constraints**: `C_rel: (F) ≥ C*0.5/2`
- **Dependencies**:
  - Parent: sulfonic acid (36)

**Group 9: Perfluoroalkylether sulfonic acids (PFESAs)**
- **ID**: 9
- **SMARTS1**: `[#6$([#6][#16](=O)(=O)([OH]))]` (sulfonic acid)
- **SMARTS2**: `[#6$([#6]O[#6])]` (ether)
- **Path**: `per` (perfluorinated path)
- **Constraints**: `C_rel: (F) ≥ (C-0.5)*0.5/2`
- **Dependencies**:
  - Parents: sulfonic acid (36), ether (31)

**Group 10: Polyfluoroalkylether sulfonic acid (PolyFESAs)**
- **ID**: 10
- **SMARTS1**: `[#6$([#6][#16](=O)(=O)([OH]))]` (sulfonic acid)
- **SMARTS2**: `[#6$([#6]O[#6])]` (ether)
- **Path**: `poly` (polyfluorinated path)
- **Constraints**: `C_rel: (F,H,Cl,Br,I) ≥ (C-1)*0.5/2`
- **Dependencies**:
  - Parents: Polyfluoroalkyl sulfonic acid (7), sulfonic acid (36), ether (31)

#### 1.3 Other Acids
**Group 11: Perfluoroalkyl sulfinic acids (PFSiAs)**
- **ID**: 11
- **SMARTS1**: `[#6$([#6][#16X3](=O)[OH1])]` (sulfinic acid)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F) ≥ (C-0.5)*0.5/2`
- **Dependencies**:
  - Parent: sulfinic acid (38)

**Group 12: Perfluoroalkyl phosphonic acids (PFPAs)**
- **ID**: 12
- **SMARTS1**: `[#6$([#6][#15](=[#8])([#8])[#8])]` (phosphonic acid)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F) ≥ (C-0.5)*0.5/2`
- **Dependencies**:
  - Parent: phosphonic acid (39)

**Group 13: Perfluoroalkyl phosphinic acids (PFPiAs)**
- **ID**: 13
- **SMARTS1**: `[#6X4$([#6][#15](=[#8])([#8])[#1,#6])]` (phosphinic acid)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F) ≥ (C-0.5)*0.5/2`
- **Dependencies**:
  - Parent: phosphinic acid (40)

---

### 2. PERFLUOROALKYL ACIDS PRECURSORS

#### 2.1 Alcohols
**Group 14: Perfluoroalkyl alcohols**
- **ID**: 14
- **SMARTS1**: `[#6X4$([#6;!$([#6]=O)][OH1])]` (alcohol, not on carbonyl)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F,O) ≥ (C-1)*0.5/2`
- **Dependencies**:
  - Parent: alcohol (29)

**Group 15: fluorotelomer alcohols**
- **ID**: 15
- **SMARTS1**: `[C$([#6][C$([CH2]O[H]),...])]` (long SMARTS for telomer)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: None
- **Dependencies**:
  - Parent: alcohol (29)

#### 2.2 Ethers
**Group 16: Perfluoropolyethers (PFPEs)**
- **ID**: 16
- **SMARTS1**: `[#6$([#6]O[#6])]` (first ether)
- **SMARTS2**: `[#6$([#6]O[#6])]` (second ether)
- **Path**: `per` (perfluorinated path)
- **Constraints**: `O ≥ 2`, `C_rel: (F) ≥ (C-1)*0.5/2`
- **Dependencies**:
  - Parent: ether (31)
  - Related: Hydrofluoroethers (17)

**Group 17: Hydrofluoroethers (HFEs)**
- **ID**: 17
- **SMARTS1**: `[#6][#8][#6]` (ether)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `O ≥ 1`, `only: C,F,O,H`, `C_rel: (F,H) ≥ (C-1)*0.5/2`
- **Dependencies**:
  - Parent: ether (31)
  - Related: Perfluoropolyethers (16)

#### 2.3 Alkenes
**Group 18: Perfluoroalkene**
- **ID**: 18
- **SMARTS1**: `[#6X3$([#6]=[#6X3])]` (alkene)
- **SMARTS2**: None
- **Path**: `per`
- **Constraints**: `C_rel: (F) ≥ C*0.5/2`
- **Dependencies**:
  - Parent: alkene (49)

**Group 19: Hydrofluoroolefins (HFOs)**
- **ID**: 19
- **SMARTS1**: `[#6X3$([#6]=[#6])]` (alkene)
- **SMARTS2**: `[#6X4$([#6]F)]` (C-F bond)
- **Path**: None
- **Constraints**: `F ≥ 1`, `only: C,F,H`
- **Dependencies**:
  - Parent: alkene (49)
  - Related: ethene (41)

#### 2.4 Alkanes
**Group 20: Hydrofluorocarbons (HFCs)**
- **ID**: 20
- **SMARTS1**: None
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `F ≥ 1`, `only: C,F,H`, `C_rel: (F,H) ≥ (C-1)*0.5/2`
- **Dependencies**: None (formula-based only)

**Group 21: Semi-fluorinated alkanes (SFAs)**
- **ID**: 21
- **SMARTS1**: `[#6X4$([#6](F)(F))]` (CF2 group)
- **SMARTS2**: `[#6X4$([#6][H,I,Br,Cl])]` (non-perfluorinated carbon)
- **Path**: `poly` (polyfluorinated path)
- **Constraints**: `F ≥ 1`, `C_rel: (Cl,I,Br,F,H) ≥ (C-1)*0.5/2`
- **Dependencies**:
  - Parent: alkane (48)

#### 2.5 Halides & Others
**Group 25: Perfluoroalkyl iodides (PFAIs)**
- **ID**: 25
- **SMARTS1**: None
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `I ≥ 1`, `C_rel: (F,I) ≥ (C-1)*0.5/2`
- **Dependencies**:
  - Parent: iodide (42)

**Group 26: Perfluoroalkane sulfonyl fluorides (PASFs)**
- **ID**: 26
- **SMARTS1**: `[#16$([#16](=O)(=O)[F])]` (sulfonyl fluoride)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F) ≥ (C-1)*0.5/2`
- **Dependencies**: None

**Group 27: Perfluoroalkyl ketones**
- **ID**: 27
- **SMARTS1**: `[#6X4$([#6]C(=O)[#6])]` (ketone)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F) ≥ C*0.5/2`
- **Dependencies**:
  - Parent: ketone (30)

**Group 28: Semi-fluoroalkyl ketones**
- **ID**: 28
- **SMARTS1**: `[#6X4$([#6]C(=O)[#6])]` (ketone)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `O ≥ 1`, `F ≥ 3`, `C_rel: (Cl,I,Br,F,H) ≥ (C-1)*0.5/2` (includes O in add_atoms)
- **Dependencies**:
  - Parent: ketone (30)

---

### 3. OTHER PFAS

**Group 22: Side-chain fluorinated aromatics**
- **ID**: 22
- **SMARTS1**: `[#6X4$([#6]a)]` (aliphatic carbon attached to aromatic)
- **SMARTS2**: None
- **Path**: `per`
- **Constraints**: None
- **Dependencies**:
  - Parent: Side-chain aromatics (51)
  - Related: azole (44), azine (45), benzodioxole (46)

**Group 23: Perfluoroalkane**
- **ID**: 23
- **SMARTS1**: None
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F) ≥ (C-1)*0.5/2`
- **Dependencies**: None (formula-based only)

**Group 24: Perfluoroalkyl-tert-amines**
- **ID**: 24
- **SMARTS1**: `[#6$([#6][#7D3])]` (tertiary amine)
- **SMARTS2**: None
- **Path**: None
- **Constraints**: `C_rel: (F) ≥ (C-1.5)*0.5/2`
- **Dependencies**:
  - Parent: amine (47)

---

## GENERIC PFAS GROUPS (IDs 29-51)

### 4. FUNCTIONAL GROUPS

#### 4.1 Oxygen-Containing
**Group 29: alcohol**
- **ID**: 29
- **SMARTS1**: `[#6$([#6;!$([#6]=O)][OH1])]` (OH not on carbonyl)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: Perfluoroalkyl alcohols (14), fluorotelomer alcohols (15)

**Group 30: ketone**
- **ID**: 30
- **SMARTS1**: `[#6$([#6][#6](=O)[#6])]` (C-CO-C)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: Perfluoroalkyl ketones (27), Semi-fluoroalkyl ketones (28)

**Group 31: ether**
- **ID**: 31
- **SMARTS1**: `[#6$([#6]([#8X2H0][#6;!$(C=O)])[#1,!$(C=O)]);#6!$([#6](=O))]` (C-O-C, not ester)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: 
  - Hydrofluoroethers (17)
  - Perfluoropolyethers (16)
  - Perfluoroalkylether carboxylic acids (4)
  - Polyfluoroalkylether carboxylic acid (5)
  - Perfluoroalkylether sulfonic acids (9)
  - Polyfluoroalkylether sulfonic acid (10)
  - benzodioxole (46)

**Group 32: ester**
- **ID**: 32
- **SMARTS1**: `[#6X4$([#6][#6](=O)O[#6,...]),#6X4$([#6]O[#6](=O)[#6,...])]` (ester)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: None

**Group 33: carboxylic acid**
- **ID**: 33
- **SMARTS1**: `[$(*[#6](=O)[OH1])]` (COOH)
- **SMARTS2**: None
- **Constraints**: None
- **Children**:
  - Perfluoroalkyl carboxylic acids (1)
  - Polyfluoroalkyl carboxylic acid (2)
  - Perfluoroalkyl dicarboxylic acids (3)

**Group 34: amide**
- **ID**: 34
- **SMARTS1**: `[$([#6][#6]([#7H0,#7H1,#7H2])=[#8]),$([#6]~[#7H0,#7H1,#7H2][#6X3]=[#8])]` (amide)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: None

**Group 35: acyl halide**
- **ID**: 35
- **SMARTS1**: `[#6$([#6][#6](=O)[#9,#17,#35,#53])]` (COX where X=halogen)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: None

#### 4.2 Sulfur-Containing
**Group 36: sulfonic acid**
- **ID**: 36
- **SMARTS1**: `[$(*[#16](=O)(=O)[O,OH1])]` (SO3H)
- **SMARTS2**: None
- **Constraints**: None
- **Children**:
  - Perfluoroalkyl sulfonic acids (6)
  - Polyfluoroalkyl sulfonic acid (7)
  - Perfluoroalkyl disulfonic acids (8)
  - Perfluoroalkylether sulfonic acids (9)
  - Polyfluoroalkylether sulfonic acid (10)

**Group 37: sulfenic acid**
- **ID**: 37
- **SMARTS1**: `[$(*[#16X2][OH1])]` (SOH)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: None

**Group 38: sulfinic acid**
- **ID**: 38
- **SMARTS1**: `[$(*[#16X3](=O)[O])]` (SO2H)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: Perfluoroalkyl sulfinic acids (11)

**Group 43: sulfonamide**
- **ID**: 43
- **SMARTS1**: `[$(*[#16]([#7])(=[#8])=[#8]),$(*[#7][#16](=[#8])=[#8])]` (SO2NH)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: None
- **Parent**: amine (47)

#### 4.3 Phosphorus-Containing
**Group 39: phosphonic acid**
- **ID**: 39
- **SMARTS1**: `[#6$([#6][#15]([#8])(=[#8])[#8])]` (PO(OH)2)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: Perfluoroalkyl phosphonic acids (12)

**Group 40: phosphinic acid**
- **ID**: 40
- **SMARTS1**: `[#6$([#6][#15](=[#8])([#8])[#1;!$([#8])])]` (P(O)(OH)H)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: Perfluoroalkyl phosphinic acids (13)

#### 4.4 Nitrogen-Containing
**Group 47: amine**
- **ID**: 47
- **SMARTS1**: `[$([#6;!$([#6]=O)][N;!$(*[#6]=O)])]` (C-N, not amide)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: 
  - Perfluoroalkyl-tert-amines (24)
  - sulfonamide (43)

#### 4.5 Halogen-Containing
**Group 42: iodide**
- **ID**: 42
- **SMARTS1**: `[#6$([#6][#6](F)[#53])]` (C-I near fluorine)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: Perfluoroalkyl iodides (25)

#### 4.6 Unsaturated
**Group 41: ethene**
- **ID**: 41
- **SMARTS1**: `[#6$([#6X3]([F,#6$([#6][#6]F),#6$([#6][#6][#6]F)])=[#6X3][F,...])]` (fluorinated ethene)
- **SMARTS2**: None
- **Constraints**: `F ≥ 2`
- **Children**: None
- **Parent**: alkene (49)

**Group 49: alkene**
- **ID**: 49
- **SMARTS1**: `[#6$([#6][#6X3]=[#6])]` (C=C)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: 
  - Perfluoroalkene (18)
  - Hydrofluoroolefins (19)
  - ethene (41)

**Group 50: alkyne**
- **ID**: 50
- **SMARTS1**: `[#6$([#6][#6]#[#6])]` (C≡C)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: None

#### 4.7 Saturated & Aromatic
**Group 48: alkane**
- **ID**: 48
- **SMARTS1**: `[C$([CX4](F)(F)([F,!$(*[#9])])[CX4](F)([F,C,O,S,N,Si,P,B,Se,Te])[F,...])]` (fluorinated alkane)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: Semi-fluorinated alkanes (21)

**Group 44: azole**
- **ID**: 44
- **SMARTS1**: `[$(*[#7r5,$([#6r5,#8r5,#16r5]:[#7r5]),$(...:[#7r5])])]` (5-membered N-heteroaromatic)
- **SMARTS2**: None
- **Constraints**: None
- **Definition**: Five-membered heterocyclic ring with ≥1 nitrogen
- **Parents**: Side-chain aromatics (51), Side-chain fluorinated aromatics (22)

**Group 45: azine**
- **ID**: 45
- **SMARTS1**: `[$(*[#7r6,$([#6r6,#8r6,#16r6]:[#7r6]),$(...:[#7r6])])]` (6-membered N-heteroaromatic)
- **SMARTS2**: None
- **Constraints**: None
- **Definition**: Six-membered heterocyclic ring with ≥1 nitrogen (e.g., pyridine)
- **Parents**: Side-chain aromatics (51), Side-chain fluorinated aromatics (22)

**Group 46: benzodioxole**
- **ID**: 46
- **SMARTS1**: `[$([#6]1[#6][#6]([#6]2[#6]([#6]1)[#8][#6][#8]2))]` (benzodioxole ring)
- **SMARTS2**: None
- **Constraints**: `F ≥ 2`
- **Definition**: Fluorinated benzodioxole
- **Parents**: ether (31), Side-chain aromatics (51), Side-chain fluorinated aromatics (22)

**Group 51: Side-chain aromatics**
- **ID**: 51
- **SMARTS1**: `[#6X4$([#6]a)]` (aliphatic carbon attached to aromatic)
- **SMARTS2**: None
- **Constraints**: None
- **Children**: 
  - Side-chain fluorinated aromatics (22)
  - azole (44)
  - azine (45)
  - benzodioxole (46)

---

## HIERARCHICAL RELATIONSHIPS

### Parent → Child Dependencies

```
Generic Groups → OECD Groups

carboxylic acid (33)
├── Perfluoroalkyl carboxylic acids (1)
│   └── Perfluoroalkylether carboxylic acids (4)
├── Polyfluoroalkyl carboxylic acid (2)
│   └── Polyfluoroalkylether carboxylic acid (5)
└── Perfluoroalkyl dicarboxylic acids (3)

sulfonic acid (36)
├── Perfluoroalkyl sulfonic acids (6)
├── Polyfluoroalkyl sulfonic acid (7)
│   └── Polyfluoroalkylether sulfonic acid (10)
├── Perfluoroalkyl disulfonic acids (8)
└── Perfluoroalkylether sulfonic acids (9)

ether (31)
├── Hydrofluoroethers (17)
│   └── Perfluoropolyethers (16) [related]
├── Perfluoropolyethers (16)
├── benzodioxole (46)
└── [connects to PFECA(4), PolyFECA(5), PFESA(9), PolyFESA(10)]

alcohol (29)
├── Perfluoroalkyl alcohols (14)
└── fluorotelomer alcohols (15)

ketone (30)
├── Perfluoroalkyl ketones (27)
└── Semi-fluoroalkyl ketones (28)

alkene (49)
├── Perfluoroalkene (18)
├── Hydrofluoroolefins (19)
└── ethene (41)

alkane (48)
└── Semi-fluorinated alkanes (21)

amine (47)
├── Perfluoroalkyl-tert-amines (24)
└── sulfonamide (43)

iodide (42)
└── Perfluoroalkyl iodides (25)

sulfinic acid (38)
└── Perfluoroalkyl sulfinic acids (11)

phosphonic acid (39)
└── Perfluoroalkyl phosphonic acids (12)

phosphinic acid (40)
└── Perfluoroalkyl phosphinic acids (13)

Side-chain aromatics (51)
├── Side-chain fluorinated aromatics (22)
│   ├── azole (44)
│   ├── azine (45)
│   └── benzodioxole (46)
├── azole (44)
├── azine (45)
└── benzodioxole (46)
```

---

## KEY CONCEPTS

### Path Types
- **`per`** (perfluorinated): Path between SMARTS1 and SMARTS2 must be fully fluorinated
- **`poly`** (polyfluorinated): Path between SMARTS1 and SMARTS2 can contain H, Cl, Br, I in addition to F

### Constraint Types
- **`rel`** (relative): Formula ratio constraints (e.g., F atoms ≥ some fraction of C atoms)
- **`gte`** (greater than or equal): Minimum count (e.g., F ≥ 1, O ≥ 2)
- **`lte`** (less than or equal): Maximum count
- **`eq`** (equal): Exact count
- **`only`**: Only these elements allowed

### Fluorination Levels
- **Perfluorinated**: All available positions occupied by fluorine (typically C_rel: F ≥ ~0.5C)
- **Polyfluorinated**: Contains both F and H/Cl/Br/I (typically C_rel: (F,H,X) ≥ ~0.5C)
- **Semi-fluorinated**: Mix of perfluorinated and non-fluorinated regions

---

## SUMMARY STATISTICS

- **Total OECD Groups**: 26 (IDs 1-28, missing 20,23)
- **Total Generic Groups**: 23 (IDs 29-51)
- **Groups with SMARTS1 only**: 35
- **Groups with SMARTS1 + SMARTS2**: 12
- **Formula-only groups**: 4 (IDs 20, 23, 25, 26 partially)
- **Groups requiring paths**: 12 (using `per` or `poly`)

