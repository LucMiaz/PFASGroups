"""
Command-line utility to parse a CSV file of SMILES and output PFAS group matches.

Usage:
    python parse_pfas_csv.py input.csv output.csv

The input CSV must have a column named 'smiles'.
"""
import sys
import pandas as pd
from PFASgroups.core import parse_PFAS_groups
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python parse_pfas_csv.py input.csv output.csv")
        sys.exit(1)
    infile = sys.argv[1]
    outfile = sys.argv[2]
    df = pd.read_csv(infile)
    results = []
    for idx, row in df.iterrows():
        smiles = row['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results.append("")
            continue
        formula = CalcMolFormula(mol)
        matches = parse_PFAS_groups(mol, formula)
        group_names = ",".join([pf.name for pf, n, m in matches])
        results.append(group_names)
    df['pfas_groups'] = results
    df.to_csv(outfile, index=False)
