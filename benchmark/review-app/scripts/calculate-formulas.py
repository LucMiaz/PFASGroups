#!/usr/bin/env python
"""
Batch calculate molecular formulas for all molecules in the database.
This script reads SMILES from the database, calculates formulas using RDKit,
and updates the database with the results.
"""

import sqlite3
import sys
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("Error: RDKit is not installed. Please install it with:")
    print("  conda install -c conda-forge rdkit")
    sys.exit(1)

def calculate_formula(smiles):
    """Calculate molecular formula from SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
        return formula
    except Exception as e:
        print(f"Error calculating formula for {smiles}: {e}")
        return None

def main():
    # Path to the database
    db_path = Path(__file__).parent.parent / "database" / "pfas_benchmark.db"
    
    if not db_path.exists():
        print(f"Error: Database not found at {db_path}")
        print("Please run the import script first to create the database.")
        print(f"  cd {Path(__file__).parent.parent}")
        print("  node scripts/import-benchmark-data.js")
        sys.exit(1)
    
    print(f"Connecting to database: {db_path}")
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    # Get all molecules without formulas
    cursor.execute("""
        SELECT id, smiles 
        FROM molecules 
        WHERE molecular_formula IS NULL
        ORDER BY id
    """)
    
    molecules = cursor.fetchall()
    total = len(molecules)
    
    if total == 0:
        print("No molecules need formula calculation. All done!")
        conn.close()
        return
    
    print(f"Found {total} molecules without formulas")
    print("Calculating formulas...")
    
    updated = 0
    failed = 0
    
    for i, (mol_id, smiles) in enumerate(molecules, 1):
        if i % 100 == 0:
            print(f"Progress: {i}/{total} ({i*100//total}%)")
        
        formula = calculate_formula(smiles)
        
        if formula:
            cursor.execute("""
                UPDATE molecules 
                SET molecular_formula = ? 
                WHERE id = ?
            """, (formula, mol_id))
            updated += 1
        else:
            failed += 1
    
    # Commit changes
    conn.commit()
    conn.close()
    
    print("\n" + "="*60)
    print("Formula calculation complete!")
    print(f"  Successfully calculated: {updated}")
    print(f"  Failed: {failed}")
    print(f"  Total processed: {total}")
    print("="*60)

if __name__ == "__main__":
    main()
