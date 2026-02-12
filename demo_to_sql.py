"""
Standalone demonstration of the to_sql functionality for PFASGroups.

This script demonstrates the new to_sql methods added to both MoleculeResult 
and ResultsModel classes. It creates mock data to show the functionality 
without requiring a full PFASgroups environment.
"""

import sqlite3
import tempfile
import os
from pathlib import Path


def demonstrate_to_sql_functionality():
    """
    Demonstrate the to_sql functionality with mock data.
    
    This shows what the new methods do:
    1. MoleculeResult.to_sql() - exports a single molecule's PFAS group data
    2. ResultsModel.to_sql() - exports multiple molecules' data in batch
    """
    
    print("=" * 70)
    print("PFASGroups to_sql() Functionality Demonstration")
    print("=" * 70)
    print()
    
    # Create a temporary database
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "pfas_demo.db"
        
        print(f"Creating demonstration database at: {db_path}")
        print()
        
        # Mock data structure matching what PFASgroups produces
        mock_results = [
            {
                "smiles": "C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O",
                "matches": [
                    {
                        "match_id": "G1",
                        "id": 1,
                        "group_name": "Perfluoroalkyl carboxylic acid",
                        "type": "PFASgroup",
                        "components": [
                            {"component": [0, 1, 2, 3, 4, 5], "SMARTS": "Perfluoroalkyl"},
                            {"component": [6, 7, 8], "SMARTS": "Carboxylic acid"},
                        ],
                    }
                ],
            },
            {
                "smiles": "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)S(=O)(=O)O",
                "matches": [
                    {
                        "match_id": "G2",
                        "id": 2,
                        "group_name": "Perfluoroalkyl sulfonic acid",
                        "type": "PFASgroup",
                        "components": [
                            {"component": [0, 1, 2, 3, 4, 5, 6, 7], "SMARTS": "Perfluoroalkyl"},
                            {"component": [8, 9, 10, 11], "SMARTS": "Sulfonic acid"},
                        ],
                    }
                ],
            },
        ]
        
        # Create database and populate with mock data
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        
        # Create tables
        cursor.execute("""
            CREATE TABLE components (
                smiles TEXT,
                group_id INTEGER,
                group_name TEXT,
                smarts_label TEXT,
                component_atoms TEXT
            )
        """)
        
        cursor.execute("""
            CREATE TABLE pfas_groups_in_compound (
                smiles TEXT,
                group_id INTEGER,
                group_name TEXT,
                match_count INTEGER
            )
        """)
        
        # Insert mock data
        for mol in mock_results:
            smiles = mol["smiles"]
            
            for match in mol["matches"]:
                if match["type"] != "PFASgroup":
                    continue
                    
                # Insert group summary
                cursor.execute("""
                    INSERT INTO pfas_groups_in_compound 
                    (smiles, group_id, group_name, match_count)
                    VALUES (?, ?, ?, ?)
                """, (smiles, match["id"], match["group_name"], 1))
                
                # Insert component details
                for comp in match["components"]:
                    atoms_str = ','.join(map(str, comp["component"]))
                    cursor.execute("""
                        INSERT INTO components 
                        (smiles, group_id, group_name, smarts_label, component_atoms)
                        VALUES (?, ?, ?, ?, ?)
                    """, (smiles, match["id"], match["group_name"], 
                          comp["SMARTS"], atoms_str))
        
        conn.commit()
        
        # Display results
        print("Database created successfully!")
        print()
        print("-" * 70)
        print("PFAS Groups in Compound Table:")
        print("-" * 70)
        
        cursor.execute("SELECT * FROM pfas_groups_in_compound")
        rows = cursor.fetchall()
        print(f"{'SMILES':<50} {'Group ID':<10} {'Group Name':<30} {'Count':<5}")
        print("-" * 70)
        for row in rows:
            smiles, gid, gname, count = row
            print(f"{smiles:<50} {gid:<10} {gname:<30} {count:<5}")
        
        print()
        print("-" * 70)
        print("Components Table:")
        print("-" * 70)
        
        cursor.execute("SELECT * FROM components")
        rows = cursor.fetchall()
        print(f"{'SMILES':<50} {'SMARTS':<25} {'Atoms':<15}")
        print("-" * 70)
        for row in rows:
            smiles, gid, gname, smarts, atoms = row
            print(f"{smiles:<50} {smarts:<25} {atoms:<15}")
        
        conn.close()
        
        print()
        print("=" * 70)
        print("Demonstration Complete!")
        print("=" * 70)
        print()
        print("The new to_sql() methods added to PFASGroups:")
        print()
        print("1. MoleculeResult.to_sql()")
        print("   - Exports a single molecule's PFAS group matches to SQL")
        print("   - Supports SQLite (filename) and PostgreSQL (connection params)")
        print("   - Parameters:")
        print("     * filename: Path to SQLite database")
        print("     * dbname, user, password, host, port: For PostgreSQL/MySQL")
        print("     * components_table: Name for components table (default: 'components')")
        print("     * groups_table: Name for groups table (default: 'pfas_groups_in_compound')")
        print("     * if_exists: 'append', 'replace', or 'fail'")
        print()
        print("2. ResultsModel.to_sql()")
        print("   - Efficiently exports all molecules in batch")
        print("   - Same parameters as MoleculeResult.to_sql()")
        print("   - Optimized for bulk operations")
        print()
        print("Environment variable support:")
        print("   - DB_USER: Default database username")
        print("   - DB_PASSWORD: Default database password")
        print("   - DB_HOST: Default host (defaults to 'localhost')")
        print("   - DB_PORT: Default port (defaults to 5432)")
        print()
        print("Example usage:")
        print()
        print("  # SQLite")
        print("  results = parse_smiles(smiles_list)")
        print("  results.to_sql(filename='pfas_results.db')")
        print()
        print("  # PostgreSQL")
        print("  results.to_sql(")
        print("      dbname='pfas_db',")
        print("      user='myuser',")
        print("      password='mypass',  # or set DB_PASSWORD env var")
        print("      host='localhost',")
        print("      port=5432")
        print("  )")
        print()


if __name__ == "__main__":
    demonstrate_to_sql_functionality()
