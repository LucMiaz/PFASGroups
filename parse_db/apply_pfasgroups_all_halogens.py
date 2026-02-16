#!/usr/bin/env python3
"""Apply PFASGroups parsing to all halogenated compounds in the clinventory database.

This script processes compounds containing F, Cl, Br, or I using PFASGroups,
optionally substituting other halogens with F to leverage the existing PFAS
classification infrastructure.

The script queries the clinventory database for compounds with halogen elements
and processes them using PFASGroups parse_mol function.
"""

import argparse
import datetime as dt
import json as json_module
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import psycopg2
from psycopg2.extras import execute_batch


@dataclass
class DBConfig:
    dbname: str
    user: str
    password: str
    host: str
    port: int


def _default_pfasgroups_path() -> Path:
    # Script is in PFASGroups root, so use current directory
    return Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply PFASGroups parsing to all halogenated compounds in clinventory database"
    )
    parser.add_argument("--db-name", default=os.environ.get("DATABASE_NAME", "clinventory"))
    parser.add_argument("--db-user", default=os.environ.get("DATABASE_USER", "django"))
    parser.add_argument("--db-password", default=os.environ.get("DJANGO_USER", ""))
    parser.add_argument("--db-host", default=os.environ.get("DATABASE_HOST", "localhost"))
    parser.add_argument("--db-port", type=int, default=int(os.environ.get("DATABASE_PORT", "5432")))

    parser.add_argument(
        "--pfasgroups-path",
        type=Path,
        default=_default_pfasgroups_path(),
        help="Path to local PFASGroups repository (default: current directory)",
    )
    
    parser.add_argument(
        "--table-name",
        default="molecules",
        help="Database table name to query (default: molecules)",
    )
    
    parser.add_argument(
        "--id-column",
        default="id",
        help="ID column name (default: id)",
    )
    
    parser.add_argument(
        "--smiles-column",
        default="smiles",
        help="SMILES column name (default: smiles)",
    )
    
    parser.add_argument(
        "--formula-column",
        default="formula",
        help="Formula column name (default: formula)",
    )
    
    parser.add_argument(
        "--use-halogen-columns",
        action="store_true",
        help="Use numeric halogen columns (F, Cl, Br, I) instead of regex on formula (for databases like zeropmdbwp2)",
    )

    parser.add_argument(
        "--run-id",
        default=f"pfasgroups_halogens_{dt.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')}",
        help="Identifier for this run",
    )

    parser.add_argument("--limit", type=int, default=None, help="Max number of compounds to process")
    parser.add_argument("--offset", type=int, default=0)
    
    parser.add_argument(
        "--halogens",
        default="F,Cl,Br,I",
        help="Comma-separated list of halogens to include (default: F,Cl,Br,I)",
    )
    
    parser.add_argument(
        "--substitute-halogens",
        action="store_true",
        help="Substitute Cl/Br/I with F for parsing (to use fluorine-based SMARTS)",
    )
    
    parser.add_argument(
        "--output-json",
        type=Path,
        help="Path to save results as JSON file",
    )
    
    parser.add_argument(
        "--save-to-db",
        action="store_true",
        help="Save results to database table pfasgroups_results",
    )
    
    parser.add_argument(
        "--results-table",
        default="pfasgroups_results",
        help="Database table name for storing results (default: pfasgroups_results)",
    )
    
    parser.add_argument("--batch-size", type=int, default=100, help="Processing batch size")
    
    return parser.parse_args()


def build_db_config(args: argparse.Namespace) -> DBConfig:
    return DBConfig(
        dbname=args.db_name,
        user=args.db_user,
        password=args.db_password,
        host=args.db_host,
        port=args.db_port,
    )


def import_pfasgroups(pfasgroups_path: Path):
    """Import PFASGroups and suppress RDKit warnings."""
    if not pfasgroups_path.exists():
        raise FileNotFoundError(f"PFASGroups path does not exist: {pfasgroups_path}")

    sys.path.insert(0, str(pfasgroups_path))

    # Import and suppress RDKit warnings
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    
    from PFASgroups import parse_mol, get_PFASGroups  # type: ignore
    from PFASgroups.results_model import ResultsModel  # type: ignore

    return parse_mol, get_PFASGroups, ResultsModel


def initialize_database_tables(conn, pfas_groups_list):
    """Initialize the three required tables: pfasgroups, pfasgroups_in_molecules, components_in_molecules."""
    
    with conn.cursor() as cur:
        # Table 1: pfasgroups - PFAS group definitions
        cur.execute("""
        CREATE TABLE IF NOT EXISTS pfasgroups (
            id INTEGER PRIMARY KEY,
            name TEXT NOT NULL,
            alias TEXT,
            category TEXT,
            pathway_type TEXT,
            smarts_primary TEXT,
            smarts_secondary TEXT,
            created_at TIMESTAMP DEFAULT NOW()
        );
        """)
        
        # Table 2: pfasgroups_in_molecules - PFAS group matches in molecules
        cur.execute("""
        CREATE TABLE IF NOT EXISTS pfasgroups_in_molecules (
            id SERIAL PRIMARY KEY,
            smiles TEXT NOT NULL,
            molecule_id INTEGER,
            group_id INTEGER,
            group_name TEXT,
            match_count INTEGER,
            created_at TIMESTAMP DEFAULT NOW()
        );
        
        CREATE INDEX IF NOT EXISTS idx_pfasgroups_in_molecules_smiles 
            ON pfasgroups_in_molecules(smiles);
        CREATE INDEX IF NOT EXISTS idx_pfasgroups_in_molecules_molecule_id 
            ON pfasgroups_in_molecules(molecule_id);
        CREATE INDEX IF NOT EXISTS idx_pfasgroups_in_molecules_group_id 
            ON pfasgroups_in_molecules(group_id);
        """)
        
        # Table 3: components_in_molecules - Component details
        cur.execute("""
        CREATE TABLE IF NOT EXISTS components_in_molecules (
            id SERIAL PRIMARY KEY,
            smiles TEXT NOT NULL,
            molecule_id INTEGER,
            group_id INTEGER,
            group_name TEXT,
            smarts_label TEXT,
            component_atoms TEXT,
            created_at TIMESTAMP DEFAULT NOW()
        );
        
        CREATE INDEX IF NOT EXISTS idx_components_in_molecules_smiles 
            ON components_in_molecules(smiles);
        CREATE INDEX IF NOT EXISTS idx_components_in_molecules_molecule_id 
            ON components_in_molecules(molecule_id);
        CREATE INDEX IF NOT EXISTS idx_components_in_molecules_group_id 
            ON components_in_molecules(group_id);
        """)
    
    conn.commit()
    
    # Populate pfasgroups table with group definitions
    with conn.cursor() as cur:
        # Check if table is empty
        cur.execute("SELECT COUNT(*) FROM pfasgroups")
        count = cur.fetchone()[0]
        
        if count == 0:
            print(f"Populating pfasgroups table with {len(pfas_groups_list)} group definitions...")
            for group in pfas_groups_list:
                cur.execute("""
                INSERT INTO pfasgroups (id, name, alias, category, pathway_type, smarts_primary, smarts_secondary)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (id) DO NOTHING
                """, (
                    group.id,
                    group.name,
                    getattr(group, 'alias', None),
                    getattr(group, 'category', None),
                    getattr(group, 'pathway_type', None),
                    str(group.smarts) if hasattr(group, 'smarts') else None,
                    str(getattr(group, 'smarts2', None)) if hasattr(group, 'smarts2') else None,
                ))
            conn.commit()
            print(f"Populated pfasgroups table with {len(pfas_groups_list)} groups")
        else:
            print(f"pfasgroups table already contains {count} groups")


def add_molecule_id_to_results(results_model, compound_id_map):
    """Add molecule_id to each result based on SMILES."""
    # ResultsModel.to_sql() uses smiles as identifier, but we want to also track molecule_id
    # We'll add molecule_id after the fact using UPDATE statements
    return compound_id_map


def fetch_halogenated_compounds(
    conn,
    halogens: List[str],
    table_name: str,
    id_column: str,
    smiles_column: str,
    formula_column: str,
    limit: Optional[int],
    offset: int,
    use_halogen_columns: bool = False,
) -> List[Tuple[int, str, str]]:
    """Fetch compounds containing specified halogens.
    
    Args:
        use_halogen_columns: If True, use numeric columns (F, Cl, Br, I).
                           If False, use regex pattern matching on formula column.
    """
    
    # Build WHERE clause for halogens
    if use_halogen_columns:
        # For databases with numeric halogen columns (like zeropmdbwp2)
        halogen_conditions = []
        for halogen in halogens:
            halogen_conditions.append(f"c.{halogen} > 0")
        where_clause = " OR ".join(halogen_conditions)
    else:
        # For databases with formula text column (like clinventory)
        # Use regex to match any of the halogens in the formula
        # Match halogen symbol followed by optional digit(s)
        halogen_pattern = "|".join([f"{h}[0-9]*" for h in halogens])
        where_clause = f"c.{formula_column} ~ '({halogen_pattern})'"
    
    query = f"""
    SELECT c.{id_column}, c.{smiles_column}, c.{formula_column}
    FROM {table_name} c
    WHERE c.{smiles_column} IS NOT NULL 
      AND c.{smiles_column} <> ''
      AND ({where_clause})
    ORDER BY c.{id_column}
    """
    
    params = []
    
    if limit is not None:
        query += "\nLIMIT %s"
        params.append(limit)
    
    query += "\nOFFSET %s"
    params.append(offset)
    
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()


def substitute_halogens_in_smiles(smiles: str) -> Tuple[str, Dict[str, int]]:
    """Replace Cl, Br, I with F and track substitutions."""
    from rdkit import Chem
    
    substitutions = {'Cl': 0, 'Br': 0, 'I': 0}
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles, substitutions
        
        # Count and replace halogens
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in ['Cl', 'Br', 'I']:
                substitutions[symbol] += 1
                atom.SetAtomicNum(9)  # Replace with F (atomic number 9)
        
        # Generate new SMILES with fluorine substitutions
        substituted_smiles = Chem.MolToSmiles(mol)
        return substituted_smiles, substitutions
        
    except Exception as e:
        print(f"Warning: Failed to substitute halogens in {smiles}: {e}")
        return smiles, substitutions


def parse_one_compound(
    compound_id: int,
    smiles: str,
    formula: str,
    parse_mol_func,
    substitute: bool,
) -> Dict:
    """Parse a single compound with PFASGroups."""
    from rdkit import Chem
    
    t0 = time.perf_counter()
    parse_time = 0.0
    
    try:
        # Optionally substitute halogens
        target_smiles = smiles
        if substitute:
            target_smiles, subs = substitute_halogens_in_smiles(smiles)
        
        # Parse molecule
        mol = Chem.MolFromSmiles(target_smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {target_smiles}")
        
        # Use PFASGroups parse_mol
        # Returns a MoleculeResult object (dict-like) with 'matches' list
        parsed_result = parse_mol_func(
            mol,
            bycomponent=True,
            include_PFAS_definitions=False
        )
        
        parse_time = time.perf_counter() - t0
        
        # Return the MoleculeResult directly with metadata
        return {
            'compound_id': compound_id,
            'smiles': smiles,
            'parse_time': parse_time,
            'status': 'success',
            'result': parsed_result,
        }
        
    except Exception as e:
        parse_time = time.perf_counter() - t0
        return {
            'compound_id': compound_id,
            'smiles': smiles,
            'parse_time': parse_time,
            'status': 'failed',
            'error': str(e),
            'result': None,
        }


def main() -> int:
    args = parse_args()
    db_config = build_db_config(args)
    
    # Parse halogen list
    halogens = [h.strip() for h in args.halogens.split(',')]
    print(f"Processing compounds with halogens: {', '.join(halogens)}")
    print(f"Database: {db_config.dbname} (table: {args.table_name})")
    
    # Import PFASGroups
    parse_mol_func, get_pfas_groups, ResultsModel = import_pfasgroups(args.pfasgroups_path)
    
    # Load PFAS groups for reference
    pfas_groups = get_pfas_groups()
    print(f"Loaded {len(pfas_groups)} PFAS group definitions")
    
    # Connect to database
    print(f"Connecting to PostgreSQL {db_config.dbname}@{db_config.host}:{db_config.port} as {db_config.user}")
    with psycopg2.connect(
        dbname=db_config.dbname,
        user=db_config.user,
        password=db_config.password,
        host=db_config.host,
        port=db_config.port,
    ) as conn:
        # Initialize database tables if saving to database
        if args.save_to_db:
            print("Initializing database tables (pfasgroups, pfasgroups_in_molecules, components_in_molecules)...")
            initialize_database_tables(conn, pfas_groups)
            print("Database tables ready")
        # Fetch compounds
        compounds = fetch_halogenated_compounds(
            conn,
            halogens=halogens,
            table_name=args.table_name,
            id_column=args.id_column,
            smiles_column=args.smiles_column,
            formula_column=args.formula_column,
            limit=args.limit,
            offset=args.offset,
            use_halogen_columns=args.use_halogen_columns,
        )
        
        total = len(compounds)
        print(f"Compounds to process: {total}")
        if total == 0:
            return 0
        
        all_results = []
        molecule_results_list = []  # For ResultsModel
        compound_id_map = {}  # Map SMILES to molecule_id
        success = 0
        failed = 0
        with_groups = 0
        
        t_run0 = time.perf_counter()
        
        for idx, (compound_id, smiles, formula) in enumerate(compounds, start=1):
            result = parse_one_compound(
                compound_id=compound_id,
                smiles=smiles,
                formula=formula,
                parse_mol_func=parse_mol_func,
                substitute=args.substitute_halogens,
            )
            all_results.append(result)
            
            if result['status'] == 'success':
                success += 1
                if result['result'] and result['result'].get('matches'):
                    with_groups += 1
                    molecule_results_list.append(result['result'])
                    compound_id_map[smiles] = compound_id
            else:
                failed += 1
            
            if idx % args.batch_size == 0:
                print(f"Processed {idx}/{total} (success={success}, failed={failed}, with_groups={with_groups})")
        
        duration = time.perf_counter() - t_run0
        
        # Print summary
        print("\nRun complete")
        print(f"  run_id: {args.run_id}")
        print(f"  processed: {total}")
        print(f"  success: {success}")
        print(f"  failed: {failed}")
        print(f"  compounds_with_groups: {with_groups}")
        print(f"  wall_time_seconds: {duration:.3f}")
        
        # Calculate statistics
        total_matches = 0
        if molecule_results_list:
            for mol_res in molecule_results_list:
                total_matches += len(mol_res.get('matches', []))
        
        if with_groups > 0:
            avg_groups = total_matches / with_groups
            print(f"  total_pfas_group_matches: {total_matches}")
            print(f"  avg_matches_per_molecule: {avg_groups:.2f}")
        
        # Save to database if requested
        if args.save_to_db and molecule_results_list:
            print(f"\nSaving results to database using ResultsModel.to_sql()...")
            
            # Create ResultsModel from results
            results_model = ResultsModel(molecule_results_list)
            
            # Build PostgreSQL connection string
            conn_string = f"postgresql://{db_config.user}:{db_config.password}@{db_config.host}:{db_config.port}/{db_config.dbname}"
            
            # Save to database using to_sql method
            results_model.to_sql(
                conn=conn_string,
                components_table="components_in_molecules",
                groups_table="pfasgroups_in_molecules",
                if_exists="append"
            )
            
            # Update molecule_id column based on SMILES mapping
            print("Updating molecule_id columns...")
            with conn.cursor() as cur:
                for smiles, mol_id in compound_id_map.items():
                    cur.execute(
                        "UPDATE pfasgroups_in_molecules SET molecule_id = %s WHERE smiles = %s AND molecule_id IS NULL",
                        (mol_id, smiles)
                    )
                    cur.execute(
                        "UPDATE components_in_molecules SET molecule_id = %s WHERE smiles = %s AND molecule_id IS NULL",
                        (mol_id, smiles)
                    )
                conn.commit()
            
            print(f"Saved {with_groups} molecules ({total_matches} matches) to database")
        
        # Optionally save to JSON
        if args.output_json:
            import json
            
            output_data = {
                'run_id': args.run_id,
                'timestamp': dt.datetime.utcnow().isoformat(),
                'config': {
                    'database': db_config.dbname,
                    'table': args.table_name,
                    'halogens': halogens,
                    'substitute_halogens': args.substitute_halogens,
                    'total_compounds': total,
                },
                'summary': {
                    'success': success,
                    'failed': failed,
                    'with_groups': with_groups,
                    'duration_seconds': duration,
                },
                'results': [
                    {
                        'compound_id': r['compound_id'],
                        'smiles': r['smiles'],
                        'status': r['status'],
                        'parse_time': r['parse_time'],
                        'num_matches': len(r['result'].get('matches', [])) if r.get('result') else 0,
                    }
                    for r in all_results
                ],
            }
            
            with open(args.output_json, 'w') as f:
                json.dump(output_data, f, indent=2)
            
            print(f"\nResults saved to: {args.output_json}")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
