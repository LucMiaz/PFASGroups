from sqlalchemy import create_engine,text
from getpass import getpass
from pubchem import dtxsid_to_mols
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
user = input("Enter username for database connection: ")
password = getpass(f"Enter password for database connection for {user}: ")
engine = create_engine(f'mariadb+mariadbconnector://{user}:{password}@localhost/toxcast')

def add_smiles(batch_size=10):
    # Fetch all chid/dtxsid pairs upfront so OFFSET drift can't happen.
    with engine.connect() as conn:
        rows = conn.execute(
            text("SELECT chid, dsstox_substance_id FROM chemical WHERE smiles IS NULL")
        ).fetchall()

    for i in tqdm(range(0, len(rows), batch_size)):
        batch = rows[i : i + batch_size]
        dtxsid_to_chid = {row[1]: row[0] for row in batch}
        dtxsids = list(dtxsid_to_chid.keys())

        mols_dict, failed = dtxsid_to_mols(dtxsids)
        if failed:
            tqdm.write(f"Failed: {failed}")

        with engine.connect() as conn:
            for dtxsid, mol_list in mols_dict.items():
                if not mol_list:
                    continue
                mol = mol_list[0]  # max_cids=1 so there is exactly one
                if mol is None:
                    continue
                conn.execute(
                    text(
                        "UPDATE chemical SET smiles=:smiles, inchi=:inchi,"
                        " inchikey=:inchikey, formula=:formula"
                        " WHERE dsstox_substance_id=:dtxsid"
                    ),
                    {
                        "smiles":    Chem.MolToSmiles(mol),
                        "inchi":     Chem.MolToInchi(mol),
                        "inchikey":  Chem.MolToInchiKey(mol),
                        "formula":   CalcMolFormula(mol),
                        "dtxsid":    dtxsid,
                    },
                )
            conn.commit()  # persist the batch

if __name__ == "__main__":
    add_smiles()