from sqlalchemy import create_engine,text
from getpass import getpass
from pubchem import dtxsid_to_smiles
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
user = input("Enter username for database connection: ")
password = getpass(f"Enter password for database connection for {user}: ")
engine = create_engine(f'mariadb+mariadbconnector://{user}:{password}@localhost/toxcast')

def add_smiles(batch_size=1000):
    with engine.connect() as conn:
        total = conn.execute(text("SELECT COUNT(*) FROM chemical")).scalar()
    with tqdm(total = total) as pbar:
        for i in range(0, total, batch_size):
            with engine.connect() as conn:
                result = conn.execute(text("SELECT chid, dsstox_substance_id FROM chemical LIMIT :limit OFFSET :offset"), {"limit": batch_size, "offset": i}).fetchall()
                dtxsids = [row[1] for row in result]
            smiles_dict = dtxsid_to_smiles(dtxsids, unique=True, join_smiles=False)
            with engine.connect() as conn:
                for dtxsid, smiles in smiles_dict.items():
                    mol = Chem.MolFromSmiles(smiles)
                    inchi = Chem.MolToInchi(mol) if mol is not None else None
                    inchikey = Chem.MolToInchiKey(mol) if mol is not None else None
                    formula = CalcMolFormula(mol) if mol is not None else None
                    conn.execute(text("UPDATE chemical SET smiles = :smiles, inchi = :inchi, inchikey = :inchikey, formula = :formula WHERE dsstox_substance_id = :dtxsid"), {"smiles": smiles, "inchi": inchi, "inchikey": inchikey, "formula": formula, "dtxsid": dtxsid})
            pbar.update(batch_size)

if __name__ == "__main__":
    add_smiles()