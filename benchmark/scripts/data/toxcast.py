"""
Populate the `smiles`, `inchi`, `inchikey` and `formula` columns of the
`chemical` table in a local ToxCast MariaDB database.

For each chemical that still lacks a SMILES string, the script looks up its
DTXSID identifier via the PubChem REST API (pubchem.dtxsid_to_mols) and
writes the RDKit-canonicalised structure back to the database.
"""

from sqlalchemy import create_engine, text
from getpass import getpass
from pubchem import dtxsid_to_mols
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# ── Database connection ───────────────────────────────────────────────────────
# Credentials are entered interactively so they are never stored in the script.
user     = input("Enter username for database connection: ")
password = getpass(f"Enter password for database connection for {user}: ")
engine   = create_engine(
    f"mariadb+mariadbconnector://{user}:{password}@localhost/toxcast"
)


def add_smiles(batch_size: int = 10) -> None:
    """Fetch structure data from PubChem and write it to the ToxCast DB.

    All rows in `chemical` where `smiles IS NULL` are resolved in batches of
    `batch_size` to avoid overloading the PubChem API.  Each batch is
    committed as a single transaction so partial progress is preserved if the
    script is interrupted.

    Parameters
    ----------
    batch_size:
        Number of DTXSIDs sent to PubChem per API call.
    """
    # Retrieve every (chid, dtxsid) pair that still needs a structure.
    # Fetching the full list upfront avoids OFFSET drift caused by rows being
    # updated while iterating with LIMIT/OFFSET.
    with engine.connect() as conn:
        rows = conn.execute(
            text("SELECT chid, dsstox_substance_id FROM chemical WHERE smiles IS NULL")
        ).fetchall()

    for i in tqdm(range(0, len(rows), batch_size)):
        batch = rows[i : i + batch_size]

        # Build a mapping so we can look up chid from dtxsid after the API call.
        dtxsid_to_chid = {row[1]: row[0] for row in batch}
        dtxsids = list(dtxsid_to_chid.keys())

        # Query PubChem for RDKit Mol objects keyed by DTXSID.
        # max_cids=1 is assumed, so each value is a one-element list.
        mols_dict, failed = dtxsid_to_mols(dtxsids)
        if failed:
            tqdm.write(f"Failed: {failed}")

        with engine.connect() as conn:
            for dtxsid, mol_list in mols_dict.items():
                if not mol_list:
                    continue
                mol = mol_list[0]  # one mol per DTXSID (max_cids=1)
                if mol is None:
                    continue

                # Write all four structure fields in a single UPDATE statement.
                conn.execute(
                    text(
                        "UPDATE chemical SET smiles=:smiles, inchi=:inchi,"
                        " inchikey=:inchikey, formula=:formula"
                        " WHERE dsstox_substance_id=:dtxsid"
                    ),
                    {
                        "smiles":   Chem.MolToSmiles(mol),
                        "inchi":    Chem.MolToInchi(mol),
                        "inchikey": Chem.MolToInchiKey(mol),
                        "formula":  CalcMolFormula(mol),
                        "dtxsid":   dtxsid,
                    },
                )
            conn.commit()  # commit the whole batch at once


if __name__ == "__main__":
    add_smiles()
