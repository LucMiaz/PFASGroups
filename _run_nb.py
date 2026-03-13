import nbformat
from nbclient import NotebookClient
from pathlib import Path

nb_path = Path(r"c:\Users\luc\git\PFASGroups\benchmark\fingerprint_structure_analysis.ipynb")
nb = nbformat.read(nb_path, as_version=4)
client = NotebookClient(nb, timeout=600, kernel_name="python3", allow_errors=False)
client.execute()
nbformat.write(nb, nb_path)
print("EXECUTION_OK")
code_cells = [c for c in nb.cells if c.get('cell_type') == 'code']
executed = sum(1 for c in code_cells if c.get('execution_count') is not None)
print(f"Executed code cells: {executed}/{len(code_cells)}")
