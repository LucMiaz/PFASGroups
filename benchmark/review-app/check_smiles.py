import sqlite3

conn = sqlite3.connect('database/pfas_benchmark.db')
cursor = conn.cursor()

# Get samples
cursor.execute('SELECT smiles FROM molecules LIMIT 20')
samples = [row[0] for row in cursor.fetchall()]

print("Sample SMILES:")
for i, s in enumerate(samples[:10], 1):
    print(f"{i}. {s}")

# Check for spaces
cursor.execute("SELECT smiles FROM molecules WHERE instr(smiles, ' ') > 0 LIMIT 5")
with_spaces = cursor.fetchall()

if with_spaces:
    print("\n\nSMILES with spaces/annotations:")
    for row in with_spaces:
        print(f"  '{row[0]}'")
else:
    print("\n\nNo SMILES with spaces found")

conn.close()
