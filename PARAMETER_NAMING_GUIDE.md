# Parameter Naming Guide - December 27, 2025

## Improved Parameter Names

The return values from parsing functions now use more descriptive names:

### Old vs New Parameter Names

| Old Name | New Name | Type | Description |
|----------|----------|------|-------------|
| `n_matches` | `match_count` | int | Number of times the group pattern was matched |
| `n_CFchain` | `chain_lengths` | list[int] | List of carbon-fluorine chain lengths found |
| `chains` | `matched_chains` | list[dict] | Detailed information about each matched chain |

### Return Value Structure

All parsing functions (`parse_smiles`, `parse_mol`, `parse_groups_in_mol`) return tuples with this structure:

```python
(PFASGroup, match_count, chain_lengths, matched_chains)
```

**Example:**
```python
from PFASgroups import parse_smiles

result = parse_smiles('FC(F)(F)C(F)(F)C(=O)O')

for group, match_count, chain_lengths, matched_chains in result:
    print(f"Group: {group.name}")
    print(f"  Matched {match_count} times")
    print(f"  Chain lengths: {chain_lengths}")
    print(f"  Number of chains: {len(matched_chains)}")
```

### CSV Output Format

When using `--format csv` with the CLI, the columns are:

| Column Name | Description |
|-------------|-------------|
| `smiles` | Input SMILES string |
| `group_id` | PFAS group ID number |
| `group_name` | Human-readable group name |
| `match_count` | Number of pattern matches |
| `chain_lengths` | Comma-separated list of chain lengths |
| `num_chains` | Total number of chains found |

**Example CSV Output:**
```bash
pfasgroups parse "FC(F)(F)C(F)(F)C(=O)O" --format csv

# Output:
smiles,group_id,group_name,match_count,chain_lengths,num_chains
FC(F)(F)C(F)(F)C(=O)O,1,Perfluoroalkyl carboxylic acids,2,2,1
FC(F)(F)C(F)(F)C(=O)O,2,Polyfluoroalkyl carboxylic acid,2,2,1
FC(F)(F)C(F)(F)C(=O)O,33,carboxylic acid,1,2,4
FC(F)(F)C(F)(F)C(=O)O,48,alkane,1,2,4
```

### JSON Output Format

JSON output now uses the new parameter names:

```json
[
  {
    "smiles": "FC(F)(F)C(F)(F)C(=O)O",
    "groups": [
      {
        "name": "Perfluoroalkyl carboxylic acids",
        "id": 1,
        "match_count": 2,
        "chain_lengths": [2],
        "num_chains": 1
      }
    ]
  }
]
```

### CLI Examples

**Parse with JSON output (default):**
```bash
pfasgroups parse "FC(F)(F)C(F)(F)C(=O)O"
pfasgroups parse "FC(F)(F)C(F)(F)C(=O)O" --pretty
```

**Parse with CSV output:**
```bash
pfasgroups parse "FC(F)(F)C(F)(F)C(=O)O" --format csv
pfasgroups parse --input compounds.txt --output results.csv --format csv
```

**Fingerprint with CSV output:**
```bash
pfasgroups fingerprint "FC(F)(F)C(F)(F)C(=O)O" --output-format csv
pfasgroups fingerprint --input compounds.txt -f dict --output-format csv -o fps.csv
```

## Migration from Old Names

If you have existing code using the old parameter names, update as follows:

**Old Code:**
```python
for group, n_matches, n_CFchain, chains in results:
    print(f"Found {n_matches} matches")
    print(f"Chain info: {n_CFchain}")
    print(f"Chains: {len(chains)}")
```

**New Code:**
```python
for group, match_count, chain_lengths, matched_chains in results:
    print(f"Found {match_count} matches")
    print(f"Chain lengths: {chain_lengths}")
    print(f"Chains: {len(matched_chains)}")
```

## Rationale

The new names are more descriptive and follow common naming conventions:

- **`match_count`**: Clearer than `n_matches` - explicitly indicates a count
- **`chain_lengths`**: More descriptive than `n_CFchain` - indicates it's a list of lengths
- **`matched_chains`**: More explicit than `chains` - indicates these are the matched chains

These names improve code readability and make the API more intuitive for new users.
