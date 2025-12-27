# PFASgroups

A Python module for parsing and visualizing PFAS (Per- and Polyfluoroalkyl Substances) groups.

## Features
- Parse SMILES strings to identify PFAS groups
- Generate PFAS fingerprints for machine learning applications
- Assign and visualize PFAS groupings
- **NEW:** Command-line interface for easy access
- **NEW:** Support for custom PFAS group definitions
- **NEW:** Flexible configuration via files, environment variables, or API
- Includes test cases and example data

## Installation

Clone the repository and install dependencies:

```sh
pip install -e .
```

After installation, the `pfasgroups` command will be available in your terminal.

## Quick Start

### Python API

```python
from PFASgroups import parse_pfas, generate_pfas_fingerprint

# Parse PFAS structures
smiles_list = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]
results = parse_pfas(smiles_list)

# Generate fingerprints
fingerprints, group_info = generate_pfas_fingerprint(smiles_list)
```

### Command Line

```bash
# Parse SMILES strings
pfasgroups parse "C(C(F)(F)F)F" "FC(F)(F)C(F)(F)C(=O)O"

# Generate fingerprints
pfasgroups fingerprint "C(C(F)(F)F)F" --format dict

# List available PFAS groups
pfasgroups list-groups
```

## Custom Configuration

Use custom pathtype definitions and PFAS groups:

```bash
# Via command line
pfasgroups parse --groups-file my_custom_groups.json "C(C(F)(F)F)F"

# Via environment variable
export PFASGROUPS_GROUPS="/path/to/custom_groups.json"

# Via Python API
from PFASgroups import set_default_config
set_default_config(groups_file='custom_groups.json')
```

## Documentation

- **[USER_GUIDE.md](USER_GUIDE.md)** - Complete documentation with examples
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - Quick reference for common tasks

## Usage Examples

See [USER_GUIDE.md](USER_GUIDE.md) for comprehensive examples including:
- Basic PFAS parsing and analysis
- Fingerprint generation for machine learning
- Custom configuration files
- Batch processing
- Integration with pandas and scikit-learn

## Licence
<a rel="license" href="http://creativecommons.org/licenses/by-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nd/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nd/4.0/">Creative Commons Attribution-NoDerivatives 4.0 International License</a>.

Contact me in case you want an exception to the No Derivatives term.

## Acknowledgments
This project is part of the [ZeroPM project](https://zeropm.eu/) (WP2) and has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 101036756. This work was developed at the [Department of Environmental Science](https://aces.su.se) at Stockholm University.<br />


<img alt="EU logo" src="https://zeropm.eu/wp-content/uploads/2021/12/flag_yellow_low.jpg" width=100/>     <a rel='zeropm_web' href="https://zeropm.eu/"/><img alt="zeropm logo" src="https://zeropm.eu/wp-content/uploads/2022/01/ZeroPM-logo.png" width=250 /></a><a rel='zeropm_web' href="https://su.se/"/><img alt="zeropm logo" src="https://eu01web.zoom.us/account/branding/p/5065401a-9915-4baa-9c16-665dcd743470.png" width=200 /></a>

[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
    
