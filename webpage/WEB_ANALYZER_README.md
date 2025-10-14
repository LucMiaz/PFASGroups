# 🧪 PFAS Group Analyzer - Web Application

A standalone HTML/JavaScript tool for identifying PFAS (Per- and Polyfluoroalkyl Substances) groups and analyzing molecular structures directly in your browser.

## ✨ Features

- 🔬 **Identify PFAS Groups** - Automatically detects 15+ different PFAS chemical classes
- 🧬 **Multiple Input Formats** - Supports SMILES and InChI notation
- 🎨 **Visual Structure Display** - Renders 2D molecular structures
- 📊 **Formula Analysis** - Shows molecular formula and atom counts
- 🚀 **No Installation Required** - Runs entirely in the browser
- 🔒 **Privacy-Friendly** - All processing happens locally, no data sent to servers

## 🎯 Quick Start

1. Open `pfas_analyzer.html` in any modern web browser
2. Enter a SMILES or InChI string (or click an example)
3. Click "Analyze" to identify PFAS groups
4. View the results with structure visualization

## 📋 Detected PFAS Groups

The analyzer can identify these PFAS chemical classes:

### Acids
- **PFCAs** - Perfluoroalkyl carboxylic acids
- **PolyFCAs** - Polyfluoroalkyl carboxylic acids
- **PFSAs** - Perfluoroalkyl sulfonic acids
- **PolyFSAs** - Polyfluoroalkyl sulfonic acids
- **PFECAs** - Perfluoroalkylether carboxylic acids
- **PFESAs** - Perfluoroalkylether sulfonic acids

### Alcohols & Ethers
- **Perfluoroalkyl alcohols**
- **Fluorotelomer alcohols** (FTOHs)
- **Hydrofluoroethers** (HFEs)

### Other Classes
- **Hydrofluorocarbons** (HFCs)
- **Perfluoroalkanes**
- **Perfluoroalkane sulfonyl fluorides** (PASFs)

## 💡 Example Molecules

### PFOA (Perfluorooctanoic acid)
```
SMILES: C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O
Detects: PFCAs, PolyFCAs
```

### PFOS (Perfluorooctanesulfonic acid)
```
SMILES: C(C(C(C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)(F)F
Detects: PFSAs, PolyFSAs
```

### GenX (HFPO-DA)
```
SMILES: C(C(C(OC(C(F)(F)F)(F)F)(F)F)(F)F)(=O)O
Detects: PFECAs, PolyFECAs
```

### 6:2 FTOH
```
SMILES: C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)(F)CCO
Detects: Fluorotelomer alcohols
```

## 🔧 Technical Details

### Technology Stack
- **RDKit.js** - WebAssembly chemistry library for molecular processing
- **Pure JavaScript** - No frameworks, vanilla JS for simplicity
- **Modern CSS** - Responsive design with gradients and animations

### Algorithm
1. **Parse Input** - Convert SMILES/InChI to molecular structure
2. **Count Atoms** - Extract formula including implicit hydrogens
3. **SMARTS Matching** - Check for functional group patterns
4. **Constraint Validation** - Verify formula-based rules
5. **Display Results** - Show matched groups with visualization

### Formula Constraints

The analyzer uses several types of constraints:

- **rel (relational)**: Element count based on other elements
  - Example: `C = F/2 + 0.5` (carbon count equals half of fluorines plus 0.5)
  
- **gte (greater than or equal)**: Minimum element count
  - Example: `F >= 1` (at least one fluorine)
  
- **only**: Allowed elements
  - Example: `["C", "F", "H"]` (only carbon, fluorine, and hydrogen)

## 🐛 Bug Fixes

### Version 1.1 (Current)
- ✅ **Fixed hydrogen counting** - Now correctly includes implicit hydrogens
- ✅ **Improved PFAS detection** - All groups now properly identified
- ✅ **Added more groups** - Expanded from 8 to 15+ PFAS classes
- ✅ **Better UI** - Added atom count breakdown and debug info

See `BUG_FIX_SUMMARY.md` for details on the hydrogen counting fix.

## 🧪 Testing

### Browser Console
Open browser developer tools (F12) to see:
- Atom counting details
- Constraint validation results
- SMARTS matching information

### Test Files
- `test_rdkit.html` - Simple test for atom counting
- `manual_check.py` - Python script for expected values
- See `TESTING_GUIDE.md` for comprehensive testing instructions

## 📊 Accuracy

The JavaScript implementation has been validated against the Python PFASgroups library:

- ✅ Molecular formula extraction matches Python
- ✅ SMARTS pattern matching consistent with RDKit
- ✅ Constraint checking produces identical results
- ⚠️ Path-finding simplified (full homologue series analysis not implemented)

## 🌐 Browser Compatibility

Works in all modern browsers that support:
- WebAssembly (for RDKit.js)
- ES6 JavaScript
- CSS Grid and Flexbox

Tested on:
- Chrome 90+
- Firefox 88+
- Edge 90+
- Safari 14+

## 📄 Files

- `pfas_analyzer.html` - Main application (standalone, single file)
- `test_rdkit.html` - Testing tool
- `TESTING_GUIDE.md` - Comprehensive testing documentation
- `BUG_FIX_SUMMARY.md` - Details on bug fixes
- `manual_check.py` - Expected results validation

## 🤝 Based On

This web application is based on the [PFASgroups Python library](https://github.com/), adapted to run entirely in the browser using RDKit.js.

## 📝 License

[Same as parent PFASgroups library]

## 🙏 Acknowledgments

- RDKit.js team for the WebAssembly port
- PFASgroups developers for the original Python implementation
- PFAS research community for chemical definitions and constraints

## 💬 Support

For issues or questions:
1. Check the `TESTING_GUIDE.md` for troubleshooting
2. Review console logs (F12) for debugging information
3. Verify input format (SMILES/InChI)

## 🔮 Future Enhancements

Potential additions:
- [ ] Full path-finding for homologue series
- [ ] Chain length calculation
- [ ] Batch processing from file
- [ ] Export results to CSV/JSON
- [ ] Additional PFAS group definitions
- [ ] Structure-activity relationship predictions

---

**Version:** 1.1  
**Last Updated:** 2025  
**Status:** ✅ Bug Fixed - Hydrogen counting corrected
