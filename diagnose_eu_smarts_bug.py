#!/usr/bin/env python3

from rdkit import Chem

mol = Chem.MolFromSmiles('FC(F)(F)N')
print('Molecule:', Chem.MolToSmiles(mol), '(trifluoromethylamine)')
print('='*60)

# Test the ACTUAL EU pattern as written in the JSON
eu_pattern_cf3 = "[#6X4!$([#6][#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]);#6X4!$([#6][#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]);#6X4!$([#6H1]);#6X4!$([#6][#17,#35,#53])](F)(F)F"

print('\nPattern 1 (EU CF3 - as written):')
try:
    patt = Chem.MolFromSmarts(eu_pattern_cf3)
    if patt:
        result = mol.HasSubstructMatch(patt)
        print(f'Result: {result}')
        if result:
            print('  ^^ BUG: Should be False!')
    else:
        print('FAILED TO COMPILE')
except Exception as e:
    print(f'ERROR: {e}')

# Now test with CORRECT SMARTS syntax (no repeating #6X4)
corrected_pattern = "[#6X4!$([#6][#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])!$([#6][#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])!$([#6H1])!$([#6][#17,#35,#53])](F)(F)F"

print('\nPattern 2 (Corrected - no repeated #6X4):')
try:
    patt = Chem.MolFromSmarts(corrected_pattern)
    if patt:
        result = mol.HasSubstructMatch(patt)
        print(f'Result: {result}')
        if not result:
            print('  ^^ CORRECT: Returns False as expected!')
    else:
        print('FAILED TO COMPILE')
except Exception as e:
    print(f'ERROR: {e}')

print('\n' + '='*60)
print('DIAGNOSIS:')
print('='*60)
print('The EU SMARTS uses semicolons (;) to chain conditions:')
print('  [#6X4!$(...);#6X4!$(...);#6X4!$(...)...](F)(F)F')
print('')
print('In SMARTS, when chaining with semicolons in brackets,')
print('you should NOT repeat the atom specification.')
print('') 
print('WRONG:  [#6X4!$(...);#6X4!$(...)](F)(F)F')
print('RIGHT:  [#6X4!$(...)!$(...)](F)(F)F')
print('OR:     [#6X4!$(...);!$(...)](F)(F)F')