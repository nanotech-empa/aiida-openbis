from rdkit import Chem
from collections import Counter

cdxml_file = "/home/jovyan/work/aiida-openbis/Python_Scripts/test.cdxml"
mol_supplier = Chem.MolsFromCDXMLFile(cdxml_file)

for mol in mol_supplier:
    smiles = Chem.MolToSmiles(mol, kekuleSmiles = True, canonical = True)
    # print("SMILES: ", smiles)

    # Count the atoms in the molecule
    atom_counts = Counter([atom.GetSymbol() for atom in mol.GetAtoms()])

    # Build the molecular formula
    formula = ""
    for atom_symbol, count in atom_counts.items():
        formula += atom_symbol
        if count > 1:
            formula += str(count)

    print(formula)

# BrC1=CC=C(C2=CC=C(Br)C=C2C3=C4C(C5=CC=CC=C5)=C(C6=CC=CC=C6)C(C7=CC=CC=C7)=C3C8=CC=CC=C8)C4=C1

