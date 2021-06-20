import openbabel, pybel
from pybel import *
from rdkit import Chem
from xyz2mol import *

## Use babel to read Smiles, and convert it to xyz, then convert it to SMILEs again
#smiles = 'C(OC1=C(OC)C=C(C(O)C)C=[C.]1)(O)C'
##smiles="c1(C(C(CO)O)O)cc(c(c(OC)c1)O)-c2c(c(OC)cc(C(C(CO)O)O)c2)O"

def smiles_convert(smiles):
    mol = readstring("smi", smiles)
    mol.make3D(forcefield='mmff94',steps=100)
    mol.write(format='xyz', filename='test.xyz', overwrite=True)
   
    charged_fragments = True
    quick = True
    atomicNumList, charge, xyz_coordinates = read_xyz_file('test.xyz')
    mol = xyz2mol(atomicNumList, charge, xyz_coordinates, charged_fragments, quick)
    # Canonical hack
    new_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
    m = Chem.MolFromSmiles(new_smiles)
    new_smiles = Chem.MolToSmiles(m, isomericSmiles=False)
   
    return new_smiles