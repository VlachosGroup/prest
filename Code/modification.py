#!/usr/bin/env python

import openbabel, pybel
from pybel import *
from rdkit import Chem
from  xyz2mol import *
from chemml.chem import Molecule

#smiles = 'COc1cc(CCC)cc(C)c1OC(OC)C(O)c1cc(OC)c(O)c(OC)c1'
#smiles = 'c1(cc(c(c(OC)c1)OC([CH]O)O)OC)CC(CO)O'
#smiles = 'c1(cc(c(c(OC)c1)OC(cO)O)OC)CC(CO)O'
#smiles = 'COc1cc(CC(O)CO)cc(OC)c1OC(O)[CH-]O'
#smiles = 'COc1cc(CC(O)CO)cc(OC)c1OC(O)CO'
#smiles = '[{Pt}]124[{Pt}]35C1(C23[{Pt}]45)C'

smiles = '[pt]124[pt]35c1(c23[pt]45)C'

mol = readstring("smi", smiles)
mol.make3D(forcefield='mmff94',steps=1)
mol.write(format='xyz', filename='babel.xyz', overwrite=True)


def smiles_convert(smiles):
    mol = readstring("smi", smiles)
    mol.make3D(forcefield='mmff94',steps=100)
    mol.write(format='xyz', filename='babel.xyz', overwrite=True)
    charged_fragments = True
    quick = True
    atomicNumList, charge, xyz_coordinates = read_xyz_file('babel.xyz')
    mol = xyz2mol(atomicNumList, charge, xyz_coordinates, charged_fragments, quick)
    # Canonical hack
    new_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
    m = Chem.MolFromSmiles(new_smiles)
    new_smiles = Chem.MolToSmiles(m, isomericSmiles=False)
    print(new_smiles)
    return new_smiles

def get_xyz_from_chemml(smiles):
    
    mol = Molecule(smiles, input_type='smiles')
    mol.hydrogens('add')
    mol.to_xyz(optimizer='MMFF', mmffVariant='MMFF94s', maxIters=300)
    
    xyz = mol.xyz.geometry
    symbols=mol.xyz.atomic_symbols
    N_num=len(symbols)
    
    xyz_out = open('chemml.xyz', 'w')
    xyz_out.write(str(N_num)+'\n')
    xyz_out.write(mol.smiles + '\n')
    for num, ele in enumerate(xyz):
        xyz_out.write('%s\t%s\n' %(symbols[num][0], ' '.join([str(i) for i in ele])))
    
    xyz_out.close()

get_xyz_from_chemml(smiles)
