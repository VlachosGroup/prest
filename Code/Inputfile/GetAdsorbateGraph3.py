import os
import networkx as nx
from RuleGeneration.Output import adsorbate
from RuleGeneration.SmilesRelated import FindAllAdsorbateGraphsOfLengthN
from RuleGeneration.ReactionRule import *
from rdkit import Chem

import numpy as np
from rdkit.Chem.rdchem import RWMol

#Process published reactions
f = open('ReactionRules.txt','w')
#f.write('Reaction Rules Generated\n')
f.close()
file = open("PublishedReactions.txt","r")
reactioncount = 0
InputFolderPath = '.\\AllFiles\\'
SurfaceAtomSymbols = ['Pt']
LatticePattern = ['[{Pt3}]12[{Pt4}][{Pt2}]1[{Pt1}]2']
for line in file:
    reactant_molecules = []
    mole = 0
    adj_mat = 0
    mol_formula = 0
    reactioncount = reactioncount + 1
    reaction = line.split(" >> ")
    reactantside = reaction[0]
    productside = reaction[1]
    reactants = reactantside.split(" + ")
    products = productside.split(" + ")
    print(len(reactants), len(products))
    outputfile = open('ReactionRules.txt','a')
    print('rule R'+str(reactioncount)+'{', file=outputfile)
    for i in range(len(reactants)):
        if i <= 1:
            print('\treactant r'+str(i+1)+'{', file=outputfile)
            outputfile.close()
            filepath = InputFolderPath + reactants[i] + '\\CONTCAR'
            mole, adj_mat, mol_formula = adsorbate.LoadByCovalentRadius(filepath, SurfaceAtomSymbols = SurfaceAtomSymbols)
#            for atoms in mole.RdkitMol.GetAtoms():
#                print(atoms.GetProp('Label'))
            reactant_molecules.append(mole)
            outputfile = open('ReactionRules.txt','a')
            print('Hydrogens:', np.int(mole.NumofH))
            print(mole.OAI.__len__())
            if mole.OAI.__len__() == 0:
                print('\t\t'+SurfaceAtomSymbols[0]+'1 labeled p5 {! connected to H with any bond}', file=outputfile)
                print('\t\t'+SurfaceAtomSymbols[0]+'2 labeled p6 single bond to p5', file=outputfile)
                print('\t\t'+SurfaceAtomSymbols[0]+'3 labeled p7 single bond to p6', file=outputfile)
                print('\t\t'+SurfaceAtomSymbols[0]+'4 labeled p8 single bond to p7', file=outputfile)
            print (FindAllAdsorbateGraphsOfLengthN(mole.RdkitMol,SurfaceAtomSymbols,1))
            print('\t}', file=outputfile)
        else: #generally the third reactant is a duplicate and generally an empty site
            print('\treactant r'+str(i+1)+' duplicates r'+str(i)+'(p5 => q1, p6 => q2, p7 => q3, p8 => q4)', file=outputfile)
            
    
    
    
    print('\tconstraints{', file=outputfile)
    for i in range(len(reactant_molecules)):
        NumofHeavyAtoms = reactant_molecules[i].NumofC + reactant_molecules[i].NumofO + reactant_molecules[i].NumofN
        print('\tr'+str(i+1)+'.size = '+str(NumofHeavyAtoms), file=outputfile),
    print('\t}', file=outputfile)
    
    print('}\n', file=outputfile)
    outputfile.close()


#AllMolecules = []
#MoleculeNames = []
#for datum in os.listdir(InputFolderPath):
#    #print(datum)
#    #datum does not capture the molecule correctly
#    #better to save molecule name based on connectivity and SMILES string
#    MoleculeNames.append(datum)
#    fpath = InputFolderPath + datum + '\\CONTCAR'
#    mole, adj_mat, mol_formula = adsorbate.LoadByCovalentRadius(fpath, SurfaceAtomSymbols = SurfaceAtomSymbols)
#    AllMolecules.append(mol_formula)
#    print FindAllAdsorbateGraphsOfLengthN(mole.RdkitMol,SurfaceAtomSymbols,1)
#    print(adj_mat)
#    
################################ USER INPUT ###################################
#InputPath = '..\\Example\\Input\\'

################################ USER INPUT ###################################
# Get SMILES
#for datum in os.listdir(InputPath):
#    fpath = InputPath + datum + '\\CONTCAR'
#    # Convert atomic coordinates (VASP) to molecular graph.
#    mole, adj_mat = adsorbate.LoadByCovalentRadius(fpath, SurfaceAtomSymbols = SurfaceAtomSymbols)
#    # Mine adsorbate graphs and print
#    print FindAllAdsorbateGraphsOfLengthN(mole.RdkitMol,SurfaceAtomSymbols,1)
#    print(adj_mat)

#InputPath_R = '.\\REACTANT\\'
#InputPath_P = '.\\PRODUCT\\'
#
#mole_r = ' '
#mole_p = ' '
#reactant_atoms_count = 0
#product_atoms_count = 0
#
#for datum_r in os.listdir(InputPath_R):
#    fpath = InputPath_R + datum_r + '\\CONTCAR'
#    mole_r, adj_mat_r, mol_formula_r = adsorbate.LoadByCovalentRadius(fpath, SurfaceAtomSymbols = SurfaceAtomSymbols)
#    print (FindAllAdsorbateGraphsOfLengthN(mole_r.RdkitMol,SurfaceAtomSymbols,1))
##    print(adj_mat_r)
##    print(adj_mat_r.shape[0])
#    reactant_atoms_count += adj_mat_r.shape[0]
#    
#print('Done with reactant')
#print(reactant_atoms_count)
#    
#for datum_p in os.listdir(InputPath_P):
#    fpath = InputPath_P + datum_p + '\\CONTCAR'
#    mole_p, adj_mat_p, mol_formula_p = adsorbate.LoadByCovalentRadius(fpath, SurfaceAtomSymbols = SurfaceAtomSymbols)
#    print (FindAllAdsorbateGraphsOfLengthN(mole_p.RdkitMol,SurfaceAtomSymbols,1))
##    print(adj_mat_p, adj_mat_p.shape[0])
#    product_atoms_count += adj_mat_p.shape[0]
#    
#print('Done with product')
#print(product_atoms_count)

#Specify reaction manually


#m = Chem.MolFromSmiles('CCC')
#patt1 = Chem.MolFromSmiles('CC')
#patt2 = Chem.MolFromSmiles('CO')
#print(m.HasSubstructMatch(patt1))
#print(m.HasSubstructMatch(patt2))

#m = Chem.MolFromSmiles('[H2]=[C]1[Pt]~[Pt]1')
#p = Chem.MolFromSmiles('[H][C]12[Pt]3~[Pt]1~[Pt]~32')
#m = Chem.MolFromSmiles('[C]')
#p = Chem.MolFromSmiles('[C]')
#print(m.HasSubstructMatch(p))


