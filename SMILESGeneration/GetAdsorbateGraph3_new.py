import os
import networkx as nx
from RuleGeneration.Output import adsorbate
#from RuleGeneration.Output import adsorbate
from RuleGeneration.ReactionRule import reactionrule
from RuleGeneration.SmilesRelated import FindAllAdsorbateGraphsOfLengthN
from RuleGeneration.ReactionRule import *
from rdkit import Chem
#from rdkit.Chem import MCS
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog

import numpy as np
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem import rdFMCS

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

try:
  basestring
except NameError:
  basestring = str

#Generating files that will be used as outputs
f = open('Inputfile/AllFiles/ReactionRules.txt','w')
db = open('Inputfile/DatabaseEntries.json','w')
sf = open('Inputfile/SmilesForReactions.txt','w')
sf.close()
db.close()
f.close()

#Input file for generating reaction mechanism
file = open("Inputfile/PublishedReactions.txt","r")
reactioncount = 0
#Specify path to DFT files
InputFolderPath = '.\\Inputfile\\AllFiles\\'
#assume metal symbol in case of metal catalysis
SurfaceAtomSymbols = ['Pt']
#Assume a lattice pattern already
#LatticePattern = ['[{Pt3}]12[{Pt4}][{Pt2}]1[{Pt1}]2']
LatticePattern = ['[{Pt}]1[{Pt}][{Pt}]1']

#Assume each line is a reaction
for line in file:
    
    outputfile = open('Inputfile/AllFiles/ReactionRules.txt','a')
    smilesfile = open('Inputfile/SmilesForReactions.txt','a')
    # Considering each line in Published reactions is an elementary reaction
    # and is atom balanced.
    print('Rule R'+str(reactioncount+1)+'{',file=outputfile)
#    outputjsonfile = open('DatabaseEntries.json','a')
#    print('{', end="", flush=True, file=outputjsonfile)
    reactant_molecules = list()
    product_molecules = list()
    adsorbate_molecules_reac = list()
    adsorbate_molecules_prod = list()
    mole = 0
    adj_mat = 0
    mol_formula = 0
    reactioncount = reactioncount + 1
    #Split reaction into reactant and products
    reaction = line.split(" >> ")
    reactantside = reaction[0]
    productside = reaction[1]
    #Split reactants 
    reactants = reactantside.split(" + ")
    productside_trim = productside.split("\n")
    productside_trim = productside_trim[0]
    #Split products
    products = productside_trim.split(" + ")
    #print(len(reactants), len(products))
    #print(reactants)
    #print(products)
    outputfile.close()
    smilesfile.close()
    #reading in the contcar file again and again
    #should read species only once
#    print('rule R'+str(reactioncount)+'{', file=outputfile)
    for i in range(len(reactants)):
        filepath = InputFolderPath + reactants[i] + '\\CONTCAR'
        #print(str(filepath))
        mole, adj_mat, mol_formula = adsorbate.LoadByCovalentRadius(str(filepath), SurfaceAtomSymbols = SurfaceAtomSymbols, Reactant_mol=True, mol_no=i+1)
#        if i <= 1:
#            print('\treactant r'+str(i+1)+'{', file=outputfile)
#            outputfile.close()
#            filepath = InputFolderPath + reactants[i] + '\\CONTCAR'
#            for atoms in mole.RdkitMol.GetAtoms():
#                print(atoms.GetProp('Label'))
        print('Molecular formula', mol_formula)
        reactant_molecules.append(mole.RdkitMol)
        adsorbate_molecules_reac.append(mole)
        assert isinstance(mole,adsorbate)
#            outputfile = open('AllFiles/ReactionRules.txt','a')
#        print('Hydrogens:', np.int(mole.NumofH))
#        print(mole.OAI.__len__())
#        if mole.OAI.__len__() == 0:
#                print('\t\t'+SurfaceAtomSymbols[0]+'1 labeled p5 {! connected to H with any bond}', file=outputfile)
#                print('\t\t'+SurfaceAtomSymbols[0]+'2 labeled p6 single bond to p5', file=outputfile)
#                print('\t\t'+SurfaceAtomSymbols[0]+'3 labeled p7 single bond to p6', file=outputfile)
#                print('\t\t'+SurfaceAtomSymbols[0]+'4 labeled p8 single bond to p7', file=outputfile)
#        print (FindAllAdsorbateGraphsOfLengthN(mole.RdkitMol,SurfaceAtomSymbols,1))
#            print('\t}', file=outputfile)
#    else: #generally the third reactant is a duplicate and generally an empty site
#            print('\treactant r'+str(i+1)+' duplicates r'+str(i)+'(p5 => q1, p6 => q2, p7 => q3, p8 => q4)', file=outputfile)
#     for i in range(len(products)):
#    print('\tconstraints{', file=outputfile)
        
    outputfile = open('Inputfile/AllFiles/ReactionRules.txt','a')
    smilesfile = open('Inputfile/SmilesForReactions.txt','a')
    reactantstring = list()
    for i in range(len(reactant_molecules)):
        assert isinstance(reactant_molecules[i],Chem.Mol)
        assert isinstance(adsorbate_molecules_reac[i],adsorbate)
        NumofHeavyAtoms = adsorbate_molecules_reac[i].NumofC + adsorbate_molecules_reac[i].NumofO + adsorbate_molecules_reac[i].NumofN
        NumofNonHeavyAtoms = adsorbate_molecules_reac[i].NumofH
        NumnofPtAtoms = adsorbate_molecules_reac[i].NumofPt
        print('Heavy: ',NumofHeavyAtoms,' NonHeavy: ',NumofNonHeavyAtoms,' CatalystAtoms: ',NumnofPtAtoms)
        smiles = Chem.MolToSmiles(reactant_molecules[i])
        print('Reactant here ',smiles)
        reactantstring.append(smiles)
        #print('\tr'+str(i+1)+'.size = '+str(NumofHeavyAtoms+NumnofPtAtoms), file=outputfile),
    
    #Combining reactant side smile strings
    s = '.'
    reactantstring = s.join(reactantstring)
    print('\t}', file=outputfile)
        
    print('\tconstraints{',file=outputfile)
    
    print('\t}',file=outputfile)
    

    for i in range(len(products)):
        filepath = InputFolderPath + products[i] + '\\CONTCAR'
        print(filepath)
        mole, adj_mat, mol_formula = adsorbate.LoadByCovalentRadius(filepath, SurfaceAtomSymbols = SurfaceAtomSymbols,          Reactant_mol=False, mol_no=i+1)
        product_molecules.append(mole.RdkitMol)
        adsorbate_molecules_prod.append(mole)
        
    productstring = list()
    for i in range(len(product_molecules)):
        assert isinstance(product_molecules[i],Chem.Mol)
        assert isinstance(adsorbate_molecules_prod[i],adsorbate)
        NumofHeavyAtoms = adsorbate_molecules_prod[i].NumofC + adsorbate_molecules_prod[i].NumofO + adsorbate_molecules_prod[i].NumofN
        NumofNonHeavyAtoms = adsorbate_molecules_prod[i].NumofH
        NumnofPtAtoms = adsorbate_molecules_prod[i].NumofPt
        print('H: ',NumofHeavyAtoms,' NH: ',NumofNonHeavyAtoms,' Pt: ',NumnofPtAtoms)
        smiles = Chem.MolToSmiles(product_molecules[i])
        print(smiles)
        productstring.append(smiles)
        #print('\tr'+str(i+1)+'.size = '+str(NumofHeavyAtoms+NumnofPtAtoms), file=outputfile),
    
    #Combining product side smile strings
    s = '.'
    productstring = s.join(productstring)
    
    #Combining reaction smile strings
    reactionstring = list()
    s= '>>'
    reactionstring.append(reactantstring)
    reactionstring.append(productstring)
    reactionstring = s.join(reactionstring)
    print(reactionstring, file=smilesfile)
    
    RRreturn = reactionrule.GenerateReactionRule(adsorbate_molecules_reac,adsorbate_molecules_prod,SurfaceAtomSymbols)

    reac_mol = Chem.Mol()
    reac_mol = Chem.RWMol(reac_mol)
    prod_mol = Chem.Mol()
    prod_mol = Chem.RWMol(prod_mol)
#    reac_mol = 0
#    prod_mol = 0
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    ps.sanitize = False
    
    print('Out Here')
    
    for i in range(len(adsorbate_molecules_reac)):
        reac = adsorbate_molecules_reac[i]
        assert isinstance(adsorbate_molecules_reac[i],adsorbate)
        if len(reac.OAI) > 0:
            #Main species to look at
            print('reactant here')
            reac_smiles = Chem.MolToSmiles(reac.RdkitMol)
            reac_mol = Chem.MolFromSmiles(reac_smiles,ps)
            break
            
    for i in range(len(adsorbate_molecules_prod)):
        prod = adsorbate_molecules_prod[i]
        assert isinstance(adsorbate_molecules_prod[i],adsorbate)
        if len(prod.OAI) > 0:
            print('product here')
            #Main species to look at
            prod_smiles = Chem.MolToSmiles(prod.RdkitMol)
            prod_mol = Chem.MolFromSmiles(prod_smiles,ps)
            break
            
    print(Chem.MolToSmiles(reac_mol),Chem.MolToSmiles(prod_mol))
    mols = [reac_mol,prod_mol]
    maxcomsub = rdFMCS.FindMCS(mols)
    print(maxcomsub.smartsString)
    Chem.MolFromSmarts(maxcomsub.smartsString)
    print('}',file=outputfile)

reactionrule.printdict()
outputfile.close() 
smilesfile.close()

mol1 = Chem.MolFromSmiles('[H][H]',ps)
mol2 = Chem.MolFromSmiles('[H]~[Pt]',ps)
mols = [mol1, mol2]
res1 = rdFMCS.FindMCS(mols) 
res2 = rdFMCS.FindMCS(mols,atomCompare=rdFMCS.AtomCompare.CompareElements)
print(res1.numAtoms,res1.numBonds,res1.smartsString)
print(res2.numAtoms,res2.numBonds,res2.smartsString)


print('here')
m1 = Chem.MolFromSmiles('[H]C([H])([H])[H]',ps)
m2 = Chem.MolFromSmiles('[H][C]([H])([H])~[Pt]',ps)
ms = [m1,m2]
r1 = rdFMCS.FindMCS(ms) 
r2 = rdFMCS.FindMCS(ms,atomCompare=rdFMCS.AtomCompare.CompareElements,bondCompare=rdFMCS.BondCompare.CompareOrder)
m3 = Chem.MolFromSmarts(r1.smartsString)
m4 = Chem.MolFromSmarts(r2.smartsString)
print(r1.numAtoms,r1.numBonds,r1.smartsString)
print(r2.numAtoms,r2.numBonds,r2.smartsString)
bis1 = m1.GetSubstructMatch(m3)
bis2 = m1.GetSubstructMatch(m4)
print('bis1: ',bis1)
print('bis2: ',bis2)

mol1 = Chem.MolFromSmiles("Cc1ccccc1")
mol2 = Chem.MolFromSmiles( "CCc1ccccc1" )
mol3 = Chem.MolFromSmiles( "Oc1ccccc1" )
mol4 = Chem.MolFromSmiles( "COc1ccccc1" )
Draw.MolsToGridImage([mol1,mol2,mol3,mol4])

print('done here')
#mols = [m1, m2] 
##bis = m1.GetSubstructMatches(m2)
##print('bis: ',bis)
##fps = [FingerprintMols.FingerprintMol(x) for x in mols]
##print('sim: ',DataStructs.FingerprintSimilarity(fps[0],fps[1]))
##
###print(DataStructs.FingerprintSimilarity(m1,m2))
##print(m1.HasSubstructMatch(m2))
#
#    
#    print('}\n', file=outputfile)
#    outputfile.close()
    
#    if len(reactants) == 1:
#        print('"reaction":"unimolecular"', file=outputjsonfile)
#    elif len(reactants) == 2:
#        print('"reaction":"bimolecular"', file=outputjsonfile)
#    else:
#        print('"reaction":"trimolecular"', file=outputjsonfile)
#    
#    print('}', file=outputjsonfile)

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


