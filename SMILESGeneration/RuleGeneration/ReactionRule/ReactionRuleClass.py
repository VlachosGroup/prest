# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 15:38:31 2019

@author: ugupta
"""
from rdkit import Chem
from RuleGeneration.Output import adsorbate
import json

reactionrule_dict = []
molecule_type = []


class reactionrule(object):
    
    def __init__(self,system,surface,active_site,reaction,reactive_fragment, \
                 elements,molecule_size_low,molecule_size_high,molecule_type, \
                 reactants,constraints,transformations,RING_rule):
        self.system = system
        self.surface = surface
        self.active_site = active_site
        self.reaction = reaction
        self.reactive_fragment = reactive_fragment
        self.elements = elements
        self.molecule_size_low = molecule_size_low
        self.molecule_size_high = molecule_size_high
        self.molecule_type = molecule_type
        self.reactants = reactants
        self.constraints = constraints
        self.transformations = transformations
        self.RING_rule = RING_rule
      
    @classmethod
    def GenerateReactionRule(cls, reactant, product,SurfaceAtomSymbols):
        #assert isinstance(reactant,__builtins__.list)
        #assert isinstance(product,__builtins__.list)
        
        reactionrule_one = dict()
#        reactionrule_dict = dict()
        if isinstance(SurfaceAtomSymbols,str):
            SurfaceAtomSymbols = [SurfaceAtomSymbols]
        else:
            assert isinstance(SurfaceAtomSymbols,list)
        
        if len(SurfaceAtomSymbols) > 0:
            reactionrule_one["system"] = "heterogeneous"
        else:
            reactionrule_one["system"] = "homogeneous"
            reactionrule_one['surface'] = 'none'
            reactionrule_one['active_site'] = 'H+'
        
        if SurfaceAtomSymbols[0] == 'Pt':
            reactionrule_one['surface'] = 'metal'
            reactionrule_one['active_site'] = 'Pt'
        else:
            reactionrule_one['surface'] = 'acid'
            reactionrule_one['active_site'] = 'bronsted'
            
        if len(reactant) == 1:
            reactionrule_one['reaction'] = 'unimolecular'
        elif len(reactant) == 2:
            reactionrule_one['reaction'] = 'bimolecular'
        else:
            reactionrule_one['reaction'] = 'trimolecular'

        nC = 0
        nH = 0
        nO = 0
        nN = 0
        nPt = 0
        elements = []
#        for i in range(len(reactant)):
        s1 = set(elements)
        reac = reactant[0]
        assert isinstance(reactant[0],adsorbate)
        if len(reac.OAI) > 0:
            #Main species to look at
            nC = reac.NumofC
            nH = reac.NumofH
            nO = reac.NumofO
            nN = reac.NumofN
            nPt = reac.NumofPt
            if nH > 0:
                item = 'H'
                if item not in s1:
                    s1.add(item)
                    elements.append(item)
            if nC > 0:
                item = 'C'
                if item not in s1:
                    s1.add(item)
                    elements.append(item)
            if nO > 0:
                item = 'O'
                if item not in s1:
                    s1.add(item)
                    elements.append(item)
            if nN > 0:
                item = 'N'
                if item not in s1:
                    s1.add(item)
                    elements.append(item)
            if nPt > 0:
                item = 'Pt'
                if item not in s1:
                    s1.add(item)
                    elements.append(item)
                    
            reac_smiles = Chem.MolToSmiles(reac.RdkitMol)
            #reac_mol = Chem.MolFromSmiles(reac_smiles,ps)
        
        TotalNHAtoms = nC + nO + nN + nPt
        Adsorbed = False
        if nPt > 0:
            Adsorbed = True
        seen = set(molecule_type)
        
        if Adsorbed == False:
            nH_temp = 2*nC + 2
            if nH == (2*nC + 2):
                item = 'paraffin'
                if item not in seen:
                    seen.add(item)
                    molecule_type.append(item)
#                    molecule_type.append('paraffin')
            elif nH == 2*nC or nH == (2*nC - 2):
                item = 'olefin'
                if item not in seen:
                    seen.add(item)
                    molecule_type.append(item)
#                    molecule_type.append('olefin')
            elif nH == (2*nC - 2) and nC >= 6:
                molecule_type.append('aromatic')
        
        reactionrule_one['molecule_type'] = molecule_type           
            
        TotalNHAtoms = nC + nO + nN + nPt
        reactionrule_one['molecule_size_low'] = TotalNHAtoms
        reactionrule_one['molecule_size_high'] = TotalNHAtoms + 1
        reactionrule_one['elements'] = elements
        reactionrule_one['reactive_fragment'] = 'C-H'
        reactionrule_one["reactants"] = [[{"C":{"c1":"main"}},{"C":{"c2":"single"}},{"H":{"h1":"single"}}],[{"Pt":{"p1":"main"}}],[{"Pt":{"q1":"duplicate"}}]]
        reactionrule_one["constraints"] = [{"size":2},{"size":4}]
        reactionrule_one["transformations"] = [{"break":[["c1","h1"]]},{"form":[["c1","p1"],["h1","q1"]]}]
        
        
        
        
        reactionrule_dict.append(reactionrule_one)
        
        print('In here!')
        
        return True

    def printdict():
        with open("DatabaseEntries.json", "w") as write_file:
            write_file.write('[')
#            entry = ','.join(map(str,reactionrule_dict))
#            write_file.write(entry)
            for entry in reactionrule_dict[:-1]:
                json.dump(entry, write_file)
                write_file.write(',\n')
            json.dump(reactionrule_dict[-1], write_file)
            write_file.write(']')
        