from ase.io import read
import numpy as np
from rdkit import Chem
from ase import Atoms as ase_Atoms
import networkx as nx

try:
  basestring
except NameError:
  basestring = str

class adsorbate(object):
    """
    This is an adsorbate graph class that converts atomic coordinates to rdkit 
    molecular graph object, Mol. Use "LoadByCovalentRadius" to initialize.
    
    Class Variables
    soan: selected organic atomic number. These atoms are considered adosbates
    rcov: covalent radius. Info available in wikipedia.
    
    Class Attributes
    ASEAtoms:                   ASE Atoms object.
    RdkitMol:                   Rdkit Mol object.
    SurfaceAtomSymbols:         List of symbols of surface atoms.
    ASEAtomIndex2RdKitAtomIndex: Index mapping from ASE atoms to Rdkit Mol
    RdKitAtomIndex2ASEAtomIndex: Index mapping from Rdkit Mol to ASE Atoms.
    """
    # selected organic atomic numbers
    soan = [1,6,7,8] #H, C, N, O 
    # atomic number -> covalent radius
    rcov = {1:0.31, 6:0.76, 7:0.71, 8:0.66, 26:1.26, 27:1.21, 28:1.21, 29:1.21,\
            44:1.16, 45:1.21, 46:1.26, 47:1.46, 75:1.21, 77:1.21 ,78:1.21, 79:1.21}
    
    NumofPt = 0    
    NumofC = 0
    NumofH = 0
    NumofO = 0
    NumofN = 0
    OAN = list()
    OAI = list()
    AdsorbedSpecies = False
    ReactantGraph = nx.Graph()
    
    def __init__(self,ASEAtoms,RdkitMol,SurfaceAtomSymbols, \
                 ASEAtomIndex2RdKitAtomIndex, RdKitAtomIndex2ASEAtomIndex, NumofC, NumofH, NumofO, NumofN, NumofPt):
        
        assert isinstance(ASEAtoms,ase_Atoms)
        assert isinstance(RdkitMol,Chem.Mol)
        assert isinstance(ASEAtomIndex2RdKitAtomIndex,dict)
        assert isinstance(RdKitAtomIndex2ASEAtomIndex,dict)
        if isinstance(SurfaceAtomSymbols,str):
            SurfaceAtomSymbols = [SurfaceAtomSymbols]
        else:
            assert isinstance(SurfaceAtomSymbols,list)
        self.ASEAtoms = ASEAtoms
        self.RdkitMol = RdkitMol
        self.SurfaceAtomSymbols = SurfaceAtomSymbols
        self.ASEAtomIndex2RdKitAtomIndex = ASEAtomIndex2RdKitAtomIndex
        self.RdKitAtomIndex2ASEAtomIndex = RdKitAtomIndex2ASEAtomIndex
        self.NumofC = NumofC
        self.NumofH = NumofH
        self.NumofO = NumofO
        self.NumofN = NumofN
        self.NumofPt = NumofPt
    
    @classmethod
    def LoadByCovalentRadius(cls,CoordinateFPath, SurfaceAtomSymbols, Reactant_mol, mol_no, \
        rfacup = 1.35,rfacdown = 0.6, z_vector = 2):
        """ 
        This function reads file using ASE read, and construts molecular graph
        in rdkit object, Mol. See manuscript for overall algorithm.
        
        
        Input List
        CoordinateFPath:    path to ASE readable coordinate file.
        SurfaceAtomSymbols: List of atomic symbols of surface atoms.
        rfacup:             Upper percentage limit for determining connectivity.
        rfacdown:           Lower percentage limit for determining connectivity.
        z_vector:           index of cell basis vector that is orthogonal to surface.
        
        Output List
        adsorbate class
        """
        #print(CoordinateFPath)
        # initialize
        ASEAtomIndex2RdKitAtomIndex = dict()
        RdKitAtomIndex2ASEAtomIndex = dict()
        Adsorbed = None #to confirm if species is adsorbed or not
        
        #Make sure surface atom symbols be a list
        if isinstance(SurfaceAtomSymbols,str):
            SurfaceAtomSymbols = [SurfaceAtomSymbols]
        else:
            assert isinstance(SurfaceAtomSymbols,list)
        # load POSCAR and AseAtoms is the ase object
        AseAtoms = read(CoordinateFPath)
        # if none given for surface layer z coordinate, average the top layer atomic coordinate
        _, SurfaceAtomIndex = DetermineSurfaceLayerZ(AseAtoms, SurfaceAtomSymbols, ZVecIndex = z_vector)

        # (p)eriodic (b)oundary (c)ondition(s)
        PBCs = [[0,0,0]]
        if AseAtoms.pbc[0]:
            temp = np.add(PBCs,[1,0,0])
            temp = np.concatenate((temp,np.add(PBCs,[-1,0,0])))
            PBCs = np.concatenate((PBCs,temp))
        if AseAtoms.pbc[1]:
            temp = np.add(PBCs,[0,1,0])
            temp = np.concatenate((temp,np.add(PBCs,[0,-1,0])))
            PBCs = np.concatenate((PBCs,temp))
        if AseAtoms.pbc[2]:
            temp = np.add(PBCs,[0,0,1])
            temp = np.concatenate((temp,np.add(PBCs,[0,0,-1])))
            PBCs = np.concatenate((PBCs,temp))
                    
        # Get organic atoms from the DFT calculations (their index and atomic number)
        ans = AseAtoms.get_atomic_numbers() # (a)tomic (n)umber(s)
        
        oai = list() #organic atom index in the atoms object
        oan = list() #organic atomic number
        for i in range(0,AseAtoms.__len__()): #check if len(AseAtoms) works
            if ans[i] in cls.soan: #only organic atoms specified here
                oai.append(i)
                oan.append(ans[i]) #saved atomic index and atomic number
                
        cls.OAI = oai
        cls.OAN = oan
         
        print('Total number of organic atoms:', len(oai))
        
            
        #print(oan)
        #print(oai.__len__())
        # Determine connectivity of the organic atoms
        adj_mat = np.zeros((len(oai),len(oai))) # adjacency matrix
        for i in range(0,len(oai)):
            for j in range(i+1,len(oai)):
                if cls._DetermineConnectivity(AseAtoms,oai[i],oai[j],PBCs,rfacup,rfacdown):
                    adj_mat[i,j] = 1
                    
        print(adj_mat)
        
        mol_formula = np.zeros(4)
        for j in range(0,len(oai)):
            if oan[j] == 6:
                mol_formula[0] += 1
            elif oan[j] == 1:
                mol_formula[1] += 1
            elif oan[j] == 8:
                mol_formula[2] += 1
            elif oan[j] == 7:
                mol_formula[3] += 1
        
        #print("Molecular formula: ", "C", mol_formula[0],"H", mol_formula[1],"O", mol_formula[2], "N", mol_formula[3])
        cls.NumofH = mol_formula[1]

        # construct mol object
        RdkitMol = Chem.Mol()
        RdkitMol = Chem.RWMol(RdkitMol)
        #reviewed till here!!!
        
        outputfile1 = open('Inputfile/AllFiles/ReactionRules.txt','a')
         
        NumofPt = 0    
        NumofC = 0
        NumofH = 0
        NumofO = 0
        NumofN = 0
        firstTime = True
        ## add atom
        ### organic atoms
        if Reactant_mol and len(oai) > 0 :
            print('\treactant r'+str(mol_no)+'{', file = outputfile1)
        for i in range(0,len(oan)):
            #print(oan[i])
            
            atom = Chem.Atom(np.int(oan[i])) #find the atom corresponding to the atomic number
            atom.SetNoImplicit(True) # this allows molecule to have radical atoms
            atom.SetBoolProp('Adsorbed',False)
            if atom.GetAtomicNum() == 6:
                NumofC += 1
                atom.SetProp('Label',atom.GetSymbol().lower() + str(NumofC))
            elif atom.GetAtomicNum() == 1:
                NumofH += 1
                atom.SetProp('Label',atom.GetSymbol().lower() + str(NumofH))
            elif atom.GetAtomicNum() == 8:
                NumofO += 1
                atom.SetProp('Label',atom.GetSymbol().lower() + str(NumofO))
            elif atom.GetAtomicNum() == 7:
                NumofN += 1
                atom.SetProp('Label',atom.GetSymbol().lower() + str(NumofN))
            
            
            #Need something better here
            if Reactant_mol :
                
                if firstTime :
                    print('\t\t',atom.GetSymbol(), 'labeled',atom.GetProp('Label'), file = outputfile1)
                    firstTime = False
                else:
                    for j in range(0,i):
                        if adj_mat[i,j] > 0 or adj_mat[j,i] > 0:
                            print('\t\t',atom.GetSymbol(), 'labeled',atom.GetProp('Label'),'single bond to',RdkitMol.GetAtomWithIdx(j).GetProp('Label'), file = outputfile1)
                            break
            RdkitMol.AddAtom(atom)
            ASEAtomIndex2RdKitAtomIndex[oai[i]] = i
            RdKitAtomIndex2ASEAtomIndex[i] = oai[i]
#            print(atom.GetProp('Label'))
            #print(oai[i], i)
       
            #print(index, i)
        if Reactant_mol and len(oai) > 0 :
            print('\t\t}', file= outputfile1)
        
        ## add bond
        ### between organic atoms
        for i in range(0,oai.__len__()):
            for j in range(i+1,oai.__len__()):
                if adj_mat[i,j] == 1:
                    RdkitMol.AddBond(i,j,order=Chem.rdchem.BondType.SINGLE) #add all bonds in the connections
                    
        ### Dont need to add zero bonds between all surface atoms
        ### between surface atoms
#        for i in range(0,len(SurfaceAtomIndex)):
#            for j in range(i+1,len(SurfaceAtomIndex)):
#                if cls._DetermineConnectivity(AseAtoms,SurfaceAtomIndex[i],SurfaceAtomIndex[j],PBCs,rfacup,rfacdown):
#                    RdkitMol.AddBond(ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[i]],ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[j]],order=Chem.rdchem.BondType.ZERO)
#                    
        ## assign radicals
        Chem.AssignRadicals(RdkitMol)#This takes into account the valency of each atom within RdkitMol
        
        ## set smilesSymbol
        for atom in RdkitMol.GetAtoms():
            if atom.GetSymbol() in ['C','O'] and atom.GetNumRadicalElectrons() == 0:
#                atom.SetProp("smilesSymbol",'[' + atom.GetSymbol() + str(atom.GetNumRadicalElectrons())+ ']')
                atom.SetProp("AdsorbateSiteType","NoBondWithSite")
            elif atom.GetNumRadicalElectrons() > 0:
#                atom.SetProp("smilesSymbol",atom.GetSymbol() + str(atom.GetNumRadicalElectrons()))
                if atom.GetNumRadicalElectrons() == 1:
                    atom.SetProp("AdsorbateSiteType","Atop")
                elif atom.GetNumRadicalElectrons() == 2:
                    atom.SetProp("AdsorbateSiteType","Bridge")
                else:
                    atom.SetProp("AdsorbateSiteType","Hollow")

            
        # Find surface binding atom. This is done by finding all the radical atoms
        rai_rdkit = list() # radical atom index for rdkit mol
        rai_ase = list() # radical atom index for rdkit ase atoms object
        for atom in RdkitMol.GetAtoms():
            if atom.GetNumRadicalElectrons() > 0:
                rai_rdkit.append(atom.GetIdx())
                rai_ase.append(oai[atom.GetIdx()])
                
        lsa_rdkit = list() #linked surface atom index for rdkit mol
        lsa_ase = list() #linked surface atom index for ase mol
        
#         ### surface atoms
#        for index in SurfaceAtomIndex:
#            atom = Chem.Atom(AseAtoms[index].symbol)
#            atom.SetBoolProp('SurfaceAtom',True)
#            atom.SetBoolProp('Occupied',False)
#            i = RdkitMol.AddAtom(atom)
#            ASEAtomIndex2RdKitAtomIndex[index] = i
#            RdKitAtomIndex2ASEAtomIndex[i] = index
        
        #To add a lattice pattern for the empty site
        emptysite_index = list()
        
        if Reactant_mol and len(oai) == 0:
            print('\treactant r'+str(mol_no)+'{', file = outputfile1)
            print('Empty site found')
            sym = SurfaceAtomSymbols[0]
            latticepatternsize = 3 #### This needs to be updated
            while (latticepatternsize > 0):
                atom = Chem.Atom(sym) #pick symbol from surfaceatoms
                atom.SetBoolProp('SurfaceAtom',True)
                atom.SetBoolProp('Occupied',False)
                k = RdkitMol.AddAtom(atom)
                emptysite_index.append(k)
                NumofPt += 1
                latticepatternsize -= 1
            
            for i in range(0,len(emptysite_index)):
                for j in range(i+1,len(emptysite_index)):
                    if (i > 0 and j == i+1) or i == 0:
#                        print(emptysite_index[i],emptysite_index[j])
                        RdkitMol.AddBond(emptysite_index[i],emptysite_index[j],order=Chem.rdchem.BondType.ZERO)
            
            print('\t\t}', file = outputfile1)            
        # Surface connectivity
        for i in range(0,len(rai_ase)):
            for j in range(0,len(SurfaceAtomIndex)):
                if cls._DetermineConnectivity(AseAtoms,rai_ase[i],SurfaceAtomIndex[j],PBCs,rfacup,rfacdown):
                    atom = Chem.Atom(AseAtoms[SurfaceAtomIndex[j]].symbol)
                    NumofPt += 1
                    atom.SetBoolProp('SurfaceAtom',True)
                    atom.SetBoolProp('Occupied',False)
                    k = RdkitMol.AddAtom(atom)
                    ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[j]] = k
                    RdKitAtomIndex2ASEAtomIndex[k] = SurfaceAtomIndex[j]
                    RdkitMol.AddBond(rai_rdkit[i],ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[j]],order=Chem.rdchem.BondType.ZERO)
                    RdkitMol.GetAtomWithIdx(ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[j]]).SetBoolProp('Occupied',True)
                    RdkitMol.GetAtomWithIdx(rai_rdkit[i]).SetBoolProp('Adsorbed',True)
                    lsa_rdkit.append(ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[j]])
                    lsa_ase.append(SurfaceAtomIndex[j])
                    
        for i in range(0,len(lsa_ase)):
            for j in range(i+1,len(lsa_ase)):
                if cls._DetermineConnectivity(AseAtoms,lsa_ase[i],lsa_ase[j],PBCs,rfacup,rfacdown):
                    RdkitMol.AddBond(lsa_rdkit[i],lsa_rdkit[j],order=Chem.rdchem.BondType.ZERO)
                    
        
#        # assign binding site.
#        for i in range(0,len(rai_rdkit)):
#            a = RdkitMol.GetAtomWithIdx(rai_rdkit[i])
#            nsurf = 0
#            for neighbor_atom in a.GetNeighbors(): #this tells what atoms this atom is bonded to
#                if neighbor_atom.GetSymbol() in SurfaceAtomSymbols:
#                    nsurf += 1
#            a.SetProp("smilesSymbol",a.GetProp("smilesSymbol") + '_' + str(nsurf) + 'fold')
                    
        if NumofPt > 0:
            Adsorbed = True
        else:
            Adsorbed = False
                
        if Adsorbed:
            if NumofPt == 1:
                KnownPtIndices = []
                for atom in RdkitMol.GetAtoms():
                    if atom.GetSymbol() == 'Pt':
                        index = atom.GetIdx()
                        KnownPtIndices.append(index)
                index1 = KnownPtIndices[0]
                NumofPt += 2
                atom = Chem.Atom('Pt')
                atom.SetBoolProp('SurfaceAtom',True)
                atom.SetBoolProp('Occupied',False)
                k1 = RdkitMol.AddAtom(atom)
                RdkitMol.AddBond(index1,k1,order=Chem.rdchem.BondType.ZERO)
                atom = Chem.Atom('Pt')
                atom.SetBoolProp('SurfaceAtom',True)
                atom.SetBoolProp('Occupied',False)
                k2 = RdkitMol.AddAtom(atom)
                RdkitMol.AddBond(index1,k2,order=Chem.rdchem.BondType.ZERO)
                RdkitMol.AddBond(k1,k2,order=Chem.rdchem.BondType.ZERO)
                #atom = Chem.Atom('Pt')
                #atom.SetBoolProp('SurfaceAtom',True)
                #atom.SetBoolProp('Occupied',False)
                #k3 = RdkitMol.AddAtom(atom)
                #RdkitMol.AddBond(k2,k3,order=Chem.rdchem.BondType.ZERO)
                #RdkitMol.AddBond(k1,k3,order=Chem.rdchem.BondType.ZERO)                        
            elif NumofPt == 2:
                KnownPtIndices = []
                for atom in RdkitMol.GetAtoms():
                    if atom.GetSymbol() == 'Pt':
                        index = atom.GetIdx()
                        KnownPtIndices.append(index)
                index1 = KnownPtIndices[0]
                index2 = KnownPtIndices[1]
                NumofPt += 1
                atom = Chem.Atom('Pt')
                atom.SetBoolProp('SurfaceAtom',True)
                atom.SetBoolProp('Occupied',False)
                k1 = RdkitMol.AddAtom(atom)
                RdkitMol.AddBond(index1,k1,order=Chem.rdchem.BondType.ZERO)
                RdkitMol.AddBond(index2,k1,order=Chem.rdchem.BondType.ZERO)
                #atom = Chem.Atom('Pt')
                #atom.SetBoolProp('SurfaceAtom',True)
                #atom.SetBoolProp('Occupied',False)
                #k2 = RdkitMol.AddAtom(atom)
                #RdkitMol.AddBond(index2,k2,order=Chem.rdchem.BondType.ZERO)
                #RdkitMol.AddBond(k1,k2,order=Chem.rdchem.BondType.ZERO)                
            elif NumofPt == 3:
                KnownPtIndices = []
                for atom in RdkitMol.GetAtoms():
                    if atom.GetSymbol() == 'Pt':
                        index = atom.GetIdx()
                        KnownPtIndices.append(index)
                index1 = KnownPtIndices[0]
                index2 = KnownPtIndices[1]
                #atom = Chem.Atom('Pt')
                #NumofPt += 0
                #atom.SetBoolProp('SurfaceAtom',True)
                #atom.SetBoolProp('Occupied',False)
                #k = RdkitMol.AddAtom(atom)
                #RdkitMol.AddBond(index1,k,order=Chem.rdchem.BondType.ZERO)
                #RdkitMol.AddBond(index2,k,order=Chem.rdchem.BondType.ZERO)
                
                    
        if NumofPt < 4 and NumofPt > 0:
            print('Species has less than 4 Pt: ',Chem.MolToSmiles(RdkitMol),NumofPt)
                
        outputfile1.close()
        
        adsorbate = cls(AseAtoms,RdkitMol,SurfaceAtomSymbols, \
                 ASEAtomIndex2RdKitAtomIndex, RdKitAtomIndex2ASEAtomIndex, NumofC, NumofH, NumofO, NumofN, NumofPt)
        
#        smi = Chem.MolToSmiles(RdkitMol)
#        print(smi) 
        return adsorbate,adj_mat, mol_formula

    @classmethod
    def _DetermineConnectivity(cls,AseAtoms,i,j,PBCs,rfacup,rfacdown):
        """
        Determine connectivity between atom i and j. See equation (1) in the 
        manuscript.
        
        Input List
        ASEAtoms:           ASE atoms containing adsorbate/surface system
        PBCs:               Periodic Boundary Conditions. e.g., (1,0,0) means 
                            cell repeats in first basis vector but not others.
        rfacup:             upper tolerance factor
        rfacdown:           lower tolerance factor
        
        Output List
        Bool:               True if connected, false if not.
        """
        xyz1 = AseAtoms[i].position
        # loop over periodic cells
        for PBC in PBCs:
            xyz2 = AseAtoms[j].position + np.dot(PBC,AseAtoms.cell)
            # Criteria:
            # TolFaclower * ideal_distance < distance < TolFacupper * ideal_distance 
            # ideal ideal_distance = Rcov(Atom1) + Rcov(Atom2)
            d = np.linalg.norm(xyz1-xyz2) # distance
            i_d = cls.rcov[AseAtoms[i].number] + cls.rcov[AseAtoms[j].number] # ideal distance
            if d <= i_d*rfacup and d >= i_d*rfacdown:
                return True
        return False
    
    
def DetermineSurfaceLayerZ(ASEAtoms, SurfaceAtomSymbols, ZVecIndex = 2, ztol = 0.5):
    """
    Find top layer surface atom z coordinates by averaging
    atoms within ztol (angstrom) of the top most atoms are selected for averaging
    
    Input List
    ASEAtoms:           ASE atoms containing adsorbate/surface system.
    SurfaceAtomSymbols: Symbol of surface atoms.
    ZVecIndex:          index of cell basis vector that is orthogonal to surface.
    ztol:               Atoms within ztol(angstrom) of the top most atoms are selected as 
                        surface atoms.
    Output List
    SurfaceLayerZ:      z coordinate of surface layer.
    SurfaceAtomIndex:   Index of surface atoms.
    """
    assert isinstance(ASEAtoms,ase_Atoms)
    # get highest surface atom coordinate
    zmax = 0
    zs = ASEAtoms.get_scaled_positions()[:,2] #scale all z coordinates
    #Go over all surface atoms and if z coordinate of atom is greater
    # than zmax, replace the zmax
    for i in range(0,len(ASEAtoms)):
        if ASEAtoms[i].symbol in SurfaceAtomSymbols and zmax < zs[i]:
            zmax = zs[i]
            
    # determine z coordinate. average out top layer
    ztol = ztol/np.linalg.norm(ASEAtoms.cell[2,:])
    SurfaceAtomIndex = list()
    SurfZs = list()
    for i in range(0,len(ASEAtoms)):
        if ASEAtoms[i].symbol in SurfaceAtomSymbols and zmax - ztol < zs[i]:
            SurfZs.append(zs[i]) #scaled position
            SurfaceAtomIndex.append(i) #surfaceatom index
    if len(SurfZs) > 0:
        SurfaceLayerZ = np.array(SurfZs).mean() #estimate mean of all atoms
    else:
        SurfaceLayerZ = 0 # initialise SurfacelayerZ with 0 if empty slab or gas-phase species.

    return SurfaceLayerZ, SurfaceAtomIndex