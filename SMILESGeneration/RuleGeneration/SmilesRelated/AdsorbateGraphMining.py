from rdkit import Chem
from rdkit.Chem import RWMol
from collections import defaultdict

try:
  basestring
except NameError:
  basestring = str

def _GetSMILES(mol,idxlist):
    tmol = mol.__copy__() #(t)emporary
    tmol = RWMol(tmol)
    for AtomIdx in range(tmol.GetNumAtoms()-1,-1,-1):
        if AtomIdx not in idxlist:
            tmol.RemoveAtom(AtomIdx)
    #print(Chem.MolToSmiles(tmol))
    return Chem.MolToSmiles(tmol)

def LumpH(molecule):
    """
    Lump hydrogen atoms as a single atom. Note that Si, Al, Mg, Na are used as 
    pseudoatoms. However, this does not affect printing SMILES, as smilesSymbol
    are appropriately set.
    """
    molecule = Chem.RWMol(molecule)
    Hidx = list()
    for i in range(0,molecule.GetNumAtoms()):
        atom = molecule.GetAtomWithIdx(i)
        if atom.GetSymbol() != 'H':
            NumH = 0
            for neighbor_atom in atom.GetNeighbors():
                if neighbor_atom.GetSymbol() == 'H':
                    NumH += 1
                    Hidx.append(neighbor_atom.GetIdx())
            if NumH == 4:
                a = Chem.Atom('Si')
                a.SetProp('smilesSymbol','H4')
                idx = molecule.AddAtom(a)
                molecule.AddBond(atom.GetIdx(),idx,Chem.rdchem.BondType.QUADRUPLE)
                molecule.GetAtomWithIdx(idx).SetNoImplicit(True)
            elif NumH == 3:
                a = Chem.Atom('Al')
                a.SetProp('smilesSymbol','H3')
                idx = molecule.AddAtom(a)
                molecule.AddBond(atom.GetIdx(),idx,Chem.rdchem.BondType.TRIPLE)
                molecule.GetAtomWithIdx(idx).SetNoImplicit(True)
            elif NumH == 2:
                a = Chem.Atom('Mg')
                a.SetProp('smilesSymbol','H2')
                idx = molecule.AddAtom(a)
                molecule.AddBond(atom.GetIdx(),idx,Chem.rdchem.BondType.DOUBLE)
                molecule.GetAtomWithIdx(idx).SetNoImplicit(True)
            elif NumH == 1:
                a = Chem.Atom('Na')
                a.SetProp('smilesSymbol','H')
                idx = molecule.AddAtom(a)
                molecule.AddBond(atom.GetIdx(),idx,Chem.rdchem.BondType.SINGLE)
                molecule.GetAtomWithIdx(idx).SetNoImplicit(True)
    Hidx.sort(reverse=True)
    for i in Hidx:
        molecule.RemoveAtom(i)
    return molecule

def GetAdsorbateGraphsOfRadius(mol,SurfaceAtomSymbols,radius):
    """
    Adsorbate Graph mining Tool. Molecular graph of the surface and adsorbate 
    are inputed. Graphs radially increases similar to Morgan Fingerprint.
    
    Input List
    mol:                RdKit Mol object of surface and adsorbate
    SurfaceAtomSymbols: Surface atom symbols
    radius:             Desired radius
    
    Output List
    SMILES_count:       dictionary of SMILES -> count
    """
    # check if mol is rdkit mol    
    # each smiles contains indexes, which are later converted to mol.
    assert isinstance(mol,Chem.rdchem.Mol)
    assert radius >= 1
    # Go through the molecule, add non-surface atoms to relevant atom list
    NSAL = list() # (N)on(S)urface (A)tom (L)ist
    SAL = list() # (S)urface (A)tom (L)ist
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in SurfaceAtomSymbols:
            NSAL.append(atom.GetIdx())
        else:
            SAL.append(atom.GetIdx())
    # Get finger prints
    smiless = list() # (F)inger(p)rint(s). It's in atom index, and will be
                    # converted to SMILES later. This is done to remove duplicates
    ## Radius = 1 case
    for idx in NSAL:
        smiless.append([idx])
    ## Radius > 1 cases
    for RA in NSAL:
        OAL = list() # (O)ld (A)tom (L)ist
        NAL = list() # (N)ew (A)tom (L)ist
        NAL.append(RA)
        for i in range(1,radius):
            ACL = list(NAL) # (A)tom to (C)heck (L)ist
            NAL = list()
            for AC in ACL:
                OAL.append(AC)
                for NeighborAtom in mol.GetAtomWithIdx(AC).GetNeighbors():
                    NeighborAtomIdx = NeighborAtom.GetIdx()
                    if NeighborAtomIdx not in ACL+OAL+NAL and \
                        NeighborAtom.GetSymbol() not in SurfaceAtomSymbols:
                        NAL.append(NeighborAtomIdx)
            smiles = sorted(OAL+NAL)
            if smiles not in smiless:
                smiless.append(smiles)
    
    # Add surface atom that is neighbor to SMILES
    for smiles in smiless:
        RSAL = list() # (R)elevant (S)urface (A)tom (L)ist
        for SA in SAL:
            for NeighborAtom in mol.GetAtomWithIdx(SA).GetNeighbors():
                NeighborAtomIdx = NeighborAtom.GetIdx()
                if NeighborAtomIdx in smiles:
                    RSAL.append(SA)
        smiles += RSAL
    # convert atom indexes to SMILES
    SMILES_count = defaultdict(int)
    for smiles in smiless:
        SMILES_count[_GetSMILES(mol,smiles)] += 1
    return SMILES_count
    
def FindAllAdsorbateGraphsOfLengthN(mol1,SurfaceAtomSymbols,NMaxEdge,NMinEdge=0,valence_based = False, debug = False):
    """
    Adsorbate Graph mining Tool. Molecular graph of the adsorbate and surface
    are inputed. Graphs of any shape with lengths specified are obtained
    
    Input List
    mol1:               RdKit Mol object of surface and adsorbate
    SurfaceAtomSymbols: Surface atom symbols
    NMaxEdge:           Maximum number of edges. if -1, find everything
    NMinEdge:           Minimum number of edges
    valence_based:      Valence based graphs obtained
    debug:              Debugging function
    
    Output List
    SMILES_count:       dictionary of SMILES -> count
    """
    # AdsorbateSMILES of length N is found. Not radial.
    # if NMaxEdge = -1, find everything
    
    # mol2 is copied, then surface atoms are removed, and the enumeration is performed
    # check if mol is rdkit mol    
    assert isinstance(mol1,Chem.rdchem.Mol)
    print(mol1)
#    mol3 = Chem.AddHs(mol1)
    mol1 = LumpH(mol1)
    mol2 = mol1.__copy__() # mol without surface atoms
    mol2 = Chem.RWMol(mol2)
    # Go through the molecule, add non-surface atoms to relevant atom list
    Mol2toMol1Map = dict()
    NSAL = list() # (N)on(S)urface (A)tom (L)ist
    SAL = list() # (S)urface (A)tom (L)ist
    mol2i = 0
    for mol1i in range (0,mol1.GetNumAtoms()):
        atom = mol1.GetAtomWithIdx(mol1i)
        if atom.GetSymbol() not in SurfaceAtomSymbols:
            NSAL.append(mol1i)
            Mol2toMol1Map[mol2i] = mol1i
            mol2i += 1
        else:
            SAL.append(mol1i)
            mol2.RemoveAtom(mol2i)
    # check N
    if NMaxEdge >= mol2.GetNumBonds():
        NMaxEdge = mol2.GetNumBonds()
    if NMinEdge >= mol2.GetNumBonds():
        NMinEdge = mol2.GetNumBonds()
    if NMaxEdge == -1:
        NMaxEdge = mol2.GetNumBonds()
    # (S)ubgraph (B)ond (I)ndex (L)ist(s)
    SBILs = list()
    if NMaxEdge > 0:
        SBILs = Chem.FindAllSubgraphsOfLengthMToN(mol2,NMinEdge,NMaxEdge)
    # valence_based
    if valence_based:
        bonds = mol1.GetBonds()
        # remove surface atom - surface atom bond
        for i in range(len(bonds)-1,-1,-1):
            if bonds[i].GetBeginAtom().GetSymbol() in SurfaceAtomSymbols and \
                bonds[i].GetEndAtom().GetSymbol() in SurfaceAtomSymbols:
                mol1.RemoveBond(bonds[i].GetBeginAtomIdx(),bonds[i].GetEndAtomIdx())
        atoms = mol1.GetAtoms()
        # Separate binding sites
        for i in range(len(atoms)-1,-1,-1):
            if atoms[i].GetSymbol() in SurfaceAtomSymbols:
                neighbors = atoms[i].GetNeighbors()
                neighboridx = list()
                for neighbor_atom in neighbors:
                    neighboridx.append(neighbor_atom.GetIdx())
                    
                if len(neighboridx) > 1: # it has connection to organic atoms
                    for j in range(1,len(neighboridx)):
                        # add valence
                        mol1.RemoveBond(atoms[i].GetIdx(),neighboridx[j])
                        new_surf_atom = Chem.Atom(SurfaceAtomSymbols)
                        new_surf_atom_idx = mol1.AddAtom(new_surf_atom)
                        mol1.AddBond(new_surf_atom_idx,neighboridx[j],Chem.rdchem.BondType.ZERO)
                    
    # debug
    if debug:
        print(Chem.MolToSmiles(mol1))
    # Get finger prints
    smiless = list() # (F)inger(p)rint(s). It's in atom index, and will be
                    # converted to SMILES later. This is done to remove duplicates
    ## length = 0 case
    if NMinEdge ==0:
        for idx in NSAL:
            smiless.append([idx])
    ## Radius > 1 cases
    # add organic atoms
    for SBIL in SBILs:
        for SBI in SBIL:
            AL = set() # (A)tom (L)ist
            for BI in SBI:
                B = mol2.GetBondWithIdx(BI)
                AL.add(Mol2toMol1Map[B.GetBeginAtomIdx()])
                AL.add(Mol2toMol1Map[B.GetEndAtomIdx()])
            smiless.append(list(AL))
    
    # Add surface atom that is neighbor to SMILES
    # add surface atoms
    for smiles in smiless:
        RSAL = list() # (R)elevant (S)urface (A)tom (L)ist
        for SA in SAL:
            for NeighborAtom in mol1.GetAtomWithIdx(SA).GetNeighbors():
                NeighborAtomIdx = NeighborAtom.GetIdx()
                if NeighborAtomIdx in smiles:
                    RSAL.append(SA)
        smiles += RSAL
    # convert atom indexes to SMILES
    SMILES_count = defaultdict(int)
    for smiles in smiless:
        SMILES_count[_GetSMILES(mol1,smiles)] += 1
    return SMILES_count
