from ase.io import read
import numpy as np
from rdkit import Chem
from ase import Atoms as ase_Atoms
import networkx as nx

inp = str('CONTCAR')
ASEAT = read(inp)