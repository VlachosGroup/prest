from ase.io import read
from ase.visualize import view

atoms = read('CHCH3/CONTCAR')
view(atoms)