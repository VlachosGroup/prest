from ase import io
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io.trajectory import Trajectory
from ase.calculators.vasp import Vasp
from ase.optimize import QuasiNewton
import numpy as np
from py_box3.ase.set_calc import print_vasp_param

atoms = io.read('qn.traj', index=-1)

'''
Freezing the bottom four layers
'''
mask = [atom.z < 10. for atom in atoms] #List of booleans indicating if the atoms are in the bottom two layers
c = FixAtoms(mask = mask) #Freeze all atoms in the bottom two layers
atoms.set_constraint(c)

'''
Setting up the VASP calculator
'''
vasp_settings = {'xc': 'PBE',       #Functional to use
                 'kpts': (6, 6, 1),  #K Point mesh
                 'encut': 400.,      #Energy cutoff
                 'ismear': -1,       #Fermi smearing
                 'sigma': 0.1,       #Width of smearing in eV
                 'ediff': 1e-6,      #Global break condition for the electronic selfconsistency loop
                 'prec': 'normal',   #Precision. Sets various other parameters
                 'lcharg': False,    #Whether LCHARG is written
                 'lwave': False,     #Whether LWAVE is written
                 'nelmin': 4,        #Minimum number of electronic selfconsistency steps
                 'nelmdl': -4,       #Number of non-selfconsistent steps at the beginning
                 'ncore': 4,         #Number of bands that are treated in parallel
                 'algo': 'fast',     #Electronic minimization algorithm
                 'lreal': 'auto',    #Determines whether projection operators are evaluated in real or reciprocal space
                 'ispin': 1}         #Whether to use spin or non-spin polarized calculations
calc = Vasp(**vasp_settings)
print_vasp_param(calc) #Showcasing py_box3.ase.print_vasp_param
atoms.set_calculator(calc)

io.write('InitialGeom.traj',atoms)
Traj=Trajectory('qn.traj','a',atoms)

'''
Setting up the optimizer and run calculation
'''
dyn = QuasiNewton(atoms,trajectory=Traj,logfile='qn.log',restart='qn.pckl')
dyn.run(fmax = 0.02) #fmax is the force criterion

'''
Output energy, atomic coordinates, distances, angle of bond
'''
print('Energy = {} eV'.format(atoms.get_potential_energy()))