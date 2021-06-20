from ase import Atoms
from ase.build import molecule
from ase.visualize import view
from ase.calculators.vasp import Vasp
from ase.optimize import LBFGS
from itertools import combinations, permutations
import numpy as np
from py_box3.ase.set_calc import print_vasp_param

'''
Creating the methane molecule
'''
cell = np.array([20., 20., 20.]) #A, Box has to be large due to periodic boundary conditions!

#Building the water molecule using ase.build
C2H6 = molecule('C2H6', cell = cell)

#Move the molecule to the center of the box
C2H6.center()

'''
Viewing the water molecule
'''
view(C2H6)

'''
Setting up the VASP calculator
'''
vasp_settings = {'xc': 'PBE',        #Functional to use
                 'kpts': (1, 1, 1),  #K Point mesh
                 'encut': 400.,      #Energy cutoff
                 'ismear': -1,        #Gaussian smearing
                 'sigma': 0.1,       #Width of smearing in eV
                 'ediff': 1e-6,      #Global break condition for the electronic selfconsistency loop
                 'prec': 'normal',   #Precision. Sets various other parameters
                 'lcharg': False,    #Whether LCHARG is written
                 'lwave': False,     #Whether LWAVE is written
                 'nelmin': 4,        #Minimum number of electronic selfconsistency steps
                 'nelmdl': -4,       #Number of non-selfconsistent steps at the beginning
                 'npar': 4,          #Number of bands that are treated in parallel
                 'algo': 'fast',     #Electronic minimization algorithm
                 'lreal': 'auto',    #Determines whether projection operators are evaluated in real or reciprocal space
                 'ispin': 1}         #Whether to use spin or non-spin polarized calculations
calc = Vasp(**vasp_settings)
print_vasp_param(calc) #Showcasing py_box3.ase.print_vasp_param
C2H6.set_calculator(calc)

'''
Setting up the optimizer and run calculation
'''
dyn = LBFGS(C2H6)
dyn.run(fmax = 0.02)

'''
Output energy
'''

print('Energy = {} eV'.format(H2O.get_potential_energy()))
