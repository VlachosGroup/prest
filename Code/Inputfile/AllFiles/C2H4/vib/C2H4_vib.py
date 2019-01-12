#!/usr/bin/env python

from ase.vibrations import Vibrations
from ase import io
from ase import Atoms
from ase.calculators.vasp import Vasp
import numpy as np
from py_box3.ase.set_calc import print_vasp_param

atoms =  io.read('POSCAR')

# Create vibration calculator
indice=[0,1,2,3,4,5]

vasp_settings = {'xc': 'PBE',       #Functional to use
                 'kpts': (1, 1, 1),  #K Point mesh
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

vib = Vibrations(atoms,indices=indice)
vib.run()
#vib.read()
vib.summary()

# Make trajectory files to visualize normal modes:
for mode in range(3*len(indice)):
    vib.write_mode(mode)
