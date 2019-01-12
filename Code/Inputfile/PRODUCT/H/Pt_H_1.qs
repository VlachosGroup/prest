#!/bin//bin/bash
#$ -cwd
#$ -pe mpi 20
#$ -l exclusive=1
#$ -o Pt_H_1.stdout
#$ -e Pt_H_1.stderr
#$ -l m_mem_free=2G
#$ -l h_cpu=36:00:00
#$ -m beas
#$ -M tavanes@udel.edu

source /etc/profile.d/valet.sh
vpkg_rollback all
vpkg_require openmpi/1.10.2-intel64-2016 ase/3.15.0:python3

export VASP_COMMAND="mpiexec -n 20 /home/work/ccei_biomass/programs/vasp5.4/vasp.5.4.1/build/std/vasp"
export VASP_PP_PATH=/home/work/ccei_biomass/programs/vasp_psp/v54/

python3 Pt_H_1.py >& Pt_H_1.out
