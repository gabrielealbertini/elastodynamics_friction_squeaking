#!/usr/bin/env python
from __future__ import print_function
import math

# DEFRIG

#from FunctionCollection import getBestPartition

file_str = dict()
file_str["input"] = """user parameters [

     simulation_name = {SIM_NAME}

     # FOLDERS
     output_folder       = {OUT_DIR}
     paraview_folder     = {OUT_DIR}
     restart_load_folder = {RESTART_DIR}
     restart_dump_folder = {RESTART_DIR}static_equilibrium/

     # MESH
     spatial_dimension = 2
     mesh              = {MESH_DIR}/block_2d_L{LGTH}_{HGHT}_msh{MESH}_quad4.msh
     antisym_setup     = {ASYMB} # does not matter for analytical

     # BOUNDARY CONDITIONS
     top_traction   = [ {TAUI:1.2e} ,-{SIGI:1.2e} , 0. ]
     bot_traction   = [-{TAUI:1.2e} , {SIGI:1.2e} , 0. ]
     left_traction  = [     0 , -{TAUI:1.2e} , 0. ]
     right_traction = [     0 ,  {TAUI:1.2e} , 0. ]

     #block_x_bcs = topleftblock_bottom_left # ,-separated
     #block_y_bcs = # ,-separated
     #block_z_bcs = topleftblock_bottom_back # ,-separated

     maximal_iteration = 2
     precision         = 1e-5
]

material elastic [
     name = slider
     rho  = {RHO}
     nu   = {NU:1.2f}
     E    = {E:1.2e}
     Plane_Stress = {PSTRESS}
]


"""


file_str["sub"] = """#!/bin/sh
#PBS -N {SIM_NAME}
#PBS -S /bin/sh
#PBS -V
#PBS -j oe
#PBS -m bea
#PBS -M ga288@cornell.edu
#PBS -l select={NB_NODES}:ncpus={CPUPN}:mpiprocs={MPIPPN}
#PBS -l walltime={HRS}:00:00

source ~/module-list.sh

cd $PBS_O_WORKDIR

mpirun -np {MPIP} ./{METHOD} {INPUT_FILE} 2>&1 | tee $PBS_JOBID.progress
"""

# -------------------------------------------------------------
# parameters
pdict = dict()

# ubwonko
pdict["OUT_DIR"]     = 'raw-output-data/'
pdict["RESTART_DIR"] = 'restart/'
pdict["MESH_DIR"]    = 'meshes'
pdict["CPUPN"]  = 48 #cpu_per_node

pdict['SIGI'] = 1e3
pdict['TAUI'] = 0e1


if 1: # green rubber
    pdict["E"]    = 1.6e6
    pdict["NU"]   = 0.45
    pdict["RHO"]  = 1050
    pdict["PSTRESS"] = 0

else:
    pdict["E"]    = 2e6
    pdict["NU"]   = 0.45#3#5
    pdict["RHO"]  = 1000
    pdict["PSTRESS"] = 0

if 0: # green rubber prony series
    E_inf = 1.196*1e6
    pdict["E"]    = E_inf # + np.sum(Ev)
    
    pdict["NU"]   = 0.45
    pdict["RHO"]  = 1134
    pdict["PSTRESS"] = 0

pdict["LGTH"] = 0.04
pdict["HTOP"] = 0.013

pdict["MESH"] = 1600 
pbc = False

# 'full_mpi_static_solution' , 'analytic_static_solution'
pdict["METHOD"]   = 'analytic_static_solution' # 2d works for now only with analytic solution
pdict["HRS"] = 1

# tagging or scotch
pdict["PART_MET"] = 'tagging' # 'scotch' , 'tagging'

if pdict["METHOD"] == "full_mpi_static_solution":
    if pdict["PART_MET"] == 'tagging': # tagging
        pdict["MPIPPN"]   = 10 # mpi_processes_per_node
        pdict['XPROC']    = 1  # has to be =1 if pbc aside from that it can be anything
        pdict['YPROC']    = 3  # should be an odd number
        # -------------------------------------------------------
        if pbc and not pdict['XPROC'] == 1:
            print("XPROC has to be zero")
            pdict['XPROC'] = 1
        pdict['MPIP'] = pdict['XPROC'] * pdict['YPROC']
        pdict['NB_NODES'] = int(math.ceil(pdict['MPIP'] / float(pdict['MPIPPN'])))
    elif pdict["PART_MET"] == 'scotch': # scotch
        pdict["MPIPPN"]   = 1 #mpi_processes_per_node
        pdict["NB_NODES"] = 1
        # -------------------------------------------------------
        pdict["MPIP"] = pdict["NB_NODES"] * pdict["MPIPPN"] # mpi process
        pdict['XPROC'] = 1 # no effect, just to make this script work
        pdict['YPROC'] = 1 # no effect, just to make this script work
elif pdict["METHOD"] == 'analytic_static_solution':
    # -------------------------------------------------------
    # method is not parallel
    pdict["MPIPPN"]   = 1 # mpi_processes_per_node
    pdict["NB_NODES"] = 1
    pdict["MPIP"] = pdict["NB_NODES"] * pdict["MPIPPN"] # mpi process
    pdict['XPROC'] = 1 # no effect, just to make this script work
    pdict['YPROC'] = 1 # no effect, just to make this script work


if __name__ == "__main__":
    model_top_and_bot = False
    series_parameter = "TAUI"
    series = [pdict['TAUI']]
    
    nb_sim = len(series)

    counter_start = 1
    for m in range(nb_sim):
        


        pdict[series_parameter] = series[m]

        if model_top_and_bot == False : #asymmetric simulation
            pdict["ASYMB"] = "true"
            pdict["HGHT"]  = "Ht{HTOP}".format(**pdict)
            pdict["MATT"]  = "_TE{E:1.2e}nu{NU:1.2f}rho{RHO}pstress{PSTRESS}".format(**pdict)
            pdict["MATB"]  = ""
        else:
            pdict["ASYMB"] = "false"
            pdict["HGHT"]  = "Ht{HTOP}b{HBOT}".format(**pdict)
            #pdict["MATL"]  = "_LE{E_LEFT:1.2e}nu{NU_LEFT:1.2f}rho{RHO_LEFT}".format(**pdict)
            #pdict["MATR"]  = "_RE{E_RIGHT:1.2e}nu{NU_RIGHT:1.2f}rho{RHO_RIGHT}".format(**pdict)

        if pbc: # periodic boundary conditions in x direction
            pdict['LPBC'] = 'pbc'
            pdict['PBC'] = '#'
            pdict['NOTPBC'] = ''
        else:
            pdict['LPBC'] = ''
            pdict['PBC'] = ''
            pdict['NOTPBC'] = '#'

        sim_name = 'equi_pure_shear_2d{MATT}{MATB}_taui{TAUI:1.2e}sigi{SIGI:1.2e}_L{LGTH}{HGHT}_elem{MESH}'.format(**pdict)

        pdict["SIM_NAME"] = sim_name

        #sub_file   = "static_2d_{}.sub".format(m + counter_start)
        #input_file = "static_2d_{}.in".format(m + counter_start)
        sub_file   = "{}.sub".format(sim_name)
        input_file = "{}.in".format(sim_name)
        pdict["INPUT_FILE"] = input_file

        

        
        with open(input_file, "w") as inp:
            print(input_file)
            print(file_str["input"].format(**pdict),file=inp)
        with open(sub_file,"w") as sub:
            print(file_str["sub"].format(**pdict),file=sub)
