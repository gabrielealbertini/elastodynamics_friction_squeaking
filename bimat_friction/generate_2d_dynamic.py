#!/usr/bin/env python
from __future__ import print_function
import math
import numpy as np

from solidmechanics import smd
from solidmechanics import LinearElasticMaterial as matlaw
import os


def read_counter():
    filename = "counter.txt"

    # Check if file exists
    if os.path.exists(filename):
        with open(filename, "r") as f:
            counter = int(f.read().strip())
    else:
        counter = 1
        with open(filename, "w") as f:
            f.write(str(counter))

    print(f"Counter: {counter}")
    return counter

def write_counter(counter):
    filename = "counter.txt"

    with open(filename, "w") as f:
        f.write(str(counter))
    
# DEFRIG

file_str = dict()
file_str["input"] = """user parameters [
     # NAMES

     simulation_name      = {SIM_NAME}
     load_simulation_name = {LOAD_SIM_NAME}

     # FOLDERS
     output_folder       = {OUT_DIR}
     paraview_folder     = paraview/
     restart_load_folder = {RESTART_DIR}static_equilibrium/
     restart_dump_folder = {RESTART_DIR}

     # BOUNDARY CONDITIONS
     {PBC}top_traction   = [ {TAUI:1.2e} , -{SIGI:1.2e} ]
     {PBC}bot_traction   = [-{TAUI:1.2e} ,  {SIGI:1.2e} ]
     {PBC}left_traction  = [     0 , -{TAUI:1.2e} ]
     {PBC}right_traction = [     0 ,  {TAUI:1.2e} ]

     {IS_SHEARVEL}shear_velocity = {SHEARVEL:1.4f}
     {IS_SHEARVEL}shear_acceleration_time = {INITACCEL:1.2e}
     {IS_SHEARVEL}keep_accelerating = {KEEPACCEL:1.4f}

     {NOTPBC}block_x_bcs = {BLKX}
     {NOTPBC}block_y_bcs = {BLKY}
     {ZSYMC}block_z_bcs  = {ZSYMB}
     
     # INTERFACE
     vertical_normals = true
     pretend_is_in_contact = {PRET_CONTACT}

     # NUCLEATION
     nuc_center    = {NUCC:1.1e}
     nuc_half_size = {NUCHS:1.1e}
     nuc_w_zone    = {NUCWZS:1.2e}
     nuc_speed     = {NUCV:1.2e}

{IHET}

     # MESH
     spatial_dimension = 2
     mesh              = {MESH_DIR}/block_2d_L{LGTH}_{HGHT}_msh{MESH}_quad4.msh
     defrig_setup      = true
     
     # TIME
     simulation_time  = {DUR}
     time_step_factor = {TSF:1.2f}

     # DUMP
     #dump_heights          = [0.0075] # <-------------------------------------------------------------------------
     dump_height_tolerance = 5e-8  

     # DUMP AT A PRECISE TIME
     # precise_dump_time = 6.51586e-05

     # FULL DUMP
     #full_dump_f = 0.0001

     # FREQUENCY THRESHOLD
     rel_slip_for_highf = 0.1
          
     # INTERFACE TEXT DUMPER 
     interface_dump_highf = {DFREQ:1.2e}
     interface_dump_lowf  = {DFREQ:1.2e} #1e-9 # <-------------------------------------------------------------------------
     interface_dump_fields = {DFLDS}
     
     # AT DISTANCE TEXT DUMPER
     at_distance_dump_highf = {DFREQ:1.2e}
     at_distance_dump_lowf  = {DFREQ:1.2e} #1e-9 # <-------------------------------------------------------------------------

     # GLOBAL TEXT DUMPER
     global_dump_highf = {DFREQ:1.2e}
     global_dump_lowf  = {DFREQ:1.2e} #1e-9 # <-------------------------------------------------------------------------
     
     # FULL PARAVIEW DUMPER  # value cannot be given in 1e6 format
     dump_paraview = false
     nb_max_paraview_dumps = 1000
     paraview_dump_f    = 0.0001
]

{FRICLAW}

{MAT}

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

source /home/ga288/module-list-ubwonko.sh

cd $PBS_O_WORKDIR

mpirun -np {MPIP} ./mixed_mode_fracture {INPUT_FILE} 2>&1 | tee $PBS_JOBID.progress

mpirun -np 6 --allow-run-as-root ./mixed_mode_fracture {INPUT_FILE} 2>&1 | tee {SIM_NAME}.progress
"""

# -------------------------------------------------------------
# parameters
pdict = dict()

# ubwonko
pdict["OUT_DIR"]     = 'raw-output-data/'
pdict["RESTART_DIR"] = 'restart/'
pdict["MESH_DIR"]    = 'meshes'
pdict["CPUPN"]  = 48 #cpu_per_node
pdict["MPIPPN"] = pdict["CPUPN"] # mpi_processes_per_node

# load
pdict["SIGI"] = 1e3 # <-------------------------------------------------------------!!!!!!!
pdict["TAUI"] = 0e1

# material
is_elastic = False
is_green_rubber = True # zhermack elite double 32

is_gen_maxwell = False
is_test_gen_maxwell = False
is_high_low_freq = False#True

if is_green_rubber:
    pdict["E"]    = 1.6e6
    pdict["NU"]   = 0.45
    pdict["RHO"]  = 1050
    pdict["PSTRESS"] = 0

else:
    pdict["E"]    = 2e6
    pdict["NU"]   = 0.45#3# 0.45
    pdict["RHO"]  = 1000
    pdict["PSTRESS"] = 0

mat = matlaw({
    smd.E   : pdict["E"],
    smd.nu  : pdict["NU"],
    smd.rho : pdict["RHO"],
    smd.pstress : pdict["PSTRESS"],
})
#mat.complete()
mat.complete()
mat.complete_wave_speeds()
mat_cd = mat[smd.cp]
mat_cs = mat[smd.cs]
print(mat)

# elastic
if is_elastic:
    pdict["MAT"] = """
material elastic [
     name = slider
     rho  = {RHO}
     nu   = {NU:1.2f}
     E    = {E:1.2e}
     Plane_Stress = {PSTRESS}

]""".format(**pdict)
    
elif is_gen_maxwell:
    if is_test_gen_maxwell:
        # use params of option -1
        
        pdict["E"]    = 1.6e6
        pdict["NU"]   = 0.45
        pdict["RHO"]  = 1050
        pdict["PSTRESS"] = 0
        #E_inf = 1.6e6-6.4e5
        pdict['EV0'] = 640000.0
        pdict['EV1'] = 0.0
        pdict['EV2'] = 0.0
        pdict['EV3'] = 0.0
        pdict['EV4'] = 0.0
        pdict['EV5'] = 0.0
        pdict['EV6'] = 0.0
        pdict['EV7'] = 0.0
        pdict['EV8'] = 0.0
        pdict['EV9'] = 0.0
        
        pdict['ETA0'] = 64.0
        pdict['ETA1'] = 1.0
        pdict['ETA2'] = 1.0
        pdict['ETA3'] = 1.0
        pdict['ETA4'] = 1.0
        pdict['ETA5'] = 1.0
        pdict['ETA6'] = 1.0
        pdict['ETA7'] = 1.0
        pdict['ETA8'] = 1.0
        pdict['ETA9'] = 1.0
    elif is_high_low_freq:        # use params of option -1 and 3
        

        Ev = [640000.0,1040000.0]
        eta = [64.0,1.04]
        
        pdict["NU"]   = 0.45
        pdict["RHO"]  = 1050
        pdict["PSTRESS"] = 0
        E_inf = 1.6e6
        pdict["E"]    = E_inf +np.sum(Ev)
        pdict['EV0'] = Ev[0]
        pdict['EV1'] = Ev[1]
        pdict['EV2'] = 0.0
        pdict['EV3'] = 0.0
        pdict['EV4'] = 0.0
        pdict['EV5'] = 0.0
        pdict['EV6'] = 0.0
        pdict['EV7'] = 0.0
        pdict['EV8'] = 0.0
        pdict['EV9'] = 0.0
        
        pdict['ETA0'] = eta[0]
        pdict['ETA1'] = eta[1]
        pdict['ETA2'] = 1.0
        pdict['ETA3'] = 1.0
        pdict['ETA4'] = 1.0
        pdict['ETA5'] = 1.0
        pdict['ETA6'] = 1.0
        pdict['ETA7'] = 1.0
        pdict['ETA8'] = 1.0
        pdict['ETA9'] = 1.0
    else:
        # use params of elite double 32
        HFD=False
        if HFD:
            # Prony series coefficients of Elite Double 32 [Haque et al., 2018]
            # more high freq damping
            Ev = np.array([0.477*1.75, 0.328*1.5, 0.221*1.5, 0.157, 0.113, 0.127, 0.114, 0.077, 0.059, 0.043])*1e6
            tau = np.array([4.99e-8, 4.83e-7, 3.55e-6, 1.72e-5, 6.04e-5,
                         2.28e-4, 1.33e-3, 9.07e-3, 6.42e-2, 4.89e-1])
            
        else:
            # Prony series coefficients of Elite Double 32 [Haque et al., 2018]
            Ev = np.array([0.477, 0.328, 0.221, 0.157, 0.113, 0.127, 0.114, 0.077, 0.059, 0.043])*1e6
            tau = np.array([4.99e-8, 4.83e-7, 3.55e-6, 1.72e-5, 6.04e-5,
                         2.28e-4, 1.33e-3, 9.07e-3, 6.42e-2, 4.89e-1])
        eta = tau*Ev
        E_inf = 1.196*1e6
        pdict["E"]    = E_inf + np.sum(Ev)
        
        pdict["NU"]   = 0.45
        pdict["RHO"]  = 1134
        pdict["PSTRESS"] = 0
        
        pdict['EV0'] = Ev[0]
        pdict['EV1'] = Ev[1]
        pdict['EV2'] = Ev[2]
        pdict['EV3'] = Ev[3]
        pdict['EV4'] = Ev[4]
        pdict['EV5'] = Ev[5]
        pdict['EV6'] = Ev[6]
        pdict['EV7'] = Ev[7]
        pdict['EV8'] = Ev[8]
        pdict['EV9'] = Ev[9]
        
        pdict['ETA0'] = eta[0]
        pdict['ETA1'] = eta[1]
        pdict['ETA2'] = eta[2]
        pdict['ETA3'] = eta[3]
        pdict['ETA4'] = eta[4]
        pdict['ETA5'] = eta[5]
        pdict['ETA6'] = eta[6]
        pdict['ETA7'] = eta[7]
        pdict['ETA8'] = eta[8]
        pdict['ETA9'] = eta[9]

    pdict["MAT"] = """
material gen_maxwell_deviatoric [
     name = slider
     rho  = {RHO}
     nu   = {NU:1.2f}
     E    = {E:1.2e}
     Ev0   = {EV0:1.2e}
     Ev1   = {EV1:1.2e}
     Ev2   = {EV2:1.2e}
     Ev3   = {EV3:1.2e}
     Ev4   = {EV4:1.2e}
     Ev5   = {EV5:1.2e}
     Ev6   = {EV6:1.2e}
     Ev7   = {EV7:1.2e}
     Ev8   = {EV8:1.2e}
     Ev9   = {EV9:1.2e}
     Eta0  = {ETA0:1.2e}
     Eta1  = {ETA1:1.2e}
     Eta2  = {ETA2:1.2e}
     Eta3  = {ETA3:1.2e}
     Eta4  = {ETA4:1.2e}
     Eta5  = {ETA5:1.2e}
     Eta6  = {ETA6:1.2e}
     Eta7  = {ETA7:1.2e}
     Eta8  = {ETA8:1.2e}
     Eta9  = {ETA9:1.2e}
     Plane_Stress = {PSTRESS}        
 ]""".format(**pdict)

    pdict['E']=E_inf
    
else: # sls dev
    if is_green_rubber:

        if 1: # option -1 peak damping at f=1.5e3
            # E 1600000.0
            pdict['EV'] = 640000.0
            # tau 1e-05
            pdict['ETA'] = 64.0
        if 1: # option 0 peak damping at f=3e3
            # E 1600000.0
            pdict['EV'] = 640000.0
            # tau 1e-05
            pdict['ETA'] = 32.0

        if 0: # option 1 peak damping at f=1.5e4
            # E 1600000.0
            pdict['EV'] = 640000.0
            # tau 1e-05
            pdict['ETA'] = 6.4

        if 0: # option 2 peak damping at f=3e4
            # E 1600000.0
            pdict['EV'] =  720000.0
            # tau 5e-06
            pdict['ETA'] = 3.6
            
        if 0: # option 3 peak damping at f=1.5e5
            # E 1600000.0
            pdict['EV'] =  1040000.0
            # tau 1e-06
            pdict['ETA'] = 1.04
    
    else:
        pdict["ETA"] = 1e2
        pdict["EV"]  = 0.4*pdict["E"]
    pdict["MAT"] = """
material sls_deviatoric [
     name = slider
     rho  = {RHO}
     nu   = {NU:1.2f}
     E    = {E:1.2e}
     Ev   = {EV:1.2e}
     Eta  = {ETA:1.2e}
     Plane_Stress = {PSTRESS}

]""".format(**pdict)



pdict["LGTH"] = 0.04 
pdict["HTOP"] = 0.02
pdict["MESH"] = 400

pbc = False

pdict['DUR'] = 0.06#4 * float(pdict["LGTH"]) / float(mat_cs) # < -------------------------
pdict['TSF'] = 0.1

pdict["NUCC"]   = 0.02
pdict["NUCHS"]  = 0.005 
pdict["NUCWZS"] = 0.001
pdict["NUCV"]   = 10.0

# 

# friction
is_fric=True#False
if is_fric:
    pdict["MUS"] = 2.1 #0.94 #0.74 #1.06
    pdict["MUK"] = 2.1 #1.90
    pdict["DC"]  = 1e-6 #2.24e-6 #1.4e-6
    pdict["TSTAR"] = 1e-4#5e-4#4#5e-5
    pdict["FLAW"] = "lswnhs{MUS:1.2f}k{MUK:1.2f}d{DC:1.2e}".format(**pdict)
    pdict["FRICLAWtmp"] = """
friction linear_slip_weakening simplified_prakash_clifton [ # rubin_ampuero [ no_regularisation [ 
     mu_s = {MUS:1.2f}
     mu_k = {MUK:1.2f}
     d_c  = {DC:1.2e}
     t_star = {TSTAR:1.2e}
]
"""#.format(**pdict)

    pdict["IHETS"] = 0.04
    pdict["IHETT"] = 0.038
    pdict["IHETDC"] = pdict["DC"]
    pdict["IHETMUS"] = 0.0
    pdict["IHETMUK"] = 0.0
    pdict["IHET"]="""
     # INTERFACE HET
     interface_heterogeneity      = true
     heterogeneity_start_position = {IHETS:1.2e}
     heterogeneity_trans_position = {IHETT:1.2e}
     heterogeneity_d_c  = {IHETDC:1.2e}
     heterogeneity_mu_s = {IHETMUS:1.2f}
     heterogeneity_mu_k = {IHETMUK:1.2f}"""#.format(**pdict)

else: # modeI
    pdict["G_C"] = 0.01e6
    pdict["TAU_C"] = 1e6 

    Gc = pdict["G_C"]
    E = pdict["E"]
    sigy = pdict["TAU_C"]
    Xc = Gc*E/(2*np.pi*sigy**2)
    Dx = pdict['LGTH']/pdict['MESH']
    print('Xc =',Xc,'Xc/Dx =' ,Xc/Dx)
    pdict["IHET"]=""
    pdict["FRICLAW"] = """
friction linear_normal_cohesive [
     G_c   = {G_C:1.2e}
     tau_c = {TAU_C:1.2e}
]
"""#.format(**pdict)
    
if False:
    pdict["G_C"]   = 1.12 # 1.12
    pdict["TAU_C"] = 4.7e6 #1e6 4.8e6 #5.3e6
    pdict["TAU_R"] = 3.7e6 #0   #3.2e6 #3.7e6
    pdict["FLAW"] = "gc{G_C:1.2f}tc{TAU_C:1.2e}tr{TAU_R:1.2e}".format(**pdict)
    pdict["FRICLAW"] = """
friction linear_cohesive [
     G_c   = {G_C:1.2e}
     tau_c = {TAU_C:1.2e}
     tau_r = {TAU_R:1.2e}
]
"""#.format(**pdict)
    
shear_vel=True
if shear_vel:
    fct =0.25
    pdict["SHEARVEL"] = fct*0.01#0.04#0.01 #m/s
    pdict["INITACCEL"] = fct*1.e-3 #1e-1
    pdict["KEEPACCEL"] = 0.0

    pdict["IS_SHEARVEL"]=''
else:
    pdict["IS_SHEARVEL"]='#'
    pdict["IHET"] = ""
if is_fric==False:
    pdict["IHET"] = "" #false

pdict["DFREQ"] = float(pdict["LGTH"]) / float(mat_cs) / 100 #float(pdict["MESH"])
pdict["DFLDS"] = 'displacement,velocity,friction_traction,contact_pressure # is_in_contact,is_sticking, # binary_shear_info,binary_normal_info,cohesive_shear,cohesive_normal,shear_gap_rate,shear_gap,normal_gap_rate,normal_gap,displacement,velocity,mu_k,mu_s,is_in_contact,is_sticking,friction_traction,contact_pressure,slip_velocity,slip,cumulative_slip' # <--------

pdict["HRS"] = 48
pdict["NB_NODES"] = 1 # <------------------------------------------------------------------------
pdict["MPIP"] = pdict["NB_NODES"] * pdict["MPIPPN"] # mpi process

#pretend_is_in_contact
pdict["PRET_CONTACT"] = "false"#"true"#

if __name__ == "__main__":

    if 1:
        series_parameter = "TSTAR"# "SHEARVEL"
        series =[1e-3]

        aspect_ratio = 8.0
        H= 0.005
        L=H*aspect_ratio
        pdict["HTOP"]=H
        pdict["LGTH"]=L
        pdict["IHETS"]= L-H/5.0
        pdict["IHETT"]= L-H/10.0
        pdict["DFREQ"]=pdict["DFREQ"]*H/0.02
        #series_parameter_1 = "IHETS"
        #series_1 = 0.04- H/5.0
        #series = [0.005]#0.02,0.013,0.01,0.005]#[1e-5, 2e-5, 3e-5, 4e-5, 5e-5] #[0.009]#[0.0075,0.005,0.0025] # 0.01
        pdict["MESH"]=int(800/2.0*aspect_ratio)
    else:
        series_parameter = "TSTAR"
        series =[1e-5]#,5e-5,1e-4,5e-4,1e-3]
    
    nb_sim = len(series)

    counter_start = read_counter() + 1
    
    
    for m in range(nb_sim):
        #if series[m]< 0.019:

            #pdict["MESH"] = 400*2
            
        counter = m + counter_start
        print(counter)
        bname ="dynamic_2d_{:0>3d}".format(counter)
        sub_file = bname+'.sub'
        input_file = bname+'.in'
        pdict["INPUT_FILE"] = input_file
        
        
        pdict[series_parameter] = series[m]
        #pdict[series_parameter_1] = series_1[m]
        print(series_parameter, series[m])
        #print(series_parameter_1, series_1[m])
        

        pdict["FRICLAW"] = pdict["FRICLAWtmp"].format(**pdict)
        pdict["IHET"]= pdict["IHET"].format(**pdict)
        
        if True: # defrig interface
            
            pdict["HGHT"] = "Ht{HTOP}".format(**pdict)
            pdict["MATT"]  = "_TE{E:1.2e}nu{NU:1.2f}rho{RHO}pstress{PSTRESS}".format(**pdict)
            pdict["MATB"]  = ""
            pdict["BLKX"] = "slider_top,slider_left,slider_right"
            pdict["BLKY"] = "slider_top,slider_left,slider_right"
            pdict["ZSYMB"] = "slider_back"


        """
        if pdict['HBOT'] == 0: #asymmetric simulation
            pdict["ASYMB"] = "true"
            pdict["HGHT"] = "Ht{HTOP}".format(**pdict)
            pdict["WTH"]  = "Wt{WTOP}".format(**pdict)
            pdict["MATT"]  = "_TE{E_TOP:1.2e}nu{NU_TOP:1.2f}rho{RHO_TOP}".format(**pdict)
            pdict["MATB"]  = ""
            pdict["ZSYMB"] = "slider_back"
            pdict["BLKX"] = "slider_top,slider_left,slider_right"
            pdict["BLKY"] = "slider_top,slider_left,slider_right"
        else:
            pdict["ASYMB"] = "false"
            pdict["HGHT"] = "Ht{HTOP}b{HBOT}".format(**pdict)
            pdict["WTH"]  = "Wt{WTOP}b{WBOT}".format(**pdict)
            pdict["MATT"]  = "_TE{E_TOP:1.2e}nu{NU_TOP:1.2f}rho{RHO_TOP}".format(**pdict)
            pdict["MATB"]  = "_BE{E_BOT:1.2e}nu{NU_BOT:1.2f}rho{RHO_BOT}".format(**pdict)
            pdict["ZSYMB"] = "slider_back,base_back"
            pdict["BLKX"] = "slider_top,slider_left,slider_right,base_bottom,base_left,base_right"
            pdict["BLKY"] = "slider_top,slider_left,slider_right,base_bottom,base_left,base_right"
        """
        if pbc: # periodic boundary conditions in x direction
            pdict['LPBC'] = 'pbc'
            pdict['PBC'] = '#'
            pdict['NOTPBC'] = ''
        else:
            pdict['LPBC'] = ''
            pdict['PBC'] = ''
            pdict['NOTPBC'] = '#'


        zsym=False

        if zsym: # symmetric w.r.t. to z
            pdict['ZSYM'] = 'zsym'
            pdict['ZSYMC'] = ''
        else:
            pdict['ZSYM'] = ''
            pdict['ZSYMC'] = '#'
        
        load_name = 'equi_pure_shear_2d{MATT}{MATB}_taui{TAUI:1.2e}sigi{SIGI:1.2e}_L{LGTH}{HGHT}_elem{MESH}'.format(**pdict)
        pdict["LOAD_SIM_NAME"] = load_name

        
        #sim_name = "dyn_3d{MATT}{MATB}_{FLAW}_taui{TAUI:1.2e}sigi{SIGI:1.2e}_nucC{NUCC:1.1f}HS{NUCHS:1.1e}WZ{NUCWZS:1.2e}V{NUCV:1.0f}_L{LGTH}{LPBC}{HGHT}{WTH}{ZSYM}msh{MESH}".format(**pdict)

        #sim_name="dynamic_2d_{}".format(counter)

        pdict["SIM_NAME"] = bname
        
        with open(input_file, "w") as inp:
            print(file_str["input"].format(**pdict),file=inp)
        with open(sub_file,"w") as sub:
            print(file_str["sub"].format(**pdict),file=sub)
        print(input_file)
        write_counter(counter)
