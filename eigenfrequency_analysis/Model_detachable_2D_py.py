# -*- coding: mbcs -*-
#Author:  Adel Djellouli
#Date:    Sept. 2025
#Purpose: Abaqus script to create a 2D plane strain block with a vertical partition,
#         preloaded with pressure and shear displacement, then perform eigenfrequency analysis.

import os
import sys
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from regionToolset import Region

# Python 2.7-safe job name
H_mm = float(os.environ.get('H_MM'))
partition_fraction = 50   #percentage of H
jobName = "Block_H%03dmm" % int(round(H_mm))
params = {}
params['L']=0.04
params['W']=0.04
params['H'] = H_mm / 1000.0  # mm -> m
params['Partition_location']=(partition_fraction/100.)*params['H']
params['Mesh_size']= max(0.016/50.,params['H']/50.)
params['Material_E']=1600000.0  # Pa
params['Material_rho']=1060.0
params['DELTAX']=params['H']/10.      # 10 mm in meters
params['P']=1.0e4            # Pa
params['numEigen']=12         # request this many eigenpairs
params['Keep_pressure_in_shear'] = True  # keep or drop pressure load in shear step

L = params['L']; W = params['W']; H = params['H']
partition_location = params['Partition_location']
mesh_size = params['Mesh_size']
Material_E = params['Material_E']
Material_rho = params['Material_rho']
DELTAX = params['DELTAX']
P = params['P']
numEigen = params['numEigen']
keep_pressure_in_shear = params['Keep_pressure_in_shear']




def printAB(string_):
    print >> sys.__stdout__,  string_


modeldb = mdb.models['Model-1']
# Geometry (same sketch)
sk = modeldb.ConstrainedSketch(name='__profile__', sheetSize=0.08)
sk.rectangle(point1=(0.0, 0.0), point2=(L, H))

# 2D part, plane strain
prt = modeldb.Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=DEFORMABLE_BODY)
prt.BaseShell(sketch=sk)
del modeldb.sketches['__profile__']

# Partition the single face with a vertical segment at x = partition_location
skp = modeldb.ConstrainedSketch(name='__partline__', sheetSize=0.08, gridSpacing=L/20.0)
skp.Line(point1=(partition_location, 0.0), point2=(partition_location, H))
f = prt.faces[0]
prt.PartitionFaceBySketch(faces=(f,), sketch=skp)
del modeldb.sketches['__partline__']

mat = modeldb.Material(name='NeoHooke')
mat.Density(table=((Material_rho,),))
#Neo-Hooke with volumetric data (filled with C10, D1)
mat.Hyperelastic(materialType=ISOTROPIC, table=((.5 * Material_E/3., 0.0),), type=NEO_HOOKE, testData=OFF, volumetricResponse=VOLUMETRIC_DATA)

# Plane strain solid section with thickness (used for mass; not deforming out-of-plane)
modeldb.HomogeneousSolidSection(name='Section-1', material='NeoHooke', thickness=0.04)
prt.Set(name='AllFaces', faces=prt.faces[:])
prt.SectionAssignment(region=prt.sets['AllFaces'], sectionName='Section-1')


#Seeding and meshing
prt.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
prt.setMeshControls(regions=prt.faces[:], elemShape=QUAD, technique=STRUCTURED)

elemTypes = (ElemType(elemCode=CPE8H, elemLibrary=STANDARD),
             ElemType(elemCode=CPE6H, elemLibrary=STANDARD))
prt.setElementType(elemTypes=elemTypes, regions=(prt.faces[:],))
prt.generateMesh()



a = modeldb.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
inst = a.Instance(name='Part-1-1', part=prt, dependent=ON)

# Edges for top/bottom
# Top edge: y = H
top_edges = inst.edges.getByBoundingBox(xMin=-1e-9, xMax=L+1e-9, yMin=H-1e-9, yMax=H+1e-9)
# Bottom edges split by partition
bot_left  = inst.edges.getByBoundingBox(xMin=-1e-9, xMax=partition_location+1e-9, yMin=-1e-9, yMax=1e-9)
bot_right = inst.edges.getByBoundingBox(xMin=partition_location-1e-9, xMax=L+1e-9, yMin=-1e-9, yMax=1e-9)
bottom_right_left = inst.edges.getByBoundingBox(xMin=-1e-9, xMax=L+1e-9, yMin=-1e-9, yMax=1e-9)
a.Set(name='TopEdge', edges=top_edges)
a.Set(name='BottomLeft', edges=bot_left)
a.Set(name='BottomRight', edges=bot_right)
a.Set(name= 'BottomRightLeft', edges=bottom_right_left)

# Surfaces (for pressure) in 2D are line surfaces built from edges:
a.Surface(name='TopSurf', side1Edges=top_edges)



modeldb.StaticStep(name='Preload_Pressure', previous='Initial', nlgeom=ON,
                     initialInc=0.01,maxInc=0.1, minInc=1e-6, maxNumInc=1000)
modeldb.StaticStep(name='Preload_Shear', previous='Preload_Pressure', nlgeom=ON,
                   initialInc=0.01,maxInc=0.1, minInc=1e-7, maxNumInc=2000)
modeldb.FrequencyStep(name='Eigen', previous='Preload_Shear', numEigen=numEigen, eigensolver=LANCZOS)

modeldb.steps['Preload_Pressure'].setValues(
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION, stabilizationMagnitude=5e-4)
modeldb.steps['Preload_Shear'].setValues(
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION, stabilizationMagnitude=5e-4)

# Pressure on the top line surface
# For plane strain with section thickness specified, using Pressure magnitude = P (Pa) is standard.
modeldb.Pressure(name='Load_TopPressure', createStepName='Preload_Pressure',
                 region=a.surfaces['TopSurf'], magnitude=P)

# Step 1: bottom roller in Y (U2=0)
modeldb.DisplacementBC(name='BC_Bottom_U2_Zero', createStepName='Preload_Pressure',
                       region=Region(edges=a.sets['BottomLeft'].edges + a.sets['BottomRight'].edges),
                       u1=UNSET, u2=0.0)

# Step 2: deactivate roller, then apply revised constraints
modeldb.boundaryConditions['BC_Bottom_U2_Zero'].deactivate('Preload_Shear')

# “stick” on the right bottom in Step 2 (tangential lock):
modeldb.DisplacementBC(name='BC_Stick_BotRight', createStepName='Preload_Shear',
                       region=Region(edges=a.sets['BottomRight'].edges),
                       u1=0.0, u2=0.0)

# Top shear in 2D: U1=DELTAX, U2=0
modeldb.DisplacementBC(name='BC_Top_Shear', createStepName='Preload_Shear',
                       region=Region(edges=a.sets['TopEdge'].edges),
                       u1=DELTAX, u2=0.0)

# Keep/drop pressure in Step 2 as before
if not keep_pressure_in_shear:
    modeldb.loads['Load_TopPressure'].deactivate('Preload_Shear')



# Deformable region for outputs (faces in 2D)
a.Set(name='AllFacesAsm', faces=inst.faces[:])
deform_region = a.sets['AllFacesAsm']

# Field outputs (U, S, E are valid for CPE*H)
for nm, step in (('F-Output-Preload1','Preload_Pressure'),
                 ('F-Output-Preload2','Preload_Shear')):
    if nm in modeldb.fieldOutputRequests: del modeldb.fieldOutputRequests[nm]
    modeldb.FieldOutputRequest(name=nm, createStepName=step,variables=('U','RF','S','E'), region=deform_region)

if 'EigenFields' in modeldb.fieldOutputRequests: del modeldb.fieldOutputRequests['EigenFields']
modeldb.FieldOutputRequest(name='EigenFields', createStepName='Eigen',variables=('U','S','E'), region=deform_region)


# Create a job
mdb.Job(name=jobName, model='Model-1', type=ANALYSIS,
        description='Prestressed eigenfrequency, H=%g mm' % H_mm,
        multiprocessingMode=DEFAULT, numCpus=4, numDomains=4)
mdb.jobs[jobName].submit()
mdb.jobs[jobName].waitForCompletion()
