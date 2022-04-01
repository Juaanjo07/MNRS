print("-------------------------")
from openseespy.opensees import *
import os
import shutil
import numpy as np
import opsvis as opsv
import matplotlib.pyplot as plt
from analisis  import *
wipe()

print("Starting MVLEM Element Gravity example")
#Crear un directorio para almacenar los datos (mejorar)
modelname="RW2_E17K"
path= "C:\\Users\\Usuario\\OneDrive\\00_EquipodeTrabajoUdeM\\01_RutinasProgramacion\\03_RutinasPython"                            #Ingrese ruta donde va a guardar eventos
folder=os.path.join(path, modelname)
if os.path.exists(folder):
    shutil.rmtree(folder)
os.makedirs(folder)

# Create ModelBuilder (with two-dimensions and 3 DOF/node)
wipe()
model('basic', '-ndm', 2, '-ndf', 3)

# --------------------------------------------
# Set geometry, nodes, boundary conditions
# --------------------------------------------

# Wall Geometry
t = 4.0;    # Wall thickness
H = 144.0;  # Wall height

# Create nodes
#    tag, X, Y
node(1, 0.0, 0.0)
#node(2, 0.0, 1.25)
#node(3, 0.0, 1.75)
#node(4, 0.0, 2.0)
#node(5, 0.0, 5.25)
#node(6, 0.0, 9.0)
#node(7, 0.0, 18.0)
#node(8, 0.0, 30.0)
#node(9, 0.0, 36.0)
#node(10, 0.0, 37.25)
#node(11, 0.0, 37.75)
#node(12, 0.0, 54.0)
#node(13, 0.0, 72.0)
#node(14, 0.0, 90.0)
#node(15, 0.0, 108.0)
#node(16, 0.0, 126.0)
node(17, 0.0, H)

# Fix supports at base of columns
#   tag, DX, DY, RZ
fix(1, 1, 1, 1)

# Set Control Node and DOF
IDctrlNode=17
IDctrlDOF=1

# ------------------------------------------------------------------------
# Define uniaxial materials
# ------------------------------------------------------------------------

# STEEL ...........................................................
# Reinforcing steel

# steel Y boundary
fyYbp=57.3 	    # fy - tension
bybp=0.0185 	# strain hardening - tension
fyYbn=63.0 	    # fy - compression
bybn=0.02 		# strain hardening - compression

# steel Y web
fyYwp=48.8   	# fy - tension
bywp=0.035   	# strain hardening - tension
fyYwn=65.0    	# fy - compression
bywn=0.02 		# strain hardening - compression

# steel misc
Es=29000.0   	# Young's modulus
R0=20.0  		# initial value of curvature parameter
a1=0.925        # curvature degradation parameter
a2=0.0015		# curvature degradation parameter

# Set MVLEM Reinforcing Ratios
rouYb= 0.029333; 	# Y boundary
rouYw= 0.003333; 	# Y web

# Build steel materials
#uniaxialMaterial('SteelMPF', matTag, fyp, fyn, E0, bp, bn, *params, a1=0.0, a2=1.0, a3=0.0, a4=1.0)
# steel Y boundary
uniaxialMaterial('SteelMPF', 1, fyYbp, fyYbn, Es, bybp, bybn,R0,a1,a2) # steel Y boundary
# steel Y web
uniaxialMaterial('SteelMPF', 2, fyYwp, fyYwn, Es, bywp, bywn,R0,a1,a2) # steel Y web

# CONCRETE ........................................................
#uniaxialMaterial('ConcreteCM', matTag, fpcc, epcc, Ec, rc, xcrn, ft, et, rt, xcrp, mon, '-GapClose', GapClose=0)
# Core concrete (confined)
fpcc=-6.9036; 	# peak compressive stress
ec0c=-0.0033;	# strain at peak compressive stress
Ecc=5091.3; 	# Young's modulus
xcrnc=1.0125;	# cracking strain - compression
rc=7.3049;		# shape parameter - compression

# Cover concrete (unconfined)
fpc= -6.2; 		# peak compressive stress
ec0= -0.0021;	# strain at peak compressive stress
ft= 0.295;		# peak tensile stress
et= 0.00008;		# strain at peak tensile stress
Ec= 4500; 		# Young's modulus
xcrnu= 1.039;	# cracking strain - compression
xcrp= 10000;		# cracking strain - tension
ru= 7;			# shape parameter - compression
rt= 1.2;			# shape parameter - tension
# Build concrete materials
# confined concrete
uniaxialMaterial('ConcreteCM', 3, fpcc, ec0c, Ecc, rc, xcrnc, ft, et, rt, xcrp, '-GapClose',1)
# unconfined concrete
uniaxialMaterial('ConcreteCM', 4, fpc, ec0, Ec, ru, xcrnu, ft, et, rt, xcrp, '-GapClose', 1)

# SHEAR ........................................................
# uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
# NOTE: large shear stiffness assigned since only flexural response
Ashweb= 192;					# Gross area of the wall cross section
G= 1875000;					# Shear Modulus
GAs= G*Ashweb; 	# Shear Stiffness

# Build shear material
uniaxialMaterial('Elastic', 5, GAs)

# ------------------------------
#  Define MVLEM elements
# ------------------------------
#element('MVLEM', eleTag, Dens, *eleNodes, m, c, '-thick', *thick, '-width', *widths, '-rho', *rho, '-matConcrete', *matConcreteTags, '-matSteel', *matSteelTags, '-matShear', matShearTag)
#Definición de parametros
thick=np.array([t,t,t,t,t,t,t,t],dtype=(float))
width=np.array([7.5, 1.5, 7.5, 7.5, 7.5, 7.5, 1.5, 7.5],dtype=(float))
rho=np.array([rouYb, 0.0, rouYw, rouYw, rouYw, rouYw, 0.0, rouYb],dtype=(float))
matconcrete=[3, 4, 4, 4, 4, 4, 4, 3]
matsteel=[1, 2, 2, 2, 2, 2, 2, 1]

element('MVLEM', 1, 0.0, 1, 17, 8, 0.4, '-thick', *thick, '-width',*width, '-rho', *rho, '-matConcrete', *matconcrete, '-matSteel', *matsteel,'-matShear', 5)


# Define gravity loads
# --------------------

#  a parameter for the axial load
P = 85.0;  # 10% of axial capacity of columns

# Create a Plain load pattern with a Linear TimeSeries
timeSeries('Linear',1)
pattern('Plain', 1, 1)
# Create nodal loads at node 17
#    nd  FX,  FY, MZ
load(IDctrlNode, 0.0, -P, 0.0)
print('------------------')
print("Axial load applied")

# ------------------------------
# End of model generation
# ------------------------------
#Reporte en los nodos
recorder('Node','-file',str(folder)+'\\Dtop.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
recorder("Node", "-file", str(folder)+"\\example.out", "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")

print('------------------')
#Importar el analisis por gravedad
gravedad()
print("Model built & gravity analysis completed.")
loadConst('-time', 0.0); # deja fija la carga de gravedad
# Print out the state of the nodes
#u17= nodeDisp(17)
#r1= nodeReaction(1)
print('--------------------------')

# ----------------------------------------------------
# Start of additional modelling for lateral loads
# ----------------------------------------------------
# Define lateral loads
# --------------------
#PushoverLoads funcion: pushover2(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8)
# o esta: pushover(Dmax,Dincr,IDctrlNode,IDctrlDOF)
F=1.0 # Reference lateral load

# Set lateral load pattern with a Linear TimeSeries
timeSeries('Linear',2)
pattern('Plain', 2000, 2)
# Create nodal loads at nodes 3 & 4
#    nd    FX  FY  MZ
load(IDctrlNode, F, 0.0, 0.0)

# characteristics of pushover analysis
Dmax=0.02*H  # maximum displacement of pushover. push to a % drift.
Dincr=0.001  # displacement increment. you want this to be small, but not too small to slow analysis
#PO2=pushover(Dmax,Dincr,IDctrlNode,IDctrlDOF)
P02=pushover2(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8)






























#Imprime el analisis del elemento
#printModel()
# 1. plot model with tag lebels
#opsv.plot_model(node_labels=1, element_labels=1,offset_nd_label=False,nodes_only=False, fmt_model='b.-')
# 2. plot supports and loads
#plt.figure()
#opsv.plot_supports_and_loads_2d()