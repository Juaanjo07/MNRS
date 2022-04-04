print("------------ Modelo W6 con 4 elementos en altura y 10 macrofibras -------------")
#------ Modelo W6 con 4 elementos en altura y 10 macrofibras -----------

from openseespy.opensees import *
import math
import numpy as np
import os
import opsvis as opsv
import shutil

print("Starting MVLEM Element Gravity example")
ModelName = 'Muro_W6-MVLEM'
path = 'C:/00_Maestría/01_ProyectoInvestigación/00_OpenSees/Rutinas_Python'
folder = os.path.join(path, ModelName)

if os.path.exists(folder):
    shutil.rmtree(folder)
os.makedirs(folder)

# Create ModelBuilder (with two-dimensions and 3 DOF/node)
wipe()
model('basic', '-ndm',2,'-ndf',3)

# --------------------------------------------
# Set geometry, nodes, boundary conditions
# --------------------------------------------

Hw = 2.4

node(1, 0.0, 0.00*Hw)
node(2, 0.0, 0.25*Hw)
node(3, 0.0, 0.50*Hw)
node(4, 0.0, 0.75*Hw)
node(5, 0.0, 1.00*Hw)

# Set Control Node and DOF
IDctrlDOF = 1
IDctrlNode = 5

fix(1, 1, 1, 1)

# ----- Definición de materiales ----
# ----- Acero de Refuerzo -----
# ---- Acero dúctil ----
Mat_Reinf_B = 1
Mat_Reinf_M = 2

FyB_T = 441900.0
by_flT = 0.00947
FyB_C = FyB_T
by_flC = by_flT
Es = 200000000.0
R0 = 20.0 
a1 = 0.92255
a2 = 0.0015
uniaxialMaterial('SteelMPF', Mat_Reinf_B, FyB_T, FyB_C, Es, by_flT, by_flC, R0, a1, a2)

# ------- Concreto --------
# ----- Concreto inconfinado ------
Mat_Uc_39 = 4

fc_uc = 39200
ec0_uc = -0.00021
Ec_uc = 3900000*math.sqrt(fc_uc/1000)
r_uc = 7
xcrn_uc = 1.039
ft_uc = 0.33*math.sqrt(fc_uc/1000)*1000
et = 0.00008
rt = 1.2
xcrp = 1000.0

uniaxialMaterial('ConcreteCM', Mat_Uc_39, -fc_uc, ec0_uc, Ec_uc, r_uc, xcrn_uc, ft_uc, et, rt, xcrp)

#---- Definición de los MinMax

Mat_MinMax_AsB = 3
Mat_MinMax_CU = 5

Max_As_B = 0.1034
Max_e_Cu = 0.003
uniaxialMaterial('MinMax', Mat_MinMax_AsB, Mat_Reinf_B, '-min', -Max_As_B, '-max', Max_As_B)
uniaxialMaterial('MinMax', Mat_MinMax_CU, Mat_Uc_39,'-min', -Max_e_Cu, '-max', Max_e_Cu)

# ----- Definir cuantias de fibras de muro -----
RowY_BI  = 0.01106
# Cuantia Elemento Transición
RowY_TB = 0.00213
# Cuantia Alma
RowY_YN = 0.00288
# Cuatia Elemento Alma con menos barras
RowY_YC = 0.00192
# Cuantia
RowY_BD =0.02580

# ----- Muro W6 -----
tf = 0.1
Lf = 0.35
tw = 0.1
L = 2.5
Lw = L - tf
Lbc = 0.15
Lmf = (Lw-tw-2*Lbc)/6
m = 10

# ------- Material de Corte ----
Mat_ElasticG = 5

AgWall = Lw*tw+Lf*tf
G = Ec_uc/2.4
GAs = G * AgWall

uniaxialMaterial('Elastic', Mat_ElasticG, GAs)

# ----------- Definición de Elementos MVLEM -----------
Nodes = [1, 2, 3, 4, 5]

thick = [Lf, tw, tw, tw, tw, tw, tw, tw, tw, tw]
width = np.array([tf, Lbc, Lmf, Lmf, Lmf, Lmf, Lmf, Lmf, Lbc, tw],dtype=(float))
rho = np.array([RowY_BI, RowY_TB , RowY_YN, RowY_YN, RowY_YN, RowY_YC, RowY_YN, RowY_YN, RowY_TB, RowY_BD],dtype=(float))
matconcrete = [Mat_Uc_39, Mat_Uc_39, Mat_Uc_39, Mat_Uc_39, Mat_Uc_39, Mat_Uc_39, Mat_Uc_39, Mat_Uc_39, Mat_Uc_39, Mat_Uc_39]
matsteel = [Mat_MinMax_AsB, Mat_MinMax_AsB, Mat_MinMax_AsB, Mat_MinMax_AsB, Mat_MinMax_AsB, Mat_MinMax_AsB, Mat_MinMax_AsB, Mat_MinMax_AsB, Mat_MinMax_AsB, Mat_MinMax_AsB]

for i in range(1, len(Nodes)-1):
    element('MVLEM', i, 0.0, i, i+1, m, 0.4, '-thick', *thick, '-width', *width, '-rho', *rho, '-matConcrete', *matconcrete, '-matSteel', *matsteel,'-matShear', 5)

# Define gravity loads
# --------------------

#  a parameter for the axial load
P = 470.0;  # 10% of axial capacity of columns

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

recorder('Node','-file',str(folder)+'\\Dtop.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
recorder("Node", "-file", str(folder)+"\\example.out", "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")
print('------------------')

#Importar el analisis por gravedad
gravedad()
print("Model built & gravity analysis completed.")
loadConst('-time', 0.0); # deja fija la carga de gravedad

# ----------------------------------------------------
# Start of additional modelling for lateral loads
# ----------------------------------------------------
# Define lateral loads
# --------------------
#PushoverLoads funcion: pushover2(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8)
# o esta: pushover(Dmax,Dincr,IDctrlNode,IDctrlDOF)
F = 1.0 # Reference lateral load

# Set lateral load pattern with a Linear TimeSeries
timeSeries('Linear',2)
pattern('Plain', 2000, 2)
# Create nodal loads at nodes 3 & 4
#    nd    FX  FY  MZ
load(IDctrlNode, F, 0.0, 0.0)

# characteristics of pushover analysis
Dmax=0.01*H  # maximum displacement of pushover. push to a % drift.
Dincr=0.001  # displacement increment. you want this to be small, but not too small to slow analysis
#PO2=pushover(Dmax,Dincr,IDctrlNode,IDctrlDOF)
P02=pushover2(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8)


