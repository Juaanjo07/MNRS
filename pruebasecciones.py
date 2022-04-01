# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 16:59:25 2022

@author: Orlando
"""
from openseespy.opensees import *
# hay que llamar al opsvis para poder plotear la visualización
import opsvis as opsv
# esta es una librería estándar para plotear otras cosasy poder crear figuras.
import matplotlib.pyplot as plt
#ploteador de python
import utilidades as ut


wipe()

# seccion
# Definición de material 
E = 24000000.0 # 24GPa pasados a KPa para hacerlo consistente
#uniaxialMaterial('Elastic',1,E)
fc = 28000.0
ec = 2*fc/E
fcu = 0.2*fc
ecu = 0.006
 
k=1.3
fcc=28000.0*k
ecc= 2*fcc/E
fucc=0.2*fcc
eucc=0.02
 
Fy=420000.0
Es=210000000.0
uniaxialMaterial('Concrete01', 2, fc, ec, fcu, ecu)
uniaxialMaterial('Concrete01', 1, fcc, ecc, fucc, eucc)
uniaxialMaterial('Steel01', 3, Fy, Es, 0.01)

# Seccion

Bcol = 0.30
Hcol = Bcol
c = 0.05  # recubrimiento 

# creación de la sección de fibra
y1col = Hcol/2.0
z1col = Bcol/2.0

y2col = 0.5*(Hcol-2*c)/3.0

nFibZ = 1
nFibZcore= 10
nFib = 20
nFibCover, nFibCore = 3, 16
As4 = 0.000127

sec30x30 = 1

col30x30 = [['section', 'Fiber', sec30x30, '-GJ', 1.0e6],
             ['patch', 'rect', 1, nFibCore, nFibZcore, c-y1col, c-z1col, y1col-c, z1col-c],
             ['patch', 'rect', 2, nFib, nFibZ, -y1col, -z1col, y1col, c-z1col],
             ['patch', 'rect', 2, nFib, nFibZ, -y1col, z1col-c, y1col, z1col],
             ['patch', 'rect', 2, nFibCover, nFibZ, -y1col, c-z1col, c-y1col, z1col-c],
             ['patch', 'rect', 2, nFibCover, nFibZ, y1col-c, c-z1col, y1col, z1col-c],
             ['layer', 'straight', 3, 4, As4, y1col-c, z1col-c, y1col-c, c-z1col],
             ['layer', 'straight', 3, 2, As4, y2col, z1col-c, y2col, c-z1col],
             ['layer', 'straight', 3, 2, As4, -y2col, z1col-c, -y2col, c-z1col],
             ['layer', 'straight', 3, 4, As4, c-y1col, z1col-c, c-y1col, c-z1col]]


# matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
# opsv.plot_fiber_section(col30x30, matcolor=matcolor)
# plt.axis('equal')
# opsv.fib_sec_list_to_cmds(col30x30)

#sección 2----------------------------------------------------------------------
Bcol = 0.30
Hcol = 0.30
c = 0.05  # recubrimiento 

# creación de la sección de fibra
y1col = Hcol/2.0
z1col = Bcol/2.0

y2col = 0.5*(Hcol-2*c)/3.0

nFibZ = 1
nFibZcore= 10
nFib = 20
nFibCover, nFibCore = 3, 16
As4 = 0.000127
As5 = 0.0002
sec30x40 = 2

vig30x40 = [['section', 'Fiber', sec30x40, '-GJ', 1.0e6],
             ['patch', 'rect', 1, nFibCore, nFibZcore, c-y1col, c-z1col, y1col-c, z1col-c],
             ['patch', 'rect', 2, nFib, nFibZ, -y1col, -z1col, y1col, c-z1col],
             ['patch', 'rect', 2, nFib, nFibZ, -y1col, z1col-c, y1col, z1col],
             ['patch', 'rect', 2, nFibCover, nFibZ, -y1col, c-z1col, c-y1col, z1col-c],
             ['patch', 'rect', 2, nFibCover, nFibZ, y1col-c, c-z1col, y1col, z1col-c],
             ['layer', 'straight', 3, 3, As5, y1col-c, z1col-c, y1col-c, c-z1col],
             ['layer', 'straight', 3, 3, As5, c-y1col, z1col-c, c-y1col, c-z1col]]

# opsv.plot_fiber_section(vig30x40, matcolor=matcolor)
# plt.axis('equal')
# opsv.fib_sec_list_to_cmds(vig30x40)

# ut.MomentCurvature(sec30x40, -20.0, 0.15)


uniaxialMaterial('Concrete02', 5, -fcc, -ecc, -fucc, -eucc)
disp = [0, 10.0, -10.0, 20.0, -20.0, 50.0, -50.0, 60.0, -120.0, 120.0]

uniaxialMaterial('MinMax', 1000, 3, '-min', -0.02, '-max', 0.02)

ut.testMaterial(1000,disp)

# disp = [10, -10, 20.0, -20.0, 30.0, -30.0]
# ut.testMaterial(3,disp)


