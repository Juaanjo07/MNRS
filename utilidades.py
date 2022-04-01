# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 12:50:00 2022

@author: orlandoaram
"""
from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np

def MomentCurvature(secTag, axialLoad, maxK, numIncr=100):
    # Script tomado de la librería de OpenSeespy de la web
    # secTag es el tag de la sección
    # axialLoad es la carga axial de la sección
    # maxK es la curvatura
    # numIncr es el número de incrementos
    # Define two nodes at (0,0)
    node(1, 0.0, 0.0)
    node(2, 0.0, 0.0)

    # Fix all degrees of freedom except axial and bending
    fix(1, 1, 1, 1)
    fix(2, 0, 1, 0)
    
    # Define element
    #                             tag ndI ndJ  secTag
    element('zeroLengthSection',  1,   1,   2,  secTag)

    # Define constant axial load
    timeSeries('Constant', 1)
    pattern('Plain', 1, 1)
    load(2, axialLoad, 0.0, 0.0)

    # Define analysis parameters
    integrator('LoadControl', 0.0)
    system('SparseGeneral', '-piv')
    test('NormUnbalance', 1e-9, 10)
    numberer('Plain')
    constraints('Plain')
    algorithm('Newton')
    analysis('Static')

    # Do one analysis for constant axial load
    analyze(1)
    loadConst('-time',0.0)

    # Define reference moment
    timeSeries('Linear', 2)
    pattern('Plain',2, 2)
    load(2, 0.0, 0.0, 1.0)

    # Compute curvature increment
    dK = maxK / numIncr

    # Use displacement control at node 2 for section analysis
    integrator('DisplacementControl', 2,3,dK,1,dK,dK)
    
    M = [0]
    curv = [0]
    
    # Do the section analysis
    for i in range(numIncr):
        analyze(1)
        curv.append(nodeDisp(2,3))
        M.append(getTime())
    plt.figure()
    plt.plot(curv,M)
    plt.xlabel('Curvatura')
    plt.ylabel('Momento (kN-m)')
    
    return M,curv

def testMaterial(matTag,displ):
    # wipe()
    
    model('basic','-ndm',2,'-ndf',3)
    # h = getNodeTags()

    node(1,0.0,0.0)
    node(2,0.0,0.0)
    
    fix(1,1,1,1)
    fix(2,1,1,0)
    
    controlnode = 2
    element('zeroLength',1,1,2,'-mat',matTag,'-dir',6)
    
    recorder('Node','-file','MPhi.out','-time','-node',2,'-dof',3,'disp')
    recorder('Element','-file','Moment.out','-time','-ele',1,'force')
    
    ratio = 1/1000
    
    timeSeries('Linear',1)
    pattern('Plain',1,1)
    load(2,0.0,0.0,1.0)
    
    constraints('Plain')
    numberer('Plain')
    system('BandGeneral')
    test('EnergyIncr',1e-6,1000)
    algorithm('Newton')
    
    currentDisp = 0.0
    Disp = [0]
    F = [0]
    nSteps = 1000
    
    for i in displ:
        Dincr = ratio*i/nSteps
        integrator('DisplacementControl',controlnode,3,Dincr)
        analysis('Static')
        
        if Dincr > 0:
            Dmax = Dincr*nSteps
            ok = 0
            while ok == 0 and currentDisp < Dmax:
                ok = analyze(1)
                currentDisp = nodeDisp(controlnode,3)
                F.append(getTime())
                Disp.append(currentDisp)
        elif Dincr < 0:
            Dmax = Dincr*nSteps
            ok = 0
            while ok == 0 and currentDisp > Dmax:
                ok = analyze(1)
                currentDisp = nodeDisp(controlnode,3)
                F.append(getTime())
                Disp.append(currentDisp)
    Fcurr = getTime()
    if ok != 0:
        print('Fallo la convergencia en ',Fcurr)
    else:
        print('Analisis completo')
    
    plt.figure()
    plt.plot(Disp,F)
    plt.xlabel('deformación unitaria (m/m)')
    plt.ylabel('esfuerzo (kPa)')
    return Disp,F
    
    
    
    
    
    
    
    
    
    
    
    