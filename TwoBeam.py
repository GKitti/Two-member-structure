import opensees as ops
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

## SI Units
mm = 0.001
kN = 1000
GPa = 1e9
kNm = kN

defaultK=0.0001

class TwoBeamStructure:
    def __init__(self, span=1, height=0.03, secD=80, secT=3, NN=5, springK=defaultK):
        self.span = span #m
        self.height = height #m
        self.secD = secD * mm
        self.secT = secT * mm
        self.springK = springK * kNm #stiffness of the joint in kNm/rad
        self.secA=0
        self.secI=0
        self.NN=NN

def CHSSection(secID,matID,D,t,nc,nr):
    # create a circular hollow cross-section and calcukate its basic cross-section properties
    R = D/2
    r = R - t
    A = math.pi * (R ** 2 - r ** 2)
    I = math.pi * (R ** 4 - r ** 4) / 4
    ops.section('Fiber',secID)
    ops.patch('circ',matID,nc,nr,0.,0.,r,R,0.,360.)
    return A,I

def Numerical(TBS, K=defaultK):

    TBS.springK=K * kNm;

    ops.wipe()
    ops.model('BasicBuilder','-ndm',2,'-ndf',3)

    display=1

    # 1. CREATE MODEL
    # material
    matID=100
    Es=210.*GPa
    ops.uniaxialMaterial('Elastic',matID,Es)

    # cross-section
    CHS=1
    TBS.secA, TBS.secI = CHSSection(CHS,matID,TBS.secD,TBS.secT,8,1)

    # create nodes
    nn=TBS.NN
    nodeNb=nn*2+1 #total number of nodes
    nodeIDA=np.zeros(nodeNb,int) #node IDs
    lA=np.zeros(nodeNb) # y coordinates of nodes
    hA=np.zeros(nodeNb) # x coordinates of nodes

    for i in range(nn+1):
        nodeIDA[i]=100+i
        lA[i]=TBS.span/(nodeNb-1)*i
        hA[i]=TBS.height/(nodeNb-1)*i*2
    for i in range(nn):
        nodeIDA[nodeNb-1-i]=200+i
        lA[nodeNb-1-i]=TBS.span-(TBS.span/(nodeNb-1)*i)
        hA[nodeNb-1-i]=TBS.height/(nodeNb-1)*i*2
    nodeIDA[nn]=0

    for i in range(nodeNb):
        ops.node(int(nodeIDA[i]),lA[i],hA[i])

    # plt.plot(lA, hA,'ro')

    # create zeroLength element nodes
    ops.node(1,TBS.span/2.0,TBS.height)
    ops.node(2,TBS.span/2.0,TBS.height)

    # end supports
    ops.fix(100,1,1,0)
    ops.fix(200,1,1,0)

    # transformations
    CorotTR=1
    ops.geomTransf('Corotational',CorotTR)

    # integration points
    gauss=9
    ops.beamIntegration('Lobatto',1,CHS,gauss)

    # create elements
    for i in range(nn - 1):
        ops.element('forceBeamColumn', 1000 + i, 100 + i, 100 + i + 1, CorotTR, 1)
    ops.element('forceBeamColumn', 1000 + nn - 1, 100 + nn - 1, 1, CorotTR, 1)
    for i in range(nn - 1):
        ops.element('forceBeamColumn', 2000 + i, 200 + i, 200 + i + 1, CorotTR, 1)
    ops.element('forceBeamColumn', 2000 + nn - 1, 200 + nn - 1, 2, CorotTR, 1)

    # zeroLength elements for the hinges
    ops.uniaxialMaterial('Elastic', 9, TBS.springK) # define the characteristics of the hinge
    ops.element('zeroLength', 1, 0, 1, '-mat', 9, '-dir', 6)
    ops.equalDOF(0, 1, 1, 2)
    ops.element('zeroLength', 2, 0, 2, '-mat', 9, '-dir', 6)
    ops.equalDOF(0, 2, 1, 2)

    #2. ANALYSIS
    # dislacement-control to find post-critical behaviour
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain',1,1)
    ops.load(0,0.,-1.0*kN,0.)

    IDctrlNode = 0
    IDctrlDOF = 2
    DispMax = TBS.height * -2.
    Steps=50
    DispIncr = DispMax / Steps
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('EnergyIncr', 1.e-10, 100)
    ops.algorithm('NewtonLineSearch')

    ops.integrator('DisplacementControl', IDctrlNode, IDctrlDOF, DispIncr)

    ops.analysis('Static')

    NDisp = np.zeros(Steps+1)
    EForce = np.zeros(Steps+1)

    for i in range(1,Steps+1):
        ok = 0
        ok = ops.analyze(1)
        if ok == 0:
            # print('analysis completed successfully')
            NDisp[i] = - ops.nodeDisp(0, 2)/mm  #mm displacement of the mid-node
            EForce[i] = ops.eleResponse(1000, 'globalForce')[1] / kN * 2 #kN total point load P

        else:
            ok = 1
            print('analysis failed to converge')

    # set finish [clock milliseconds]
    # set runTime [expr ($finish-$start)/1000.0]

    #if ok == 0:
    #    print('analysis completed successfully in ... seconds')

    return NDisp,EForce

def Analytical(TBS, K=defaultK):
    steps=100

    K=K*kNm

    L = TBS.span
    H = TBS.height
    EA = TBS.secA * 210 * GPa
    EI = TBS.secI * 210 * GPa
    tT0 = H / (L / 2)
    T0 = math.atan(tT0)
    cT0 = math.cos(T0)
    sT0 = math.sin(T0)

    Ppin = np.zeros(steps)
    Psr = np.zeros(steps)
    Pr = np.zeros(steps)
    PpinI = np.zeros(steps)
    PrI = np.zeros(steps)
    dA = np.zeros(steps)

    for i in range(steps):
        T=T0-2*T0/(steps-1)*i
        cT = math.cos(T)
        sT = math.sin(T)
        tT = math.tan(T)
        dA[i] = (tT0-tT) * L / 2
        #'exact'
        Ppin[i] = 2 * EA * sT * (1- cT0 / cT)
        Psr[i]    = 4 * K / L * cT ** 2 * (T0 - T)
        Psr[i]    = Ppin[i] + Psr[i]/ (1+(K*L*cT**2)/(6*EI*cT0**3))
        Pr[i] = Ppin[i] + math.pi ** 4 * EI / (4 * L ** 2) * cT0 ** 5 * (tT0 - tT)
        #Pr[i]   = Ppin[i] + math.pi ** 4 * EI / (8* L ** 2) * cT**6 * (tT0-tT) *(2*cT-3*tT*sT+5*sT*tT0)
        #first order
        PpinI[i] = EA * T*(T**2-T0**2)
        PrI[i] = PpinI[i] + 24 * EI / (L ** 2) * (T0 - T)
        # output in mm and kN
        dA[i] = dA[i] / mm;
        Ppin[i] = Ppin[i]/ kN;
        Psr[i]  = Psr[i] / kN;
        Pr[i]   = Pr[i]  / kN;
        PpinI[i]= PpinI[i] / kN;
        PrI[i]  = PrI[i] / kN;
    return Ppin, Psr, Pr, dA, PpinI, PrI
