import operator
import functools as ft 
import matplotlib.pyplot as plt
import numpy as np
from models import sch, vadModel
from tools.functions import runkut4, runkut4_ECG, ECG_McSharry, equacoes_ECG
from tools.plot import runPlot, runCSV
from tools.printTool import initPrint
from tqdm import tqdm
from init import Data
from config.Cmd import Cmd  
from tools.read import readFromFile


initPrint('Preparing Simulation with ventricular assist device...\nThis can take a while.')

#Init sch model
model = sch.Simaan()
model.setECG()
model.eletricalParameters()

#Init Vad
#vad = vadModel.VAD(model, deltaPVAD=Cmd.args.pvad[0])
vad = vadModel.VAD(model, 
  deltaPVAD=readFromFile('deltaPVAD')
)


vad.x.extend(model.x)

ii = 0
ejectmode = 'Sync'

for i in tqdm(range(model.n-1), desc="Processing", ncols=100, colour='#FE8E86'):

    #ramp function
    model.rampFunction(i, vad)
    vad.valveManager(i, model)
    vad.modeManager(i, ii, ejectmode)
      
    # Calculate de External Bladder pressure from the Drive Pressure]
    A2 = -1/(vad.Rd*vad.Cd)
    B2 = -A2

    vad.Pex = runkut4(model.passo,A2,vad.Pex,B2,vad.Pd)
    vad.Pex2 = vad.Pex + np.dot(vad.xdot[2],0.15)
    vad.Ri = vad.Ri0 + np.exp(-0.25*vad.x[5])
    vad.Ro = 0.00015*np.absolute(vad.Qo) + 0.05

    xo, xi, zi, zo, ri, ro, rri, rro = vad.returnValues()

    delta = (model.Dm/model.Rm)  + (model.Da/model.Ra)
    a11 = xi*ro - rri
    a12 = ri - xi*rro
    a13 = (xi*zo - zi)/vad.Cp
    a16 = zi*(model.Ea[i] + model.Emin)
    a21 = ro - xo*rri
    a22 = xo*ri - rro
    a23 = (zo - xo*zi)/vad.Cp
    a26 = xo*zi*(model.Ea[i] + model.Emin)
    a46 = (model.Da/model.Ra)*(model.Ea[i] + model.Emin)
    a66 =  -delta*(model.Ea[i] + model.Emin)
    a86 = (model.Dm/model.Rm)*(model.Ea[i] + model.Emin)
    a88 = -((1/model.Rs) + model.Dm/model.Rm)

    #A e B são matrizes variantes no tempo
    A = np.zeros((8,8), dtype=float)
    
    A[0] = np.true_divide([a11, a12, a13, -xi*zo, 0, a16, 0,0],(1-xi*xo))
    A[1] = np.true_divide([a21, a22, a23, -zo, 0, a26, 0,0],(1-xi*xo))
    A[2] = [1, -1, 0, 0, 0, 0, 0, 0]
    A[3] = np.true_divide([0,1,0, -model.Da/model.Ra, -1, a46, 0, 0],model.Cao)
    A[4] = np.true_divide([0,0,0, 1, -model.Rc, 0, -1, 0],model.Ls)
    A[5] = [-1,0,0, model.Da/model.Ra, 0, a66, 0, model.Dm/model.Rm]
    A[6] = np.true_divide([0,0,0, 0, 1, 0, -1/model.Rs, 1/model.Rs],model.Cs)
    A[7] = np.true_divide([0,0,0, 0, 0, a86, 1/model.Rs, a88],model.Cae);

    #Matriz B, 5X1
    B = np.zeros((8,4), dtype=float)

    B[0,0] = -A[0][2]
    B[1,0] = -A[1][2]
    B[0,1] = A[0][2]*vad.Cp
    B[1,1] = A[1][2]*vad.Cp

    for i_index in range(8):
      B[i_index,2] = -A[i_index][5]
      B[i_index,3] = -A[i_index][5]*(model.Emin*model.dV)/(model.Ea[i]+model.Emin)
    
    vad.x = runkut4(model.passo,A,vad.x,B,[vad.Vd_vad, vad.Pex2, model.Vd, 1])
    vad.xdot = np.dot(A,vad.x) + np.dot(B,[vad.Vd_vad, vad.Pex2, model.Vd, 1])
    vad.gammad[i+1] = vad.gammad[i] + model.passo/model.tc

    if vad.gammad[i+1] >= 1:
      ii = ii + 1
      vad.gammaHandler(i+1, ii)

    #Reference Change
    vad.refereceChange(i)
    vad.att(model,i+1)

    model.Pve[i + 1] = (model.Vve[i + 1] - model.Vd)*(model.Ea[i + 1] + model.Emin) - model.Emin*model.dV
    vad.Pee[i + 1] = vad.Pe[ii]

runPlot(model, vad)
runCSV(model, vad)