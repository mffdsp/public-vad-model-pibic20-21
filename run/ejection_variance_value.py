import operator
import functools as ft 
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from models import sch, vadModel
import os
from tools.functions import runkut4, runkut4_ECG, ECG_McSharry, equacoes_ECG
from tools.printTool import initPrint


initPrint('Preparing tests...')

ejection_variance = list(np.arange(0,1,0.1))
ejection_variance_values = [(0,0)]

for j in tqdm(range(len(ejection_variance)-1), desc="Processing" , ncols=100, colour='#149414'):
    model = sch.Simaan(end_t=10)
    model.setECG()
    model.eletricalParameters()
    #model.setECG()
    vad = vadModel.VAD(model, deltaPVAD=ejection_variance[j])
    vad.x.extend(model.x)

    print(" -> test for ejection_variance =", ejection_variance[j])
    
    ii = 0
    for i in range(model.n-1):
        #ramp function
        if model.Pae[i] >= model.Pve[i]:
            model.Dm = 1
        else:
            if model.Dm == 1 and model.Da == 0:
                vad.deltaS[i] = 1
                model.Dm = 0

        if model.Pve[i] >= model.Pao[i]:
            model.Da = 1
        else:
            model.Da = 0


        if vad.Di == 0: #IF VALVE IS CLOSE
            if model.Pve[i] >= vad.Pt:                
                vad.Di = 1
                vad.Do = 0
                vad.delay_in = i
        else:
            if vad.Qi <= 0:
                vad.Di = 0
                vad.Do = 1

        #Identical Logic used for outlet valve
        if vad.Do == 0:
            if vad.Pt >= model.Pao[i]:
                vad.Do = 1
                vad.Di = 0
                vad.delay_out = i
        else:
            if vad.Qo <= 0: 
                vad.Do = 0
                vad.Di = 1

        
        ejectmode = 'Sync'

        if ejectmode == 'fill2empty':
            if vad.mode == 0:
                if vad.x >= vad.max_vol:
                    vad.mode = 1
                    vad.Pd = model.Pas[i]
                    vad.t_eject_c = i + vad.t_eject
                else:
                    vad.Pd = vad.Pf
            else:
                if i == vad.t_eject_c:
                    vad.mode = 0
                    vad.Pd = vad.Pf
                else:
                    vad.Pd = model.Pas[i]
        elif ejectmode == 'Sync':
            if vad.mode == 0:
                if vad.deltaPVAD[i] == 1:
                    vad.mode = 1
                    vad.Pd = vad.Pe[ii]
                    vad.t_eject_c = i + vad.t_eject
            else:
                if i == vad.t_eject_c:
                    vad.mode = 0
                    vad.Pd = vad.Pf
        
        # Calculate de External Bladder pressure from the Drive Pressure]
        A2 = -1/(vad.Rd*vad.Cd)
        B2 = -A2

        vad.Pex = runkut4(model.passo,A2,vad.Pex,B2,vad.Pd)


        vad.Pex2 = vad.Pex + np.dot(vad.xdot[2],0.15)
        vad.Ri = vad.Ri0 + np.exp(-0.25*vad.x[5])
        vad.Ro = 0.00015*np.absolute(vad.Qo) + 0.05

        xo = np.dot(vad.Do, vad.Lp)/(vad.Lp+vad.Lo)
        xi = np.dot(vad.Di, vad.Lp)/(vad.Lp+vad.Li)
        zi = vad.Di/(vad.Lp + vad.Li)
        zo = vad.Do / (vad.Lp+vad.Lo)
        ri = np.dot(vad.Di, vad.Rp)/(vad.Lp+vad.Li)
        ro = (vad.Do*vad.Rp)/(vad.Lp+vad.Lo)
        rri = (np.dot(vad.Do,vad.Rp)+vad.Ri)/(np.dot(vad.Do, vad.Lp)+vad.Li)
        rro = (np.dot(vad.Di,vad.Rp)+vad.Ro)/(np.dot(vad.Di, vad.Lp)+vad.Lo)

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
        #Matriz A, 5X5
        
        # x = [ x1 x2 x3 x4 x5 ] 
        A = np.zeros((8,8), dtype=float)
        #A = [[0],[0],[0],[0], [0], [0], [0], [0]]
        
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
            vad.gammad[i+1] = 0
            vad.deltaPVAD[i+1] = 1
            
            #Controller
            ii = ii + 1

            vad.mPao = np.mean(vad.aux)

            vad.aux = []


            vad.e[ii] = vad.mPao_ref - vad.mPao
            vad.ierror = vad.ierror + vad.e[ii]
            vad.de = vad.e[ii] - vad.e[ii-1]

            if ii == 0:
                vad.Pe[ii]  = (vad.b1/vad.a0)*vad.e[ii]
            elif ii == 1:
                vad.Pe[ii] = (-vad.a1/vad.a0)*vad.Pe[ii-1] + (vad.b0/vad.a0)*vad.e[ii] + (vad.b1/vad.a0)*vad.e[ii-1]
            else:
                vad.Pe[ii] = (-vad.a1/vad.a0)*vad.Pe[ii-1] + (-vad.a2/vad.a0) * vad.Pe[ii-2] + (vad.b0/vad.a0) * vad.e[ii] + (vad.b1/vad.a0)*vad.e[ii-1] + (vad.b2/vad.a0) * vad.e[ii-2]

            if vad.Pe[ii] > 280:
                vad.Pe[ii] = 280
            elif vad.Pe[ii]  < 100:
                vad.Pe[ii] = 100
            else:
                vad.Pe[ii] = vad.Pe[ii]

            #Reference Change
            if i>= 19999 and i < 39999:
                vad.mPao_ref = 110


        vad.mPao_v[i+1] = vad.mPao
        vad.mPao_ref_v[i+1] = vad.mPao_ref

        vad.Qi = vad.x[0]
        vad.Qo = vad.x[1]
        vad.Vvad = vad.x[2]
        model.Pao[i+1] = vad.x[3]
        model.Qa[i+1] = vad.x[4]
        model.Vve[i+1] = vad.x[5]
        model.Pas[i+1] = vad.x[6]
        model.Pae[i+1] = vad.x[7]
        
        vad.aux.append(model.Pao[i+1])

        model.Pve[i+1] = (model.Vve[i+1] - model.Vd)*(model.Ea[i+1] + model.Emin) - model.Emin*model.dV

        vad.Pee[i+1] = vad.Pe[ii]

        vad.Pp = (vad.Vvad - vad.Vd_vad)/vad.Cp
        vad.Px = vad.Pp + vad.Pex2
        vad.Pt = vad.Px + vad.Rp*(vad.Qi - vad.Qo) + vad.Lp* (vad.xdot[0] - vad.xdot[1])
        vad.Pdrive = vad.Pex
        vad.Pexternal = vad.Pex2

        #Atualizar as variáveis de estado do próximo ciclo
        #model.att(i + 1,vad);
        if(np.absolute(vad.mPao_v[i+1] - vad.mPao_ref_v[i+1]) <= 1 and vad.mPao_ref>=110):
            ejection_variance_values.append((round(model.T[i],3),round(ejection_variance[j],2)))
            break    

ejection_variance_values = sorted(ejection_variance_values, key=lambda tupla: tupla[0])
np.savetxt("./testResults/ejection_results_110.csv", ejection_variance_values, delimiter=",", fmt='%.3f')