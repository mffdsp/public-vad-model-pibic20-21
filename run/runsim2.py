import operator
import functools as ft 
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from models import simaanModel
from tools.printTool import initPrint
from tools.functions import runkut4, runkut4_ECG, ECG_McSharry, equacoes_ECG
import os

initPrint('Preparing Simulation without ventricular assist device...')
model = simaanModel.Simaan()

for i in tqdm(range(model.n-1), desc="Processing", ncols=100, colour='#FE8E86'):
    
    #ramp function
    if model.Pae[i] >= model.Pve[i]:
      model.Dm = 1
    else:
      model.Dm = 0

    if model.Pve[i] >= model.Pao[i]:
      model.Da = 1
    else:
      model.Da = 0
    
    #A e B são matrizes variantes no tempo
    #Matriz A, 5X5
    a13 = (model.Da)/(model.Ra) * (model.E[i])
    a33 = -(((model.Dm)/(model.Rm)) + ((model.Da)/(model.Ra)))*model.E[i]
    a53 = (model.Dm)/(model.Rm) * model.E[i]
    a55 = -((1/(model.Rs))+((model.Dm)/(model.Rm)))

    # x = [ x1 x2 x3 x4 x5 ] 
    A = [[0], [0], [0], [0], [0]];
    A[0] = [-((model.Da)/(model.Ra)), -1, a13, 0, 0]/np.array(model.Cao)
    A[1] = [1, -(model.Rc), 0, -1, 0]/np.array(model.Ls)
    A[2] = [(model.Da)/(model.Ra), 0, a33, 0, (model.Dm)/(model.Rm)]
    A[3] = [0, 1, 0, -1/(model.Rs), 1/(model.Rs)]/np.array(model.Cs);
    A[4] = [0, 0, a53, 1/(model.Rs), a55]/np.array(model.Cae);
   

    #Matriz B, 5X1
    B = [-(((model.Da)/(model.Ra))*model.E[i]*model.Vo)/model.Cao,
         0,
         ((model.Dm)/(model.Rm) + (model.Da)/(model.Ra))*model.E[i]*model.Vo,
         0,
         ( -((model.Dm)/(model.Rm))*model.E[i]*model.Vo)/model.Cae]

    model.x =  runkut4(model.passo, A, model.x, B)
    
    
    #Atualizar as variáveis de estado do próximo ciclo
    model.Pao[i+1] =  model.x[0]
    model.Qa[i+1]  =  model.x[1]
    model.Vve[i+1] =  model.x[2]
    model.Pas[i+1] =  model.x[3]
    model.Pae[i+1] =  model.x[4]
    model.Pve[i+1] = model.E[i+1] * (model.Vve[i+1] - model.Vo)

  
plt.plot(model.Pao)
plt.show()