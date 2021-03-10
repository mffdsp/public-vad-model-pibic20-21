import numpy as np
import scipy as sp
import sys
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.integrate import quad
from tools.functions import runkut4, runkut4_ECG, ECG_McSharry, equacoes_ECG
from tools.read import readFromFile

class Simaan():
    def __init__(self):
        #Simulation Time
        self.start_t = readFromFile('start_t')
        self.passo = readFromFile('step')
        self.end_t = readFromFile('end_t')
 
        #Uses the already created Time Scale
        self.T = list(np.arange(self.start_t,self.end_t,self.passo)) #
        self.n = len(self.T)
        
        #Cardiovascular System
        self.HR = readFromFile('HR')
        self.Emax = readFromFile('Emax') #Amplitude da função elastância
        self.Emin = readFromFile('Emin')
        self.tc = 60/self.HR #Intervalo de tempo referente à duração de um ciclo cardíaco
        self.t_max = 0.2 + 0.15*self.tc; #Tempo máximo da duração de um ciclo cardíaco
        
        self.En = self.elastance(self.T) #Elastância normalizada
        self.E = (self.Emax-self.Emin)*self.En + self.Emin
        self.Ea = self.En*self.Emax
        
        #Cardiovascular System Model Parameters (from Simaan2009)
        
        #Valores de resistores,fontes e bobina 
        
        #Alocando
        self.alocate()
 
        #Initial Conditions
        self.Pao[0] = readFromFile('Pao')
        self.Qa[0]  = readFromFile('Qa')
        self.Vve[0] = readFromFile('Vve')
        self.Pas[0] = readFromFile('Pas')
        self.Pae[0] = readFromFile('Pae')
 
        #x = [x1 x2 x3 x4 x5]
        self.x = np.transpose([self.Pao[0], self.Qa[0], self.Vve[0], self.Pas[0], self.Pae[0]])
 
        #Estado inicial dos diodos
        self.Dm = 0
        self.Da = 0

        #ECG
        self.TH = []
        self.signal_P = 0
        self.signal_Q = 0
        self.signal_R = 0
        self.signal_S = 0
        self.signal_T = 0
        self.bi = [] #Largura da equação Gaussiana
        self.w = 0 #É a velocidade angular da trajetória a medida que se move em torno do ciclo limite

        self.xECG = []
        self.yECG = []
        self.zECG = []
    
    def alocate(self):
        self.Pao = np.zeros_like(self.T)#Pressão na aorta
        self.Qa  = np.zeros_like(self.T) #Fluxo na aorta
        self.Vve = np.zeros_like(self.T) #Volume no ventrículo esquerdo
        self.Pas = np.zeros_like(self.T) #Pressão na aorta sistêmica
        self.Pae = np.zeros_like(self.T) #Pressão no átrio esquerdo
        self.Pve = np.zeros_like(self.T) #Pressão no ventrículo esquerdo
        self.Dm_ = np.zeros_like(self.T) 
        self.Da_ = np.zeros_like(self.T) 

    def elastance(self, t):
        tn = np.asarray(t)%self.tc/self.t_max;
        En = 1.55 * np.power(np.asarray(tn)/.7, 1.9) / (1 + np.power(np.asarray(tn)/.7, 1.9)) / (1 + np.power(np.asarray(tn)/1.17, 21.9))
        return En

    def setECG(self, signal_P=1.2, signal_Q=-5, signal_R=25, signal_S=-7.5, signal_T=0.75):
        self.TH = np.array([-1/3, -1/12, 0, 1/12, 1/2])*np.pi
        self.signal_P = signal_P
        self.signal_Q = signal_Q
        self.signal_R = signal_R
        self.signal_S = signal_S
        self.signal_T = signal_T
        self.bi = [0.25, 0.10, 0.50, 0.10, 0.40]
        self.w = (2*np.pi)/self.tc #É a velocidade angular da trajetória a medida que se move em torno do ciclo limite

        [self.xECG, self.yECG, self.zECG] = ECG_McSharry(self.passo, self.n, self.TH, self.signal_P, self.signal_Q, self.signal_R, self.signal_S, self.signal_T, self.bi, self.w)
        
 
    def att(self, index, vad):
       self.x =  runkut4(model.passo, A, model.x, B,i)

       #CVS VARIABLES
    
       self.Pao[index] =  self.x[0]
       self.Qa[index]  =  self.x[1]
       self.Vve[index] =  self.x[2]
       self.Pas[index] =  self.x[3]
       self.Pae[index] =  self.x[4]
       self.Pve[index] = self.E[index] * (self.Vve[index] - self.Vo)

       #VAD VARIABLES
       vad.Ri[index] = vad.Ri[index] + np.exp(-0.25*self.Vve[index])
    
    def eletricalParameters(self, Rs=1, Rm=0.005, Ra=0.0060, Rc=0.0398, Cae=4.4, Cs=1.33, Cao=0.08, Ls=0.0005, Vo=15, Vd=5):
        self.Rs  = Rs #Resistência sistêmica vascular
        self.Rm  = Rm#Resistência da válvula mitral
        self.Ra  = Ra #Resistência da válvula aorta
        self.Rc  = Rc #Characteristics Resistance
        self.Cae = Cae #Elastância no átrio esquerdo
        self.Cs  = Cs #Elastância sistêmica
        self.Cao = Cao #Elastância na aorta
        self.Ls  = Ls #Inertância sanguinea na aorta
        self.Vo = Vo #Pressão Inicial
        self.Vd = Vd
        self.dV = self.Vo - self.Vd
        self.Pve[0] = self.E[0]* (self.Vve[0] - self.Vo)
    
    def rampFunction(self, index, vad):
        if self.Pae[index] >= self.Pve[index]:
            self.Dm = 1
        else:
            if self.Dm == 1 and self.Da == 0:
                vad.deltaS[index] = 1
            self.Dm = 0
        if self.Pve[index] >= self.Pao[index]:
            self.Da = 1
        else:
            self.Da = 0

