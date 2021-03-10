import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.integrate import quad
from tools.functions import runkut4, runkut4_ECG, ECG_McSharry, equacoes_ECG

class Simaan():
    def __init__(self):
        #Simulation Time
        self.start_t = 0 
        self.passo = 0.0001
        self.end_t = 4

        #Uses the already created Time Scale
        self.T = list(np.arange(self.start_t,self.end_t,self.passo)) #
        self.n = len(self.T)
        
        #Cardiovascular System
        self.HR = 60
        self.Emax = 2.5 #Amplitude da função elastância
        self.Emin = 0.06
        self.tc = 60/self.HR #Intervalo de tempo referente à duração de um ciclo cardíaco
        self.t_max = 0.2 + 0.15*self.tc; #Tempo máximo da duração de um ciclo cardíaco
      
        self.E = self.elastance(self.T) #Elastância normalizada
        

        #Cardiovascular System Model Parameters (from Simaan2009)

        #Valores de resistores,fontes e bobina 
        self.Rs  = 1.0000 #Resistência sistêmica vascular
        self.Rm  = 0.0050 #Resistência da válvula mitral
        self.Cae = 4.4000 #Elastância no átrio esquerdo
        self.Ra  = 0.0010 #Resistência da válvula aorta
        self.Rc  = 0.0398 #Characteristics Resistance
        self.Cs  = 1.3300 #Elastância sistêmica
        self.Cao = 0.0800 #Elastância na aorta
        self.Ls  = 0.0005 #Inertância sanguinea na aorta

        self.Vo = 10 #Pressão Inicial

        #PreAllocating
        self.Pao = np.zeros_like(self.T) #Pressão na aorta
        self.Qa  = np.zeros_like(self.T) #Fluxo na aorta
        self.Vve = np.zeros_like(self.T) #Volume no ventrículo esquerdo
        self.Pas = np.zeros_like(self.T) #Pressão na aorta sistêmica
        self.Pae = np.zeros_like(self.T) #Pressão no átrio esquerdo
        self.Pve = np.zeros_like(self.T) #Pressão no ventrículo esquerdo
        self.Dm_ = np.zeros_like(self.T) 
        self.Da_ = np.zeros_like(self.T) 


        #Initial Conditions
        self.Pao[0] = 90
        self.Qa[0]  = 0
        self.Vve[0] = 140 
        self.Pas[0] = 90
        self.Pae[0] = 10

        self.Pve[0] = self.E[0]* (self.Vve[0] - self.Vo) 

        #x = [x1 x2 x3 x4 x5]
        self.x = np.transpose([self.Pao[0], self.Qa[0], self.Vve[0], self.Pas[0], self.Pae[0]])

        #Initial States of diodes
        self.Dm = 0
        self.Da = 0

    def elastance(self, t):

        tn = np.asarray(t)%self.tc/self.t_max;
        En = 1.55 * np.power(np.asarray(tn)/.7, 1.9) / (1 + np.power(np.asarray(tn)/.7, 1.9)) / (1 + np.power(np.asarray(tn)/1.17, 21.9))
        return (self.Emax-self.Emin)*En + self.Emin

