import operator
import functools as ft 
import matplotlib.pyplot as plt
import numpy as np
from tools.functions import runkut4, runkut4_ECG, ECG_McSharry, equacoes_ECG

class VAD():
  def __init__(self, model, deltaPVAD=0.4):
    self.max_vol = 100
    self.Vd_vad = 107
    self.Vd_vad_v = list(np.zeros_like(model.T))
    self.Vd_vad_v[0] = self.Vd_vad
   
    self.V_vad = 110

    self.V_total = 370 #Volume total de sangue

    #Valores humanos para a cânula
    self.DAV_cannula_length_out = 32
    self.DAV_cannula_length_in = 24
    self.DAV_cannula_area = (0.9**2)*np.pi
    self.DAV_rho_blood = 1

    self.DAV_length = 10 #cm
    self.DAV_area = 5 #cm

    #self.Ri = list(np.zeros_like(model.T))
    #self.Ri[0] = np.exp(-0.25*model.Vve[0])
    self.Ri = 0.15
    self.Ri0 = 0.15
    self.Ro = 0.05

    self.Li = 1.75*(9/4)*(self.DAV_rho_blood*self.DAV_cannula_length_in)/(1359.5*self.DAV_cannula_area)
    self.Lo = (self.DAV_rho_blood*self.DAV_cannula_length_out)/(1359.5*self.DAV_cannula_area)

    self.Lp = (9/4)*(self.DAV_rho_blood*self.DAV_length)/(1359.5*self.DAV_area) #A ausência desse indutor não altera o comportamento do DAV


    self.Rp = 0.05
    self.Cp = 2.0
    self.Rd = 0.01
    self.Cd = 4.0


    self.ts_DAV = 0.50 #Tempo de duração da sistole do DAV
    self.t_eject = (model.tc*self.ts_DAV)/model.passo # Ejection time for VAD in fill-to-empty operation
    self.t_eject_c = self.t_eject

    self.Pf = 0

    self.Qi = 0
    self.Qo = 0

    self.deltaS = list(np.zeros_like(model.T))
    self.deltaPVAD = list(np.zeros_like(model.T))
    self.gammad = list(np.zeros_like(model.T))


    self.x = [self.Qi, self.Qo, self.V_vad]
    self.xdot = [0,0,0,0,0,0,0,0]

    self.deltaPVAD = list(np.zeros_like(model.T))
    self.deltaPVAD[0] = deltaPVAD
    self.deltaS[0] = 1
    self.gammad[0] = 1 - self.deltaPVAD[0]

    self.Pt = 0
    self.Px = 0
    self.Pex = 0
    self.Pex2 = 0
    self.Pd = 0
    self.Pee = list(np.zeros_like(model.T))
    self.Pe = list(np.zeros_like(model.T))
    self.Pe[0] = 90
    self.Pee[0] = self.Pe[0]
    self.Pp = (self.V_vad - self.Vd_vad)/self.Cp
    self.Pdrive = 0
    self.Pexternal = self.Pd

    self.Di = 0
    self.Do = 0

    self.open = 1

    self.Vr = model.Cae*model.Pae[0]
    self.Vp = model.Cs*model.Pas[0]
    self.Vao = model.Cao*model.Pao[0]

    self.Pp = (self.V_vad - self.Vd_vad)/self.Cp

    model.Vve[0] = self.V_total - self.V_vad - self.Vr - self.Vp - self.Vao
    model.Pve[0] = (model.Vve[0] - model.Vd)*(model.Ea[0] + model.Emin) - model.Emin*model.dV

    self.delay_in  = -100
    self.delay_out = -100

    self.mode = 0

    #Reference
    self.mPao_ref = 90
    self.mPao_v = list(np.zeros_like(model.T))
    self.mPao_v[0] = self.mPao_ref

    self.mPao_ref_v = list(np.zeros_like(model.T))
    self.mPao_ref_v[0] = self.mPao_ref

    self.aux = []
    self.mPao = 80

    self.ip = 0
    self.ig = 0

    self.ea = 0
    self.e = list(np.zeros_like(model.T))
    self.e[0] = 0
    self.ierror = 0
    self.ierror_max = 100
    self.Kp = 1
    self.Ki = 4

    self.Kd = 0.5
    self.N = 8


    self.b0 = self.Kp*(1+model.n) + self.Ki*(1+model.n) + self.Kd*model.n
    self.b1 = -(self.Kp*(2+model.n) + self.Ki + 2*self.Kd*model.n)
    self.b2 = self.Kp + self.Kd*model.n
    self.a0 = 1 + model.n
    self.a1 = -(2+model.n)
    self.a2 = 1

  def att(self, model, index):

    self.mPao_v[index] = self.mPao
    self.mPao_ref_v[index] = self.mPao_ref

    self.Qi, self.Qo, self.Vvad, model.Pao[index], model.Qa[index], model.Vve[index], model.Pas[index], model.Pae[index] = self.x
    
    self.aux.append(model.Pao[index])

    self.Pp = (self.Vvad - self.Vd_vad)/self.Cp
    self.Px = self.Pp + self.Pex2
    self.Pt = self.Px + self.Rp*(self.Qi - self.Qo) + self.Lp* (self.xdot[0] - self.xdot[1])
    self.Pdrive = self.Pex
    self.Pexternal = self.Pex2
  
  def refereceChange(self, index, mPao_ref=110):
    if index>= 19999 and index < 39999:
        self.mPao_ref = mPao_ref

  def valveManager(self, index, model):
    if self.Di == 0: #IF VALVE IS CLOSE
      if model.Pve[index] >= self.Pt:                
        self.Di = 1
        self.Do = 0
        self.delay_in = index
    else:
      if self.Qi <= 0:
        self.Di = 0
        self.Do = 1

    #Identical Logic used for outlet valve
    if self.Do == 0:
      if self.Pt >= model.Pao[index]:
        self.Do = 1
        self.Di = 0
        self.delay_out = index
    else:
      if self.Qo <= 0: 
        self.Do = 0
        self.Di = 1

  def modeManager(self, index, index0, ejectmode):
    
    if ejectmode == 'fill2empty':
      if self.mode == 0:
        if self.x >= self.max_vol:
            self.mode = 1
            self.Pd = model.Pas[index]
            self.t_eject_c = index + self.t_eject
        else:
            self.Pd = self.Pf
      else:
        if index == self.t_eject_c:
            self.mode = 0
            self.Pd = self.Pf
        else:
            self.Pd = model.Pas[index]
    elif ejectmode == 'Sync':
      if self.mode == 0:
        if self.deltaPVAD[index] == 1:
            self.mode = 1
            self.Pd = self.Pe[index0]
            self.t_eject_c = index + self.t_eject
      else:
        if index == self.t_eject_c:
            self.mode = 0
            self.Pd = self.Pf
  def returnValues(self):
    return np.dot(self.Do, self.Lp)/(self.Lp+self.Lo), np.dot(self.Di, self.Lp)/(self.Lp+self.Li), self.Di/(self.Lp + self.Li), self.Do / (self.Lp+self.Lo), np.dot(self.Di, self.Rp)/(self.Lp+self.Li), (self.Do*self.Rp)/(self.Lp+self.Lo), (np.dot(self.Do,self.Rp)+self.Ri)/(np.dot(self.Do, self.Lp)+self.Li), (np.dot(self.Di,self.Rp)+self.Ro)/(np.dot(self.Di, self.Lp)+self.Lo)
  
  def gammaHandler(self, index, index0):
    self.gammad[index] = 0
    self.deltaPVAD[index] = 1
    
    #Controller
    self.mPao = np.mean(self.aux)
    self.aux = []
    self.e[index0] = self.mPao_ref - self.mPao
    self.ierror = self.ierror + self.e[index0]
    self.de = self.e[index0] - self.e[index0-1]

    if index0 == 0:
      self.Pe[index0]  = (self.b1/self.a0)*self.e[index0]
    elif index0 == 1:
      self.Pe[index0] = (-self.a1/self.a0)*self.Pe[index0-1] + (self.b0/self.a0)*self.e[index0] + (self.b1/self.a0)*self.e[index0-1]
    else:
      self.Pe[index0] = (-self.a1/self.a0)*self.Pe[index0-1] + (-self.a2/self.a0) * self.Pe[index0-2] + (self.b0/self.a0) * self.e[index0] + (self.b1/self.a0)*self.e[index0-1] + (self.b2/self.a0) * self.e[index0-2]

    if self.Pe[index0] > 280:
      self.Pe[index0] = 280
    elif self.Pe[index0]  < 100:
      self.Pe[index0] = 100
    else:
      self.Pe[index0] = self.Pe[index0]

