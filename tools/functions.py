import operator
import functools as ft 
import matplotlib.pyplot as plt
import numpy as np

def ECG_McSharry(Ts,N,TH,aiP,aiQ,aiR,aiS,aiT,bi,w):
    T = np.arange(N)

    xECG = list(np.zeros_like(T))
    yECG = list(np.zeros_like(T))
    zECG = list(np.zeros_like(T))

    xECG[0] = -1 
    
    x = [xECG[0], yECG[0], zECG[0]]

    for i in range(N-1):
      x = runkut4_ECG(Ts,x,TH,aiP,aiQ,aiR,aiS,aiT,bi,w)
      xECG[i+1] = x[0]
      yECG[i+1] = x[1]
      zECG[i+1] = x[2]
      
    return [xECG, yECG, zECG]

    
def runkut4_ECG(passo,x,TH,aiP,aiQ,aiR,aiS,aiT,bi,w):

    x = np.array(x)

    xdot = equacoes_ECG(x,TH,aiP,aiQ,aiR,aiS,aiT,bi,w)
    kx1 = passo*np.array(xdot)

    x1 = x + 0.5*np.array(kx1)

    xdot = equacoes_ECG(x1,TH,aiP,aiQ,aiR,aiS,aiT,bi,w)
    kx2 = passo*np.array(xdot)
    
    x1 = x + 0.5*np.array(kx2)
    xdot = equacoes_ECG(x1,TH,aiP,aiQ,aiR,aiS,aiT,bi,w)
    kx3 = passo*np.array(xdot)
 
    x1 = x + kx3
    xdot = equacoes_ECG(x1,TH,aiP,aiQ,aiR,aiS,aiT,bi,w)
    kx4 = passo*np.array(xdot)
 
    value = np.asarray(x + np.array((kx1 + 2*np.array(kx2) + 2*np.array(kx3) + kx4))/6)

    return value

def equacoes_ECG(x,TH,aiP,aiQ,aiR,aiS,aiT,bi,w):
    k = 0

    alpha = 1-np.sqrt(np.power(x[0],2) + np.power(x[1],2))
    th = np.arctan2(x[1], x[0])
    ai = [aiP, aiQ, aiR, aiS, aiT]
    xdot = [[0],[0],[0]]

    for i in range(5):
      k = k - (ai[i]*(th-TH[i])*np.exp( - (np.power(th-TH[i],2)/ (2*np.power(bi[i],2)))))
      

    xdot[0] = alpha*np.array(x[0]) - w*np.array(x[1])
    xdot[1] = w*np.array(x[0]) + alpha*np.array(x[1])
    xdot[2] = -1*np.array(x[2]) + k
  
    return xdot

def runkut4(passo,A,x,B,u=1):
    
    xdot = np.dot(A, x) + np.dot(B,u)
    kx1 = passo*xdot
  
    x1 = x + 0.5*kx1
    xdot = np.dot(A, x1) + np.dot(B,u)
    kx2 = passo*xdot
    
    x1 = x + 0.5*kx2
    xdot = np.dot(A, x1) + np.dot(B,u)
    kx3 = passo*xdot
 
    x1 = x + kx3
    xdot = np.dot(A, x1) + np.dot(B,u)
    kx4 = passo*xdot
 
    value = x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6
 
    return value