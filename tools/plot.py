import operator
import functools as ft 
import matplotlib.pyplot as plt
import numpy as np
import os
_path = os.path.abspath(__file__) 
path = './plotResults/'

#PLOT VAD
#plt.plot(model.T, vad.mPao_ref_v)
def runPlot(model, vad):

    print("\rPlotting...", end='')
    plt.figure(figsize=(20,5))
    plt.xlabel('Time(s)')
    plt.ylabel('Pressure(mmHg)')

    plt.plot(model.T, vad.mPao_v, label='Aortic Pressure')
    plt.plot(model.T, vad.mPao_ref_v, label='Aortic Pressure Reference')
    plt.title('Mean Aortic Pressure')
    plt.legend()
    plt.savefig('./plotResults/fig1.png', bbox_inches='tight')

    plt.figure(figsize=(20,5))
    plt.xlabel('Time(s)')
    plt.ylabel('Pressure(mmHg)')

    plt.plot(model.T, model.Pao, label='Aortic Pressure')
    plt.plot(model.T, np.dot(vad.deltaS,np.max(model.Pao)), label='DeltaS')
    plt.plot(model.T, np.dot(vad.deltaPVAD,np.max(model.Pao)), label='DeltaPVAD')
    plt.title('Aortic, Left Ventricule and Left Atrial Pressure')
    plt.legend()
    plt.savefig('./plotResults/fig2.png', bbox_inches='tight')

    plt.figure(figsize=(20,5))
    plt.xlabel('Time(s)')
    plt.ylabel('Volume(ml)')

    plt.plot(model.T, vad.Pee, label='Ejection Pressure')
    plt.title('Pee')
    plt.legend()
    plt.savefig('./plotResults/fig3.png', bbox_inches='tight')

    #plt.plot(model.Vve, vad.Ri)
    #plt.title('Modelo de colapso ventricular; Resistência de entrada do DAV em função do volume ventricular')
    #plt.show()

    #plt.plot(model.Pao)
    #plt.title('Modelo de colapso ventricular; Resistência de entrada do DAV em função do volume ventricular')
    #plt.show()
    #PLOT

    plt.figure(figsize=(20,5))
    plt.plot(model.T, model.Pve)
    plt.plot(model.T, model.Pao)
    plt.title('Pressão no ventriculo esquerdo e aorta')
    plt.savefig('./plotResults/fig4.png', bbox_inches='tight')

    plt.figure(figsize=(20,5))
    plt.plot(model.Vve)
    plt.title('Volume no ventrículo esquerdo')
    plt.savefig(path + 'fig5.png', bbox_inches='tight')

    plt.figure(figsize=(20,5))
    plt.plot(model.Qa)
    plt.title('Fluxo na aorta')
    plt.savefig('./plotResults/fig6.png', bbox_inches='tight')

    plt.figure(figsize=(20,5))
    plt.plot(model.T, model.zECG)
    plt.title('zECG')
    plt.savefig('./plotResults/fig7.png', bbox_inches='tight')
    print('\rPlot is ready at ' + path)


def runCSV(model, vad):
    import csv  
    with open(path + 'output.csv', 'w', encoding='utf-8', newline=''    ) as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['Time (ms)', 'Pao', 'Qa', 'Vve', 'Pas', 'Pae', 'Pve'])
        
        for i in range(len(model.Pao)):
            writer.writerow(np.around([(i + 1)/10 , model.Pao[i],model.Qa[i], model.Vve[i], model.Pas[i], model.Pae[i], model.Pve[i]], 3))