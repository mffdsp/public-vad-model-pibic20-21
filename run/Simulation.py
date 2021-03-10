import os
class Simulation(object):
    @staticmethod
    def startSimulation(mode='vad', inputPath='../input.txt'):
        if mode=='vad':
            os.system('python .\\run\\runsim.py ' + inputPath)
        elif mode=='simaan':
            os.system('python .\\run\\runsim2.py ' + inputPath)
        elif mode=='test':
            os.system('python .\\run\\ejection_variance_value.py ' + inputPath)