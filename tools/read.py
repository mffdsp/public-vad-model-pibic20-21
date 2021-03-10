from config.Cmd import Cmd 
import sys

def readFromFile(st):
    with open(Cmd.args.input) as f:
        try:
            values = f.readlines()
            for i in range(len(values)):
                atualStr = values[i]
                if st in atualStr:
                    returnValue = atualStr.split("=", 1)[1]
                    return float(returnValue.rstrip())
            raise Exception("INPUT ERROR: searching from " + st + "in input:\n" + str(values))
            
        except (Exception, BaseException) as e:
            print(e)
            sys.exit('')
