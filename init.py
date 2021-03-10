from config.Cmd import Cmd
import os
import sys 
import run

args = Cmd.args

Data = '<Dat020k02ma0...>'
if __name__ == "__main__":
    try:
        if args.mode == "simaan":
            from run import runsim2
        elif args.mode == "test": 
            from run import ejection_variance_value
        else:
            from run import runsim
    except (Exception, BaseException) as e:
        print(e)
        sys.exit(1)    