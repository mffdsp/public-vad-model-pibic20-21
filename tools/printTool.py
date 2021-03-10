def printTool(size=100, actual=0):

    parseValue = (actual+2)*100/size
    print('\rProcessing', parseValue, '%', end="")
    if parseValue % 2 == 0:
        loadChar = ' \\'
    else:
        loadChar = ' /'
    print(loadChar, end="")

def initPrint(str):
    import os

    os.system('cls' if os.name == 'nt' else 'clear')
    print(str)