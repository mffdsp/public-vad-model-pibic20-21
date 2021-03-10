import sys
import argparse

class Cmd:
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mode", choices=['vad', 'simaan', 'test'], default='mode_vad', help="")
    parser.add_argument("-endt", metavar='N', type=float, nargs=1, default=[10], help='an integer for the end time')
    parser.add_argument("-pvad", metavar='N', type=float, nargs=1, default=[0.34], help='an integer for the inicial pression')
    parser.add_argument("input", nargs='?', help="The file containing parameters")

    args = parser.parse_args()