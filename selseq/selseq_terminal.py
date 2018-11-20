#!/usr/bin/env python3

from selseq_high import *
import sys
import argparse

parser = argparse.ArgumentParser()

sys.argv = ['Argparse_learn.py','--test_virus','--directory','./']

parser.add_argument('--test_virus', action='store_true', help='Test for virus sequence')
parser.add_argument('--directory', help='directory where all files are located')
#parser.add_argument('', help='')
#parser.add_argument('', help='')
#parser.add_argument('', help='')

print (parser.parse_args())


if parser.parse_args().directory != None:
    direc = os.path.abspath(parser.parse_args().directory)
    HOME_DIRECTORY = direc + '/'
    print(direc)

if parser.parse_args().test_virus:
    ln = selseq_launcher()
    #ln.selectively()
    ln.total()
    print('end high')

