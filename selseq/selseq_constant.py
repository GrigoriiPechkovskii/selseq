#!/usr/bin/env python
import os
import sys


NAME_EXP = 'LIN_Experemen1_'
TABLE_CSV = 'group_table_virus.csv'
QUERY_SEQ = 'GCA_InfluenzaA1.faa'
GROUP_TAG = ['Group']
group_for_selection = ['Group']

#if sys.platform == 'linux':
	#Linux
HOME_DIRECTORY  = os.getcwd() + '/'

HOME_DIRECTORY = '/home/grigoriipechkovskii/Desktop/selseq/test_virus/'
TABLE_CSV = HOME_DIRECTORY + TABLE_CSV

REDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'redata/'
ALNDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'alndata/'

#Windows
#HOME_DIRECTORY  = os.getcwd() + '\\'
#REDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'redata\\'
#ALNDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'alndata\\'

HOME_FILES = os.listdir(HOME_DIRECTORY)
