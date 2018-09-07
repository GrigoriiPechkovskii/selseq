#!/usr/bin/env python3
import os
import sys


NAME_EXP = 'LIN_Experemen1_'
TABLE_CSV = 'group_table_virus.csv'
QUERY_SEQ = 'InfluenzaA1.faa'
GROUP_TAG = ['Group']
group_for_selection = ['Group']
PERCENT_RANGE  = 5
PERSENT_THRESHHOLD = 70

if sys.platform == 'linux':
	#Linux
	HOME_DIRECTORY  = os.getcwd() + '/'

	HOME_DIRECTORY = '/home/grigoriipechkovskii/Desktop/selseq/test_virus/'	

	REDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'redata/'
	ALNDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'alndata/'
	TABLE_CSV = HOME_DIRECTORY + TABLE_CSV

if sys.platform == 'win32':
	#Windows
	HOME_DIRECTORY  = os.getcwd() + '\\'
	HOME_DIRECTORY = 'C:\\Users\\Grin\\Desktop\\selseq\\test_virus\\'
	REDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'redata\\'
	ALNDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'alndata\\'
	TABLE_CSV = HOME_DIRECTORY + TABLE_CSV

HOME_FILES = os.listdir(HOME_DIRECTORY)
