#!/usr/bin/env python3
import os
import sys


NAME_EXP = 'LIN_Experemen1_'
TABLE_CSV = 'group_table_virus.csv'
QUERY_SEQ = 'InfluenzaA1.faa'
QUERY_SEQ_LIST = ['tbl_InfluenzaA1.csv','tbl_InfluenzaA3.csv']
GROUP_TAG = ['Group1','Group2']
group_for_selection = ['Group1','Group2']
PERCENT_RANGE  = 5
PERSENT_THRESHHOLD = 70

if sys.platform == 'linux':
	#Linux
	HOME_DIRECTORY  = os.getcwd() + '/'

	HOME_DIRECTORY = '/home/grigoriipechkovskii/Desktop/selseq/test_virus/'	

	REDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'redata/'
	ALNDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'alndata/'
	TABLE_CSV = HOME_DIRECTORY + TABLE_CSV
	slash = '/'

if sys.platform == 'win32':
	#Windows
	HOME_DIRECTORY  = os.getcwd() + '\\'
	HOME_DIRECTORY = 'C:\\Users\\Grin\\Desktop\\selseq\\test_virus\\'
	REDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'redata\\'
	ALNDATA_DIRECTORY = HOME_DIRECTORY + NAME_EXP +'alndata\\'
	TABLE_CSV = HOME_DIRECTORY + TABLE_CSV
	slash = '\\'
HOME_FILES = os.listdir(HOME_DIRECTORY)

ASSEMBLE_FILES = [file for file in HOME_FILES if file.endswith('.faa')]