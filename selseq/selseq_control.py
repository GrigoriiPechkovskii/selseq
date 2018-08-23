#by Pechkovskii Grigorii grigorii.pechkovskii@gmail.com
'''
Module for testing proselseq
'''
import time
import os
from selseq_constant import *
time_block_start = time.time()
time_initial = time.time()

def initial_timecheck_line():
    report_line = '{0: <30}{1: <15}{2: <12}'.format('name','block_min','initial_min') + '\n'
    print(report_line,end='')
    with open(HOME_DIRECTORY + 'time_report.txt','a') as time_report:
        time_report.write(report_line)

def timecheck(name = 'name_check'):
    '''TimeCheck for testing'''
    global time_block_start
    check_block_sec = round(time.time() - time_block_start,3)
    check_block_min = round(check_block_sec / 60,3)
    check_initial_sec = round(time.time() - time_initial,3)
    check_initial_min = round(check_initial_sec / 60,3)       
    report_line = '{0:.<30}{1:.<15}{2:.<5}'.format(name,check_block_min,check_initial_min) + '\n'
    print(report_line,end='')
    with open(HOME_DIRECTORY + 'time_report.txt','a') as time_report:
        time_report.write(report_line)    
    time_block_start = time.time()


import shutil
def del_dir(HOME_DIRECTORY,redata_dir,alndata_dir):
    files = os.listdir(HOME_DIRECTORY)    
    if os.access(redata_dir,os.F_OK):
        shutil.rmtree(redata_dir)
    if os.access(alndata_dir,os.F_OK):
        shutil.rmtree(alndata_dir)
    for file in files:
        if 'aln_' in file:
            shutil.rmtree(HOME_DIRECTORY + file)
            
def check_file_in_table(file_name_faa):
    #для проверки таблицы        
    genom_prok = open(TABLE_CSV)
    assemble_table = []
    for i in genom_prok:
        if '#' not in i:
            assemble_table += [i.split(',')[0].strip()]#индекс

    for file in file_name_faa:
        if file.strip() not in assemble_table:
            print('not in genom_prok',file)
    for assemble in assemble_table:
        if assemble.strip() not in file_name_faa:
            print('not in files',assemble)    
    genom_prok.close()

#тест происходит от rename file

class Control():
    w = 1
    def __init__(self,redata_dir,alndata_dir):
        print('start test')
        self.redata_dir = redata_dir
        self.alndata_dir = alndata_dir
        self.startseq_count = 0
        self.startseq_list =[]
        self.endseq_count = 0
        self.endseq_list = []
        self.startseq_not_in_end = 0
        self.output = []
        
    def _checkStartSeq(self):        
        rename_opened = open(self.redata_dir + 'joint_file')      
        for rename_line in rename_opened:
            if '>' in rename_line:
                self.startseq_count += 1
                self.startseq_list += [rename_line]
        rename_opened.close()

    def _loopEndSeq(self):
        self._checkStartSeq()
        files = os.listdir(self.alndata_dir)
        for aln_file in files:
            aln_opened = open(self.alndata_dir + aln_file)
            for aln_line in aln_opened:
                if '>' in aln_line:
                    self.endseq_count += 1
                    self.endseq_list += [aln_line]                    
                if ('>' in aln_line) and (aln_line not in self.startseq_list):
                    self.startseq_not_in_end += 1 #
                    print('startseq_not_in_end = ', aln_line)
        aln_opened.close()
                    
    def _loopEndSeq2(self):
        self._loopEndSeq()        
        for startseq in self.startseq_list:
            if self.endseq_list.count(startseq) != 1:                    
                    self.output += [str(self.endseq_list.count(startseq)) + 'X' + startseq]
    def group_sum(self):
        pass
    
    def doControl(self):
        self._loopEndSeq2()
        [print(i) for i in self.output]        
        
        print('{0:.<30}{1:.>5}\n{2:.<30}{3:.>5}\n{4:.<30}{5:.>5}\n{6:.<30}{7:.>5}'
              .format('startseq_count', str(self.startseq_count),
              'endseq_count',str(self.endseq_count),
              'startseq_not_in_end', str(self.startseq_not_in_end),
              'difference',str(self.startseq_count - self.endseq_count)))  