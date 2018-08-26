print('start identity_percent')
from selseq_main import find_tag, SequenceFasta
import os
from selseq_constant import *
import matplotlib.pyplot as plt
import pandas as pd

directory_files = os.listdir(ALNDATA_DIRECTORY)

def decorator_aln():
    def wrap_aln():
        pass

def calculate_identity_percent(aln_file,itself=True):
    '''input file the aligned sequence
       output pd.Series identity_percent
       itself - parametr for calculate identity percent the same sequence
    '''
    aln_file = SequenceFasta(aln_file)
    aln_file.seq_process(strip=False)
    data_persent = pd.Series()    
    
    if itself and len(aln_file.seq_lst) == 1:
        data_persent[find_tag('seq_id',aln_file.name_lst[0])+'and'+find_tag('seq_id',aln_file.name_lst[0])] = 100
    else:
        n=0
        identical = 0
        for seq_id_1 in range(len(aln_file.seq_lst)):            
            n += 1
            for seq_id_2 in range(n,len(aln_file.seq_lst)):                
                for character1, character2 in zip(aln_file.seq_lst[seq_id_1],aln_file.seq_lst[seq_id_2]):
                    if character1 == character2:
                        identical +=1
                data_persent[find_tag('seq_id',aln_file.name_lst[seq_id_1])+'and'+find_tag('seq_id',aln_file.name_lst[seq_id_2])] = identical / len(aln_file.seq_lst[seq_id_1]) * 100
                identical = 0
    return data_persent

def enumeration_aln(directory_files,fun):
    data = pd.Series()
    for file in directory_files:
        if file.endswith('.aln'):
            #print(fun(ALNDATA_DIRECTORY + file))
            data = data.append(calculate_identity_percent(ALNDATA_DIRECTORY + file))
    return data

data = enumeration_aln(directory_files,calculate_identity_percent)
'''
data = pd.Series()
for file in directory_files:
    if file.endswith('.aln'):
        #print (file)
        print(calculate_identity_percent('C:\\Users\\Grin\\Desktop\\selseq\\test_virus\\LIN_Experemen1_alndata\\' + file))
        data = data.append(calculate_identity_percent('C:\\Users\\Grin\\Desktop\\selseq\\test_virus\\LIN_Experemen1_alndata\\' + file))
'''

#plt.hist(data.values, bins=100, alpha=1,color='blue',edgecolor='black')
plt.show()
print('end indentity_persent')

