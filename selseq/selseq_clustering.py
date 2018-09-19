print('start identity_percent')
#from selseq_main import *
#from selseq_constant import *

#from selseq_extra import *

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
import re
import subprocess
import shutil
from selseq_main import *

def clustering_kmeans_aln(aln_file,itself=True):
    '''input file the aligned sequence
       output clustering by kmeans files       
    '''
    aln_file = calculate_identity_percent(aln_file,itself=True)

    if any((aln_file.identity_matrix<60).any()):                
        kmeans = KMeans(n_clusters=2)
        kmeans.fit(aln_file.identity_matrix)
        y_kmeans = kmeans.predict(aln_file.identity_matrix)      
                            
        for kmeans_index in range(len(y_kmeans)):
            name_aln_file_0 = aln_file.file_dir[0:-4] + '_0'
            name_aln_file_1 = aln_file.file_dir[0:-4] + '_1'
            if y_kmeans[kmeans_index] == 0:
                with open(name_aln_file_0,'a') as aln_clustered:                        
                    aln_clustered.write(aln_file.name_lst[kmeans_index] + aln_file.seq_lst[kmeans_index].replace('-','').replace('\n','') + '\n')
            
            if y_kmeans[kmeans_index] == 1:
                with open(name_aln_file_1,'a') as aln_clustered:                        
                    aln_clustered.write(aln_file.name_lst[kmeans_index] + aln_file.seq_lst[kmeans_index].replace('-','').replace('\n','') + '\n')
        
        subprocess.call('muscle ' + '-in ' +name_aln_file_0 + ' -out ' + name_aln_file_0 + '.aln 2>' + HOME_DIRECTORY + '111',shell = True)
        subprocess.call('muscle ' + '-in ' +name_aln_file_1 + ' -out ' + name_aln_file_1 + '.aln 2>' + HOME_DIRECTORY + '111',shell = True)

        clustering_kmeans_aln(name_aln_file_0 + '.aln',itself=True)
        clustering_kmeans_aln(name_aln_file_1 + '.aln',itself=True)

        os.remove(name_aln_file_0)
        os.remove(name_aln_file_1)
        os.remove(aln_file.file_dir)
    
    else:        
        return aln_file


def calculate_identity_percent(aln_file,itself=True):
    '''input file the aligned sequence
       output SequenceFasta with identity_percent and identity_matrix
       itself - parametr for calculate identity percent the alone sequence
    '''
    aln_file = SequenceFasta(aln_file)
    aln_file.seq_process(strip=False)
    data_persent = pd.Series()
    identity_matrix = pd.DataFrame()    

    if itself and len(aln_file.seq_lst) == 1:        
        data_persent[find_tag('seq_id',aln_file.name_lst[0])+'and'+find_tag('seq_id',aln_file.name_lst[0])] = 110
        aln_file.data_persent = data_persent
        identity_matrix = pd.DataFrame([])
        aln_file.identity_matrix = identity_matrix
        return aln_file

    else:
        name_lst_seq_id = []
        for name_seq in aln_file.name_lst:
            name_lst_seq_id.append(find_tag('seq_id',name_seq))
        array_100 = np.zeros((len(aln_file.name_lst), len(aln_file.name_lst))) +100
        identity_matrix = pd.DataFrame(array_100,columns=name_lst_seq_id,index=name_lst_seq_id)
        
        n=0
        identical = 0
        for seq_id_1 in range(len(aln_file.seq_lst)):            
            n += 1
            for seq_id_2 in range(n,len(aln_file.seq_lst)):                
                for character1, character2 in zip(aln_file.seq_lst[seq_id_1],aln_file.seq_lst[seq_id_2]):
                    if character1 == character2:
                        identical +=1
                seq_1 = find_tag('seq_id',aln_file.name_lst[seq_id_1])
                seq_2 = find_tag('seq_id',aln_file.name_lst[seq_id_2])
                persent_identical = identical / len(aln_file.seq_lst[seq_id_1]) * 100
                data_persent[seq_1+'and'+seq_2] = persent_identical
                identity_matrix[seq_1][seq_2] = persent_identical
                identity_matrix[seq_2][seq_1] = persent_identical
                identical = 0
        aln_file.data_persent = data_persent
        aln_file.identity_matrix = identity_matrix

    return aln_file


def clustering_aln(directory):
    directory_files = os.listdir(directory)
    for file in directory_files:
        if file.endswith('.aln'):            
            clustering_kmeans_aln(ALNDATA_DIRECTORY + file,itself=True)

def enumeration_identity_percent(directory):
    '''Just for plot'''    
    data_persent_for_plot = pd.Series()    
    directory_files = os.listdir(directory)
    for file in directory_files:
        if file.endswith('.aln'):            
            aln_file = calculate_identity_percent(ALNDATA_DIRECTORY + file,itself=True)            
            data_persent_for_plot = data_persent_for_plot.append(aln_file.data_persent)
    return data_persent_for_plot
            

print('end indentity_persent')

