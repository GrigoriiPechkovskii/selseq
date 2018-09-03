print('start identity_percent')
from selseq_main import find_tag, SequenceFasta
import os
from selseq_constant import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import re
import subprocess
import shutil
#from selseq_extra import * 


def calculate_identity_percent(aln_file,itself=True):
    '''input file the aligned sequence
       output pd.Series identity_percent
       itself - parametr for calculate identity percent the alone sequence
    '''  
    aln_file = SequenceFasta(aln_file)
    aln_file.seq_process(strip=False)
    data_persent = pd.Series()
    identity_matrix = pd.DataFrame()    
    
    if itself and len(aln_file.seq_lst) == 1:        
        data_persent[find_tag('seq_id',aln_file.name_lst[0])+'and'+find_tag('seq_id',aln_file.name_lst[0])] = 100
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
                data_persent[seq_1+'and'+seq_2] = identical / len(aln_file.seq_lst[seq_id_1]) * 100
                identity_matrix[seq_1][seq_2] = identical / len(aln_file.seq_lst[seq_id_1]) * 100
                identity_matrix[seq_2][seq_1] = identical / len(aln_file.seq_lst[seq_id_1]) * 100
                identical = 0
        aln_file.data_persent = data_persent
        aln_file.identity_matrix = identity_matrix       
        
        return aln_file        


def clustering_aln(directory):    
    data_persent_for_plot = pd.Series() 
    print('i am worked')
    directory_files = os.listdir(directory)
    for file in directory_files:
        if file.endswith('.aln'):            
            aln_file = calculate_identity_percent(ALNDATA_DIRECTORY + file,itself=True)
            data_persent_for_plot = data_persent_for_plot.append(aln_file.data_persent)            
            if aln_file.identity_matrix.shape != (0, 0):                               
                print(aln_file.file_dir)
                while any((aln_file.identity_matrix<60).any()):
                    print('loop')
                    
                    kmeans = KMeans(n_clusters=2)
                    kmeans.fit(aln_file.identity_matrix)
                    y_kmeans = kmeans.predict(aln_file.identity_matrix)                    
                    muscle_set = set()
                    for kmeans_index in range(len(y_kmeans)):
                        if y_kmeans[kmeans_index] == 0:
                            with open(aln_file.file_dir[0:-4] + '_0','a') as aln_clustered: #del .aln
                                muscle_set.add(aln_file.file_dir[0:-4] + '_0')
                                aln_file.file_dir[0:-4] + '_0'
                                aln_clustered.write(aln_file.name_lst[kmeans_index] + aln_file.seq_lst[kmeans_index].replace('-','').replace('\n','') + '\n')
                        
                        elif y_kmeans[kmeans_index] == 1:                            
                            with open(aln_file.file_dir[0:-4] + '_1','a') as aln_clustered: #del .aln
                                muscle_set.add(aln_file.file_dir[0:-4] + '_1')
                                aln_clustered.write(aln_file.name_lst[kmeans_index] + aln_file.seq_lst[kmeans_index].replace('-','').replace('\n','') + '\n')
                    
                    os.remove(aln_file.file_dir)
                    
                    print(muscle_set)
                    for seq_file in muscle_set:
                        subprocess.call('muscle ' + '-in ' +seq_file + ' -out ' + seq_file + '.aln 2>' + HOME_DIRECTORY + '111',shell = True)
                        os.remove(seq_file)    
                        aln_file = calculate_identity_percent(seq_file + '.aln',itself=True)
                else:
                    print('else===',aln_file.file_dir)               
                                        
    return data_persent_for_plot
                    

#data_persent_for_plot = clustering_aln(ALNDATA_DIRECTORY)
#data = enumeration_aln(directory_files,calculate_identity_percent)

#plt.hist(data.values, bins=100, alpha=1,color='blue',edgecolor='black')

#plt.show()
print('end indentity_persent')

