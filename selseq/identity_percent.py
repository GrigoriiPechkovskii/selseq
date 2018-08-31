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
#from selseq_extra import * 

directory_files = os.listdir(ALNDATA_DIRECTORY)

def decorator_aln():
    def wrap_aln():
        pass

def calculate_identity_percent(aln_file,itself=True):
    '''input file the aligned sequence
       output pd.Series identity_percent
       itself - parametr for calculate identity percent the alone sequence
    '''
    global identity_matrix
    global data_persent

    aln_file = SequenceFasta(aln_file)
    aln_file.seq_process(strip=False)
    data_persent = pd.Series()
    identity_matrix = pd.DataFrame()

       
    
    if itself and len(aln_file.seq_lst) == 1:        
        data_persent[find_tag('seq_id',aln_file.name_lst[0])+'and'+find_tag('seq_id',aln_file.name_lst[0])] = 100
        identity_matrix = 0
        return identity_matrix

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
        
        #print(identity_matrix)
        return aln_file
        #return data_persent

def enumeration_aln(directory_files,fun):
    data = pd.Series()
    for file in directory_files:
        if file.endswith('.aln'):
            #print(fun(ALNDATA_DIRECTORY + file))
            data = data.append(calculate_identity_percent(ALNDATA_DIRECTORY + file))
    return data

def enumeration_aln2(directory_files):    
    for file in directory_files:
        if file.endswith('.aln'):
            
            aln_file = calculate_identity_percent(ALNDATA_DIRECTORY + file,itself=False)
            print(aln_file.file_dir)
            if identity_matrix.shape != (1,1):
                #print(identity_matrix)           
                #print(aln_file.file_dir)
                print(aln_file.data_persent.sort_values(ascending=False,inplace=False))                

                if any(aln_file.data_persent<60):
                    kmeans = KMeans(n_clusters=2)
                    kmeans.fit(identity_matrix)
                    y_kmeans = kmeans.predict(identity_matrix)

                    print(identity_matrix.index)  
                    print(y_kmeans)
                    muscle_list = []
                    for kmeans_index in range(len(y_kmeans)):
                        if y_kmeans[kmeans_index] == 0:
                            with open(aln_file.file_dir[0:-4] + '_0','a') as aln_clustered: #del .aln
                                muscle_list.append(aln_file.file_dir[0:-4] + '_0')
                                aln_clustered.write(aln_file.name_lst[kmeans_index] + aln_file.seq_lst[kmeans_index].replace('-','').replace('\n','') + '\n')
                        
                        elif y_kmeans[kmeans_index] == 1:                            
                            with open(aln_file.file_dir[0:-4] + '_1','a') as aln_clustered: #del .aln
                                muscle_list.append(aln_file.file_dir[0:-4] + '_1')
                                aln_clustered.write(aln_file.name_lst[kmeans_index] + aln_file.seq_lst[kmeans_index].replace('-','').replace('\n','') + '\n')
                        
                    for seq_file in muscle_list:
                        subprocess.call('muscle ' + '-in ' +seq_file + ' -out ' + seq_file + '.aln 2>' + HOME_DIRECTORY + '111',shell = True)    
                    #print(aln_file.seq_lst[kmeans_index])
enumeration_aln2(directory_files)
#data = enumeration_aln(directory_files,calculate_identity_percent)

#plt.hist(data.values, bins=100, alpha=1,color='blue',edgecolor='black')

'''kmeans = KMeans(n_clusters=2)
kmeans.fit(identity_matrix)
y_kmeans = kmeans.predict(identity_matrix)'''


plt.show()
print('end indentity_persent')

