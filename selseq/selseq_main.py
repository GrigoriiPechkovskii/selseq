#!/grin/bin/env python3
#by Grigorii Pechkovskii, grigorii.pechkovskii@gmail.com

print('start programm')

import os
import subprocess
import math
import numpy as np
import re
from selseq_constant import *

def make_tags(tags_names,tags_contents,string):
	'''Adds tags to string and return changed string with tags'''
	if type(tags_names) is str and type(tags_contents) is str:
		tags_name_open,tags_name_close ='<' + tags_names + '>','</' + tags_names + '>'		
		string = string + tags_name_open + tags_contents + tags_name_close
		return string
	else:		
		for tags_name, tags_content in zip(tags_names,tags_contents):
			tags_name_open,tags_name_close ='<' + tags_name + '>','</' + tags_name + '>'		
			string = string + tags_name_open + tags_content + tags_name_close
		return string

def find_tag(tags_name,string):
    '''Find tags in string and return tag contents'''
    tags_name_open ='<' + tags_name + '>'
    tags_name_close = '</' + tags_name + '>'
    patern = tags_name_open + '(.*)' + tags_name_close
    find_tag = re.search(patern,string) #смотри на сырые строки их нету
    return find_tag.group(1)

class SequenceFasta():
    '''Processes a fasta file'''
    def __init__(self,file_dir):
        self.file_dir=file_dir
        self.name_lst = []
        self.seq_lst = []
    def just_named(self):
        '''Works only with sequence names'''
        with open(self.file_dir) as file_opened:
            for line in file_opened:
                if '>' in line:
                    self.name_lst.append(line)
    def seq_process(self,strip=False):
        seq_tmp = ''
        with open(self.file_dir) as file_opened:            
            for line in file_opened:
                if '>' in line:
                    self.seq_lst.append(seq_tmp)
                    seq_tmp = ''
                    if strip:
                        self.name_lst.append(line.rstrip())
                    else:
                        self.name_lst.append(line)	
                else:
                    if strip:
                        seq_tmp += line.rstrip()
                    else:
                        seq_tmp += line
        self.seq_lst.append(seq_tmp)
        del self.seq_lst[0]
        self.seq_len = len(self.seq_lst)


class Cluster():
    '''Return dict with matching the name of the genome and its properties (cluster)'''
    def __init__(self,CSV,file_name,clusters = GROUP_TAG):
        self.CSC =CSV
        self.file_name = file_name
        self.clusters = clusters
        self.csvprep = []
        self.table_head = []
        self.clusters_indexing = []
        self.cluster_dic = {}        
        self.group_dict = {}        
        self._CsvPrep()

    def _CsvPrep(self):
        table_opened = open (self.CSC)
        for line in table_opened:
            self.csvprep += [[i.strip() for i in line.strip().split(',')]]
        table_opened.close()

    def _TableHead(self):
        for line in self.csvprep:
            if '#' in line[0]:
                self.table_head += line

    def _ClusterIndexing(self):
        self._TableHead()
        for clusters_name in self.clusters:
            self.clusters_indexing += [self.table_head.index(clusters_name)]
            
    def _Groups(self):
        group_set = set()
        for index in self.clusters_indexing:            
            for line in self.csvprep:                
                if line[index] not in self.table_head:
                    group_set.add(self.table_head[index]+'='+line[index])                    
            self.group_dict[self.table_head[index]] = group_set            
            group_set = set()
        
    def ClusterDic(self):
        self._ClusterIndexing()
        self._Groups()
        cluster = []
        for line in self.csvprep:
            for name in self.file_name:
                if name in line:
                    for index in self.clusters_indexing:                        
                        cluster += [self.table_head[index] + '=' + line[index].strip('"')] #то что пойдет в имя белка
                    self.cluster_dic[name] = cluster
                    cluster = []
        return self.cluster_dic

def make_file_name_faa(files):
	file_name_faa  = [] 
	for file in files:
	    if '.faa' in file:
	        file_name_faa += [file[0:-4]]
	return file_name_faa


def UnpackCluster(lst):
    val_str = ''
    for val in lst:
        val_str += ',' + val#Note splitter
    return val_str.strip(',')

def rename_fasta(files,cluster_dic):     
    num = 0
    for file_name in files:
        if '.faa' in file_name:
            file_opened = open (HOME_DIRECTORY + file_name, 'r' )
            file_reopened = open (REDATA_DIRECTORY + file_name, 'a' )
            for file_line in file_opened:
                if '>' in file_line:
                    num += 1
                    string = '>'
                    file_line = make_tags(['seq_id','assemble','main','cluster_tag'],['seq_num_' + str(num),file_name[0:-4],file_line.strip('>\n'),UnpackCluster(cluster_dic[file_name[0:-4]])],string) + '\n' 
                    #file_line2 ='>' + 'seq_num_' + str(num) + '<'+ file_name[0:-4]+ '<' + file_line.strip('>\n') + '&&&' + UnpackCluster(cluster.cluster_dic[file_name[0:-4]])+'\n'  #внимание на разделитель
                    file_reopened.write(file_line)                          

                else:
                   file_reopened.write(file_line)
            file_opened.close()
            file_reopened.close()


def joining_files(files,REDATA_DIRECTORY):
    for file_name in files:
        if '.faa' in file_name: #!!!!!
            file_opened = open ( REDATA_DIRECTORY + file_name,'r' )
            file_reopened = open (REDATA_DIRECTORY + 'joint_file', 'a' )
            var_name_file = file_opened.read()
            file_reopened.write(var_name_file)
            file_opened.close()
            file_reopened.close()

def blast(DB,QUERY_SEQ,OUT_BD,OUT_TBL):

    if sys.platform == 'linux':
        subprocess.call('/home/grigoriipechkovskii/bio/ncbi-blast-2.7.1+/bin/makeblastdb -in '+ REDATA_DIRECTORY + DB + ' -dbtype prot -out ' + REDATA_DIRECTORY + OUT_BD + ' >' + HOME_DIRECTORY + '111',stdout=subprocess.PIPE,shell=True)
        subprocess.call('/home/grigoriipechkovskii/bio/ncbi-blast-2.7.1+/bin/blastp -db '+ REDATA_DIRECTORY + OUT_BD + ' -query '+ REDATA_DIRECTORY + QUERY_SEQ + ' -out '+ REDATA_DIRECTORY + OUT_TBL + ' -outfmt 10 -evalue 0.001 2>' + HOME_DIRECTORY + '111', shell=True)
    
    if sys.platform == 'win32':    
        subprocess.call('makeblastdb -in '+ REDATA_DIRECTORY + DB + ' -dbtype prot -out ' + REDATA_DIRECTORY + OUT_BD +' >' + HOME_DIRECTORY + '111',stdout=subprocess.PIPE,shell=True)
        subprocess.call('blastp -db '+ REDATA_DIRECTORY + OUT_BD + ' -query '+ REDATA_DIRECTORY + QUERY_SEQ + ' -out '+ REDATA_DIRECTORY + OUT_TBL +' -outfmt 10 -evalue 0.001 2>' + HOME_DIRECTORY + '111', shell=True)

def parsing_balst_table(tbl):
    
    files = os.listdir(ALNDATA_DIRECTORY)
    base_used = list()
    table_opened = open(REDATA_DIRECTORY + tbl, 'r')
    for table_line in table_opened: 
        table_list = table_line.split(',')
        if float(table_list[2])>80:
            
            id_list_ref =find_tag('seq_id',table_list[0]) #список полное название референта
            id_list_base =find_tag('seq_id',table_list[1]) #список полное название найденного белка
            assemble_base = find_tag('assemble',table_list[1])            

            if table_list[1] not in base_used:                
                omics = SequenceFasta(REDATA_DIRECTORY + assemble_base + '.faa')
                omics.seq_process()                                                
                base_used += [table_list[1]]
                
                for prot in range(0,len(omics.name_lst)):                   
                    if table_list[1] in omics.name_lst[prot]:                       
                        protein_rec = open (REDATA_DIRECTORY + id_list_ref+ '.faa', 'a') 
                        protein_rec.write(omics.name_lst[prot] + omics.seq_lst[prot])
                        protein_rec.close()  
    table_opened.close()

print('start programm tail')

def make_tail(REDATA_DIRECTORY):

    files = os.listdir(REDATA_DIRECTORY)
    set_aln = set()

    for file_aln in files:  
        if '.faa' and 'seq_num' in file_aln :
            file_aln_opened = open(REDATA_DIRECTORY + file_aln)
            for line_file_aln in file_aln_opened:
                if '>' in line_file_aln:
                    set_aln.add(line_file_aln.strip()) 
    file_aln_opened.close()


    files = os.listdir(REDATA_DIRECTORY)

    #tail = open(REDATA_DIRECTORY + 'tail', 'a')
    prot_save = str()
    name_save = str()
    tail = open(REDATA_DIRECTORY + 'tail', 'a')
    for file_GCF in files:
        if 'GCA' in file_GCF:

            omics = SequenceFasta(REDATA_DIRECTORY + file_GCF)
            omics.seq_process()
            for index in range(0,len(omics.name_lst)):
               if omics.name_lst[index].strip() not in set_aln:
                    tail.write(omics.name_lst[index] + omics.seq_lst[index])

    tail.close()

def make_muscle(REDATA_DIRECTORY):

    files = os.listdir(REDATA_DIRECTORY)
    for protein_file in files:
        if 'seq_num' in protein_file:
            subprocess.call('muscle ' + '-in ' +REDATA_DIRECTORY + protein_file + ' -out ' + ALNDATA_DIRECTORY + protein_file[0:-4] + '.aln 2>' + HOME_DIRECTORY + '111',shell = True)

def into(direct):
    home_files = os.listdir(HOME_DIRECTORY)
    prot_count = 0
    matrix_genom = []
    files = os.listdir(direct)
        
    for file_genom in home_files:
        if 'GC' in file_genom:
            matrix_genom += [file_genom[0:-4]]

    into = open(direct + 'into.csv', 'w')
    into.write(',' + ','.join(matrix_genom))
    into.close()    
    for file_protein in files:
            
            if 'seq_num' in file_protein:
                
                proteom_opened = open(direct + file_protein)
                proteom = proteom_opened.read()
                proteom_opened.close()
                into = open (direct + 'into.csv', 'a')
                into.write('\n'+file_protein)
                for genom in matrix_genom:
                    if genom in proteom:
                        into.write(',' + str(proteom.count(genom)))
                    else:
                        into.write(',' + '0')
                into.close()
                
print('start grouping')

class Selection():

    def __init__(self,group=dict):
        
        self.files = os.listdir(ALNDATA_DIRECTORY)
        self.group = group
        self.seq_num_list = []
        self.seq_direct()
        
    def seq_num_prep(self,file):
        seq_opened = open(ALNDATA_DIRECTORY + file)
        self.seq_num_list = seq_opened.read().split('>')
        self.seq_num_list.pop(0)
        self.seq_num_list = ['>'+i for i in self.seq_num_list]
        seq_opened.close()
        
    def seq_direct(self):        
        for file in self.files:
            if 'seq_num' in file:
                #self.seq_num_prep(file)

                omics = SequenceFasta(ALNDATA_DIRECTORY + file)
                omics.seq_process()

                for index in range(0,len(omics.name_lst)):                    
                    omics.seq_lst[index] = omics.seq_lst[index].replace('-','').replace('\n','') + '\n'                    
                    for key in self.group:
                        if not os.access(HOME_DIRECTORY + NAME_EXP+'aln_'+ key + '/',os.F_OK):
                            os.mkdir(HOME_DIRECTORY + NAME_EXP+'aln_'+ key + '/') 
                        for val in self.group[key]:
                            if not os.access(HOME_DIRECTORY +NAME_EXP+'aln_'+ key + '/' + 'aln_' + val+ '/',os.F_OK):
                                os.mkdir(HOME_DIRECTORY +NAME_EXP+'aln_' + key + '/' + 'aln_' + val+ '/')                                           
                        
                            if val in omics.name_lst[index]:                                                    
                                seq_group_opened = open(HOME_DIRECTORY +NAME_EXP+'aln_'+ key + '/' +  'aln_' +val + '/' + file,'a')
                                seq_group_opened.write(omics.name_lst[index] + omics.seq_lst[index])                                
                                seq_group_opened.close()
                    
def make_group_dict_for_selection(group_for_selection,group_dict):
    global a
    group_dict_for_selection = {}
    for gr in group_for_selection:        
        for dic_items in group_dict.items():            
            if gr in dic_items:
                group_dict_for_selection.update(dict.fromkeys([dic_items[0]],dic_items[1]))
                
    return group_dict_for_selection

def grouping_HOME_DIRECTORY(HOME_DIRECTORY):                    
	group_HOME_DIRECTORY = {}
	HOME_DIRECTORY_lev1 = os.listdir(HOME_DIRECTORY)
	for file_lev1 in HOME_DIRECTORY_lev1:    
	    if 'aln_' in file_lev1:        
	        dir_lst = []
	        HOME_DIRECTORY_lev2 = os.listdir(HOME_DIRECTORY + file_lev1)
	        for file_lev2 in HOME_DIRECTORY_lev2:
	            aln_group_dir = HOME_DIRECTORY + file_lev1 + '/' + file_lev2 + '/'
	            if os.path.isdir(aln_group_dir):
	                dir_lst += [aln_group_dir]
	                                       
	        group_HOME_DIRECTORY[HOME_DIRECTORY + file_lev1 + '/'] = dir_lst
	return group_HOME_DIRECTORY

def make_group_muscle(group_HOME_DIRECTORY):
    for group_lst,subgroup_lst in group_HOME_DIRECTORY.items():
        for subgroup in subgroup_lst:
            files = os.listdir(subgroup)
            for file in files:
                if 'seq_num' in file:
                    subprocess.call('muscle ' + '-in ' +subgroup + file + ' -out ' + subgroup + file[0:-4] + '.aln 2>' + HOME_DIRECTORY + '111',shell = True)
                    #os.remove(subgroup + file)

def entropy_calculate(subgroup):
    HOME_DIRECTORY  = os.getcwd()
    files = os.listdir(subgroup)    

    for file_aln in files:
        if 'seq_num' in file_aln:
            
            matrix_aln_name = []
            matrix_aln_prot = []
            matrix_aln_prot2 = str()
        
            file_aln_opened = open (subgroup + file_aln)
            for line_file_aln in file_aln_opened:
                if '>' in line_file_aln:
                    matrix_aln_name += [line_file_aln.strip()]
                    matrix_aln_prot += [matrix_aln_prot2]
                    matrix_aln_prot2 = str()                
                else:
                    matrix_aln_prot2 += line_file_aln.strip()               

            matrix_aln_prot += [matrix_aln_prot2]
            matrix_aln_prot.pop(0)       
            file_aln_opened.close()

            #для нахождение номера рефереса
            ref_number = int()

############+++++++++
            #расчет энтропии
            entropy_list = list()
            column_matrix = list()
            for amin_num in range(0,len(matrix_aln_prot[0])):        
                    column_prot = str() 
                    for prot_num in range(0,len(matrix_aln_prot)):
                        column_prot += matrix_aln_prot[prot_num][amin_num]
                    column_matrix += [column_prot]
            
        
            for protein in column_matrix:
                set_amino = str()
                entropy = int()
            
                for amino in protein:
                    if amino in set_amino:
                        continue
                    else:
                        entropy += protein.count(amino)/len(protein)* math.log2(protein.count(amino)/len(protein))
                        set_amino += amino
                entropy_list += [round(abs(entropy),3)]            
            report_entropy_seq = open (subgroup + 'report_entropy_seq.txt', 'a')
            report_entropy_seq.write (file_aln + ',' + str(entropy_list).strip('[,]').replace(' ', '') + '\n')
            report_entropy_seq.close()
        
            mean_entropy = np.mean(entropy_list)
            sum_entropy = np.sum(entropy_list)
        
            report_entropy_stat = open (subgroup + 'report_entropy_stat.txt', 'a')
            report_entropy_stat.write (file_aln + ',' + str(round(mean_entropy,5)) + ',' + str(round(sum_entropy,5)) + '\n')
            report_entropy_stat.close()
        

def common (direct):
    into_list = list()
    into_list_proteom = list()
    intotbl_opened = open(direct + 'into.csv')
    #открывает таблицу вхождения белков и зписывает их в список
    prot_sum = 0
    for into_line in intotbl_opened:
        if 'GCA' in into_line:
            into_list_proteom += into_line.strip().split(',')
        else:
            into_list += [into_line.strip().split(',')]
            prot_sum += 1
    del into_list_proteom[0]  
    
    sum_list = list()

    for proteom_1 in range (1, len(into_list_proteom)+1):    
        sum_count = int()
        for inti_str in into_list:
            sum_count += int(inti_str[proteom_1])
    
        sum_list.append(sum_count)
                
    common = str()
    common_line = str()

    for i in into_list_proteom:
        common_line += ',' + i
    common_line += '\n'
    
    for proteom_1 in range (1, len(into_list_proteom)+1):
        common_line += into_list_proteom[proteom_1 - 1]
        for proteom_2 in range (1, len(into_list_proteom)+1):        
            common_count = int()
            for inti_str in into_list:
            
                common_count += min(int(inti_str[proteom_1]),int(inti_str[proteom_2]))
            common_line += ',' + str( round(common_count) ) #Убери деление будут числа а не проценты

        
               
        common_line +=  '\n'
    
    report_common2 = open (direct + 'report_common2.csv','w')
    report_common2.write(common_line.strip())
    report_common2.close()



if __name__ == '__main__':
    print('start __main__')

    