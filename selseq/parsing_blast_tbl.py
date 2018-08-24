import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
from selseq_constant import *
from selseq_main import find_tag, SequenceFasta
print('start tbl_pars_pandas')

#df = pd.read_csv('tbl.csv',names = ['fir','sec','per','4444','5555','6666','7777','8888','9999','10000','111111','122222'])

full_seq_threshold = pd.Series()

def seq_choice(seq_store_in,percent_range = 3):
	'''Input seq_store_in is pd.Series all sbjct sequence of one query sequence 
	Output seq_threshold is pd.Series what chosen
	>>> seq_store_in = pd.Series({'<seq_id>seq_num_9</seq_id><assemble>GCA_2<...': 100,'<seq_id>seq_num_3</seq_id><assemble>GCA_3</assemble>...' : 100,'<seq_id>seq_num_8</seq_id><assemble>GCA_2</assemble>...' : 100})
	>>> seq_choice(seq_store_in)
	{'seq_num_3': 'GCA_3', 'seq_num_8': 'GCA_22', 'seq_num_9': 'GCA_2'}
	'''	
	global full_seq_threshold
	seq_store_in.sort_values(ascending=False,inplace=True)	
	seq_threshold= seq_store_in[seq_store_in>80]	
	seq_store_in = seq_store_in.drop(seq_threshold.index)
	persent_last = 80
	
	for name,percent in seq_store_in.items():
		if persent_last - percent <= percent_range:
			persent_last = percent
			seq_threshold[name] = percent
	
	seq_id_assemble_base = {}
	seq_id_persent = {}

	for name in seq_threshold.index:
		id_base = find_tag('seq_id',name)
		assemble_base = find_tag('assemble',name)
		seq_id_assemble_base[id_base] = assemble_base
		seq_id_persent[id_base] = seq_threshold[name]	 
	return seq_id_assemble_base,seq_id_persent 

base_used = []

def division(seq_id_assemble_base):
	for seq_id, assemble in seq_id_assemble_base.items():
		if seq_id not in base_used:	
			base_used.append(seq_id)
			assemble_fasta = SequenceFasta(REDATA_DIRECTORY + assemble + '.faa')
			assemble_fasta.seq_process()
			for prot_number in range(0,len(assemble_fasta.name_lst)):
				if seq_id in assemble_fasta.name_lst[prot_number]:
					protein_rec = open (REDATA_DIRECTORY + id_list_ref+ '.faa', 'a') 
					protein_rec.write(assemble_fasta.name_lst[prot_number] + assemble_fasta.seq_lst[prot_number])
					protein_rec.close()	
plt_dic = {}
def division2(seq_choiced):
	'''The function is intended for distribution of sequences from the blast table
	and for create plot data.
	Input two dictionary in turple
	First seq_id and asseble
	{'seq_num_3': 'GCA_3', 'seq_num_8': 'GCA_2', 'seq_num_9': 'GCA_2'}
	Second seq_id and persent
	{'seq_num_3': '100', 'seq_num_8': '85', 'seq_num_9': '90'}       
	'''
	global plt_dic
	
	for seq_id, assemble in seq_choiced[0].items():
		#print (id_list_ref)
		if seq_id not in base_used:	
			base_used.append(seq_id)
			plt_dic[seq_id + assemble] = seq_choiced[1][seq_id]
			assemble_fasta = SequenceFasta(REDATA_DIRECTORY + assemble + '.faa')
			assemble_fasta.seq_process()
			for prot_number in range(0,len(assemble_fasta.name_lst)):				
				if seq_id == find_tag('seq_id',assemble_fasta.name_lst[prot_number]):					
					protein_rec = open (REDATA_DIRECTORY + id_list_ref+ '.faa', 'a')					
					
					protein_rec.write(assemble_fasta.name_lst[prot_number] + assemble_fasta.seq_lst[prot_number])
					protein_rec.close()

def parsing_blast(tbl):
	global id_list_ref	
	seq_store = pd.Series()
	table_opened = open(REDATA_DIRECTORY + tbl)
	ref_seq_last = 'start'
	for table_line in table_opened:
		table_list = table_line.split(',')
		ref_seq = table_list[0]
		#id_list_ref = find_tag('seq_id',ref_seq)
		#adding in seq_store:pd.Series()
		if ref_seq==ref_seq_last or ref_seq_last == 'start':
			#Проверка если в таблице есть тот же белок но с большим процентом, возможно не обязательная
			if table_list[1] in seq_store:
				if seq_store[table_list[1]]<float(table_list[2]):
					seq_store[table_list[1]] = float(table_list[2])		
					ref_seq_last = table_list[0]
			else:
				seq_store[table_list[1]] = float(table_list[2])		
				ref_seq_last = table_list[0]

		else:
			id_list_ref = find_tag('seq_id',ref_seq_last)
			assemble_base = find_tag('assemble',table_list[1])
			#print('seq_store',seq_store)
			seq_choiced = seq_choice(seq_store)
			
			#print(seq_choiced[1])	
			division2(seq_choiced)

			seq_store = pd.Series()
			seq_store[table_list[1]] = float(table_list[2])
			ref_seq_last = table_list[0]

	id_list_ref = find_tag('seq_id',ref_seq_last)
	assemble_base = find_tag('assemble',table_list[1])
	seq_choiced = seq_choice(seq_store)	
	division2(seq_choiced)
	#print(seq_id_assemble_base)
	#print(seq_store)
	table_opened.close()	
#parsing_blast('tbl.csv')
'''
if __name__ == '__main__':
	import doctest
	doctest.testmod() 
'''


#w = seq_choice(seq_store)


print('end tbl_pars_pandas')


