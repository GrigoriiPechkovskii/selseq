'''Module for parsing blast output format 10'''

import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import os
#import re
from selseq_constant import *
from selseq_main import find_tag, SequenceFasta

print('start selseq_parsing_blast')


full_seq_threshold = pd.Series()
base_used = set()


class ParsingBlast():
    def __init__(self,*blast_tbl):

        self.blast_tbl = blast_tbl
        self.tbl = None
        self.id_query = ''
        self.choiced_full = dict()
        self.choiced = dict()
        self.sbjct_used = set()
        self.seq_id_assemble_sbjct = {}  
        self.sbjct_id_used = set() 

    def parsing_generator(self):
        '''Make a generator wich yield pd.Series one blast sbjct'''
        seq_store = pd.Series()
        
        with open(REDATA_DIRECTORY + self.tbl) as table_opened:
            query_last = 'start'
            for table_line in table_opened:
                table_list = table_line.split(',')
                query= table_list[0] #ref_seq
                sbjct = table_list[1]
                percent_tbl = table_list[2]
                #adding in seq_store:pd.Series()
                if query==query_last or query_last == 'start':
                    self.id_query = find_tag('seq_id',query)
                    #Check if the table has the same protein but with a high percentage, perhaps not mandatory
                    if sbjct in seq_store:
                        if seq_store[sbjct]<float(percent_tbl):
                            seq_store[sbjct] = float(percent_tbl)        
                            query_last = query
                    else:
                        seq_store[sbjct] = float(percent_tbl)        
                        query_last = query
                else:
                    self.id_query = find_tag('seq_id',query_last)
                    yield seq_store                    
                    seq_store = pd.Series()
                    seq_store[sbjct] = float(percent_tbl)
                    query_last = query
        yield seq_store

    def seq_choice(self,seq_store_in):
        self.choiced = dict()
        seq_store_in.sort_values(ascending=False,inplace=True)    
        seq_threshold= seq_store_in[seq_store_in>PERSENT_THRESHHOLD]    
        seq_store_in = seq_store_in.drop(seq_threshold.index)
        persent_last = PERSENT_THRESHHOLD
        
        for name,percent in seq_store_in.items():
            if persent_last - percent <= PERCENT_RANGE:
                persent_last = percent
                seq_threshold[name] = percent

        seq_id_threshold = {find_tag('seq_id',seq) for seq in set(seq_threshold.index)}
        seq_id_threshold = set(seq_threshold.index)#!!!maybe need change
        
        self.choiced[self.id_query] = seq_id_threshold
        
        self.seq_choised_union(self.choiced)

        seq_id_assemble_base = {}
        seq_id_persent = {}

        for name in seq_threshold.index:
            id_base = find_tag('seq_id',name)
            assemble_base = find_tag('assemble',name)
            seq_id_assemble_base[id_base] = assemble_base
            seq_id_persent[id_base] = seq_threshold[name]     
        return seq_id_assemble_base#,seq_id_persent

    def seq_choised_union(self,choiced):       
        #if not choiced[self.id_query].issubset(self.sbjct_used):   #Attantion maybe it not requied     
            #self.sbjct_used = base_used.union(choiced[self.id_query])
        #print('self.choiced_full===',self.choiced_full)
        if len(self.choiced_full.values()) == 0:
            self.choiced_full[self.id_query] =choiced[self.id_query]

        else:            
            flag = True
            for val_set in self.choiced_full:
                if not self.choiced_full[val_set].isdisjoint(choiced[self.id_query]):                    
                    self.choiced_full[val_set].update(choiced[self.id_query])#
                    flag = False
                
            if flag:                
                self.choiced_full[self.id_query] =choiced[self.id_query]

    def parsing_with_union(self,*blast_tbl_add):#!!!do not undrstand why var blast_tbl_add exist
        '''Parsing blast table with using union set'''        
        if len(blast_tbl_add) > 0:
            self.blast_tbl = blast_tbl_add 
        for tbl in self.blast_tbl:
            self.tbl = tbl
            for gen in self.parsing_generator():
                #print('!!!!!!',gen)
                self.seq_choice(gen)

    def distribution(self):
        for id_query in self.choiced_full:
            for sbjct in self.choiced_full[id_query]:
                sbjct_assemble = find_tag('assemble',sbjct)
                sbjct_id = find_tag('seq_id',sbjct)
                if sbjct_id not in self.sbjct_id_used:                                    
                    self.sbjct_id_used.add(sbjct_id)
                    assemble_fasta = SequenceFasta(REDATA_DIRECTORY + sbjct_assemble + '.faa')
                    assemble_fasta.seq_process()
                    for index in range(len(assemble_fasta.name_lst)):
                        if sbjct_id == find_tag('seq_id',assemble_fasta.name_lst[index]):                    
                            with open(REDATA_DIRECTORY + id_query + '.faa','a') as new_fasta:
                                    new_fasta.write(assemble_fasta.name_lst[index] + assemble_fasta.seq_lst[index])
            
if __name__ == '__main__':             
    pb = ParsingBlast(*['join_assemble1_tbl.csv','join_assemble2_tbl.csv','join_assemble3_tbl.csv'])
    pb.parsing_with_union()
    pb.distribution()
    



print('end selseq_parsing_blast')


