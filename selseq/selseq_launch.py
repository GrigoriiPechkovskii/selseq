#!/usr/bin/env python
from selseq_constant import *


from selseq_main import *
from parsing_blast_tbl import *
import selseq_control as selseq_control



print('start __main__')

#----------
selseq_control.del_dir(HOME_DIRECTORY,REDATA_DIRECTORY,ALNDATA_DIRECTORY)
#----------
selseq_control.initial_timecheck_line()

os.mkdir(REDATA_DIRECTORY)
os.mkdir(ALNDATA_DIRECTORY)


file_name_faa = make_file_name_faa(HOME_FILES)

selseq_control.check_file_in_table(file_name_faa)

cluster = Cluster(TABLE_CSV,file_name_faa,GROUP_TAG)

cluster.ClusterDic()

rename_fasta(HOME_FILES,cluster.cluster_dic)


joining_files(HOME_FILES,REDATA_DIRECTORY)

selseq_control.timecheck('beforblast')

blast('joint_file',QUERY_SEQ, 'blastdb','tbl.csv')

selseq_control.timecheck('blast')

parsing_blast('tbl.csv')
#parsing_balst_table('tbl.csv')

selseq_control.timecheck('tbl')

print('start programm tail')

make_tail(REDATA_DIRECTORY)

selseq_control.timecheck('tail')

blast('tail','tail','tailblastdb','tbltail.csv')

parsing_blast('tbltail.csv')
#parsing_balst_table('tbltail.csv')

selseq_control.timecheck('blatstltail')

make_muscle(REDATA_DIRECTORY)

selseq_control.timecheck('mucsle')

into(ALNDATA_DIRECTORY)

selseq_control.timecheck('into')

print('start grouping')

files = os.listdir(REDATA_DIRECTORY)

group_dict_for_selection = make_group_dict_for_selection(group_for_selection,cluster.group_dict)            

sel = Selection(group_dict_for_selection)

group_HOME_DIRECTORY = group_HOME_DIRECTORY(HOME_DIRECTORY)

make_group_muscle(group_HOME_DIRECTORY)                      

selseq_control.timecheck('mucsle')


#fig = plt.figure()
plt_dic = pd.Series(plt_dic)
plt.hist(plt_dic.values, bins=100, alpha=1,color='blue',edgecolor='black');
plt.show()
print(plt_dic)



#===================================================================              
for group_lst,subgroup_lst in group_HOME_DIRECTORY.items():
    for subgroup in subgroup_lst:
        into(subgroup)
        entropy_calculate(subgroup)
        common (subgroup)
        
entropy_calculate(ALNDATA_DIRECTORY)
common (ALNDATA_DIRECTORY)
selseq_control.timecheck('calculate')
#=============================================================================================

Cp = selseq_control.Control(REDATA_DIRECTORY,ALNDATA_DIRECTORY)
Cp.doControl()

selseq_control.timecheck('control')