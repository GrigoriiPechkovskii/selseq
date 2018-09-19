#!/usr/bin/env python
from selseq_constant import *
from selseq_clustering import *
from selseq_main import *
from parsing_blast_tbl import *
import selseq_control as selseq_control
from selseq_plot import *
from selseq_parsing_blast import *

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

#blast('joint_file',QUERY_SEQ, 'blastdb','tbl.csv')
#blast_total()
blast_selectively(assemble_query_files)
selseq_control.timecheck('blast')

#pb = ParsingBlast(*['join_assemble1_tbl.csv','join_assemble2_tbl.csv','join_assemble3_tbl.csv'])
pb = ParsingBlast(*['tbl_InfluenzaA1.csv','tbl_InfluenzaA3.csv'])
pb.parsing_with_union()
pb.distribution()
#parsing_blast('tbl.csv')
#parsing_balst_table('tbl.csv')

selseq_control.timecheck('Parsing_tbl')

print('start programm tail')

make_tail(REDATA_DIRECTORY)

selseq_control.timecheck('tail')

blast('tail','tail','tailblastdb','tbltail.csv')

pbt = ParsingBlast(*['tbltail.csv'])
pbt.parsing_with_union()
pbt.distribution()
#parsing_blast('tbltail.csv')
#parsing_balst_table('tbltail.csv')

selseq_control.timecheck('blatstltail')

make_muscle(REDATA_DIRECTORY)

selseq_control.timecheck('mucsle')

into(ALNDATA_DIRECTORY,name='into_before.csv')
plot_hist_frequency_into(ALNDATA_DIRECTORY,'into_before.csv',name_plot=ALNDATA_DIRECTORY.rsplit(slash,2)[-2]+'_plot_into_before.png')       
plot_hist_frequency_into_pie_chart(ALNDATA_DIRECTORY,'into_before.csv',name_plot=ALNDATA_DIRECTORY.rsplit(slash,2)[-2]+'_plot_pie_into_before.png')       

selseq_control.timecheck('into')

print('start grouping')

files = os.listdir(REDATA_DIRECTORY)

group_dict_for_selection = make_group_dict_for_selection(group_for_selection,cluster.group_dict)            

#sel = Selection(group_dict_for_selection)
#group_HOME_DIRECTORY = grouping_HOME_DIRECTORY(HOME_DIRECTORY)
#make_group_muscle(group_HOME_DIRECTORY)                      

selseq_control.timecheck('mucsle')

data_persent_for_plot_before = enumeration_identity_percent(ALNDATA_DIRECTORY)

clustering_aln(ALNDATA_DIRECTORY)

sel = Selection(group_dict_for_selection)
group_HOME_DIRECTORY = grouping_HOME_DIRECTORY(HOME_DIRECTORY)
make_group_muscle_with_plot(group_HOME_DIRECTORY)  

data_persent_for_plot_after = enumeration_identity_percent(ALNDATA_DIRECTORY)

selseq_control.timecheck('clustering_aln')



print(len(data_persent_for_plot_before))
plot_hist_frequency(data_persent_for_plot_before.values[data_persent_for_plot_before.values != 110],
                    ALNDATA_DIRECTORY,'data_persent_for_plot_before.png')

print(len(data_persent_for_plot_after))
plot_hist_frequency(data_persent_for_plot_after.values[data_persent_for_plot_after.values != 110],
                    ALNDATA_DIRECTORY,'data_persent_for_plot_after.png')

plot_hist_frequency_for_two(data_persent_for_plot_before.values[data_persent_for_plot_before.values != 110],
                    data_persent_for_plot_after.values[data_persent_for_plot_after.values != 110],
                    ALNDATA_DIRECTORY,'data_persent_for_plot_befor_and_after.png')

print(len(plt_dic))
data_persent_for_plot_blast = pd.Series(plt_dic)
plot_hist_frequency(data_persent_for_plot_blast.values[data_persent_for_plot_blast.values != 100],
                    ALNDATA_DIRECTORY,'data_persent_for_plot_blast.png')


#===================================================================              
for group_lst,subgroup_lst in group_HOME_DIRECTORY.items():
    for subgroup in subgroup_lst:
        into(subgroup)
        plot_hist_frequency_into(subgroup,'into.csv',name_plot=subgroup.rsplit('/',2)[-2]+'_plot_into.png')
        plot_hist_frequency_into_pie_chart(subgroup,'into.csv',name_plot=subgroup.rsplit(slash,2)[-2]+'_plot_pie_into.png')
        entropy_calculate(subgroup)
        common (subgroup)
into(ALNDATA_DIRECTORY)

plot_hist_frequency_into(ALNDATA_DIRECTORY,'into.csv',name_plot=ALNDATA_DIRECTORY.rsplit(slash,2)[-2]+'_plot_into.png')       

plot_hist_frequency_into_pie_chart(ALNDATA_DIRECTORY,'into.csv',name_plot=ALNDATA_DIRECTORY.rsplit(slash,2)[-2]+'_plot_pie_into.png')       

entropy_calculate(ALNDATA_DIRECTORY)
common (ALNDATA_DIRECTORY)
selseq_control.timecheck('calculate')
#=============================================================================================

Cp = selseq_control.Control(REDATA_DIRECTORY,ALNDATA_DIRECTORY)
Cp.doControl()


selseq_control.timecheck('control')