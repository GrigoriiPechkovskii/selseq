#by Grigorii Pechkovskii
'''Selseq hight class for launch programm'''
print('start high')

from selseq_constant import *
from selseq_clustering import *
from selseq_main import *
from parsing_blast_tbl import *
import selseq_control as selseq_control
from selseq_plot import *
from selseq_parsing_blast import *


class selseq_launcher():
    def __init__(self):
        pass

    def start(self):

        selseq_control.del_dir(HOME_DIRECTORY,REDATA_DIRECTORY,ALNDATA_DIRECTORY)
        selseq_control.initial_timecheck_line()
        os.mkdir(REDATA_DIRECTORY)
        os.mkdir(ALNDATA_DIRECTORY)
        selseq_control.check_file_in_table(ASSEMBLE_NAMES)

    def selectively(self):    

        self.start()

        #if TABLE_CSV exsist
        cluster = Cluster(TABLE_CSV,ASSEMBLE_NAMES,GROUP_TAG)
        cluster.ClusterDic()
        rename_fasta(HOME_FILES,cluster.cluster_dic)

        #joining_files(HOME_FILES,REDATA_DIRECTORY)
        db_for_blast(ASSEMBLE_FILES,'joint_file')

        selseq_control.timecheck('before_blast_selectively')
        
        blast_selectively(assemble_query_files)
        selseq_control.timecheck('blast')
        pb = ParsingBlast(*QUERY_SEQ_LIST)
        pb.parsing_with_union()
        pb.distribution()
        selseq_control.timecheck('Parsing_tbl')                
        make_tail(REDATA_DIRECTORY)
        selseq_control.timecheck('tail')
        blast('tail','tail','tailblastdb','tbltail.csv')        
        pb.parsing_with_union(*['tbltail.csv'])
        pb.distribution()
        selseq_control.timecheck('blatstltail')

        make_muscle(REDATA_DIRECTORY)
        selseq_control.timecheck('mucsle')

        self.into_plot_before()

        #files = os.listdir(REDATA_DIRECTORY)
        group_dict_for_selection = make_group_dict_for_selection(group_for_selection,cluster.group_dict)            
        
        self.data_persent_for_plot_before = enumeration_identity_percent(ALNDATA_DIRECTORY)

        clustering_aln(ALNDATA_DIRECTORY)

        sel = Selection(group_dict_for_selection)

        self.group_HOME_DIRECTORY = grouping_HOME_DIRECTORY(HOME_DIRECTORY)
        make_group_muscle_with_plot(self.group_HOME_DIRECTORY)
        selseq_control.timecheck('make_group_muscle_with_plot')  
        self.group_HOME_DIRECTORY[ALNDATA_DIRECTORY] = [ALNDATA_DIRECTORY]

        self.high_stat()
        self.high_plot()
        self.high_control()        

    def high_plot(self):

        data_persent_for_plot_after = enumeration_identity_percent(ALNDATA_DIRECTORY)
        selseq_control.timecheck('clustering_aln')

        print(len(self.data_persent_for_plot_before))
        plot_hist_frequency(self.data_persent_for_plot_before.values[self.data_persent_for_plot_before.values != 110],
                            ALNDATA_DIRECTORY,'data_persent_for_plot_before.png')
        print(len(data_persent_for_plot_after))
        plot_hist_frequency(data_persent_for_plot_after.values[data_persent_for_plot_after.values != 110],
                            ALNDATA_DIRECTORY,'data_persent_for_plot_after.png')
        plot_hist_frequency_for_two(self.data_persent_for_plot_before.values[self.data_persent_for_plot_before.values != 110],
                            data_persent_for_plot_after.values[data_persent_for_plot_after.values != 110],
                            ALNDATA_DIRECTORY,'data_persent_for_plot_befor_and_after.png')
        print(len(plt_dic))
        data_persent_for_plot_blast = pd.Series(plt_dic)
        plot_hist_frequency(data_persent_for_plot_blast.values[data_persent_for_plot_blast.values != 100],
                            ALNDATA_DIRECTORY,'data_persent_for_plot_blast.png')

        for subgroup_lst in self.group_HOME_DIRECTORY.values():
            for subgroup in subgroup_lst:
                plot_hist_frequency_into(subgroup,'into.csv',name_plot=subgroup.rsplit('/',2)[-2]+'_plot_into.png')
                plot_hist_frequency_into_pie_chart(subgroup,'into.csv',name_plot=subgroup.rsplit(slash,2)[-2]+'_plot_pie_into.png')

    def into_plot_before(self):

        into(ALNDATA_DIRECTORY,name='into_before.csv')
        plot_hist_frequency_into(ALNDATA_DIRECTORY,'into_before.csv',name_plot=ALNDATA_DIRECTORY.rsplit(slash,2)[-2]+'_plot_into_before.png')       
        plot_hist_frequency_into_pie_chart(ALNDATA_DIRECTORY,'into_before.csv',name_plot=ALNDATA_DIRECTORY.rsplit(slash,2)[-2]+'_plot_pie_into_before.png')       
        selseq_control.timecheck('into_before')

    def high_stat(self):
        for subgroup_lst in self.group_HOME_DIRECTORY.values():
            for subgroup in subgroup_lst:
                into(subgroup)
                entropy_calculate(subgroup)
                common (subgroup)

    def high_control(self):
        Cp = selseq_control.Control(REDATA_DIRECTORY,ALNDATA_DIRECTORY)
        Cp.doControl()
        selseq_control.timecheck('control')



    def total(self):
        self.start()

        #if TABLE_CSV exsist
        cluster = Cluster(TABLE_CSV,ASSEMBLE_NAMES,GROUP_TAG)
        cluster.ClusterDic()
        rename_fasta(HOME_FILES,cluster.cluster_dic)

        #joining_files(HOME_FILES,REDATA_DIRECTORY)
        db_for_blast(ASSEMBLE_FILES,'joint_file')

        selseq_control.timecheck('before_blast_total')
        
        blast_total()#!
        selseq_control.timecheck('blast')
        pb = ParsingBlast(*['join_assemble1_tbl.csv','join_assemble2_tbl.csv','join_assemble3_tbl.csv'])
        pb.parsing_with_union()
        pb.distribution()
        selseq_control.timecheck('Parsing_tbl') 

        #make_tail(REDATA_DIRECTORY)
        #selseq_control.timecheck('tail')
        #blast('tail','tail','tailblastdb','tbltail.csv')        
        #pb.parsing_with_union(*['tbltail.csv'])
        #pb.distribution()
        #selseq_control.timecheck('blatstltail')

        make_muscle(REDATA_DIRECTORY)
        selseq_control.timecheck('mucsle')

        self.into_plot_before()

        #files = os.listdir(REDATA_DIRECTORY)
        group_dict_for_selection = make_group_dict_for_selection(group_for_selection,cluster.group_dict)            
        
        self.data_persent_for_plot_before = enumeration_identity_percent(ALNDATA_DIRECTORY)

        clustering_aln(ALNDATA_DIRECTORY)

        sel = Selection(group_dict_for_selection)

        self.group_HOME_DIRECTORY = grouping_HOME_DIRECTORY(HOME_DIRECTORY)
        make_group_muscle_with_plot(self.group_HOME_DIRECTORY)
        selseq_control.timecheck('make_group_muscle_with_plot')  
        self.group_HOME_DIRECTORY[ALNDATA_DIRECTORY] = [ALNDATA_DIRECTORY]

        self.high_stat()
        self.high_plot()
        self.high_control()     


if __name__ == '__main__':   
    ln = selseq_launcher()
    #ln.selectively()
    ln.total()
    print('end high')