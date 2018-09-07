from selseq_constant import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

print('start selseq_plot')

def plot_hist_frequency_into(directory,into,name_plot='plot_into.png'):
    plt.style.use('ggplot')
    fig = plt.figure()
    csv_data = pd.read_csv(directory + into)
    data_sum_row = csv_data.sum(axis=1)    
    plt.title(name_plot[0:-4])
    plt.xlabel("sum_sequence")
    plt.ylabel("frequency")
    plt.hist(data_sum_row.values, bins=100, alpha=1,color='green',edgecolor='black')
    plt.savefig(directory+name_plot ,fmt='png')

def plot_hist_frequency_into_pie_chart(directory,into,name_plot='plot_into.png'):
    plt.style.use('ggplot')
    fig = plt.figure()
    csv_data = pd.read_csv(directory + into)
    data_sum_row = csv_data.sum(axis=1)
    unique_elements, counts_elements = np.unique(data_sum_row, return_counts=True)
    unique_elements = [str(elements) for elements in unique_elements]    
    plt.title(name_plot[0:-4])
    #colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']
    explode = [0.04]*len(unique_elements)
    plt.pie(counts_elements, explode=explode,labels=unique_elements,
        autopct=lambda p : '{:.2f}%  ({:,.0f})'.format(p,p * sum(counts_elements)/100), shadow=True, startangle=140)
    plt.axis('equal')
    plt.savefig(directory+name_plot ,fmt='png')    


def plot_hist_frequency(data,directory_to_save,name_plot='plot_frequency.png'):
    '''data is np.array'''
    plt.style.use('ggplot')
    fig = plt.figure()
    plt.title(name_plot[0:-4])
    plt.xlabel("persent")
    plt.ylabel("frequency")    
    plt.hist(data, bins=100, alpha=1,color='green',edgecolor='black')
    plt.savefig(directory_to_save+name_plot ,fmt='png')

def plot_hist_frequency_for_two(data,data2,directory_to_save,name_plot='plot_frequency.png'):
    '''data is np.array'''
    plt.style.use('ggplot')
    fig = plt.figure()
    plt.title(name_plot[0:-4])
    plt.xlabel("persent")
    plt.ylabel("frequency")    
    plt.hist(data, bins=100, alpha=0.5,color='red',edgecolor='black',label='1')
    plt.hist(data2, bins=100, alpha=0.5,color='blue',edgecolor='black',label='2')    
    plt.legend()
    plt.savefig(directory_to_save+name_plot ,fmt='png')

#plot_hist_frequency_into_pie_chart(ALNDATA_DIRECTORY,'into.csv',name_plot=ALNDATA_DIRECTORY.rsplit('/',2)[-2]+'_plot_into.png')
#
data = np.random.randn(1000)
#data2 = np.random.randn(1000)

#plot_hist_frequency_into(ALNDATA_DIRECTORY,'into.csv',name_plot=ALNDATA_DIRECTORY.rsplit('/',2)[-2]+'_plot_into.png')
#plot_hist_frequency_for_two(data,data2,ALNDATA_DIRECTORY,name_plot='plot_frequency.png')


#plt.show()

print('end selseq_plot')