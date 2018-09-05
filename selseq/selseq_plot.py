from selseq_constant import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

print('start selseq_plot')

def plot_hist_frequency_into(directory,into,name_plot='plot_into.png'):
    csv_data = pd.read_csv(directory + into)
    data_sum_row = csv_data.sum(axis=1)
    fig = plt.figure()
    plt.hist(data_sum_row.values, bins=100, alpha=1,color='blue',edgecolor='black')
    plt.savefig(directory+name_plot ,fmt='png')
    #plt.show()
#plot_hist_frequency_into(ALNDATA_DIRECTORY,'into.csv')

print('end selseq_plot')