#!/usr/bin/env python

import os 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib
rcparams = {}
rcparams['text.usetex'] = True
rcparams['axes.linewidth'] = 0.5
rcparams['font.family'] = 'Times New Roman'
rcparams['font.size'] = 16
matplotlib.rcParams.update(rcparams)

PATH = "/Users/yasmeenasali/Files/Marka_Lab/Run_Data/ER12_O2_LHV/" #change this on your computer
PLOT_PATH = "/Users/yasmeenasali/Files/Marka_Lab/GWHEN/Plots/O2_cwb_histograms"
dirs = os.listdir(PATH)
params = {'rho': 6, 'netCC': 7, 'netED': 8, 
          'duration': 56, 'frequency': 57, 'bandwidth': 60}
events = len(dirs)

def get_parameter(dir, parameter):
    with open(f'{PATH}/{dir}/eventDump.txt', 'r') as file:
        data = file.readlines()
        for line in data:
            for value in line.split():
                if '{}:'.format(parameter) in value:
                    print(line)

def get_all_parameters(dir, params):
    with open(f'{PATH}/{dir}/eventDump.txt', 'r') as file:
        data = file.readlines()
        vals = []
        for param in params.keys():
            idx = params[param]
            vals.append(data[idx].split()[1:])
        return np.concatenate(vals)

full_data = np.zeros([events, 9]) #where 9 is num of params (including params w double values)
i = 0
for dir in dirs:
    full_data[i] = get_all_parameters(dir, params)
    i += 1

col_names = ['Rho', 'NetCC', 'NetED', 'Duration1', 'Duration2', 
             'Frequency1', 'Frequency2', 'Bandwidth1', 'Bandwidth2']
full_data = pd.DataFrame(full_data, columns = col_names) 

bins = np.linspace(5.8,10,30)
plt.hist(full_data['Rho'], bins=bins, color = 'b', alpha=0.5)
plt.xlabel('Rho')
plt.savefig(f'{PLOT_PATH}/rho_hist.pdf')

bins = np.linspace(0.2,1,30)
plt.hist(full_data['NetCC'], bins=bins, color = 'b', alpha=0.5)
plt.xlabel('NetCC')
plt.savefig(f'{PLOT_PATH}/netCC_hist.pdf')

bins = np.linspace(0,2.5,30)
plt.hist(full_data['NetED'], bins=bins, color = 'b', alpha=0.5)
plt.xlabel('NetED')
plt.savefig(f'{PLOT_PATH}/netED_hist.pdf')

bins = np.linspace(0,2000,40)
plt.hist(full_data['Frequency1'], bins=bins, color = 'b', alpha=0.5)
plt.hist(full_data['Frequency2'], bins=bins, color = 'r', alpha=0.5)
plt.xlabel('Frequency')
plt.savefig(f'{PLOT_PATH}/frequency_hist.pdf')

bins = np.linspace(0,1.1,40)
plt.hist(full_data['Duration1'], bins=bins, color = 'b', alpha=0.5)
plt.hist(full_data['Duration2'], bins=bins, color = 'r', alpha=0.5)
plt.xlabel('Duration')
plt.savefig(f'{PLOT_PATH}/duration_hist.pdf')

bins = np.linspace(0,820,40)
plt.hist(full_data['Bandwidth1'], bins=bins, color = 'b', alpha=0.5)
plt.hist(full_data['Bandwidth2'], bins=bins, color = 'r', alpha=0.5)
plt.xlabel('Bandwidth')
plt.savefig(f'{PLOT_PATH}/bandwidth_hist.pdf')


