#!/home/yasmeen.asali/miniconda3/envs/ligo-py36/bin/python3.6

DESC = """
A plotting module 
"""

#TODO FIX ARGPARSE
import argparse 
parser = argparse.ArgumentParser(description=DESC)
parser.add_argument('-p', '--parameter', type=str, help=("Generate Subplots for a specified parameter"))
args = parser.parse_args()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import get_data

rcparams = {}
rcparams['text.usetex'] = True
rcparams['axes.linewidth'] = 0.5
rcparams['font.family'] = 'Times New Roman'
rcparams['font.size'] = 16
matplotlib.rcParams.update(rcparams)
   

def best_scale(parameter):
    if parameter == 'Odds Ratio':
        bins_dat = np.logspace(np.log10(10**-30), np.log10(0.002), 50)
        x_scale = 'log'
        y_scale = 'log'
    elif parameter == 'SNR':
        bins_dat = np.linspace(0, 50, 50)
        x_scale = 'linear'
        y_scale = 'log'
    elif parameter == 'Sky Area':
        bins_dat = np.logspace(np.log10(0.1), np.log10(36000), 50)
        x_scale = 'log'
        y_scale = 'log'
    elif parameter == 'p Terrestrial':
        bins_dat = np.logspace(np.log10(10**-10),np.log10(0.1),50)
        x_scale = 'log'
        y_scale = 'log'
    elif parameter == 'Distance':
        bins_dat = np.logspace(np.log10(0.1), np.log10(8600), 50)
        x_scale = 'log'
        y_scale = 'log'
    elif parameter == 'Neutrino Count':
        bins_dat = np.linspace(0, 22, 50) 
        x_scale = 'linear'
        y_scale = 'log'
    else:
        print("Invalid Parameter")
        exit()

    return bins_dat, x_scale, y_scale

def plot_parameters(numpy_file, RUN, parameter = 'SNR'):
    dat = get_data.downselect_npz(numpy_file, parameter = parameter)
    dat_far = get_data.downselect_npz(numpy_file, parameter = parameter, farcut = True)

    bins_dat, x_scale, y_scale = best_scale(parameter) 

    dens = False
    plt.figure(1)
    plt.hist(dat, bins=bins_dat, color='b', alpha=0.5, density=dens, label="{} Events".format(len(dat)))
    plt.hist(dat_far, bins=bins_dat, color='r', alpha=0.5, density=dens, label="{} Events".format(len(dat_far)))
    plt.xlabel(parameter)
    plt.ylabel('Number of Events')
    plt.xscale(x_scale)
    plt.yscale(y_scale)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/{RUN}/{parameter}_with_far_cut.pdf')
    plt.close(1)

 
def pipeline_subplots(numpy_file, RUN, parameter = 'SNR'):
    
    dat_gstlal, dat_pycbc, dat_spiir, dat_MBTAOnline = get_data.downselect_npz_pipeline(numpy_file, parameter = parameter)
    dat_far_gstlal, dat_far_pycbc, dat_far_spiir, dat_far_MBTAOnline = get_data.downselect_npz_pipeline(numpy_file, parameter = parameter, farcut = True)
    
    fig, axs = plt.subplots(4, 2, figsize = (8, 12))
    bins_dat, x_scale, y_scale = best_scale(parameter) 

    for (m,n), subplot in np.ndenumerate(axs):
        subplot.set_xscale(x_scale)
        subplot.set_yscale(y_scale)
    
    dens = False
        
    axs[0,0].hist(dat_gstlal, bins=bins_dat, color='b', alpha=0.5, density=dens)
    axs[0,0].set_xlabel(parameter)
    axs[0,0].set_title('gstlal ({} Events)'.format(len(dat_gstlal)))

    axs[0,1].hist(dat_far_gstlal, bins=bins_dat, color='b', alpha=0.5, density=dens)
    axs[0,1].set_xlabel(parameter)
    axs[0,1].set_title('gstlal FAR cut ({} Events)'.format(len(dat_far_gstlal)))
    
    axs[1,0].hist(dat_pycbc, bins=bins_dat, color='c', alpha=0.5, density=dens)
    axs[1,0].set_xlabel(parameter)
    axs[1,0].set_title('pycbc ({} Events)'.format(len(dat_pycbc)))
    
    axs[1,1].hist(dat_far_pycbc, bins=bins_dat, color='c', alpha=0.5, density=dens)
    axs[1,1].set_xlabel(parameter)
    axs[1,1].set_title('pycbc FAR cut ({} Events)'.format(len(dat_far_pycbc)))
    
    axs[2,0].hist(dat_spiir, bins=bins_dat, color='r', alpha=0.5, density=dens)
    axs[2,0].set_xlabel(parameter)
    axs[2,0].set_title('spiir ({} Events)'.format(len(dat_spiir)))
    
    axs[2,1].hist(dat_far_spiir, bins=bins_dat, color='r', alpha=0.5, density=dens)
    axs[2,1].set_xlabel(parameter)
    axs[2,1].set_title('spiir FAR cut ({} Events)'.format(len(dat_far_spiir)))
    
    axs[3,0].hist(dat_MBTAOnline, bins=bins_dat, color='g', alpha=0.5, density=dens)
    axs[3,0].set_xlabel(parameter)
    axs[3,0].set_title('MBTAOnline ({} Events)'.format(len(dat_MBTAOnline)))
    
    axs[3,1].hist(dat_far_MBTAOnline, bins=bins_dat, color='g', alpha=0.5, density=dens)
    axs[3,1].set_xlabel(parameter)
    axs[3,1].set_title('MBTAOnline FAR cut ({} Events)'.format(len(dat_far_MBTAOnline)))

    plt.tight_layout()
    plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/{RUN}/{parameter}_by_pipeline_zoom.pdf')
    plt.close()

def plot_odds_ratio(np_file, RUN):
    odds_ratio = get_data.downselect_npz(np_file, parameter='Odds Ratio', ptercut1=True)    
    odds_ratio_ptercut = get_data.downselect_npz(np_file, parameter='Odds Ratio', ptercut2=True)    

    cut1 = 'p Terr $< 10^{-1}$'
    cut2 = 'p Terr $> 10^{-1}$'

    bins_dat = np.logspace(np.log10(10**-30), np.log10(0.002), 50)

    plt.figure(1)
    plt.hist(odds_ratio, bins=bins_dat, color='b', alpha=0.5, density=False, label="{} Events\n({})".format(len(odds_ratio), cut1))
    plt.hist(odds_ratio_ptercut, bins=bins_dat, color='r', alpha=0.5, density=False, label="{} Events\n({})".format(len(odds_ratio_ptercut), cut2))
    plt.xlabel(r'Odds Ratio')
    plt.ylabel('Number of Events')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/{RUN}/Odds_Ratio_ptercut_10e-1_above_and_below.pdf')
    plt.close(1)

def plot_pter_high_odds(RUN):
    parameter="p Terrestrial"

    dat_gstlal, dat_pycbc, dat_spiir, dat_MBTAOnline = get_data.downselect_npz_pipeline(numpy_file, parameter = parameter, oddscut = True)

    fig, axs = plt.subplots(2, 2, figsize = (8, 8))

    bins_dat = np.logspace(np.log10(10**-15), np.log10(1), 50)
    x_scale = 'log'
    y_scale = 'log'

    for (m,n), subplot in np.ndenumerate(axs):
        subplot.set_xscale(x_scale)
        subplot.set_yscale(y_scale)

    axs[0,0].hist(dat_gstlal, bins=bins_dat, color='b', alpha=0.5)
    axs[0,0].set_xlabel(parameter)
    axs[0,0].set_title('gstlal FAR+OR cut ({} Events)'.format(len(dat_gstlal)))
    
    axs[0,1].hist(dat_pycbc, bins=bins_dat, color='c', alpha=0.5)
    axs[0,1].set_xlabel(parameter)
    axs[0,1].set_title('pycbc FAR+OR cut ({} Events)'.format(len(dat_pycbc)))
    
    axs[1,0].hist(dat_spiir, bins=bins_dat, color='r', alpha=0.5)
    axs[1,0].set_xlabel(parameter)
    axs[1,0].set_title('spiir FAR+OR cut ({} Events)'.format(len(dat_spiir)))
    
    axs[1,1].hist(dat_MBTAOnline, bins=bins_dat, color='g', alpha=0.5)
    axs[1,1].set_xlabel(parameter)
    axs[1,1].set_title('MBTAOnline FAR+OR cut ({} Events)'.format(len(dat_MBTAOnline)))

    plt.tight_layout()
    plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/bkg_{RUN}/p_Terrestiral_by_pipeline_highest_Odds_Ratios.pdf')
    

if __name__ == "__main__":
    PATH='/home/yasmeen.asali/GWHEN/O3B_subthreshold/data'
    RUN = 'O3B_subthreshold'
    #numpy_file = f'{PATH}/bkg_{RUN}/subthresh_full_results.npz' 
    numpy_file = f'{PATH}/O3B_cbc_full_data.npz' 
    
    #plot_parameters(numpy_file, RUN, parameter = 'SNR')
    #plot_parameters(numpy_file, RUN, parameter = 'Sky Area')
    #plot_parameters(numpy_file, RUN, parameter = 'p Terrestrial')
    #plot_parameters(numpy_file, RUN, parameter = 'Distance')

    #pipeline_subplots(numpy_file, RUN, parameter = 'SNR')
    #pipeline_subplots(numpy_file, RUN, parameter = 'Sky Area')
    pipeline_subplots(numpy_file, RUN, parameter = 'p Terrestrial')
    #pipeline_subplots(numpy_file, RUN, parameter = 'Distance')

    #pipeline_subplots(numpy_file, RUN, parameter = 'Odds Ratio')
    #pipeline_subplots(numpy_file, RUN, parameter = 'Neutrino Count')

    #plot_odds_ratio(numpy_file, RUN)
    
    #plot_pter_high_odds(RUN)
