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
        bins_dat = np.logspace(np.log10(10**-10),np.log10(1),50)
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
    elif parameter == 'Amplitude':
        bins_dat = 30
        x_scale, y_scale = 'linear', 'linear'
    elif parameter == 'Central Frequency':
        bins_dat = 30
        x_scale, y_scale = 'linear', 'linear'
    elif parameter == 'Bandwidth':
        bins_dat = 30
        x_scale, y_scale = 'linear', 'linear'
    elif parameter == 'Duration':
        bins_dat = 30
        x_scale, y_scale = 'linear', 'linear'
    else:
        print("Invalid Parameter")
        exit()

    return bins_dat, x_scale, y_scale

def plot_O3A_v_O3B(O3A_file, O3B_file, parameter = 'SNR'):
    dat_O3A, dat_O3B = get_data.downselect_O3A_O3B_files(O3A_file, O3B_file, parameter = parameter)

    bins_dat, x_scale, y_scale = best_scale(parameter) 

    dens = False
    plt.figure(1)
    plt.hist(dat_O3A, bins=bins_dat, color='b', alpha=0.5, density=dens, label="O3A Events\n({})".format(len(dat_O3A)))
    plt.hist(dat_O3B, bins=bins_dat, color='r', alpha=0.5, density=dens, label="O3B Events\n({})".format(len(dat_O3B)))
    plt.xlabel(parameter)
    plt.ylabel('Number of Events')
    plt.xscale(x_scale)
    plt.yscale(y_scale)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/O3_subthreshold/O3A_vs_O3B_{parameter}.pdf')
    plt.close(1)

    
    dat_far_O3A, dat_far_O3B = get_data.downselect_O3A_O3B_files(O3A_file, O3B_file, parameter = parameter, farcut=True)

    plt.figure(2)
    plt.hist(dat_far_O3A, bins=bins_dat, color='b', alpha=0.5, density=dens, label="O3A with FAR $<$ 1/day\n({})".format(len(dat_far_O3A)))
    plt.hist(dat_far_O3B, bins=bins_dat, color='r', alpha=0.5, density=dens, label="O3B with FAR $<$ 1/day\n({})".format(len(dat_far_O3B)))
    plt.xlabel(parameter)
    plt.ylabel('Number of Events')
    plt.xscale(x_scale)
    plt.yscale(y_scale)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/O3_subthreshold/O3A_vs_O3B_{parameter}_with_far_cut.pdf')
    plt.close(2)

def plot_O3A_v_O3B_cwb(O3A_file, O3B_file, parameter = 'SNR'):
    dat_O3A = get_data.downselect_npz_cwb(O3A_file, parameter = parameter)
    dat_O3B = get_data.downselect_npz_cwb(O3B_file, parameter = parameter)

    bins_dat, x_scale, y_scale = best_scale(parameter) 

    dens = False
    plt.figure(1)
    plt.hist(dat_O3A, bins=bins_dat, color='b', alpha=0.5, density=dens, label="O3A Events\n({})".format(len(dat_O3A)))
    plt.hist(dat_O3B, bins=bins_dat, color='r', alpha=0.5, density=dens, label="O3B Events\n({})".format(len(dat_O3B)))
    plt.xlabel(parameter)
    plt.ylabel('Number of Events')
    plt.xscale(x_scale)
    plt.yscale(y_scale)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/CWB/O3A_vs_O3B_{parameter}.pdf')
    plt.close(1)

    dat_far_O3A = get_data.downselect_npz_cwb(O3A_file, parameter = parameter, farcut=True)
    dat_far_O3B = get_data.downselect_npz_cwb(O3B_file, parameter = parameter, farcut=True)

    plt.figure(2)
    plt.hist(dat_far_O3A, bins=bins_dat, color='b', alpha=0.5, density=dens, label="O3A with FAR $<$ 1/day\n({})".format(len(dat_far_O3A)))
    plt.hist(dat_far_O3B, bins=bins_dat, color='r', alpha=0.5, density=dens, label="O3B with FAR $<$ 1/day\n({})".format(len(dat_far_O3B)))
    plt.xlabel(parameter)
    plt.ylabel('Number of Events')
    plt.xscale(x_scale)
    plt.yscale(y_scale)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/CWB/O3A_vs_O3B_{parameter}_with_far_cut.pdf')
    plt.close(2)


def plot_parameters(numpy_file, parameter = 'SNR', background=False):
    if background == True:
        dat = get_data.downselect_npz(numpy_file, parameter = parameter, background=True)
        dat_far = get_data.downselect_npz(numpy_file, parameter = parameter, farcut = True, background=True)
    else:
        dat = get_data.downselect_npz(numpy_file, parameter = parameter)
        dat_far = get_data.downselect_npz(numpy_file, parameter = parameter, farcut = True)

    bins_dat, x_scale, y_scale = best_scale(parameter) 

    dens = False
    plt.figure(1)
    plt.hist(dat, bins=bins_dat, color='b', alpha=0.5, density=dens, label="All Events\n({})".format(len(dat)))
    plt.hist(dat_far, bins=bins_dat, color='r', alpha=0.5, density=dens, label="FAR $<$ 1/day\n({})".format(len(dat_far)))
    plt.xlabel(parameter)
    plt.ylabel('Number of Events')
    plt.xscale(x_scale)
    plt.yscale(y_scale)
    plt.grid()
    plt.legend()
    if background == True:
        plt.title('O3A Subthreshold Background Results')
        plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/O3A_bkg_{RUN_BKG}/with_titles/{parameter}_with_far_cut_full_results.pdf')
    else:
        plt.tight_layout()
        #plt.title('O3 Subthreshold Results')
        plt.savefig(f'/home/yasmeen.asali/public_html/GWHEN/O3_subthreshold/{parameter}_with_far_cut_tight_layout.pdf')
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
    PATH_BKG ='/home/yasmeen.asali/GWHEN/O3A_subthreshold/data'
    RUN_BKG = 'Feb2020'
    numpy_file_bkg = f'{PATH_BKG}/bkg_{RUN_BKG}/subthresh_full_results.npz' 

    PATH ='/home/yasmeen.asali/GWHEN/analysis'
    numpy_file = f'{PATH}/O3_full_run_data.npz' 
    O3A_cwb = f'{PATH}/O3A_cwb_full_data.npz' 
    O3B_cwb = f'{PATH}/O3B_cwb_full_data.npz' 
    
    PATH_O3A = '/home/yasmeen.asali/GWHEN/O3A_subthreshold/data'
    O3A_file = f'{PATH_O3A}/O3A_cbc_data.npy'
    PATH_O3B = '/home/yasmeen.asali/GWHEN/O3B_subthreshold/data'
    O3B_file = f'{PATH_O3B}/O3B_cbc_full_data.npz'    

    plot_O3A_v_O3B_cwb(O3A_cwb, O3B_cwb, parameter = 'SNR')
    plot_O3A_v_O3B_cwb(O3A_cwb, O3B_cwb, parameter = 'Amplitude')
    plot_O3A_v_O3B_cwb(O3A_cwb, O3B_cwb, parameter = 'Central Frequency')
    plot_O3A_v_O3B_cwb(O3A_cwb, O3B_cwb, parameter = 'Duration')
    plot_O3A_v_O3B_cwb(O3A_cwb, O3B_cwb, parameter = 'Bandwidth')

    #plot_O3A_v_O3B(O3A_file, O3B_file, parameter = 'SNR')
    #plot_O3A_v_O3B(O3A_file, O3B_file, parameter = 'Sky Area')
    #plot_O3A_v_O3B(O3A_file, O3B_file, parameter = 'p Terrestrial')
    #plot_O3A_v_O3B(O3A_file, O3B_file, parameter = 'Distance')

    #plot_parameters(numpy_file, parameter = 'SNR')
    #plot_parameters(numpy_file, parameter = 'Sky Area')
    #plot_parameters(numpy_file, parameter = 'p Terrestrial')
    #plot_parameters(numpy_file, parameter = 'Distance')

    #plot_parameters(numpy_file, parameter = 'Odds Ratio', background = True)
    #plot_parameters(numpy_file, parameter = 'Neutrino Count', background = True)

    #pipeline_subplots(numpy_file, RUN, parameter = 'SNR')
    #pipeline_subplots(numpy_file, RUN, parameter = 'Sky Area')
    #pipeline_subplots(numpy_file, RUN, parameter = 'p Terrestrial')
    #pipeline_subplots(numpy_file, RUN, parameter = 'Distance')

    #pipeline_subplots(numpy_file, RUN, parameter = 'Odds Ratio')
    #pipeline_subplots(numpy_file, RUN, parameter = 'Neutrino Count')

    #plot_odds_ratio(numpy_file, RUN)
    
    #plot_pter_high_odds(RUN)
