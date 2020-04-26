#!/home/yasmeen.asali/miniconda3/envs/ligo-py36/bin/python3.6

DESC='''A module for generating full background results files for O3A/B subthreshold runs. 
Also contains function for downselecting data.

NOTE: Some of the functions in this script are broken

'''

import numpy as np
import json
import csv
import ligo.skymap.moc
from astropy.table import Table
from astropy.io import fits
import pandas as pd

#Note this function needs to be adapted to work for both .npy and .npz 
def generate_full_bkg(O3A_file, bkg_file, PATH, RUN):
    '''
    Read in the full O3A data from 'O3A_cbc_events.npy' and cross match with a .npz format background results file.
    '''
    superevent_dataset = np.load(O3A_file, allow_pickle=True)
    superevent_id = superevent_dataset[0]
    ifos = superevent_dataset[1]
    tag = superevent_dataset[2]
    pipeline = superevent_dataset[3]
    snr = superevent_dataset[4]
    far = superevent_dataset[5]
    skyarea = superevent_dataset[6]
    pter = superevent_dataset[7]
    dist = superevent_dataset[8]
    
    bkg_dataset = np.load(bkg_file) #requires npz format
    odds_ratio = bkg_dataset['Odds_Ratio']
    neutrino_count = bkg_dataset['Neutrino_Count']
    bkg_superevent_id = bkg_dataset['Superevent_ID']
    
    bkg_pipeline = []
    bkg_ifos = []
    bkg_tag = []
    bkg_snr = []
    bkg_far = []
    bkg_skyarea = []
    bkg_pter = []
    bkg_dist = []
    i = 0
    for element in bkg_superevent_id:
        idx = np.where([i == element for i in superevent_id])[0]
        bkg_pipeline.append(pipeline[idx][0]) 
        bkg_snr.append(float(snr[idx][0]))
        bkg_far.append(float(far[idx][0]))
        bkg_skyarea.append(float(skyarea[idx][0]))
        bkg_pter.append(float(pter[idx][0]))
        bkg_dist.append(float(dist[idx][0]))
        bkg_ifos.append(ifos[idx][0])
        bkg_tag.append(tag[idx][0])
        i += 1
        if i % 10000 == 0:
            print('Iteration number {} of {}'.format(i, len(bkg_superevent_id)))
    
    kwargs = {
            "Superevent_ID": bkg_superevent_id, 
            "Odds_Ratio": odds_ratio, 
            "Neutrino_Count": neutrino_count, 
            "Pipeline": bkg_pipeline,
            "SNR": bkg_snr, 
            "FAR": bkg_far,
            "Sky_Area": bkg_skyarea, 
            "p_Terrestrial": bkg_pter, 
            "Mean_Distance": bkg_dist, 
            "IFOs": bkg_ifos,
            "Tag": bkg_tag
        }
    
    np.savez(f'{PATH}/bkg_{RUN}/subthresh_full_results.npz', **kwargs)    

#note: keys are diff for O3B
def downselect_npz(results_file, parameter='FAR', farcut=False, ptercut1=False, ptercut2=False):
    data = np.load(results_file, allow_pickle=True)
    superevent_id = data['Superevent'] 
    #odds_ratio = data['Odds_Ratio'] 
    #neutrinos = data['Neutrino_Count'] 
    pipeline = data['Pipeline'] 
    snr = data['SNR']
    far = data['FAR'] 
    skyarea = data['Skyarea'] 
    pter = data['p_Terrestrial'] 
    dist = data['Distance'] 
    #ifos = data['IFOs']
    #tag = data['Tag']

    #set lower bound on odds ratio at 10^-30
    #odds_ratio[odds_ratio<10**-30] = 10**-30

    if farcut == True:
        #cut on FAR at 1/day
        farcut=(1.1574 * 10 ** (-5))
        idx_farcut = np.where([i < farcut for i in far])
        far = far[idx_farcut]
        snr = snr[idx_farcut]
        pipeline = pipeline[idx_farcut]
        superevent_id = superevent_id[idx_farcut]
        skyarea = skyarea[idx_farcut]
        pter = pter[idx_farcut]
        dist = dist[idx_farcut]
        #odds_ratio = odds_ratio[idx_farcut]
        #neutrinos = neutrinos[idx_farcut]
        #ifos = ifos[idx_farcut]
        #tag = tag[idx_farcut]
    
    if ptercut1 == True:
        #cut on pter
        ptercut=(10 ** (-1))
        idx_ptercut = np.where([i < ptercut for i in pter])
        far = far[idx_ptercut]
        snr = snr[idx_ptercut]
        pipeline = pipeline[idx_ptercut]
        superevent_id = superevent_id[idx_ptercut]
        skyarea = skyarea[idx_ptercut]
        pter = pter[idx_ptercut]
        dist = dist[idx_ptercut]
        #odds_ratio = odds_ratio[idx_ptercut]
        #neutrinos = neutrinos[idx_ptercut]
        #ifos = ifos[idx_ptercut]
        #tag = tag[idx_ptercut]
    
    if ptercut2 == True:
        #cut on pter
        ptercut=(10 ** (-1))
        idx_ptercut = np.where([i > ptercut for i in pter])
        far = far[idx_ptercut]
        snr = snr[idx_ptercut]
        pipeline = pipeline[idx_ptercut]
        superevent_id = superevent_id[idx_ptercut]
        skyarea = skyarea[idx_ptercut]
        pter = pter[idx_ptercut]
        dist = dist[idx_ptercut]
        #odds_ratio = odds_ratio[idx_ptercut]
        #neutrinos = neutrinos[idx_ptercut]
        #ifos = ifos[idx_ptercut]
        #tag = tag[idx_ptercut]
    
    if parameter == 'FAR':
        dat = far
    elif parameter == 'SNR':
        dat = snr
    #elif parameter == 'Odds Ratio':
    #    dat = odds_ratio
    #elif parameter == 'Neutrino Count':
    #    dat = neutrino_count
    elif parameter == 'Sky Area':
        dat = skyarea
    elif parameter == 'p Terrestrial':
        dat = pter
    elif parameter == 'Distance':
        dat = dist

    return dat

def downselect_npz_pipeline(results_file, parameter='FAR', farcut=False, oddscut=False):
    data = np.load(results_file, allow_pickle=True)
    superevent_id = data['Superevent'] 
    #odds_ratio = data['Odds_Ratio'] 
    #neutrinos = data['Neutrino_Count'] 
    pipeline = data['Pipeline'] 
    snr = data['SNR']
    far = data['FAR'] 
    skyarea = data['Skyarea'] 
    pter = data['p_Terrestrial'] 
    dist = data['Distance'] 
    #ifos = data['IFOs']
    #tag = data['Tag']

    #set lower bound on odds ratio at 10^-30
    #odds_ratio[odds_ratio<10**-30] = 10**-30

    if farcut == True:
        #cut on FAR at 1/day
        farcut=(1.1574 * 10 ** (-5))
        idx_farcut = np.where([i < farcut for i in far])
        far = far[idx_farcut]
        snr = snr[idx_farcut]
        pipeline = pipeline[idx_farcut]
        superevent_id = superevent_id[idx_farcut]
        skyarea = skyarea[idx_farcut]
        pter = pter[idx_farcut]
        dist = dist[idx_farcut]
        #odds_ratio = odds_ratio[idx_farcut]
        #neutrinos = neutrinos[idx_farcut]
        #ifos = ifos[idx_farcut]
        #tag = tag[idx_farcut]
        
    if oddscut == True:
        #cut on odds ratio (100 largest elements)
        oddscut=odds_ratio[np.argsort(odds_ratio)[-100]]
        idx_oddscut = np.where([i > oddscut for i in odds_ratio])
        far = far[idx_oddscut]
        snr = snr[idx_oddscut]
        pipeline = pipeline[idx_oddscut]
        superevent_id = superevent_id[idx_oddscut]
        skyarea = skyarea[idx_oddscut]
        pter = pter[idx_oddscut]
        dist = dist[idx_oddscut]
        #odds_ratio = odds_ratio[idx_oddscut]
        #neutrinos = neutrinos[idx_oddscut]
        #ifos = ifos[idx_oddscut]
        #tag = tag[idx_oddscut]

    idx_gstlal = np.where([i == 'gstlal' for i in pipeline])
    idx_pycbc = np.where([i == 'pycbc' for i in pipeline])
    idx_spiir = np.where([i == 'spiir' for i in pipeline])
    idx_MBTAOnline = np.where([i == 'MBTAOnline' for i in pipeline])
    
    if parameter == 'FAR':
        dat_gstlal = far[idx_gstlal]
        dat_pycbc = far[idx_pycbc]
        dat_spiir = far[idx_spiir]
        dat_MBTAOnline = far[idx_MBTAOnline]
    elif parameter == 'SNR':
        dat_gstlal = snr[idx_gstlal]
        dat_pycbc = snr[idx_pycbc]
        dat_spiir = snr[idx_spiir]
        dat_MBTAOnline = snr[idx_MBTAOnline]
    #elif parameter == 'Odds Ratio':
    #    dat_gstlal = odds_ratio[idx_gstlal]
    #    dat_pycbc = odds_ratio[idx_pycbc]
    #    dat_spiir = odds_ratio[idx_spiir]
    #    dat_MBTAOnline = odds_ratio[idx_MBTAOnline]
    #elif parameter == 'Neutrino Count':
    #    dat_gstlal = neutrinos[idx_gstlal]
    #    dat_pycbc = neutrinos[idx_pycbc]
    #    dat_spiir = neutrinos[idx_spiir]
    #    dat_MBTAOnline = neutrinos[idx_MBTAOnline]
    elif parameter == 'Sky Area':
        dat_gstlal = skyarea[idx_gstlal]
        dat_pycbc = skyarea[idx_pycbc]
        dat_spiir = skyarea[idx_spiir]
        dat_MBTAOnline = skyarea[idx_MBTAOnline]
    elif parameter == 'p Terrestrial':
        dat_gstlal = pter[idx_gstlal]
        dat_pycbc = pter[idx_pycbc]
        dat_spiir = pter[idx_spiir]
        dat_MBTAOnline = pter[idx_MBTAOnline]
    elif parameter == 'Distance':
        dat_gstlal = dist[idx_gstlal]
        dat_pycbc = dist[idx_pycbc]
        dat_spiir = dist[idx_spiir]
        dat_MBTAOnline = dist[idx_MBTAOnline]
    else:
        print("Invalid Parameter. Must be either FAR, SNR, Odds Ratio, Neutrino Count, Sky Area, p Terrestrial, or Distance")
        dat_gstlal = 0
        dat_pycbc = 0
        dat_spiir = 0
        dat_MBTAOnline = 0

    return dat_gstlal, dat_pycbc, dat_spiir, dat_MBTAOnline

if __name__ == "__main__": 
    print('Generating O3A Subthreshold Background File')
    PATH = '/home/yasmeen.asali/GWHEN/O3A_subthreshold/data'
    RUN = 'Feb2020'
    O3A_file = f'{PATH}/O3A_cbc_data.npy'
    bkg_file = f'{PATH}/bkg_{RUN}/subthresh_significance_outputs.npz'
    generate_full_bkg(O3A_file, bkg_file, PATH, RUN)
