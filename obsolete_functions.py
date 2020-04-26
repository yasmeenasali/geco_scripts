#!/home/yasmeen.asali/miniconda3/envs/ligo-py36/bin/python3.6

DESC='''A module for searching through O3 subthreshold data.

NOTE: This file contains obsolete functions to parse and plot data from O3A/B. 
Certain imports might be missing, and most of the functions likely won't run. 
'''

import numpy as np
import json
import csv
import ligo.skymap.moc
from astropy.table import Table
from astropy.io import fits
import pandas as pd

def generate_has_remnant_file(superevents_cbc):
    superevent_dataset = []
    for superevent in superevents_cbc:
        with open("superevents/{}/p_astro.json".format(superevent), "r") as read_file:
            dataP = json.load(read_file)
            pBBH = dataP['BBH']
        try:
            with open("superevents/{}/em_bright.json".format(superevent), "r") as read_file:
                dataEM = json.load(read_file)
                em = dataEM['HasRemnant']
        except FileNotFoundError:
            print('EM Bright File not Found for {}'.format(superevent)) 
        superevent_dataset.append([superevent, float(pBBH), float(em)])
    filename = 'cbc_events_has_remnant.npy'
    np.save('data/{}'.format(filename), superevent_dataset)
    return filename

def get_tag(value, dataP):
    tag = list(dataP.keys())[list(dataP.values()).index(value)]
    return tag 

def generate_IFOs_tags_file(superevents_cbc):
    superevent_dataset =[]
    for superevent in superevents_cbc:
        with open("superevents/{}/preferred_event_data.json".format(superevent), "r") as read_file:
            data = json.load(read_file)
            IFOs = data['extra_attributes']['CoincInspiral']['ifos'].split(',')
        with open("superevents/{}/p_astro.json".format(superevent), "r") as read_file:
            dataP = json.load(read_file)
            tag = get_tag(max(dataP.values()), dataP)
            if tag == 'Terrestrial':
                max2 = list(sorted(dataP.values()))[-2]
                tag = get_tag(max2, dataP)
        superevent_dataset.append([superevent, IFOs, tag])
        print('Finished {}'.format(superevent))
    filename = 'superevent_ifos_tag.npy'
    np.save(filename, superevent_dataset)
    return filename

def sort_npy_file(superevent_list, filename, col, num=False):
    dset = np.load(filename, allow_pickle=True)
    superevent_unsorted = dset[:,0]
    data_unsorted = dset[:,col]
    data = []
    for element in superevent_list:
        idx_A = np.where([i == element for i in superevent_unsorted])
        if num == True:
            data.append(float(data_unsorted[idx_A]))
        else:
            data.append(data_unsorted[idx_A])
    return data

def cross_match_O3A_subthresh_files(superevent_file, ifos_file, skyarea_file):
    superevent_dataset = np.load(superevent_file, allow_pickle=True)
    superevent_id = superevent_dataset[:,0].tolist()
    snr_str = superevent_dataset[:,1].tolist()
    far_str = superevent_dataset[:,2].tolist()
    pipeline = superevent_dataset[:,3].tolist()
    pter_str = superevent_dataset[:,4].tolist()
    dist_str = superevent_dataset[:,5].tolist()
   
    skyarea = sort_npy_file(superevent_id, skyarea_file, 4, num=True)
    ifos_list = sort_npy_file(superevent_id, ifos_file, 1)
    tag_list = sort_npy_file(superevent_id, ifos_file, 2)
    
    ifos = []
    tag = []
    snr = []
    far = []
    pter = []
    dist = []
    for i in range(0, len(ifos_list)):
        ifos.append(ifos_list[i][0])
        tag.append(str(tag_list[i][0]))
        snr.append(float(snr_str[i]))
        far.append(float(far_str[i]))
        pter.append(float(pter_str[i]))
        dist.append(float(dist_str[i]))

    master_header = ['Superevent_ID', 'IFOs', 'Tag', 'Pipeline', 'SNR', 'FAR', 'Sky_Area', 'p_Terrestrial', 'Mean_Distance']
    np.save('O3A_cbc_data.meta', master_header)
    master_array = [superevent_id, ifos, tag, pipeline, snr, far, skyarea, pter, dist]
    np.save('O3A_cbc_data.npy', master_array)

def cross_match_bkg_subthresh_files(superevent_file, coinc_neutrino_file, skyarea_file):
    '''
    NOTE: THIS FUNCTION IS LARGELY OUT OF DATE
    background data is now saved in a different format than what this function assumes
    '''
    superevent_dataset = np.load(superevent_file, allow_pickle=True)
    superevent_id = superevent_dataset[:,0]
    snr = superevent_dataset[:,1]
    far = superevent_dataset[:,2]
    pipeline = superevent_dataset[:,3]
    pter = superevent_dataset[:,4]
    dist = superevent_dataset[:,5]

    skyarea_dset = np.load(skyarea_file, allow_pickle=True)
    superevent_unsorted = skyarea_dset[:,0]
    skyarea_unsorted = skyarea_dset[:,4]
    skyarea = []
    for element in superevent_id:
        idx_A = np.where([i == element for i in superevent_unsorted])
        skyarea.append(float(skyarea_unsorted[idx_A]))
    skyarea = np.array(skyarea)
    
    lines = open(coinc_neutrino_file, 'r').readlines()
    
    odds_ratio = []
    neutrino_count = []
    bkg_superevent_id = []
    for i in range(0, len(lines)):
        sid = lines[i].split()[2][26:].split('-', 1)[0]
        bkg_superevent_id.append(sid)
        odds = float(lines[i].split()[0])
        odds_ratio.append(odds)
        count = float(lines[i].split()[1])
        neutrino_count.append(count)
    
    bkg_pipeline = []
    bkg_snr = []
    bkg_far = []
    bkg_skyarea = []
    bkg_pter = []
    bkg_dist = []
    i = 0
    for element in bkg_superevent_id:
        idx = np.where([i == element for i in superevent_id])[0]
        bkg_pipeline.append(pipeline[idx]) 
        bkg_snr.append(float(snr[idx][0]))
        bkg_far.append(float(far[idx][0]))
        bkg_skyarea.append(float(skyarea[idx][0]))
        bkg_pter.append(float(pter[idx][0]))
        bkg_dist.append(float(dist[idx][0]))
        i += 1
        if i % 10000 == 0:
            print('Iteration number {} of {}'.format(i, len(bkg_superevent_id)))
    
    dataframe = pd.DataFrame(
        {
            "Superevent_ID": bkg_superevent_id, 
            "Odds_Ratio": odds_ratio, 
            "Neutrino_Count": neutrino_count, 
            "Pipeline": bkg_pipeline,
            "SNR": bkg_snr, 
            "FAR": bkg_far,
            "Sky_Area": bkg_skyarea, 
            "p_Terrestrial": bkg_pter, 
            "Mean_Distance": bkg_dist, 
        }
    )
    
    dataframe.to_pickle("bkg_coinc_data.pkl")
    
    master_header = ['Superevent_ID', 'Odds_Ratio', 'Neutrino_Count', 'Pipeline', 'SNR', 'FAR', 'Sky_Area', 'p_Terrestrial', 'Mean_Distance']
    np.save('bkg_coinc_data.meta', master_header)
    master_array = [bkg_superevent_id, odds_ratio, neutrino_count, bkg_pipeline, bkg_snr, bkg_far, bkg_skyarea, bkg_pter, bkg_dist]
    np.save('bkg_coinc_data.npy', master_array)

def downselect_pter_O3A(farcut=False):
    superevent_dataset = np.load('superevent_snr_far_pter_dist.npy', allow_pickle=True)
    superevent_id = superevent_dataset[:,0]
    far = superevent_dataset[:,2]
    pipeline = superevent_dataset[:,3]
    pter = superevent_dataset[:,4]
   
    pter = np.array([float(i) for i in pter])
    far = np.array([float(i) for i in far])

    pter[pter<10**-10] = 10**-10
    
    if farcut == True:
        #cut on FAR
        farcut=(1.1574 * 10 ** (-5))
        idx_farcut = np.where([i < farcut for i in far])
        pipeline = pipeline[idx_farcut]
        superevent_id = superevent_id[idx_farcut]
        pter = pter[idx_farcut]
    
    idx_gstlal = np.where([i == 'gstlal' for i in pipeline])
    idx_pycbc = np.where([i == 'pycbc' for i in pipeline])
    idx_spiir = np.where([i == 'spiir' for i in pipeline])
    idx_MBTAOnline = np.where([i == 'MBTAOnline' for i in pipeline])
        
    dat_gstlal = pter[idx_gstlal]
    dat_pycbc = pter[idx_pycbc]
    dat_spiir = pter[idx_spiir]
    dat_MBTAOnline = pter[idx_MBTAOnline]

    return dat_gstlal, dat_pycbc, dat_spiir, dat_MBTAOnline

def downselect_odds_bkg(numpy_file, farcut=False):
    data = np.load(numpy_file, allow_pickle=True)
    odds_ratio = data[1]
    far = data[5]

    odds_ratio[odds_ratio<10**-30] = 10**-30

    if farcut == True:
        #cut on FAR
        farcut=(1.1574 * 10 ** (-5))
        idx_farcut = np.where([i < farcut for i in far])
        odds_ratio = odds_ratio[idx_farcut]

    return odds_ratio

def downselect_npy(numpy_file, parameter='FAR', farcut=False):
    data = np.load(numpy_file, allow_pickle=True)
    superevent_id = data[0]
    odds_ratio = data[1]
    neutrinos = data[2]
    pipeline = data[3]
    snr = data[4]
    far = data[5]
    skyarea = data[6]
    pter = data[7]
    dist = data[8]

    for i in range(0, len(superevent_id)):
        pipeline[i] = pipeline[i][0]

    odds_ratio[odds_ratio<10**-30] = 10**-30

    if farcut == True:
        #cut on FAR
        farcut=(1.1574 * 10 ** (-5))
        idx_farcut = np.where([i < farcut for i in far])
        pipeline = pipeline[idx_farcut]
        superevent_id = superevent_id[idx_farcut]
        skyarea = skyarea[idx_farcut]
        pter = pter[idx_farcut]
        dist = dist[idx_farcut]
        odds_ratio = odds_ratio[idx_farcut]
        neutrinos = neutrinos[idx_farcut]
    
    idx_gstlal = np.where([i == 'gstlal' for i in pipeline])
    idx_pycbc = np.where([i == 'pycbc' for i in pipeline])
    idx_spiir = np.where([i == 'spiir' for i in pipeline])
    idx_MBTAOnline = np.where([i == 'MBTAOnline' for i in pipeline])
    
    if parameter == 'FAR':
        dat_gstlal = far[idx_gstlal]
        dat_pycbc = far[idx_pycbc]
        dat_spiir = far[idx_spiir]
        dat_MBTAOnline = far[idx_MBTAOnline]
    elif parameter == 'Odds Ratio':
        dat_gstlal = odds_ratio[idx_gstlal]
        dat_pycbc = odds_ratio[idx_pycbc]
        dat_spiir = odds_ratio[idx_spiir]
        dat_MBTAOnline = odds_ratio[idx_MBTAOnline]
    elif parameter == 'Neutrino Count':
        dat_gstlal = neutrinos[idx_gstlal]
        dat_pycbc = neutrinos[idx_pycbc]
        dat_spiir = neutrinos[idx_spiir]
        dat_MBTAOnline = neutrinos[idx_MBTAOnline]
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
        print("Invalid Parameter. Must be either SNR, Sky Area, p Terrestrial, or Distance")
        dat_gstlal = 0
        dat_pycbc = 0
        dat_spiir = 0
        dat_MBTAOnline = 0

    return dat_gstlal, dat_pycbc, dat_spiir, dat_MBTAOnline

def downselect_dataframe(df, parameter='FAR', farcut=False):
    #TODO unfinished function 
    superevent_id = df['Superevent_ID'].values
    odds_ratio = df['Odds_Ratio'].values
    neutrinos = df['Neutrino_Count'].values
    pipeline = df['Pipeline'].values
    f_far = df['FAR'].values
    f_skyarea = df['Sky_Area'].values
    f_pter = df['p_Terrestrial'].values
    f_dist = df['Mean_Distance'].values
    
    far = np.zeros(len(superevent_id))
    skyarea = np.zeros(len(superevent_id))
    pter = np.zeros(len(superevent_id))
    dist = np.zeros(len(superevent_id))

    for i in range(0, len(superevent_id)):
        far[i] = float(f_far[i][0])
        skyarea[i] = float(f_skyarea[i][0])
        pter[i] = float(f_pter[i][0])
        dist[i] = float(f_dist[i][0])
     
    odds_ratio[odds_ratio<10**-30] = 10**-30

    if farcut == True:
        #cut on FAR
        farcut=(1.1574 * 10 ** (-5))
        idx_farcut = np.where([i < farcut for i in far])
        pipeline = pipeline[idx_farcut]
        superevent_id = superevent_id[idx_farcut]
        skyarea = skyarea[idx_farcut]
        pter = pter[idx_farcut]
        dist = dist[idx_farcut]
        odds_ratio = odds_ratio[idx_farcut]
        neutrinos = neutrinos[idx_farcut]
    
    idx_gstlal = np.where([i == 'gstlal' for i in pipeline])
    idx_pycbc = np.where([i == 'pycbc' for i in pipeline])
    idx_spiir = np.where([i == 'spiir' for i in pipeline])
    idx_MBTAOnline = np.where([i == 'MBTAOnline' for i in pipeline])
    
    if parameter == 'FAR':
        dat_gstlal = far[idx_gstlal]
        dat_pycbc = far[idx_pycbc]
        dat_spiir = far[idx_spiir]
        dat_MBTAOnline = far[idx_MBTAOnline]
    elif parameter == 'Odds Ratio':
        dat_gstlal = odds_ratio[idx_gstlal]
        dat_pycbc = odds_ratio[idx_pycbc]
        dat_spiir = odds_ratio[idx_spiir]
        dat_MBTAOnline = odds_ratio[idx_MBTAOnline]
    elif parameter == 'Neutrino Count':
        dat_gstlal = neutrinos[idx_gstlal]
        dat_pycbc = neutrinos[idx_pycbc]
        dat_spiir = neutrinos[idx_spiir]
        dat_MBTAOnline = neutrinos[idx_MBTAOnline]
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
        print("Invalid Parameter. Must be either SNR, Sky Area, p Terrestrial, or Distance")
        dat_gstlal = 0
        dat_pycbc = 0
        dat_spiir = 0
        dat_MBTAOnline = 0

    return dat_gstlal, dat_pycbc, dat_spiir, dat_MBTAOnline
    
