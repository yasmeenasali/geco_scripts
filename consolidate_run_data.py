#!/home/yasmeen.asali/miniconda3/envs/ligo-py36/bin/python3.6

DESC='''A module for searching through O3 subthreshold data. Available functions:
    get_90_skyarea: Calculate the 90 area for a given superevent path 
    test_skyarea: Test the above function (no arguments)
    get_distmean: Get the mean posterior distance for a given superevent path
    generate_data_file: Generate a numpy array file for all 03A/B CBC events with the following cols:
         Superevent ID, SNR, FAR, pipeline, p Terrestrial, Mean Posterior Distance, 
         and optionally the 90% Sky Area

This file can be run as a script. 

'''

import argparse 
parser = argparse.ArgumentParser(description=DESC, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-r', '--run', type=str, help=("Specify which run to generate file for (ex. O3A, O3B)"))
parser.add_argument('--cbc', default=False, help=("Generate CBC File"))
parser.add_argument('--cwb', default=False, help=("Generate CWB File"))
args = parser.parse_args()

import numpy as np
import json
import csv
import ligo.skymap.moc
from astropy.table import Table
from astropy.io import fits
import pandas as pd

def get_90_skyarea(PATH):
    data = Table.read(PATH)
    if data['UNIQ'].dtype == 'uint64': #dtype must be 'uint64' for ligo.skymap version < 0.1.9, 'int64' otherwise
        data['UNIQ'].dtype = 'int64'
    ligo_data_nested = Table(ligo.skymap.moc.rasterize(data))
    #(order, ipix) = ligo.skymap.moc.uniq2nest(ipix_uniq)
    #Downselect 90% region 
    num_pix = len(ligo_data_nested) #* u.pixel
    sphere = 4*np.pi #* u.sr
    conversion = sphere / num_pix
    prob_per_pix = ligo_data_nested["PROBDENSITY"] * conversion
    ligo_data_nested['PROB'] = prob_per_pix
    ipix = np.array(range(1, num_pix+1))
    ligo_data_nested['PIX'] = ipix
    sorted_data = np.sort(ligo_data_nested, order="PROB")
    sorted_data_reverse = Table(sorted_data)[::-1]
    prob = 0 
    idx = []
    while prob < .9:
        for i in range(0, len(sorted_data_reverse)):
            row = sorted_data_reverse[i]
            if (prob > .9):
                break
            prob += row["PROB"]
            idx.append(i)
    ligo_data_90 = sorted_data_reverse[idx]
    sq_deg = 41253
    skyarea = (len(ligo_data_90)/num_pix) * sq_deg
    return float(skyarea) 

def test_skyarea():
    #test get_90_skyarea function
    superevent='S190814bv'
    skymap = Table.read("superevents/{}/bayestar.multiorder.fits".format(superevent))
    skyarea = get_90_skyarea(skymap) 
    print(skyarea)

def get_distmean(PATH):
    data = fits.open(PATH)
    dist = data[1].header['DISTMEAN']
    return float(dist)

def generate_data_file(PATH, RUN, superevents_cbc, get_skyarea=False):
    superevent_data = []
    snr_data = []
    far_data = []
    pipeline_data = []
    pter_data = []
    dist_data = []
    skyarea_data = []
    for superevent in superevents_cbc:
        with open("{}/superevents/{}/preferred_event_data.json".format(PATH, superevent), "r") as read_file:
            data = json.load(read_file)
            snr = data['extra_attributes']['CoincInspiral']['snr']
            pipeline = data['pipeline']
        with open("{}/superevents/{}/superevent_data.json".format(PATH, superevent), "r") as read_file_S:
            dataS = json.load(read_file_S)
            far = dataS['far']
        with open("{}/superevents/{}/p_astro.json".format(PATH, superevent), "r") as read_file:
            dataP = json.load(read_file)
            pter = dataP['Terrestrial']
        skymap = "{}/superevents/{}/bayestar.multiorder.fits".format(PATH, superevent)
        dist = get_distmean(skymap) 
        superevent_data.append(superevent)
        snr_data.append(float(snr))
        far_data.append(float(far))
        pipeline_data.append(pipeline)
        pter_data.append(float(pter))
        dist_data.append(dist)
        if get_skyarea == True:
            skyarea = get_90_skyarea(skymap)
            skyarea_data.append(skyarea)
        print('Finished {}'.format(superevent))
    if get_skyarea == True:
        kwargs = {'Superevent': superevent_data, 'SNR': snr_data, 'FAR': far_data, 
                  'Pipeline': pipeline_data, 'p_Terrestrial': pter_data, 'Distance': dist_data, 'Skyarea': skyarea_data}
    else:
        kwargs = {'Superevent': superevent_data, 'SNR': snr_data, 'FAR': far_data, 
                  'Pipeline': pipeline_data, 'p_Terrestrial': pter_data, 'Distance': dist_data}
    np.savez(f'{RUN}_cbc_full_data', **kwargs)

def generate_cwb_data_file(PATH, RUN, superevents_cwb):
    superevent_data = []
    snr_data = []
    far_data = []
    amplitude_data = []
    bandwidth_data = []
    central_freq_data = []
    duration_data = []
    for superevent in superevents_cwb:
        try:
            with open("{}/superevents/{}/preferred_event_data.json".format(PATH, superevent), "r") as read_file:
                data = json.load(read_file)
                amplitude = data['extra_attributes']['MultiBurst']['amplitude']
                bandwidth = data['extra_attributes']['MultiBurst']['bandwidth']
                central_freq = data['extra_attributes']['MultiBurst']['central_freq']
                duration = data['extra_attributes']['MultiBurst']['duration']
                snr = data['extra_attributes']['MultiBurst']['snr']
                far = data['far']
                superevent_data.append(superevent)
                snr_data.append(snr)        
                far_data.append(far)        
                amplitude_data.append(amplitude)        
                bandwidth_data.append(bandwidth)        
                central_freq_data.append(central_freq)        
                duration_data.append(duration)        
        except FileNotFoundError as err:
            continue
    kwargs = {'Superevent': superevent_data, 'SNR': snr_data, 'FAR': far_data, 
              'Amplitude': amplitude_data, 'Bandwidth': bandwidth_data, 
              'Central_Freq': central_freq_data, 'Duration': duration_data}
    np.savez('{}_cwb_full_data'.format(RUN), **kwargs)

if __name__ == "__main__": 
    RUN = args.run
    if args.cbc:
        print(f'Generating {RUN} Subthreshold CBC Events File')
        PATH = f'/home/yasmeen.asali/GWHEN/{RUN}_subthreshold'
        superevents_cbc = np.load(f'{PATH}/data/{RUN}_cbc_superevent_ids.npy')
        generate_data_file(PATH, RUN, superevents_cbc, get_skyarea=True)
        print("Finished generating file") 
    if args.cwb:
        print(f'Generating {RUN} Subthreshold CWB Events File')
        PATH = f'/home/yasmeen.asali/GWHEN/{RUN}_subthreshold'
        if RUN == 'O3A':
            superevents_cwb = np.load(f'{PATH}/data/superevent_ID_lists/missing_event_ids_cwb.npy')
        if RUN == 'O3B':
            superevents_cwb = np.load(f'{PATH}/data/O3B_cwb_superevent_ids.npy')
        generate_cwb_data_file(PATH, RUN, superevents_cwb)
        print("Finished generating file") 
