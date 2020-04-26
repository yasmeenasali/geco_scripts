#!/home/yasmeen.asali/miniconda3/envs/ligo-py36/bin/python3.6

DESC = '''

A script for downloading data from GraceDb for O3A/O3B

'''

from os import mkdir, path 
from errno import EEXIST
from ligo.gracedb.rest import GraceDb, HTTPError
import json
import numpy as np

def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''
    try:
        mkdir(mypath)
    except OSError as exc:
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

client = GraceDb()

def download_run_superevent_ids(run = 'O3A'):
    #Download full list of O3A superevents
    if run == 'O3A':
        query = 'gpstime: 1238112018 .. 1253977218'
    #Download full list of O3B superevents
    elif run == 'O3B':
        query = 'gpstime: 1256601600 .. 1269388860'
    else:
        print('Invalid run')

    superevent_iterator = client.superevents(query) 
    superevent_ids = [superevent['superevent_id'] for superevent in superevent_iterator]
    
    #in O3A there were 8049 events, in O3B there were 7198 
    print(len(superevent_ids)) #check that all were downloaded 
    np.save('superevent_ids', superevent_ids)

def get_event_data(superevent, pref_event, data):
    with open('superevents/{}/superevent_data.json'.format(superevent), 'w') as outfile_json:
        json.dump(data, outfile_json)
    response_pe = client.event(pref_event)
    with open('superevents/{}/preferred_event_data.json'.format(superevent), 'w') as outfile_pejson:
        json.dump(response_pe.json(), outfile_pejson)

def save_skymap_pastro_files(event_id, superevent):
   try: 
       response_sky = client.files(event_id, 'bayestar.multiorder.fits') 
   except HTTPError as err_fname: 
       response_sky = client.files(event_id, 'bayestar.fits') 
   with open('superevents/{}/bayestar.multiorder.fits'.format(superevent), 'wb') as outfile_sky: 
       outfile_sky.write(response_sky.read())
   response_pastro = client.files(event_id, 'p_astro.json')
   with open('superevents/{}/p_astro.json'.format(superevent), 'wb') as outfile_pastro: 
       outfile_pastro.write(response_pastro.read())

def save_embright_file(superevent):
    response_json = client.superevent(superevent)
    data = response_json.json()
    pref_event = data['preferred_event']
    try: 
        response_em = client.files(pref_event, 'subthreshold.em_bright.json') 
    except HTTPError as err_fname: 
        response_em = client.files(pref_event, 'em_bright.json') 
    with open('superevents/{}/em_bright.json'.format(superevent), 'wb') as outfile_em: 
        outfile_em.write(response_em.read())

def get_CBC_event_files(superevent):
    response_json = client.superevent(superevent)
    data = response_json.json()
    pref_event = data['preferred_event']
    try:
        save_skymap_pastro_files(superevent, superevent)
    except HTTPError as err_pref:
        response_json = client.superevent(superevent)
        data = response_json.json()
        pref_event = data['preferred_event']
        save_skymap_pastro_files(pref_event, superevent)
    get_event_data(superevent, pref_event, data)

def get_CWB_event_files(superevent):
    response_json = client.superevent(superevent)
    data = response_json.json()
    pref_event = data['preferred_event']
    try:
        response_sky = client.files(superevent, 'cWB.fits.gz') 
    except HTTPError as err_pref:
        response_sky = client.files(pref_event, 'cWB.fits.gz') 
    with open('superevents/{}/cWB.fits.gz'.format(superevent), 'wb') as outfile_sky: 
       outfile_sky.write(response_sky.read())
    trigger_file = 'trigger_{}.txt'.format(data['t_0'])
    response_trigger = client.files(pref_event, trigger_file)
    with open('superevents/{0}/{1}'.format(superevent, trigger_file), 'wb') as outfile_trigger:
        outfile_trigger.write(response_trigger.read())     
    response_eventlog = client.files(pref_event, 'event.log') 
    with open('superevents/{}/event.log'.format(superevent), 'wb') as outfile_el:
        outfile_el.write(response_eventlog.read())     
    get_event_data(superevent, pref_event, data)

#4/21/20 downloading O3B skymaps
#downloading cbc and cwb 

all_superevent_ids = np.load('data/O3B_superevent_ids.npy') #all O3B superevents

#cwb_superevent_ids = np.load('data/missing_event_ids_cwb.npy') #all O3A cWB events
#cbc_superevent_ids = np.load('data/superevent_ids_cbc_events.npy') #all O3A CBC events

cbc_event_ids = []
cwb_event_ids = []
missing_events = []
for superevent in all_superevent_ids:
    mkdir_p('superevents/{}'.format(superevent)) 
    try:
        try:
            get_CBC_event_files(superevent)
            #save_embright_file(superevent)
            print(f"Finished CBC Event {superevent}")
            cbc_event_ids.append(superevent)
        except HTTPError as err_no_bayestar_skymap:
            get_CWB_event_files(superevent)
            print(f"Finished CWB Event {superevent}")
            cwb_event_ids.append(superevent)
    except HTTPError as err_unknown:
        print(f'Skipping event {superevent}')
        missing_events.append(superevent)

print('Total CBC events: {}'.format(len(cbc_event_ids)))
np.save('O3B_cbc_superevent_ids', cbc_event_ids) 

print('Total CWB events: {}'.format(len(cwb_event_ids)))
np.save('O3B_cwb_superevent_ids', cwb_event_ids) 

print('Skipped {} Events'.format(len(missing_events)))
np.save('O3B_missing_superevent_ids', missing_events)
