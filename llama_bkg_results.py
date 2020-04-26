#!/home/yasmeen.asali/miniconda3/envs/ligo-py36/bin/python3.6
import numpy as np 

'''
A module for parsing results from llama background run results from the significance text file.
'''

def get_superevent_ID(line):
    path = line.split()[2]
    ID = path.split('/')[2].split('-')[0]
    return ID

def read_lines(lines):
    odds_list = []
    neutrino_list = []
    superevent_list = []
    for line in lines:
        element_list = line.split()
        odds = float(element_list[0])
        neutrino_num = float(element_list[1])
        ID = get_superevent_ID(line)
        odds_list.append(odds)
        neutrino_list.append(neutrino_num)
        superevent_list.append(ID)
    return odds_list, neutrino_list, superevent_list

if __name__ == "__main__":
    filename = 'subthresh_significance_outputs.txt'
    with open(filename, 'r', newline='') as f:
        lines = [line.strip() for line in f]
        lines = lines[2:] #remove header lines
        odds_list, neutrino_list, superevent_list = read_lines(lines)
        kwargs = {'Odds_Ratio': odds_list, 'Neutrino_Count': neutrino_list, 'Superevent_ID': superevent_list}
        np.savez('subthresh_significance_outputs.npz', **kwargs)

