#!/usr/bin/env python3

import sys
import re

#takes in GFA
def load_map(fname,mode='minigraph'):
    with open(fname,'r') as fin:
        IDs = dict()
        for line in fin:
            if line[0] == 'S':
                parts = line.split()
                if mode == 'minigraph':
                    IDs[parts[4].split(':')[-1] + ':'+parts[5].split(':')[-1]+'-'+parts[3].split(':')[-1]] = parts[1]
                elif mode == 'vg': 
                    IDs[parts[4].split(':')[-1]+'['+  parts[5].split(':')[-1]+']'] = parts[1]
        return IDs

def ID_labels(fname,maps):
    with open(sys.argv[1],'r') as fin:
        labels = dict()    
        for line in fin:
            parts = line.split()
            ids = parts[2].replace('+','').split(',')
            for _id in ids:
                labels[_id] = maps[parts[1]]
        return labels
 
def relabel_ids(line,mappings):
    parts = line.split()
    st = ''
    for strand, _id in zip(*(iter(re.split('(>|<)',parts[5])[1:]),) * 2):
        st+=f'{strand}{mappings.get(_id,"REF")}'
    sys.stdout.write(line.replace(parts[5],st))

if sys.argv[1] == '--minigraph':
    ID_map = load_map(sys.argv[2],'minigraph')
elif sys.argv[1] == '--vg':
    ID_map = ID_labels(sys.argv[3],load_map(sys.argv[2],'vg'))

for line in sys.stdin:
    relabel_ids(line,ID_map)
