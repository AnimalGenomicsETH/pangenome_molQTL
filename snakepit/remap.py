#!/usr/bin/env python3

import sys
import re
#import bisect


def get_reference_node(keys,coordinate,nodes):
    idx = bisect.bisect_left(keys,coordinat)
    if idx == len(keys) or keys[idx] != coordinate:
        idx -= 1
    return nodes[keys[idx]]

#def build_node_dict():


#takes in GFA
def load_map(fname,mode='minigraph'):
    with open(fname,'r') as fin:
        IDs = dict()
        #ref_IDs = defaultdict(dict)
        for line in fin:
            if line[0] == 'S':
                parts = line.split()
                #if 'ARS' in line:
                #    continue
                if mode == 'minigraph':
                    start = parts[5].split(':')[-1]
                    #if 'ARS' in parts[4]:
                    #    chrom = parts[4].split(':')[-1].split('_')[0]
                    #    ref_IDs[chrom][start] = parts[1]
                    #else:
                    IDs[parts[4].split(':')[-1] + ':'+start+'-'+str(int(start)+int(parts[3].split(':')[-1]))] = parts[1]
                elif mode == 'vg': 
                    IDs[parts[4].split(':')[-1]+'['+  parts[5].split(':')[-1]+']'] = parts[1]
        return IDs

def ID_labels(fname,maps):
    with open(fname,'r') as fin:
        labels = dict()    
        for line in fin:
            #if 'ARS' in line:
            #    continue
            parts = line.split()
            ids = parts[2].replace('+','').split(',')
            for _id in ids:
                labels[_id] = maps[parts[1]]
        return labels
 
def relabel_ids(line,mappings):
    parts = line.split()
    st = ''
    for strand, _id in zip(*(iter(re.split('(>|<)',parts[5])[1:]),) * 2):
        st+=f'{strand}{mappings.get(_id,"XXREF")}'
    #if not st:
        #do bisecting
    sys.stdout.write(line.replace(parts[5],st or 'REF'))

if sys.argv[1] == '--minigraph':
    ID_map = load_map(sys.argv[2],'minigraph')
elif sys.argv[1] == '--vg':
    ID_map = ID_labels(sys.argv[3],load_map(sys.argv[2],'vg'))
for line in sys.stdin:
    relabel_ids(line,ID_map)
