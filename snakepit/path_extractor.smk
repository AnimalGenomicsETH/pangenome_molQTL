def read_paths(fname,paths):
    with open(fname,'r') as fin:
        for line in fin:
            parts = line.split()
            chrom = parts[0].split('_')[0]
            if parts[5][0] == '.':
                continue
            
            paths[f'{chrom}_{parts[3][2:]}->{chrom}_{parts[4][2:]}'].add(parts[5].split(':')[0].replace('s',f'{chrom}_'))

from collections import defaultdict
def extract_haplotype_walks(haplotypes):
    paths = defaultdict(set)

    for i in range(1,30):
        for haplotype in haplotypes:
            read_paths(f'ARS_run_5/{haplotype}.path.{i}.L50.bed',paths)

    with open('test.paths','w') as fout:
        for key, v in paths.items():
            fout.write(key +',' +','.join(v) + '\n')


extract_haplotype_walks(['Os1_hifiasm','Os2_hifiasm','Od1_hifiasm','Od2_hifiasm','B31_hifiasm','B32_hifiasm','B41_hifiasm','B42_hifiasm','O_hifiasm','H_hifiasm','P_hifiasm','H_clr','A_clr','S_ont'])
