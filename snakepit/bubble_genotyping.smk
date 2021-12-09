import re
from pathlib import PurePath, Path
from collections import defaultdict
from itertools import product

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base,ext='',**kwargs):
    if base == 'SV':
        base_dir = list(config['assemblies'].keys())[0] + '_run_{run}'
    elif base == 'VG':
        base_dir = list(config['assemblies'].keys())[0] + '_run_{run}_VG'
    elif base == 'PG':
        base_dir = 'pangenie'
    elif base == 'fastq':
        base_dir = '/cluster/scratch/alleonard'
    elif base == 'main':
        base_dir = ''
    elif base == 'fasta':
        base_dir = 'fasta'
    elif base == 'pbsv':
        base_dir = 'pbsv'
    else:
        raise Exception('Base not found')
    if ext and ext[0] == '.':
        return f'{base_dir}{ext}'.format_map(Default(kwargs))
    return str(PurePath(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

wildcard_constraints:
    L = r'\d+'

include: 'pangenie.smk'
include: 'pbsv.smk'

def capture_logic():
    targets = []


    for caller in ('mg',):#'vg'):
        targets.append(get_dir('VG',f'annotated.L50.{caller}.df',run='TEST'))
    targets.append(get_dir('PG','samples.all.pangenie_genotyping.vcf'))

    #for sample in config['samples']:
    #    targets.append(get_dir('VG',f'{sample}.all.L50.mg.gaf',run='TEST'))
    #    targets.append(get_dir('PG',f'{sample}.all.L50.vg.gaf',run='TEST'))
    #    targets.append(get_dir('PG',f'{sample}.all.pangenie_phasing.vcf.gz'))
    
    return targets

rule all:
    input:
        capture_logic()
        #get_dir('VG','test.all.L50.vg.gaf',run='TEST'),
        #get_dir('VG','test.all.L50.mg.gaf',run='TEST')

rule merge_minigraph:
    input:
        gfas = (get_dir('SV','{chr}.L{L}.gfa',chr=CHR) for CHR in range(1,30))
    output:
        get_dir('VG','all.L{L}.gfa')
    threads: 1
    resources:
        mem_mb = 5000
    run:
        for i,gfa in enumerate(input.gfas,1):
            shell(f"sed -e 's/s\([0-9]\+\)/{i}_\\1/g' {gfa} >> {output}")

rule merge_intersections:
    input:
        gfas = (get_dir('SV','{chr}.L{L}.join.df',chr=CHR) for CHR in range(1,30))
    output:
        get_dir('SV','all.L{L}.join.df')
    run:
        with open(output[0],'w') as fout:
            for i,fname in enumerate(input.gfas):
                with open(fname,'r') as fin:
                    first_line = fin.readline()
                    if i == 0:
                        fout.write(first_line)
                    for line in fin:
                        parts = line.split(',')
                        _id = parts[-2]
                        new_id = _id.replace('_INV','').split('_')
                        r_id = new_id[-1]+'_'+new_id[0][1:]
                        fout.write(line.replace(_id,r_id))

rule count_gaf_node_support:
    input:
        get_dir('VG','{sample}.all.L{L}.{caller}.gaf')
    output:
        get_dir('VG','{sample}.all.L{L}.{caller}.node_counts')
    shell:
        '''
        awk '{{split($6,b,/>|</); for (key in b) {{ if(b[key]~/[[:digit:]]/) print  b[key] }} }}' {input} | sort -V | uniq -c | sort -k1,1nr > {output}
        #awk '$6~/[[:digit:]]/ {{split($6,b,/>|</); for (key in b) {{if(b[key]~/[[:digit:]]/) print  b[key] }} }}' {input} | sort -V | uniq -c | sort -k1,1nr > {output}
        '''

rule annotate_variants:
    input:
        df = get_dir('SV','all.L{L}.join.df'),
        counts = (get_dir('VG','{sample}.all.L{L}.{caller}.node_counts',sample=S) for S in config['samples'])
    output:
        get_dir('VG','annotated.L{L}.{caller}.df')
    run:
        sample_counter = defaultdict(lambda : defaultdict(int))
        for sample,sample_f in zip(config['samples'],input.counts):
            with open(sample_f,'r') as fin:
                for line in fin:
                    print(line)
                    parts = line.rstrip().split()
                    sample_counter[parts[1]][sample] = int(parts[0])
        print('made it thus far')
        with open(input.df,'r') as fin_df, open(output[0],'w') as fout:
            for i,line in enumerate(fin_df):
                if i == 0:
                    fout.write(line.rstrip() +',' + ','.join(config['samples'])+'\n')
                else:
                    parts = line.split(',')
                    counts = ','.join(map(str,(sample_counter[parts[-2]][sample] for sample in config['samples'])))
                    fout.write(line.rstrip() + ','+counts+ '\n')



rule gfatools_noseq:
    input:
        get_dir('VG','all.L{L}.gfa')
    output:
        get_dir('VG','all.L{L}.noseq.gfa')
    resources:
        mem_mb = 5000,
        walltime = '30'
    shell:
        'gfatools view -S {input} > {output}'
 
### minigraph
rule minigraph_sr:
    input:
        gfa = get_dir('VG','all.L{L}.gfa'),
        fastq = get_dir('fastq','{sample}.fastq.gz') #        fastq = lambda wildcards: config['samples'][wildcards.sample]
    output:
        gaf = get_dir('VG','{sample}.all.L{L}.mg.gaf')
    threads: 16
    resources:
        mem_mb = 2000,
        walltime = '4:00'
    shell:
        '''
        minigraph -t {threads} -xsr --vc {input.gfa} {input.fastq} > {output.gaf} # | {workflow.basedir}/remap.py --minigraph {input.gfa} > {output.gaf}
        '''

### VG
rule vg_construct:
    input:
        get_dir('VG','all.L{L}.gfa')
    output:
        gbz = multiext(get_dir('VG','all.L{L}'),'.giraffe.gbz','.min','.dist','.xg','.chopped.P_lines')
    params:
        lambda wildcards, output: PurePath(output[1]).with_suffix('')
    threads: 4
    resources:
        mem_mb = 15000,
        walltime = '24:00'
    shell:
        '''
        singularity exec -B $(pwd):$(pwd) -B $TMPDIR:$TMPDIR /cluster/work/pausch/alex/images/vg_v1.36.0.sif \
        /bin/bash -c "vg autoindex --workflow giraffe --request XG \
        -g {input} -t {threads} -p {params} -T $TMPDIR; \
        vg view -g {output[3]} | awk '$1=="P"&&$2!~/ARS/' > {output[4]}"
        '''

rule vg_snarl:
    input:
        ''
    output:
        ''
    shell:
        '''
        vg snarls -T chr${i}.vg >> all.snarls
        '''

rule vg_giraffe:
    input:
        gbz = multiext(get_dir('VG','all.L{L}'),'.giraffe.gbz','.min','.dist','.chopped.P_lines'),
        gfa = get_dir('VG','all.L{L}.noseq.gfa'),
        fastq = get_dir('fastq','{sample}.fastq.gz') #lambda wildcards: config['samples'][wildcards.sample]
    output:
        gaf = get_dir('VG','{sample}.all.L{L}.vg.gaf')
    threads: 18
    resources:
        mem_mb = 3500
    params:
        lambda wildcards, input: Path(input.fastq).resolve().parent
    shell:
        '''
        singularity exec -B $(pwd):$(pwd) -B {params}:{params} /cluster/work/pausch/alex/images/vg_v1.36.0.sif \
        /bin/bash -c "vg giraffe -t {threads} -Z {input.gbz[0]} \
        -m {input.gbz[1]} -d {input.gbz[2]} -o gaf \
        -i -f {input.fastq}" > {output} #| {workflow.basedir}/remap.py --vg {input.gfa} {input.gbz[3]} | cut -f-1,5-12,15- > {output.gaf}
        #awk '$6 ~ /^[0-9]+$/'
        '''
