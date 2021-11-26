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
    elif base == 'mash':
        base_dir = 'mash'
    elif base == 'fasta':
        base_dir = 'fasta'
    else:
        raise Exception('Base not found')
    if ext and ext[0] == '.':
        return f'{base_dir}{ext}'.format_map(Default(kwargs))
    return str(PurePath(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

wildcard_constraints:
    L = r'\d+'

#include: 'pangenie.smk'

rule all:
    input:
        get_dir('VG','test.all.L50.vg.gaf',run='TEST'),
        get_dir('VG','test.all.L50.mg.gaf',run='TEST')

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
        fastq = lambda wildcards: config['samples'][wildcards.sample]
    output:
        gaf = get_dir('VG','{sample}.all.L{L}.mg.gaf')
    threads: 18
    resources:
        mem_mb = 3500
    shell:
        '''
        minigraph -t {threads} -xsr {input.gfa} {input.fastq} | {workflow.basedir}/remap.py --minigraph {input.gfa} > {output.gaf}
        #awk '$6~/>/||$6~/</ {{print $6}}'
        '''

### VG
rule vg_construct:
    input:
        get_dir('VG','all.L{L}.gfa')
    output:
        gbz = multiext(get_dir('VG','all.L{L}'),'.giraffe.gbz','.min','.dist','.chopped.xg','.chopped.P_lines')
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
        vg view -g {output[3]} | awk '/P/' > {output[4]}"
        '''

rule vg_giraffe:
    input:
        gbz = multiext(get_dir('VG','all.L{L}'),'.giraffe.gbz','.min','.dist','.chopped.P_lines'),
        gfa = get_dir('VG','all.L{L}.noseq.gfa'),
        fastq = lambda wildcards: config['samples'][wildcards.sample]
    output:
        gaf = get_dir('VG','{sample}.all.L{L}.vg.gaf')
    threads: 18
    resources:
        mem_mb = 3500
    shell:
        '''
        singularity exec -B $(pwd):$(pwd) /cluster/work/pausch/alex/images/vg_v1.36.0.sif \
        /bin/bash -c "vg giraffe -t {threads} -Z {input.gbz[0]} \
        -m {input.gbz[1]} -d {input.gbz[2]} -o gaf \
        -i -f {input.fastq} | {workflow.basedir}/remap.py --vg {input.gfa} {input.gbz[3]} | cut -f-1,5-12,15- > {output.gaf}"
        #awk '$6 ~ /^[0-9]+$/'
        '''
