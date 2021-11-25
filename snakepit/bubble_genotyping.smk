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

rule all:
    input:
        get_dir('VG','test.all.L50.vg.gaf')


rule merge_minigraph:
    input:
        gfas = (get_dir('SV','{chr}.L{L}.gfa',chr=CHR) for CHR in range(1,30))
    output:
        get_dir('VG','all.L{L}.gfa')
    run:
        for i,gfa in enumerate(input.gfas,1):
            shell("sed -e 's/s\([0-9]\+\)/{i}_\1/g' {gfa} >> {output")

#rule minigraph_sr:
#    input:
#        gfa = '',
#        sr = ''
#    output:
#        'gaf'
#    threads: 24
#    resources:
#        mem_mb = 3500
#    shell:
#        '''
#        minigraph -t {threads} -xsr {input.gfa} {input.sr} | awk '$6~/>/||$6~/</ {print $6}' > {output}
#        '''
#
#rule gfatools_noseq:
#    input:
#        ''graph
#    output:
#        ''noseq
#    shell:
#        'gfatools view -S {input} > {output}'
#
#rule genotype_gaf:
#    input:
#        gfa = ''
#        gaf = ''
#    output:
#        ''
#    run:
#        with open(input,'r') as fin:
#            for line in fin:
#                splits =  list(zip(*(iter(re.split('>|<',line)[1:]),) * 2))
#                for hits in splits:
#                    orientation = hits[0]
#                    location = hits.split(':')
#                    code = location[0][1:]
#                    start = location[1].split('-')[0]
#                    shell(f"awk '/{code}/&&/{start}/ {print $2}' >> {{output}}")



### VG
rule vg_construct:
    input:
        get_dir('VG','all.L{L}.gfa')
    output:
        gbz = multiext(get_dir('VG','all.L{L}'),'.giraffe.gbz','.min','.dist','.xg','.gfa')
    params:
        lambda wildcards, output: PurePath(output[1].with_suffix('')
    shell:
        '''
        vg autoindex --workflow giraffe --request XG -g {input.gfa} -p {params} -T $TMPDIR
        vg view -g {output[3]} > {output[4]}
        '''

rule vg_giraffe:
    input:
        gbz = multiext(get_dir('VG','all.L{L}'),'.giraffe.gbz','.min','.dist'),
        fastq = lambda wildcards: config['samples'][wildcards.sample]
    output:
        gaf = get_dir('VG','{sample}.all.L{L}.chopped.{caller}.gaf')
    shell:
        '''
        /bin/bash -c "vg giraffe -t {threads} -Z {input.gbz[0]} -m {input.gbz[1]} -d {input.gbz[2]} -i -f {input.fastq} -o gaf | cut  -f-1,5-12,15- > {output.gaf}"
        '''


rule reassociate_IDs:
    input:
        get_dir('VG','{sample}.all.L{L}.chopped.vg.gaf')
    output:
        get_dir('VG','{sample}.all.L{L}.vg.gaf')
    shell:
        ''
