

def get_sample_location(sample):
    fastqs = []
    for R in (1,2):
        fastqs.append(str(Path(f'/cluster/work/pausch/inputs/fastq/BTA/{sample}_R{R}.fastq.gz').resolve()))
    return fastqs
    
rule fastp:
    input:
        lambda wildcards: get_sample_location(wildcards.sample)
    output:
        temp(get_dir('fastq','{sample}.fastq'))
    threads: 4
    resources:
        mem_mb = 6000
    shell:
        '''
        fastp -w {threads} -i {input[0]} -I {input[1]} --stdout -g --thread {threads} --html /dev/null --json /dev/null --dont_eval_duplication > {output}
        '''

rule pangenie:
    input:
        reference = config['reference'],
        vcf = config['panel'],
        fastq = get_dir('fastq','{sample}.fastq') # '/cluster/scratch/alleonard/{sample}.fastq'#lambda wildcards: config['samples'][wildcards.sample]
    output:
        get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf')
    params:
        prefix = lambda wildcards, output: str(PurePath(output[0]).with_suffix('')).replace(f'_{wildcards.pangenie_mode}',''),
        phasing = lambda wildcards: '-p' if wildcards.pangenie_mode == 'phasing' else ''
    threads: 8
    resources:
        mem_mb = 15000
    shell:
        '''
        PanGenie -i {input.fastq} -r {input.reference} -v {input.vcf} -t {threads} -j {threads} -s {wildcards.sample} {params.phasing} -u -o {params.prefix}
        '''

rule bgzip_tabix:
    input:
        get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf')
    output:
        multiext(get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf.gz'),'','.tbi')
    threads: 2
    resources:
        mem_mb = 4000
    shell:
        '''
        bgzip --threads {threads} -c {input} > {output[0]}
        tabix -p vcf {output[0]}
        '''
