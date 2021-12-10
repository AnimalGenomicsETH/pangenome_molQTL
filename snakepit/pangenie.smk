

def get_sample_location(sample):
    fastqs = []
    for R in (1,2):
        fastqs.append(str(Path(f'/cluster/work/pausch/inputs/fastq/BTA/{sample}_R{R}.fastq.gz').resolve()))
    return fastqs
    
rule fastp:
    input:
        lambda wildcards: get_sample_location(wildcards.sample)
    output:
        temp(get_dir('fastq','{sample}.fastq.gz'))
    threads: 4
    resources:
        mem_mb = 6000
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        fastp -w {threads} -i {input[0]} -I {input[1]} --stdout -g --thread {threads} --html /dev/null --json /dev/null --dont_eval_duplication | pigz -p {threads} - > {output}
        '''

rule pangenie:
    input:
        reference = config['reference'],
        vcf = config['panel'],
        fastq = get_dir('fastq','{sample}.fastq.gz') # '/cluster/scratch/alleonard/{sample}.fastq'#lambda wildcards: config['samples'][wildcards.sample]
    output:
        get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf')
    params:
        prefix = lambda wildcards, output: str(PurePath(output[0]).with_suffix('')).replace(f'_{wildcards.pangenie_mode}',''),
        phasing = lambda wildcards: '-p' if wildcards.pangenie_mode == 'phasing' else ''
    threads: 12
    resources:
        mem_mb = 15000,
        disk_scratch = 120,
        walltime = '24:00'
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        pigz -p {threads} -c -d {input.fastq} > $TMPDIR/{wildcards.sample}.fastq
        PanGenie -i $TMPDIR/{wildcards.sample}.fastq -r {input.reference} -v {input.vcf} -t {threads} -j {threads} -s {wildcards.sample} -g {params.phasing} -o {params.prefix}
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

rule merge_pangenie:
    input:
        (get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf.gz',sample=S) for S in config['samples'])
    output:
        get_dir('PG','samples.all.pangenie_{pangenie_mode}.vcf')
    threads: 6
    resources:
        mem_mb = 4000
    shell:
        '''
        bcftools merge --threads {threads} -o {output} {input}
        '''

rule compare_pangenie:
    input:
        pangenie = get_dir('PG','{sample}.all.pangenie_genotyping.vcf.gz')
    output:
        get_dir('concordance','{sample}.genotype_concordance_summary_metrics')
    params:
        gc_out = lambda wildcards,output: PurePath(output[0]).parent / f'{wildcards.sample}'
    resources:
        mem_mb = 15000,
        walltime = '60'
    envmodules:
        'gcc/8.2.0',
        'picard/2.25.7'
    shell:
        '''
        picard GenotypeConcordance CALL_VCF={input.pangenie} TRUTH_VCF={config[reference_vcf]} CALL_SAMPLE={wildcards.sample} TRUTH_SAMPLE={wildcards.sample} O={params.gc_out}
        '''
