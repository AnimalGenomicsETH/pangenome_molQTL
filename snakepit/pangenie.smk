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

rule jellyfish_count:
    input:
        lambda wildcards: get_sample_location(wildcards.sample),
        paths = 'pangenie.2772146623768787664.path_segments.fasta'
    output:
        temp(get_dir('fastq','{sample}.jf'))
    threads: 6
    resources:
        mem_mb = 8000
    shell:
        '''
        fastp -w {threads} -i {input[0]} -I {input[1]} --stdout -g --thread {threads} --html /dev/null --json /dev/null --dont_eval_duplication | jellyfish count -L 1 -U 10000 -m 31 -s 3000000000 -p 126 -c 7 -C -t {threads} --if {input.paths} -o {output} /dev/fd/0
        '''

rule pangenie:
    input:
        reference = config['reference'],
        vcf = config['panel'],
        fastq = get_dir('fastq','{sample}.fastq.gz') # '/cluster/scratch/alleonard/{sample}.fastq'#lambda wildcards: config['samples'][wildcards.sample]
    output:
        get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf.OLD')
    params:
        prefix = lambda wildcards, output: str(PurePath(output[0]).with_suffix('')).replace(f'_{wildcards.pangenie_mode}',''),
        phasing = lambda wildcards: '-p' if wildcards.pangenie_mode == 'phasing' else ''
    threads: 12
    resources:
        mem_mb = 12500,
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

rule pangenie_cereal:
    input:
        reference = config['reference'],
        vcf = config['panel']
    output:
        'pangenie.2772146623768787664.path_segments.fasta' #how to get the hash name?
    threads: 6
    resources:
        mem_mb = 10000
    shell:
        '''
        /cluster/work/pausch/alex/software/Alex-pangenie/build/src/PanGenie -T -i /dev/null -r {input.reference} -v {input.vcf} -t {threads}
        touch {output}
        '''

rule direct_pangenie:
    input:
        reference = config['reference'],
        vcf = config['panel'],
        fastq = lambda wildcards: get_sample_location(wildcards.sample),
        paths = 'pangenie.2772146623768787664.path_segments.fasta'
    output:
        get_dir('PG','{sample}.all.pangenie_genotyping.vcf')
        #temp(multiext(get_dir('PG','{sample}.all.pangenie_'),'genotyping.vcf','histogram.histo'))
    params:
        prefix = lambda wildcards, output: str(PurePath(output[0]).with_suffix('')).replace(f'_genotyping','')
        #phasing = lambda wildcards: '-p' if wildcards.pangenie_mode == 'phasing' else ''
    threads: 6
    resources:
        mem_mb = 16000,
        disk_scratch = 75,
        walltime = '24:00'
    shell:
        '''
        fastp -w {threads} -i {input.fastq[0]} -I {input.fastq[1]} --stdout -g --thread {threads} --html /dev/null --json /dev/null --dont_eval_duplication | jellyfish count -L 1 -U 10000 -m 31 -s 3000000000 -p 126 -c 7 -C -t {threads} --if {input.paths} -o $TMPDIR/{wildcards.sample}.jf /dev/fd/0
        /cluster/work/pausch/alex/software/Alex-pangenie/PanGenie-static -i $TMPDIR/{wildcards.sample}.jf -r {input.reference} -v {input.vcf} -t {threads} -j {threads} -s {wildcards.sample} -g -o {params.prefix}
        '''

rule bgzip_tabix:
    input:
        get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf')
    output:
        temp(multiext(get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf.gz'),'','.tbi'))
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
       vcf = (get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf.gz',sample=S) for S in config['samples']),
       tbi = (get_dir('PG','{sample}.all.pangenie_{pangenie_mode}.vcf.gz.tbi',sample=S) for S in config['samples'])
    output:
        get_dir('PG','samples.all.pangenie_{pangenie_mode}.vcf.gz')
    threads: 6
    resources:
        mem_mb = 4000
    shell:
        '''
        bcftools merge --threads {threads} -o {output} {input.vcf}
        tabix -p vcf {output}
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
