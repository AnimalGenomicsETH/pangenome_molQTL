def get_sample_location(sample):
    fastqs = []
    for R in (1,2):
        fastqs.append(str(Path(f'{config["fastq"]}{sample}_R{R}.fastq.gz').resolve()))
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

rule pangenie_index:
    input:
        reference = config['reference'],
        vcf = config['panel']
    output:
        multiext('pangenie','.cereal','.path_segments.fasta')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 1
    resources:
        mem_mb = 15000
    shell:
        '''
        pangenie -B {params.prefix} -i /dev/null -r {input.reference} -v {input.vcf} -t {threads}
        '''

rule jellyfish_count:
    input:
        lambda wildcards: get_sample_location(wildcards.sample),
        pangenie_index = rules.pangenie_index.output[1]
    output:
        temp(get_dir('fastq','{sample}.jf'))
    threads: 6
    resources:
        mem_mb = 8000
    shell:
        '''
        fastp -w {threads} -i {input[0]} -I {input[1]} --stdout -g --thread {threads} --html /dev/null --json /dev/null --dont_eval_duplication | jellyfish count -L 1 -U 10000 -m 31 -s 3000000000 -p 126 -c 7 -C -t {threads} --if {input.pangenie_index} -o {output} /dev/fd/0
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

rule pangenie_genotype:
    input:
        reference = config['reference'],
        jellyfish = rules.jellyfish_count.output[0],
        pangenie_index = rules.pangenie_index.output
    output:
        get_dir('PG','{sample}.all.pangenie_genotyping.vcf')
    params:
        prefix = lambda wildcards, output: str(PurePath(output[0]).with_suffix('')).replace(f'_genotyping',''),
        index = lambda wildcards, input: PurePath(input.pangenie_index[0]).with_suffix('')
    threads: 8
    resources:
        mem_mb = 13000,
        scratch = '75G',
        walltime = '4:00'
    shell:
        '''
        pangenie -i {input.jellyfish} -r {input.reference} -v {params.index} -t {threads} -j {threads} -s {wildcards.sample} -g -o {params.prefix}
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

rule extract_SVs:
    input:
        get_dir('PG','samples.all.pangenie_{pangenie_mode}.vcf.gz')
    output:
        get_dir('PG','samples.all.pangenie_{pangenie_mode}.SVs.vcf.gz')
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools view --threads {threads} -i 'abs(ILEN)>50' -o {output} {input}
        '''

rule merge_with_population_SR:
    input:
        DV = config['DV-SR'],
        pangenie = get_dir('PG','samples.all.pangenie_{pangenie_mode}.vcf.gz')
    output:
        get_dir('PG','samples.all.pangenie_{pangenie_mode}_DV.vcf.gz')
    threads: 8
    resources:
        mem_mb = 2500,
        disk_scratch = 75,
        walltime = '24:00'
    shell:
        '''
        bcftools concat -a -D --threads {threads} -Ou {input}| \
        bcftools view -e 'FILTER="MONOALLELIC"' -Ou - | \
        bcftools norm --threads {threads} -f {config[reference]} -m -any -Ou - | \
        bcftools norm --threads {threads} -f {config[reference]} -d none -Ou - | \
        bcftools norm --threads {threads} -f {config[reference]} -m +any -Ou - | \
        bcftools sort -T $TMPDIR -Ou - | \
        bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%REF\_%ALT' -o {output} -
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
