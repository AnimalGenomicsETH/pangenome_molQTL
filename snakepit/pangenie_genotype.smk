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

rule vcfwave:
    input:
        vcf = config['panel']
    output:
        'pangenie_panel.vcfwave.vcf'
    threads: 8
    resources:
        mem_mb = 1500,
        walltime = '1h'
    shell:
        '''
        vcfwave -t {threads} --quiet {input} > {output}
        '''

rule pangenie_index:
    input:
        reference = config['reference'],
        vcf = rules.vcfwave.output #config['panel']
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

rule jellyfish_count_index:
    input:
        rules.pangenie_index.output[1]
    output:
        'pangenie.jf'
    threads: 4
    resources:
        mem_mb = 7000
    shell:
        '''
        jellyfish count -L 1 -U 10000 -m 31 -s 3000000000 -p 126 -c 7 -C -t {threads} -o {output} {input}
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
        walltime = '24h'
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        pigz -p {threads} -c -d {input.fastq} > $TMPDIR/{wildcards.sample}.fastq
        PanGenie -i {input} -r {input.reference} -v {input.vcf} -t {threads} -j {threads} -s {wildcards.sample} -g {params.phasing} -o {params.prefix}
        '''

rule pangenie_genotype:
    input:
        reference = config['reference'],
        fastq = lambda wildcards: get_sample_location(wildcards.sample),
        pangenie_index = rules.pangenie_index.output
    output:
        get_dir('PG','{sample}.all.pangenie_genotyping.vcf')
    params:
        prefix = lambda wildcards, output: str(PurePath(output[0]).with_suffix('')).replace(f'_genotyping',''),
        index = lambda wildcards, input: PurePath(input.pangenie_index[0]).with_suffix('')
    threads: 12
    resources:
        mem_mb = 8000,
        scratch = '75G',
        walltime = '4h'
    shell:
        '''
        fastp -w {threads} -i {input.fastq[0]} -I {input.fastq[1]} --stdout -g --thread {threads} --html /dev/null --json /dev/null --dont_eval_duplication | seqtk seq -A > $TMPDIR/reads.fa
        pangenie -i $TMPDIR/reads.fa -r {input.reference} -v {params.index} -t {threads} -j {threads} -s {wildcards.sample} -g -o {params.prefix}
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

rule beagle5_impute:
    input:
        '/cluster/work/pausch/alex/eQTL_GWAS/variants/DV-SR/cohort.autosomes.WGS.vcf.gz'
    output:
        '/cluster/work/pausch/alex/eQTL_GWAS/variants/DV-SR/cohort.autosomes.WGS.imputed.vcf.gz'
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        name = lambda wildcards, output: PurePath(output[0]).name
    threads: 24
    resources:
        mem_mb = 5000,
        walltime = '24h'
    shell:
        '''
        java -jar -Xss25m -Xmx95G /cluster/work/pausch/alex/software/beagle.22Jul22.46e.jar gt={input} nthreads={threads} out={params.prefix}
        mv {output[0]} $TMPDIR/{params.name}
        bcftools reheader -f {config[reference]}.fai -o {output[0]} $TMPDIR/{params.name}
        tabix -fp vcf {output[0]}
        '''


rule merge_with_population_SR:
    input:
        DV = rules.beagle5_impute.output,#config['DV-SR'],
        pangenie = get_dir('PG','samples.all.pangenie_{pangenie_mode}.vcf.gz')
    output:
        get_dir('PG','samples.all.pangenie_{pangenie_mode}_DV.vcf.gz')
    threads: 8
    resources:
        mem_mb = 2500,
        disk_scratch = 75,
        walltime = '24h'
    shell:
        '''
        bcftools concat -a -D --threads {threads} {input}| \
        grep -v "MONOALLELIC" |\
        bcftools norm --threads {threads} -f {config[reference]} -m -any -Ou - | \
        bcftools norm --threads {threads} -f {config[reference]} -d none -Ou - | \
        bcftools norm --threads {threads} -f {config[reference]} -m +any -Ou - | \
        bcftools sort -T $TMPDIR -Ou - | \
        bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' -o {output} -
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
