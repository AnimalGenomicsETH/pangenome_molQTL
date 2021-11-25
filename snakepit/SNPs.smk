rule all:
    input:
        expand('{asm}.ARS.vcf',asm=config['assemblies'])

rule minimap2_align:
    input:
        '/cluster/work/pausch/alex/assembly/SV_graph/fasta/{asm}.fasta'
    output:
        '{asm}.ARS.paf'
    threads: 8
    resources:
        mem_mb = 5000,
        walltime = '2:00'
    shell:
        'minimap2 -cx asm5 -t {threads} --cs {config[reference]} {input} > {output}'

rule paf_variants:
    input:
        '{asm}.ARS.paf'
    output:
        '{asm}.ARS.vcf'
    resources:
        mem_mb = 10000,
        walltime = '30'
    shell:
        'sort {input} -k6,6 -k8,8n | paftools.js call -f {config[reference]} - > {output}'


rule prep_vcf:
    input:
        '{asm}.ARS.vcf'
    output:
        multiext('{asm}.ARS.snp.vcf.gz','','.tbi')
    shell:
        '''
        bcftools reheader -s <(echo {wildcards.asm}) {input} | bcftools view -O z -o {output[0]} -v snps - 
        tabix -p vcf {output[0]}
        '''

rule bcftools_merge:
    input:
        expand('{asm}.ARS.snp.vcf.gz',asm=config['assemblies'])
    output:
        multiext('merged.ARS.snp.vcf.gz','','.tbi')
    threads: 2
    resources:
        mem_mb = 2000,
        walltime = '30'
    shell:
        'bcftools merge --threads 2 -O z -o {output[0]} {input}; tabix -p vcf {output[0]}'
