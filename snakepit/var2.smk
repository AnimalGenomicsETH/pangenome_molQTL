rule all:
    input:
        expand('variant_calling/{sample}.FORCE.sniffles.vcf.gz',sample=config['HiFi'])

rule sniffles_call:
    input:
        bam = '/cluster/work/pausch/alex/CCS_eQTL_cohort/ARS/{sample}.mm2.cram'
    output:
        vcf = temp('variant_calling/{sample}.sniffles.vcf.gz'),
        snf = temp('variant_calling/{sample}.sniffles.snf')
    threads: 4
    resources:
        mem_mb = 2500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.bam} --reference {config[reference]} --sample-id {wildcards.sample} --threads {threads} --vcf {output.vcf} --snf {output.snf}
        '''

rule sniffles_merge:
    input:
        snfs = expand('variant_calling/{sample}.FORCE.sniffles.snf',sample=config['HiFi'])
    output:
        vcf = 'variant_calling/samples.FORCE.sniffles.vcf.gz'
    threads: 2
    resources:
        mem_mb = 2500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.snfs} --reference {config[reference]} --threads {threads} --vcf {output.vcf}
        '''

rule sniffles_genotype:
    input:
        bam = '/cluster/work/pausch/alex/CCS_eQTL_cohort/ARS/{sample}.mm2.cram',
        vcf = '/cluster/work/pausch/alex/eQTL_GWAS/SV_ACCURACY/pangenie.vcf'
    output:
        vcf = temp('variant_calling/{sample}.FORCE.sniffles.vcf.gz'),
    threads: 4
    resources:
        mem_mb = 1500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.bam} --genotype-vcf {input.vcf} --reference {config[reference]} --sample-id {wildcards.sample} --threads {threads} --vcf {output.vcf}
        '''
