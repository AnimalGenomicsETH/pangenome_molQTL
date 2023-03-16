rule all:
    input:
        expand('variant_calling/{sample}.denovo.sniffles.vcf.gz',sample=config['HiFi']),
        'variant_calling/samples.denovo.sniffles.vcf.gz'

rule sniffles_call:
    input:
        bam = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long-read_alignments/PacBio_CCS/eQTL/{sample}.mm2.cram'
    output:
        vcf = temp('variant_calling/{sample}.denovo.sniffles.vcf.gz'),
        snf = temp('variant_calling/{sample}.denovo.sniffles.snf')
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
        snfs = expand('variant_calling/{sample}.{call}.sniffles.snf',sample=config['HiFi'],allow_missing=True)
    output:
        vcf = 'variant_calling/samples.{call}.sniffles.vcf.gz'
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
        vcf = lambda wildcards: '' if wildcards.call == 'denovo' else '/cluster/work/pausch/alex/eQTL_GWAS/SV_ACCURACY/genotyping.vcf.gz' #'/cluster/work/pausch/alex/eQTL_GWAS/SV_ACCURACY/pangenie.vcf'
    output:
        vcf = temp('variant_calling/{sample}.{call}.sniffles.vcf.gz')
    params:
    threads: 4
    resources:
        mem_mb = 1500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.bam} --genotype-vcf {input.vcf} --reference {config[reference]} --sample-id {wildcards.sample} --threads {threads} --vcf {output.vcf}
        '''
