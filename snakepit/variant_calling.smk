rule all:
    input:
        expand('variant_calling/samples.{V}.vcf.gz',V=('sniffles','bcftools')),
        'DV-SR/cohort.autosomes.WGS.vcf.gz',
        'DV-LR/cohort.autosomes.WGS.vcf.gz'

rule minimap2_align:
    input:
        lambda wildcards: config['HiFi'][wildcards.sample]
    output:
        bam = temp('alignments/{sample}.HiFi.bam'),
        csi = temp('alignments/{sample}.HiFi.bam.csi')
    threads: 12
    resources:
        mem_mb = 6000,
        disk_scratch = 100,
        walltime = '24:00'
    shell:
        '''
        minimap2 -axmap-hifi -t {threads} {config[reference]} {input} |  samtools sort - -m 3000M -@ 4 -T $TMPDIR --write-index -o {output.bam}
        '''

rule sniffles_call:
    input:
        bam = 'alignments/{sample}.HiFi.bam'
    output:
        vcf = temp('variant_calling/{sample}.sniffles.vcf.gz'),
        snf = temp('variant_calling/{sample}.sniffles.snf')
    threads: 12
    resources:
        mem_mb = 5000
    conda:
        '/cluster/work/pausch/alex/software/miniconda3/envs/sniffles'
    shell:
        '''
        sniffles --input {input.bam} --reference {config[reference]} --threads {threads} --vcf {output.vcf} --snf {output.snf}
        '''

rule sniffles_merge:
    input:
        snfs = expand('variant_calling/{sample}.sniffles.snf',sample=config['HiFi'])
    output:
        vcf = 'variant_calling/samples.sniffles.vcf.gz'
    threads: 4
    resources:
        mem_mb = 5000
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.snfs} --reference {config[reference]} --threads {threads} --vcf {output.vcf}
        '''

rule bcftools_merge:
    input:
        vcfs = expand('variant_calling/{sample}.sniffles.vcf.gz',sample=config['HiFi'])
    output:
        vcf = 'variant_calling/samples.bcftools.vcf.gz'
    threads: 4
    resources:
        mem_mb = 5000
    shell:
        '''
        bcftools merge --threads {threads} -o {output} {input.vcfs}
        tabix -fp vcf {output}
        '''

### IMPORT DEEPVARIANT MODULES
module DV_SR_calling:
    snakefile:
        '/cluster/work/pausch/alex/BSW_analysis/snakepit/deepvariant.smk'
    config: config["DV-SR-workflow"]
    prefix: 'DV-SR'

use rule * from DV_SR_calling as DV_SR_*

module DV_LR_calling:
    snakefile:
        '/cluster/work/pausch/alex/BSW_analysis/snakepit/deepvariant.smk'
    config: config["DV-LR-workflow"]
    prefix: 'DV-LR'

use rule * from DV_LR_calling as DV_LR_*

use rule deepvariant_make_examples from DV_LR_calling as DV_LR_deepvariant_make_examples with:
    input:
        ref = multiext(config['reference'],'','.fai'),
        bam = multiext('alignments/{animal}.HiFi.bam','','.csi')

## Assembly rules:
rule minimap2_align_asm:
    input:
        lambda wildcards: config['asm'][wildcards.sample]
    output:
        'alignments/{sample}.asm.paf'
    threads: 4
    resources:
        mem_mb = 10000
    shell:
        '''
        minimap2 -cx asm5 -t {threads} --cs {config[reference]} {input} > {output}
        '''

rule paftools_call:
    input:
        'alignments/{sample}.asm.paf'
    output:
        'variant_calling/{sample}.asm.paf'
    threads: 1
    resources:
        mem_mb = 2000
    shell:
        '''
        sort -k6,6 -k8,8n {input} | paftools.js call -f {config[reference]} {input} > {output.vcf}
        '''
