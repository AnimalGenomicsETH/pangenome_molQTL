rule all:
    input:
        expand('variant_calling/samples.{V}.vcf.gz',V=('sniffles',)),
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
        snfs = expand('variant_calling/{sample}.sniffles.snf',sample=config['HiFi'])
    output:
        vcf = 'variant_calling/samples.sniffles.vcf.gz'
    threads: 2
    resources:
        mem_mb = 2500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.snfs} --reference {config[reference]} --threads {threads} --vcf {output.vcf}
        '''

rule bcftools_autosomes:
    input:
        'variant_calling/samples.sniffles.vcf.gz'
    output:
        'variant_calling/samples.sniffles.autosomes.vcf'
    threads: 2
    resources:
        mem_mb = 4000
    params:
        regions = ','.join(map(str,range(1,30)))
    shell:
        '''
        bcftools view --threads {threads} -r {params.regions} -o {output} {input}
        '''

rule bcftools_split_panel:
    output:
        SV = 'variant_calling/panel.SV.vcf',
        small = 'variant_calling/panel.small.vcf.gz'
    params:
        SV_size = config['SV_size'],
        bcf = '$TMPDIR/normed.bcf'
    threads: 2
    resources:
        mem_mb = 4000,
        disk_scratch = 10
    shell:
        '''
        bcftools norm --threads {threads} -f {config[reference]} -m -any {config[panel]} -Ou {input} > {params.bcf}
        bcftools view -i 'abs(ILEN)>={params.SV_size}' -o {output.SV} {params.bcf}
        bcftools view -i 'abs(ILEN)<{params.SV_size}' -o {output.small} {params.bcf}
        tabix -fp vcf {output.small}
        '''

rule jasmine_intersect:
    input:
        read = 'variant_calling/samples.sniffles.autosomes.vcf',
        asm = 'variant_calling/panel.SV.vcf'
    output:
        'jasmine.vcf'
    params:
        _input = lambda wildcards, input: ','.join(input)
    conda:
        'jasmine'
    threads: 2
    resources:
        mem_mb = 3000,
        disk_scratch = 5
    shell:
        '''
        jasmine --comma_filelist file_list={params._input} threads={threads} out_file={output} out_dir=$TMPDIR spec_reads=0 genome_file={config[reference]} min_seq_id=.5 --pre_normalize
        '''

rule bcftools_isec:
    input:
        read = 'DV-{read}/cohort.autosomes.WGS.vcf.gz',
        asm = 'variant_calling/panel.small.vcf.gz'
    output:
        isec = 'isec_{read}_{strictness}.txt',
        _count = 'isec_{read}_{strictness}.count'
    threads: 2
    resources:
        mem_mb = 3000
    shell:
        '''
        bcftools isec --threads {threads} -n +1 -o {output.isec} -c {wildcards.strictness} {input}
        cut -f 5 {output.isec} | sort | uniq -c > {output._count}
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
