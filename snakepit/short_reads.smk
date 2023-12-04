rule all:
    input:
        'SR_SV/cohort.delly.forced.vcf.gz',
        expand('SR_SV/{sample}_{caller}',caller=('insurveyor',),sample=config['HiFi'])

rule picard_add_MQ:
    input:
        '/cluster/work/pausch/inputs/bam/BTA_eQTL/{sample}.bam'
    output:
        multiext('SR_SV/{sample}.bam','','.bai')
    envmodules:
        'gcc/11.4.0',
        'picard/3.1.1'
    threads: 1
    resources:
        mem_mb = 7500
    shell:
        '''
        picard FixMateInformation -I {input} -O {output[0]}
        samtools index {output[0]}
        '''

rule insurveyor:
    input:
        rules.picard_add_MQ.output
    output:
        _dir = directory('SR_SV/{sample}_insurveyor'),
        #vcf = 'SR_SV/{sample}.insurveyor.vcf.gz'
    threads: 4
    resources:
        mem_mb = 5000
    shell:
        '''
        mkdir -p {output._dir}
        python /cluster/work/pausch/alex/software/INSurVeyor/insurveyor.py --threads {threads} --samplename {wildcards.sample} {input[0]} {output._dir} {config[reference]}
        '''

rule delly_call_denovo:
    input:
        '/cluster/work/pausch/inputs/bam/BTA_eQTL/{sample}.bam'
    output:
        'SR_SV/{sample}.delly.denovo.bcf'
    threads: 4
    resources:
        mem_mb = 4000
    shell:
        '''
        /cluster/work/pausch/alex/software/delly/src/delly call -g {config[reference]} -o {output} {input}
        bcftools index {output}
        '''

rule delly_merge:
    input:
        expand(rules.delly_call_denovo.output,sample=config['HiFi'])
    output:
        'SR_SV/cohort.delly.denovo.bcf'
    threads: 4
    resources:
        mem_mb = 2500
    shell:
        '''
        /cluster/work/pausch/alex/software/delly/src/delly merge --precise --pass --minsize 50 --vaf 0.01 -o {output} {input}
        '''

rule delly_call_forced:
    input:
        bam = '/cluster/work/pausch/inputs/bam/BTA_eQTL/{sample}.bam',
        panel = rules.delly_merge.output
    output:
        'SR_SV/{sample}.delly.forced.vcf.gz'
    threads: 4
    resources:
        mem_mb = 4000
    shell:
        '''
        /cluster/work/pausch/alex/software/delly/src/delly call -g {config[reference]} -v {input.panel} -o $TMPDIR/sample.bcf {input.bam}
        bcftools view -o {output} $TMPDIR/sample.bcf
        tabix -fp vcf {output}
        '''

rule bcftools_merge:
    input:
        expand(rules.delly_call_forced.output,sample=config['HiFi'])
    output:
        'SR_SV/cohort.delly.forced.vcf.gz'
    shell:
        '''
        bcftools merge -m id -o {output} {input}
        tabix -fp vcf {output}
        '''
