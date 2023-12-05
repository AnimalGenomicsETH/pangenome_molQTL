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
        vcf = 'SR_SV/{sample}_insurveyor/out.pass.vcf.gz'
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
    envmodules:
        'boost'
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
    envmodules:
        'boost'
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
    envmodules:
        'boost'
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
        lambda wildcards: expand(rules.delly_call_forced.output,sample=config['HiFi']) if wildcards.caller == 'delly.forced' else (expand(rules.survtyper.output['vcf'],sample=config['HiFi']) if wildcards.caller == 'insurveyor.forced' else expand(rules.insurveyor.output['vcf'],sample=config['HiFi']))
    output:
        'SR_SV/cohort.{caller}.vcf.gz'
    shell:
        '''
        bcftools merge -m id -o {output} {input}
        tabix -fp vcf {output}
        '''

rule delly_filter:
    input:
        expand(rules.bcftools_merge.output,caller='delly.forced')
    output:
        'SR_SV/cohort.delly.filtered.vcf.gz'
    envmodules:
        'gcc/9.3.0',
        'boost/1.74.0',
        'gsl/2.7.1'
    threads: 1
    resources:
        mem_mb = 2500
    shell:
        '''
        /cluster/work/pausch/alex/software/delly/src/delly filter -f germline --altaf 0.01 --minsize 50 --pass -o $TMPDIR/sample.bcf {input}
        bcftools view -o {output} $TMPDIR/sample.bcf
        tabix -fp vcf {output}
        '''

rule survtyper:
    input:
        vcf = expand(rules.bcftools_merge.output,caller='insurveyor'),
        bam = '/cluster/work/pausch/inputs/bam/BTA_eQTL/{sample}.bam'
    output:
        _dir = directory('SR_SV/{sample}_insurveyor_forced'),
        vcf = 'SR_SV/{sample}_insurveyor_forced/genotyped.vcf.gz'
    threads: 4
    resources:
        mem_mb = 5000
    shell:
        '''
        python /cluster/work/pausch/alex/software/SurVTyper/survtyper.py --threads {threads} --samplename {wildcards.sample} {input.vcf} {input.bam} {output._dir} {config[reference]}
        '''

rule merge_callers:
    input:
        delly = rules.delly_filter.output,
        insurveyor = expand(rules.bcftools_merge.output,caller='insurveyor.forced')
    output:
        'SR_SV/SVs.vcf.gz'
    threads: 2
    shell:
        '''
        bcftools view {input.delly} -e "SVTYPE=='INS'" -o $TMPDIR/deldel.vcf.gz
        tabix -p vcf $TMPDIR/deldel.vcf.gz
        bcftools concat -a --threads 2 -D -o {output} $TMPDIR/deldel.vcf.gz {input.insurveyor}
        tabix -p vcf {output}
        '''
