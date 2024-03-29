from pathlib import PurePath

rule all:
    input:
        'SR_SV/SVs.vcf.gz',
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
        mem_mb = 15000,
        walltime = '24h'
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
        mem_mb = 7000
    shell:
        '''
        mkdir -p {output._dir}
        python /cluster/work/pausch/alex/software/INSurVeyor/insurveyor.py --threads {threads} --samplename {wildcards.sample} {input[0]} {output._dir} {config[reference]}
        '''

rule SurVClusterer:
    input:
        expand(rules.insurveyor.output['vcf'],sample=config['HiFi'])
    output:
        multiext('SR_SV/cohort.insurveyor','.sv','.vcf.gz','.vcf.gz.tbi')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 6
    resources:
        mem_mb = 4000
    shell:
        '''
        awk '$1=="-" {{print $2, "SR_SV/"$2"_insurveyor/out.pass.vcf.gz"}}' config/SR_SV.yaml > $TMPDIR/samples.fofn
        /cluster/work/pausch/alex/software/SurVClusterer/clusterer $TMPDIR/samples.fofn {config[reference]} -t {threads} --min-overlap-precise 0.95 --max-dist-precise 25 --overlap-for-ins -o {params.prefix}
        tabix -fp vcf {output[1]}
        '''

rule delly_call_denovo:
    input:
        '/cluster/work/pausch/inputs/bam/BTA_eQTL/{sample}.bam'
    output:
        'SR_SV/{sample}.delly.denovo.bcf'
    envmodules:
        'gcc/9.3.0',
        'boost/1.74.0',
        'gsl/2.6'
    threads: 4
    resources:
        mem_mb = 4000
    shell:
        '''
        /cluster/work/pausch/alex/software/delly/src/delly call -g {config[reference]} -o {output} {input}
        #bcftools index {output}
        '''

rule delly_merge:
    input:
        expand(rules.delly_call_denovo.output,sample=config['HiFi'])
    output:
        'SR_SV/cohort.delly.denovo.bcf'
    envmodules:
        'gcc/9.3.0',
        'boost/1.74.0',
        'gsl/2.6'
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
        'gcc/9.3.0',
        'boost/1.74.0',
        'gsl/2.6'
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
    resources:
        mem_mb = 10000
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
        'gsl/2.6'
    threads: 1
    resources:
        mem_mb = 2500
    shell:
        '''
        /cluster/work/pausch/alex/software/delly/src/delly filter -f germline --altaf 0.01 --minsize 50 --pass -o $TMPDIR/sample.bcf {input}
        bcftools view -o {output} $TMPDIR/sample.bcf
        tabix -fp vcf {output}
        '''

#awk '$1=="-" {print $2, "SR_SV/"$2"_insurveyor/out.pass.vcf.gz"}' ../config/SR_SV.yaml > ../insurveyor.fofn
rule survtyper:
    input:
        vcf = rules.SurVClusterer.output[1],
        bam = '/cluster/work/pausch/inputs/bam/BTA_eQTL/{sample}.bam'
    output:
        _dir = directory('SR_SV/{sample}_insurveyor_forced'),
        vcf = 'SR_SV/{sample}_insurveyor_forced/genotyped.vcf.gz'
    threads: 6
    resources:
        mem_mb = 8000,
        walltime = '4h'
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
        bcftools concat -a --threads 2 -D -o $TMPDIR/concat.vcf.gz $TMPDIR/deldel.vcf.gz {input.insurveyor}
        tabix -p vcf $TMPDIR/concat.vcf.gz
        bcftools view -e INFO/INCOMPLETE_ASSEMBLY!=0 -o {output} $TMPDIR/concat.vcf.gz 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
        tabix -p vcf {output}
        '''
