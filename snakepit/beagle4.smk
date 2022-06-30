from pathlib import PurePath
import subprocess

wildcard_constraints:
    chunk = r'\d*'

chromosomes = list(map(str,range(1,30)))

rule all:				      
    input:
        'pangenie_DV_imputed/samples.all.pangenie_genotyping_DV.imputed.phased.vcf.gz'

rule determine_scatter:
    output:
        temp('scatter.txt')
    params:
        window = config['window'],
        overlap = config['overlap']
    run:
        with open(f'{config["reference"]}.fai') as fai, open(output[0],'w') as fout:
            for line in fai:
                chrom, size = line.split()[:2]
                if chrom not in chromosomes:
                    continue
                for coord in range(1,int(size),params.window):
                    fout.write(f'{chrom}:{max(1,coord-params.overlap)}-{coord+params.window+params.overlap}\n')

checkpoint bcftools_scatter:
    input:
        vcf = 'pangenie/samples.all.pangenie_genotyping_DV.vcf.gz',
        #scatter = 'scatter.txt'
    output:
        temp(directory('pangenie_DV_scatter'))
    threads: 4
    resources:
        mem_mb = 2500
    params:
        scatter = ','.join(chromosomes),
        prefix = 'raw.chr'
    shell:
        '''
        bcftools +scatter {input.vcf} --threads {threads} -s {params.scatter} -Oz -p {params.prefix} -o {output}
        '''

rule beagle4_impute:
    input:
        lambda wildcards: f'pangenie_DV_scatter/{"raw" if wildcards.out == "phased" else "raw"}.chr{{chunk}}.vcf.gz'
    output:
        temp('pangenie_DV_scatter/{out,imputed|phased}.chr{chunk}.vcf.gz')
    threads: 8
    resources:
        mem_mb = 4000,
        walltime = '4:00'
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
	    ne = 200,
        mem = lambda wildcards, threads, resources: int(resources.mem_mb*threads/1000),
        window = config['window'],
        beagle = config['beagle5']
    shell:
        '''
        java -Xmx{params.mem}g -Xss50m -jar {params.beagle} \
        gt={input} \
        nthreads={threads} \
        out={params.prefix}
        tabix -fp vcf {output}
        '''

def empty_vcf(vcf):
    #hardcode removing the .imputed to get the original vcf filename
    return subprocess.getoutput(f"zgrep -cm1 -v '#' {vcf.replace('.imputed','')}").strip() != "0"

def aggregate_scatter(wildcards):
    checkpoint_dir = checkpoints.bcftools_scatter.get(**wildcards).output[0]
    return sorted([f'pangenie_DV_scatter/phased.chr{chunk}.vcf.gz' for chunk in glob_wildcards(PurePath(checkpoint_dir).joinpath('raw.chr{chunk,\d*}.vcf.gz')).chunk])

rule merge_vcfs:
    input:
        aggregate_scatter
    output:
        'pangenie_DV_imputed/samples.all.pangenie_genotyping_DV.imputed.phased.vcf.gz'
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools concat --threads {threads} -a -d exact -o {output} {input}
        tabix -p vcf {output}
        '''

rule beagle_phase_vcf:
    input:
        'pangenie_DV_imputed/samples.all.pangenie_genotyping_DV.imputed.vcf.gz'
    output:
        'pangenie_DV_imputed/samples.all.pangenie_genotyping_DV.imputed.phased.X.vcf.gz'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        mem = lambda wildcards, threads, resources: int(threads*resources.mem_mb/1024),
        ne = 200
    threads: 12
    resources:
        mem_mb = 10000,
        walltime = '24:00'
    shell:
        '''
        java -Xmx{params.mem}g -jar {config[beagle5]} \
        gt={input} \
        out={params.out} \
        nthreads={threads}
        '''
