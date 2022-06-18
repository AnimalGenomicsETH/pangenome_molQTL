

reference_file="/cluster/work/pausch/inputs/ref/BTA/UCD1.2/ARS-UCD1.2_Btau5.0.1Y.fa"
OUT_DIR="/cluster/work/pausch/alex/BSW_OBV_SVs/pangenie_imputed/"
IN_DIR="/cluster/work/pausch/alex/BSW_OBV_SVs/pangenie/"


#tools
BEAGLE='/cluster/work/pausch/alex/software/beagle.27Jan18.7e1.jar'


##parameters for split_vcf
#nvar=50_000
#overlap=5_000


#nvar=200_000
#overlap=5_000
nvar=60_000
overlap=5_000



##paramterrs for Beagle
window= int(nvar*1.5 + overlap) #beagles default window size is 50k variants [this is the biggest possible window -- by split.vcf]

##chromosomes=list(range (1,30)) + ["X", "Y"]
chromosomes = list(map(str,range(1,30)))
#chromosomes = [el for el in chromosomes if el !=26]
#chromosomes=["Y"]

from pathlib import PurePath

wildcard_constraints:
    chunk = r'\d*:\d*-\d*'

rule all:				      
    input:
        'pangenie_DV_imputed/samples.all.pangenie_genotyping_DV.imputed.vcf.gz'
        #expand (OUT_DIR + "CHR{chr}/PHASED/clean_beagle4.1.vcf.gz", chr=chromosomes),
        #expand (OUT_DIR + "CHR{chr}/PHASED/clean_beagle4.1.vcf.gz.tbi", chr=chromosomes)


rule determine_scatter:
    output:
        'scatter.txt'
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
        scatter = 'scatter.txt'
    output:
        directory('pangenie_DV_scatter')
    threads: 4
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools +scatter {input.vcf} --threads {threads} -S {input.scatter} -Oz -p chunk- -o {output}
        '''

rule beagle4:
    input:
        'pangenie_DV_scatter/chunk-{chunk}.vcf.gz'
    output:
        'pangenie_DV_scatter/chunk-{chunk}.imputed.vcf.gz'
    threads: 6
    resources:
        mem_mb = 4000
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
	    ne = 200,
        mem = lambda wildcards, threads, resources: int(resources.mem_mb*threads/1000),
        window = config['window'],
        beagle = config['beagle']
    shell:
        '''
        java -Xmx{params.mem}g -Xss50m -jar {params.beagle} \
        gl={input} \
        ne={params.ne} \
        gprobs=true \
        nthreads={threads} \
        out={params.prefix}
        '''

def aggregate_scatter(wildcards):
    checkpoint_dir = checkpoints.bcftools_scatter.get(**wildcards).output[0]
    return expand('pangenie_DV_scatter/chunk-{chunk}.imputed.vcf.gz',chunk=glob_wildcards(PurePath(checkpoint_dir).joinpath('chunk-{chunk,\d*:\d*-\d*}.vcf.gz')).chunk)

rule merge_vcfs:
    input:
        aggregate_scatter
    output:
        'pangenie_DV_imputed/samples.all.pangenie_genotyping_DV.imputed.vcf.gz'
    threads: 2
    resources:
        mem_mb = 2500
    script:
        '''
        bcftools concat --allow-overlaps --ligate --threads {threads} -o {output} {input}
        tabix -p vcf {output}
        '''
        


