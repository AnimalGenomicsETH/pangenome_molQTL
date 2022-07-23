from pathlib import PurePath

wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    chunk = r'\d+',
    chrom = r'\d+',
    MAF = r'\d+',
    vcf = r'(eQTL|gwas)/\S+'

rule all:
    input:
        'test/chip.17.67869210.mlma'

rule bcftools_view:
    input:
        'pangenie/samples.annotated.vcf.gz'
    output:
        'test/region.{chr}.{region}.vcf.gz'
    params:
        region = lambda wildcards: wildcards.chr + ':' + f'{int(wildcards.region)-2e6}-{int(wildcards.region)+2e6}'
    shell:
        '''
        bcftools view -O v {input} {params.region} | perl -pe "s/\s\.:/\t.\/.:/g" | bgzip -c > {output}
        '''

rule split_chromomosome:
    input:
        lambda wildcards: config[wildcards._vcf]
    output:
        'gwas/{_vcf}.{chrom}.vcf.gz'
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools view --threads {threads} -o {output} {input} {wildcards.chrom}
        tabix -fp vcf {output}
        '''

rule exclude_duplicates:
    input:
        panel = config['vcf_phased'],
        chip = config['chip']
    output:
        exclude = 'gwas/chip.exclude'
    shell:
        'grep -f <(bcftools query -l {input.panel}) <(bcftools query -l {input.chip}) > {output.exclude}'

rule beagle_impute_chip:
    input:
        panel = 'gwas/vcf_phased.{chrom}.vcf.gz',
        chip = 'gwas/chip.{chrom}.vcf.gz',
        exclude = 'gwas/chip.exclude'
    output:
        chip = 'gwas/chip.beagle.{chrom}.vcf.gz'
    threads: 12
    resources:
        mem_mb = lambda wildcards: 20000 if int(wildcards.chrom) > 2 else 25000,
        walltime = lambda wildcards: '4:00' if int(wildcards.chrom) > 2 else '24:00'
    params:
        chip = lambda wildcards, output: PurePath(output.chip).with_suffix('').with_suffix(''),
        mem = lambda wildcards, threads, resources: int(threads*resources.mem_mb/1024),
        ne = 200,
        beagle = config['beagle5']
    retries: 1
    shell:
        '''
        java -Xmx{params.mem}g -jar {params.beagle} \
        ref={input.panel} \
        gt={input.chip} \
        out={params.chip} \
        ne={params.ne} \
        nthreads={threads} \
        excludesamples={input.exclude}
        tabix -fp vcf {output}
        '''

rule plink2_convert:
    input:
        'gwas/chip.beagle.{chrom}.normed.vcf.gz'
    output:
        'gwas/chip.beagle.{chrom}.bim'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix(''),
        mem = lambda wildcards, threads, resources: int(threads*resources.mem_mb)
    threads: 4
    resources:
        mem_mb = 15000
    shell:
        '''
        plink2 --vcf {input} --chr-set 30 --make-bed --memory {params.mem} --threads {threads} --out {params.out} 
        '''

rule gcta_make_grm:
    input:
        'gwas/chip.beagle.{chrom}.bim'
    output:
        'gwas/chip.beagle.{chrom}.{MAF}.grm.bin'
    params:
        _input = lambda wildcards, input: PurePath(input[0]).with_suffix(''),
        _output = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    threads: 8
    resources:
        mem_mb = 6000
    shell:
        '''
        gcta --thread-num {threads} --autosome-num 29 --make-grm --bfile {params._input} --maf .{wildcards.MAF} --out {params._output}
        '''

rule gcta_mlma:
    input:
        bim = 'gwas/chip.beagle.{chrom}.bim',
        grm = 'gwas/chip.beagle.{chrom}.{MAF}.grm.bin',
        phenotype = 'fpr_sire.phen'#cluster/work/pausch/naveen/GWAS/GCTA/RESULT/fpr/all/sire/phenotypes.txt'
    output:
        'gwas/chip.beagle.{chrom}.{MAF}.mlma'
    params:
        bim = lambda wildcards, input: PurePath(input.bim).with_suffix(''),
        grm = lambda wildcards, input: PurePath(input.grm).with_suffix('').with_suffix(''),
        _output = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    threads: 12
    resources:
        mem_mb = 5000
    shell:
        '''
        gcta --thread-num {threads} --autosome-num 29 --mlma --bfile {params.bim} --grm {params.grm} --out {params._output} --pheno {input.phenotype}
        '''

rule plink_PCA:
    input:
        'gwas/chip.beagle.pgen'
    output:
        'gwas/chip.beagle.PCA.eigenvec'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    shell:
        '''
        plink2 --pfile samples.imputed --pca --out samples.imputed.PCA --memory 8000 --threads 2 --chr-set 30
        '''

#plink2 --pfile samples.imputed --pheno ../all_phenotypes.txt --out samples.imputed.GLM2 --glm --covar samples.imputed.PCA.eigenvec --memory 8000 --threads 2 --chr-set 30

rule gcta:
    input:
        'test/chip.{chr}.{region}.beagle.vcf.gz'
    output:
        'test/chip.{chr}.{region}.mlma'
    threads: 4
    resources:
        mem_mb = 5000
    shell:
        '''
        plink2 --vcf {input} --make-pgen --out test/{wildcards.chr}.{wildcards.region}
        gcta --mlma --pfile test/{wildcards.chr}.{wildcards.region} --pheno vzr.phen --out {output} --thread-num {threads}
        '''

rule normalise_vcf:
    input:
        lambda wildcards: '{vcf}.vcf.gz' if 'gwas' in wildcards.vcf else config['vcf_phased']# gwas/chip.beagle.{chr}.vcf.gz
        #config['vcf_phased']
    output:
        #temp
        '{vcf}.normed.vcf.gz'
        #'gwas/chip.beagle.{chr}.normed.vcf.gz'
        #'eQTL/variants.normed.vcf.gz'
    threads: 4
    resources:
        mem_mb = 2500,
        disk_scratch = 50
    shell:
        '''
        bcftools norm --threads {threads} -f {config[reference]} -m -any {input} -Ou | \
        bcftools sort -T $TMPDIR -Ou - | \
        bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' -o {output} -
        tabix -fp vcf {output}
        '''

rule exclude_MAF:
    input:
        'resources/variants.normed.vcf.gz'
    output:
        'resources/exclude_sites.{MAF}.txt'
    shell:
        '''
        bcftools query -f '%ID\n' -i 'MAF<0.{wildcards.MAF}' {input} > {output}
        '''

def get_pass(_pass,input):
    if _pass == 'permutations':
        return f'--permute {config["permutations"]}'
    elif _pass == 'conditionals':
        return f'--mapping {input.mapping}'
    elif _pass == 'nominals':
        return '--nominal 1.0'
rule qtltools_parallel:
    input:
        vcf = 'resources/variants.normed.vcf.gz',
        exclude = 'resources/exclude_sites.{MAF}.txt',
        bed = lambda wildcards: config['genes'][wildcards.qtl],
        cov = lambda wildcards: config['covariates'][wildcards.qtl],
        mapping = lambda wildcards: '{qtl}/permutations_all.{MAF}.thresholds.txt' if wildcards._pass == 'conditionals' else []
    output:
        merged = temp('{qtl}/{_pass}.{chunk}.{MAF}.txt')
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input),#f'--permute {config["permutations"]}' if wildcards._pass == 'permutations' else f'--mapping {input.mapping}',
        debug = '--silent' if 'debug' in config else '',
        grp = lambda wildcards: '--grp-best' if wildcards.qtl == 'sQTL' else ''
    threads: 1
    resources:
        mem_mb = 12500,
        walltime = lambda wildcards: '24:00' if wildcards._pass == 'permutations' else '4:00'
    shell:
        '''
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} {params._pass} {params.grp} --window {config[window]} --normal --chunk {wildcards.chunk} {config[chunks]} --out {output} {params.debug}
        '''

localrules: qtltools_gather, qtltools_postprocess

rule qtltools_gather:
    input:
        expand('{{qtl}}/{{_pass}}.{chunk}.{{MAF}}.txt',chunk=range(1,config['chunks']+1))
    output:
        '{qtl}/{_pass}.{MAF}.txt'
    resources:
        mem_mb = 3000,
        walltime = '20'
    params:
        sort_key = lambda wildcards: '-k9,9n -k10,10n' if wildcards.qtl == 'eQTL' else '-k11,11n -k12,12n'
    shell:
        '''
        sort {params.sort_key} {input} > {output}
        '''

rule qtltools_FDR:
    input:
        '{qtl}/permutations.{MAF}.txt'
    output:
        '{qtl}/permutations_all.{MAF}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    shell:
        '''
        Rscript /cluster/work/pausch/alex/software/qtltools/scripts/qtltools_runFDR_cis.R {input} 0.05 {params.out}
        '''

rule qtltools_postprocess:
    input:
        '{qtl}/conditionals.{MAF}.txt'
    output:
        '{qtl}/significant_hits.{minS}.{MAF}.fastman'
    params:
        print_key = lambda wildcards: '$9"\t"($11-$10)"\t"$10"\t"$8"\tN\teQTL\t2\t"$20"\t"$19"\t"$18' if wildcards.qtl == 'eQTL' else '$11"\t"($13-$12)"\t"$12"\t"$10"\tN\teQTL\t2\t"$22"\t"$21"\t"$20'
    shell:
        '''
        echo "CHR\tsize\tBP\tSNP\tA1\tTEST\tNMISS\tBETA\tSTAT\tP" > {output}
        awk -v L={wildcards.minS} '($11-$10)>=L {{print {params.print_key}}}' {input} >>  {output}
        '''



# plink2 --cow --vcf {input.vcf} --make-bed --out --threads {threads} --memory {params.memory}
# gcta --pfile --autosome-num 29 --maf {params.maf} --make-grm --out {params.prefix} --threads 4


## DEADCODE ##

#awk ' function abs(v) {return v < 0 ? -v : v} abs(length($4)-length($5))>100' samples.GLM | awk 'NR>1 {print $1"\t"$3"\t"$2"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11"\t"$12}'
#python -c "exec(\"IDs = {L.split()[0]:L.split()[1].rstrip() for L in open('IDs2.txt')}\nfor line in open('eQTL/conditionals.sorted.txt'): print(line.replace('.',IDs['-'.join(line.split()[8:10])],1), end='')\")" > eQTL/conditionals.ID.txt

##java -Xmx40g -jar /cluster/work/pausch/alex/software/beagle.08Feb22.fa4.jar gt={input.panel} out={params.panel} ne=200 nthreads={threads}
#/cluster/work/pausch/alex/XENA/qtltools/bin/QTLtools gwas --bed aligned_genes.bed.gz --vcf imputed_chip_unique.vcf.gz --cov /cluster/work/pausch/xena/eQTL/covariates.txt --normal --out gwas_results.txt
#sed -e '/60436342/{r temp.g3' -e 'd}' DV.eQTL.merged.no_missing.vcf | bgzip -c >  DV.eQTL.merged.no_missing.vcf2.gz;  tabix -fp vcf DV.eQTL.merged.no_missing.vcf2.gz


#awk '{print "0\t"$2"\t"$3}' /cluster/work/pausch/xena/gwas/phenotypes/from_naveen/vzr/BV/sire/phenotypes.txt > vzr.phen
#grep -f <(bcftools query -l imputed_chip.vcf.gz.vcf.gz) <(bcftools query -l chips_1.vcf.gz) > exclude_samples.txt
#zcat impute_region.vcf.gz | perl -pe "s/\s\.:/\t.\/.:/g" | bgzip -c > impute_region2.vcf.gz
#java -jar /cluster/work/pausch/alex/software/beagle.08Feb22.fa4.jar gt=impute_region2.vcf.gz out=imputed_chip.vcf.gz ne=100


#bcftools view -S ^<(bcftools query -l imputed_chipX.vcf.gz | grep : | awk '{print $0}') -o imputed_chip_unique.vcf.gz imputed_chipX.vcf.gz

#zgrep 1_94218322_INDEL imputed_chip_unique.vcf.gz | cut -f 10- | sed -r -e "s/:\S*//g" | tr '\t' '\n' | sort | uniq -c
#"java -jar /cluster/work/pausch/alex/software/beagle.08Feb22.fa4.jar ref=imputed_chip.vcf.gz.vcf.gz gt=chips_1.vcf.gz out=imputed_chipX.vcf.gz ne=100 nthreads=4 excludesamples=exclude_samples.txt
