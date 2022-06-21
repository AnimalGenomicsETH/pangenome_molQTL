from pathlib import PurePath

wildcard_constraints:
    _pass = r'permutations|conditionals',
    chunk = r'\d*',
    MAF = r'\d*'

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

rule beagle_phase_vcf:
    input:
        config['vcf']
    output:
        config['vcf_phased']
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        mem = lambda wildcards, threads, resources: int(threads*resources.mem_mb/1024)
    threads: 24
    resources:
        mem_mb = 5000,
        walltime = '24:00'
    shell:
        '''
        java -Xmx{params.mem}g -jar {config[beagle]} \
        gt={input} \
        ne=200 \
        out={params.out} \
        nthreads={threads}
        '''

rule tabix_split:
    input:
        lambda wildcards: config[wildcards.vcf]
    output:
        'gwas/{vcf}.{chr}.vcf.gz'
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        tabix -h {input} {wildcards.chr} | bgzip -c -@ {threads} > {output}
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
        panel = 'gwas/vcf_phased.{chr}.vcf.gz',
        chip = 'gwas/chip.{chr}.vcf.gz',
        #chip = config['chip'],
        exclude = 'gwas/chip.exclude'
    output:
        chip = 'gwas/chip.beagle.{chr}.vcf.gz'
    threads: 12
    resources:
        mem_mb = 25000,
        walltime = lambda wildcards: '4:00' if int(wildcards.chr) > 2 else '24:00'
    params:
        chip = lambda wildcards, output: PurePath(output.chip).with_suffix('').with_suffix(''),
        mem = lambda wildcards, threads, resources: int(threads*resources.mem_mb/1024)
    shell:
        '''
        java -Xmx{params.mem}g -jar {config[beagle]} ref={input.panel} gt={input.chip} out={params.chip} ne=200 nthreads={threads} excludesamples={input.exclude}
        tabix -fp vcf {output}
        '''

rule bcftools_norm:
    input:
        'gwas/chip.beagle.{chr}.vcf.gz'
    output:
        'gwas/chip.beagle.{chr}.normed.vcf.gz'
    threads: 4
    resources:
        mem_mb = 2500,
        disk_scratch = 50
    shell:
        '''
        bcftools norm --threads {threads} -f {config[reference]} -m -any {input} | bcftools sort -T $TMPDIR - | bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' -o {output} -
        tabix -fp vcf {output}
        '''

rule bcftools_concat:
    input:
        expand('gwas/chip.beagle.{chr}.normed.vcf.gz',chr=range(1,30))
    output:
        'gwas/chip.beagle.normed.vcf.gz'
    threads: 6
    resources:
        mem_mb = 5000,
        walltime = '24:00'
    shell:
        '''
        bcftools concat --threads {threads} -o {output} {input}
        tabix -fp vcf {output}
        '''
 
rule plink_convert:
    input:
        'gwas/chip.beagle.vcf.gz'
    output:
        'gwas/chip.beagle.pgen'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    shell:
        '''
        plink2 --vcf {input} --chr-set 30 --make-pgen --memory 8000 --threads {threads} --out {params.out} 
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

rule qtltools_prepare:
    input:
        config['vcf_phased']
    output:
        #temp
        'eQTL/variants.normed.vcf.gz'
    threads: 4
    resources:
        mem_mb = 2500,
        disk_scratch = 50
    shell:
        '''
        bcftools norm --threads {threads} -f {config[reference]} -m -any {input} | bcftools sort -T $TMPDIR - | bcftools annotate --threads {threads} --set-id '%CHROM\_%POS\_%TYPE\_%REF\_%ALT' -o {output} -
        tabix -fp vcf {output}
        '''

rule exclude_MAF:
    input:
        'eQTL/variants.normed.vcf.gz'
    output:
        'eQTL/exclude_sites.{MAF}.txt'
    shell:
        '''
        bcftools query -f '%ID\n' -i 'MAF<0.{wildcards.MAF}' {input} > {output}
        '''

rule qtltools_parallel:
    input:
        vcf = 'eQTL/variants.normed.vcf.gz', #config['vcf_phased'],
        #bed = config['bed'],
        #cov = config['cov'],
        #mapping = lambda wildcards: '{qtl}/permutations_all.thresholds.txt' if wildcards._pass == 'conditionals' else []
        bed = lambda wildcards: config['genes'][wildcards.qtl],
        cov = lambda wildcards: config['covariates'][wildcards.qtl],
        mapping = lambda wildcards: '{qtl}/permutations_all.thresholds.{MAF}.txt' if wildcards._pass == 'conditionals' else [],
        exclude = 'eQTL/exclude_sites.{MAF}.txt'
    output:
        merged = temp('{qtl}/{_pass}.{chunk}.{MAF}.txt')
        #merged = temp('eQTL/{_pass}.{chunk}.{MAF}.txt')
    params:
        _pass = lambda wildcards,input: f'--permute {config["permutations"]}' if wildcards._pass == 'permutations' else f'--mapping {input.mapping}',
        debug = '--silent' if 'debug' in config else '',
        grp = lambda wildcards: '--grp-best' if wildcards.qtl == 'sQTL' else ''
    threads: 1
    resources:
        mem_mb = 1024,
        walltime = lambda wildcards: '4:00' if wildcards._pass == 'permutations' else '30'
    shell:
        '''
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} {params._pass} {params.grp} --window {config[window]} --normal --chunk {wildcards.chunk} {config[chunks]} --out {output} {params.debug}
        '''

localrules: qtltools_gather, qtltools_postprocess

rule qtltools_gather:
    input:
        expand('eQTL/{{_pass}}.{chunk}.{{MAF}}.txt',chunk=range(1,config['chunks']+1))
    output:
        'eQTL/{_pass}.{MAF}.txt'
    resources:
        mem_mb = 3000,
        walltime = '20'
    shell:
        '''
        sort -k9,9n -k10,10n {input} > {output}
        '''

rule qtltools_FDR:
    input:
        'eQTL/permutations.{MAF}.txt'
    output:
        'eQTL/permutations_all.thresholds.{MAF}.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    shell:
        '''
        Rscript /cluster/work/pausch/alex/software/qtltools/scripts/qtltools_runFDR_cis.R {input} 0.05 {params.out}
        '''

rule qtltools_postprocess:
    input:
        'eQTL/conditionals.{MAF}.txt'
    output:
        'eQTL/significant_hits.{minS}.{MAF}.fastman'
    shell:
        '''
        echo "CHR\tsize\tBP\tSNP\tA1\tTEST\tNMISS\tBETA\tSTAT\tP" > {output}
        awk -v L={wildcards.minS} '($11-$10)>=L {{print $9"\t"($11-$10)"\t"$10"\t"$8"\tN\teQTL\t2\t"$20"\t"$19"\t"$18}}' {input} >>  {output}
        '''



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
