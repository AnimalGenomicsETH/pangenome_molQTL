from pathlib import PurePath

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

rule beagle_impute:
    input:
        panel = 'test/region.{chr}.{region}.vcf.gz',
        chip = 'test/chips.vcf.gz'
    output:
        panel = 'test/region.{chr}.{region}.beagle.vcf.gz',
        exclude = 'test/region.{chr}.{region}.exclude',
        chip = 'test/chip.{chr}.{region}.beagle.vcf.gz'
    threads: 6
    resources:
        mem_mb = 5000
    params:
        panel = lambda wildcards, output: PurePath(output.panel).with_suffix('').with_suffix(''),
        chip = lambda wildcards, output: PurePath(output.chip).with_suffix('').with_suffix('')
    shell:
        '''
        java -Xmx40g -jar /cluster/work/pausch/alex/software/beagle.08Feb22.fa4.jar gt={input.panel} out={params.panel} ne=200 nthreads={threads}
        
        grep -f <(bcftools query -l {output.panel}) <(bcftools query -l {input.chip}) > {output.exclude}
        java -Xmx40g -jar /cluster/work/pausch/alex/software/beagle.08Feb22.fa4.jar ref={output.panel} gt={input.chip} out={params.chip} ne=200 nthreads={threads} excludesamples={output.exclude}

        '''

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



rule qtltools_parallel:
    input:
        ''
    output:
        'X'
    shell:
        '''
        parallel --ungroup --jobs {threads} "QTLtools cis --silent --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --permute 100 --chunk {{}} config[chunks] --out {outputpermutations_$j\_20.txt ::: {1..100}
        cat {output.permutations} > {output.nomimal}
        '''
#/cluster/work/pausch/alex/XENA/qtltools/bin/QTLtools gwas --bed aligned_genes.bed.gz --vcf imputed_chip_unique.vcf.gz --cov /cluster/work/pausch/xena/eQTL/covariates.txt --normal --out gwas_results.txt
#sed -e '/60436342/{r temp.g3' -e 'd}' DV.eQTL.merged.no_missing.vcf | bgzip -c >  DV.eQTL.merged.no_missing.vcf2.gz;  tabix -fp vcf DV.eQTL.merged.no_missing.vcf2.gz


#awk '{print "0\t"$2"\t"$3}' /cluster/work/pausch/xena/gwas/phenotypes/from_naveen/vzr/BV/sire/phenotypes.txt > vzr.phen
#grep -f <(bcftools query -l imputed_chip.vcf.gz.vcf.gz) <(bcftools query -l chips_1.vcf.gz) > exclude_samples.txt
#zcat impute_region.vcf.gz | perl -pe "s/\s\.:/\t.\/.:/g" | bgzip -c > impute_region2.vcf.gz
#java -jar /cluster/work/pausch/alex/software/beagle.08Feb22.fa4.jar gt=impute_region2.vcf.gz out=imputed_chip.vcf.gz ne=100


#bcftools view -S ^<(bcftools query -l imputed_chipX.vcf.gz | grep : | awk '{print $0}') -o imputed_chip_unique.vcf.gz imputed_chipX.vcf.gz

#zgrep 1_94218322_INDEL imputed_chip_unique.vcf.gz | cut -f 10- | sed -r -e "s/:\S*//g" | tr '\t' '\n' | sort | uniq -c
#"java -jar /cluster/work/pausch/alex/software/beagle.08Feb22.fa4.jar ref=imputed_chip.vcf.gz.vcf.gz gt=chips_1.vcf.gz out=imputed_chipX.vcf.gz ne=100 nthreads=4 excludesamples=exclude_samples.txt
