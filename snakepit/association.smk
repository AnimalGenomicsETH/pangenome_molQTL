from pathlib import PurePath

wildcard_constraints:
    _pass = r'permutations|conditionals|nominals',
    chunk = r'\d+',
    chrom = r'\d+',
    MAF = r'\d+',
    vcf = r'(eQTL|gwas)/\S+'


rule all:
    input:
        expand('QTL/{QTL}/Testis_{variants}/conditionals.01.txt',QTL=('eQTL','sQTL'),variants=config['variants'])

localrules: concat_genes
rule concat_genes:
    input:
        lambda wildcards: expand('/cluster/work/pausch/xena/eQTL/gene_counts/{tissue_code}/QTLtools/{chromosome}_gene_counts.gz',chromosome=range(1,30),tissue_code={'Testis':'testis','Epididymis_head':'epi_h','Vas_deferens':'vas_d'}[wildcards.tissue]) if wildcards.QTL=='eQTL' else expand('/cluster/work/pausch/xena/eQTL/sQTL/{tissue_code}/removed/phenotypes/leafcutter.clus.{chromosome}.bed.gz',chromosome=range(1,30),tissue_code={'Testis':'testis','Epididymis_head':'epi_h','Vas_deferens':'vas_d'}[wildcards.tissue])
    output:
        'aligned_genes/{QTL}.{tissue}.bed.gz',
        'aligned_genes/{QTL}.{tissue}.bed.gz.tbi' 
    shell:
        '''
        zcat {input} | sort -u -k1,1n -k2,2n | bgzip -@ 2 -c > {output[0]}
        tabix -p bed {output[0]}
        '''

rule normalise_vcf:
    input:
        lambda wildcards: config['variants'][wildcards.variants]
    output:
        'QTL/{variants}.vcf.gz'
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
        rules.normalise_vcf.output
    output:
        'QTL/{variants}.exclude_sites.{MAF}.txt'
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
        return f'--nominal {config.get("nominal",0.05)}'

rule qtltools_parallel:
    input:
        vcf = rules.normalise_vcf.output,
        exclude = rules.exclude_MAF.output,
        bed = rules.concat_genes.output,
        cov = lambda wildcards: config['covariates']['eQTL'][wildcards.tissue],
        mapping = lambda wildcards: 'QTL/{QTL}/{tissue}_{variants}/permutations_all.{MAF}.thresholds.txt' if wildcards._pass == 'conditionals' else []
    output:
        merged = temp('QTL/{QTL}/{tissue}_{variants}/{_pass}.{chunk}.{MAF}.txt')
    params:
        _pass = lambda wildcards,input: get_pass(wildcards._pass,input),#f'--permute {config["permutations"]}' if wildcards._pass == 'permutations' else f'--mapping {input.mapping}',
        debug = '--silent' if 'debug' in config else '',
        grp = lambda wildcards: '--grp-best' if wildcards.QTL == 'sQTL' else ''
    threads: 1
    resources:
        mem_mb = 12500,
        walltime = lambda wildcards: '24h' if wildcards._pass == 'permutationsX' else '4h'
    shell:
        '''
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} {params._pass} {params.grp} --window {config[window]} --normal --chunk {wildcards.chunk} {config[chunks]} --out {output} {params.debug}
        '''

localrules: qtltools_gather

rule qtltools_gather:
    input:
        expand(rules.qtltools_parallel.output,chunk=range(0,config['chunks']+1),allow_missing=True)
    output:
        'QTL/{QTL}/{tissue}_{variants}/{_pass}.{MAF}.txt'
    resources:
        mem_mb = 3000,
        walltime = '20'
    params:
        sort_key = lambda wildcards: '-k9,9n -k10,10n' if wildcards.QTL == 'eQTL' else '-k11,11n -k12,12n'
    shell:
        '''
        sort {params.sort_key} {input} > {output}
        '''

rule qtltools_FDR:
    input:
        expand(rules.qtltools_gather.output,_pass='permutations',allow_missing=True)
    output:
        'QTL/{QTL}/{tissue}_{variants}/permutations_all.{MAF}.thresholds.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    envmodules:
        'gcc/8.2.0',
        'r/4.2.2'
    shell:
        '''
        Rscript /cluster/work/pausch/alex/software/qtltools/scripts/qtltools_runFDR_cis.R {input} 0.05 {params.out}
        '''
