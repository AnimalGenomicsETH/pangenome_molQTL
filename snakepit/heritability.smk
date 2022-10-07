
import os
from pathlib import PurePath

rule all:
    input:
        'GRM/main.hsq'

rule partition_IDs:
    input:
        lambda wildcards: config['vcfs'][wildcards.vcf]
    output:
        SV = 'bfiles/{vcf}.SV.ids',
        small = 'bfiles/{vcf}.small.ids'
    shell:
        '''
        bcftools query -i 'abs(ILEN)>=50' -f '%ID\n' {input} > {output.SV}
        bcftools query -e 'abs(ILEN)>=50' -f '%ID\n' {input} > {output.small}
        '''

rule plink_make_bed:
    input:
        vcf = lambda wildcards: config['vcfs'][wildcards.vcf],
        ids = 'bfiles/{vcf}.{variants}.ids'
    output:
        'bfiles/{vcf}.{chromosome}.{variants}.bim'
    params:
        bfile = lambda wildcards, output: PurePath(output[0]).with_suffix(''),
        memory = lambda wildcards, threads, resources: int(threads*resources.mem_mb)
    threads: 1
    resources:
        mem_mb = 20000,
        walltime = '15'
    shell:
        '''
        plink2 --native --vcf {input.vcf} --threads {threads} --extract {input.ids} --memory {params.memory} --maf 0.01:minor --cow --make-bed --chr {wildcards.chromosome} --out {params.bfile}
        '''

rule gcta_grm:
    input:
        expand('bfiles/{{vcf}}.{chromosome}.{{variants}}.bim',chromosome=range(1,30))
    output:
        'GRM/{vcf}.{variants}.grm.bin'
    params:
        _input = lambda wildcards, input: '\\n'.join([str(PurePath(I).with_suffix('')) for I in input]),
        _output = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    threads: 8
    resources:
        mem_mb = 6000
    shell:
        '''
        echo -e "{params._input}" > $TMPDIR/mbfile
        gcta --thread-num {threads} --autosome-num 30 --make-grm-bin --mbfile $TMPDIR/mbfile --out {params._output}
        '''

localrules: prep_covars
rule prep_covars:
    input:
        config['covariates']['eQTL']
    output:
        'GRM/covar.txt'
    shell:
        '''
        awk '{{for(i=1;i<=NF;i++)a[i][NR]=$i}} END {{for(i in a)for(j in a[i])printf"%s"(j==NR?RS:FS),a[i][j]}}' {input} | awk 'NR>1{{print "0",$0}}' > {output}
        '''

localrules: prep_phenotypes
checkpoint prep_phenotypes:
    input:
        config['genes']['eQTL']
    output:
        directory('GRM/molecular_phenotypes')
    shell:
        '''
        mkdir -p {output}
        zcat {input} | awk -v P={output} '{{if (NR==1){{for(i=7;i<=NF;i++)a[i]=$i}} else {{GENE=$4;for(i=7;i<=NF;i++)b[i]=$i;for(i in a)print 0,a[i],b[i] >> P"/"GENE".phen"}} }}'
        '''

def aggregate_phenotypes(wildcards):
      checkpoint_output = checkpoints.prep_phenotypes.get(**wildcards).output[0]
      return expand("GRM/{{vcf}}.{phenotype}.hsq", phenotype=glob_wildcards(os.path.join(checkpoint_output, "{phenotype}.phen")).phenotype)

localrules: gcta_reml
rule gcta_reml:
    input:
        grm = expand('GRM/{{vcf}}.{variants}.grm.bin',variants=('SV','small')),
        covar = rules.prep_covars.output,
        phenotype = 'GRM/molecular_phenotypes/{phenotype}.phen'
    output:
        'GRM/{vcf}.{phenotype}.hsq'
    params:
        grm = lambda wildcards, input: '\\n'.join([PurePath(I).with_suffix('').with_suffix('') for I in input.grm]),
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    shell:
        '''
        gcta --reml --mgrm-bin {params.grm} --pheno {input.phenotype} --qcovar {input.covar} --reml-no-lrt --reml-alg 2 --reml-maxit 10000 --out {params.out}
        '''
#find GRM/molecular_phenotypes/ -type f -name '*.phen' -exec parallel -I@@ -j 6 gcta --reml --grm-bin GRM/main.small --pheno @@ --qcovar GRM/covar.txt --reml-no-lrt --reml-alg 2 --reml-maxit 10000 --out {.}.small ::: {} \+ > /dev/null

localrules: gather_hsq
rule gather_hsq:
    input:
        aggregate_phenotypes
    output:
        'GRM/{vcf,main}.hsq'
    shell:
        '''
        awk '$1=="V(G)/Vp" {{print FILENAME,$2,$3}}' {input} | sed 's/\.hsq//g' > {output}
        '''
#awk '{if($1=="V(G1)/Vp"){a[FILENAME][1]=$2;a[FILENAME][2]=$3}else{if($1=="V(G2)/Vp"){a[FILENAME][3]=$2;a[FILENAME][4]=$3}else{if($3=="V(G)/Vp"){a[FILENAME][5]=$4;a[FILENAME][6]=$5}}}} END {for(F in a){print F,a[F][1],a[F][2],a[F][3],a[F][4],a[F][5],a[F][6]}}' GRM/molecular_phenotypes/*hsq | sed 's/GRM\/molecular_phenotypes\///g;s/\.hsq//g' > GRM/main.hsq


        
