
import os
from pathlib import PurePath

rule all:
    input:
        'GRM/main.hsq'

rule gcta_grm:
    input:
        expand('bfiles/{{vcf}}.{chromosome}.bim',chromosome=range(1,30))
    output:
        'GRM/{vcf}.grm.bin'
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
        grm = rules.gcta_grm.output[0],
        covar = rules.prep_covars.output,
        phenotype = 'GRM/molecular_phenotypes/{phenotype}.phen'
    output:
        'GRM/{vcf}.{phenotype}.hsq'
    params:
        grm = lambda wildcards, input: PurePath(input[0]).with_suffix('').with_suffix(''),
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    shell:
        '''
        gcta --reml --grm-bin {params.grm} --pheno {input.phenotype} --qcovar {input.covar} --reml-no-lrt --out {params.out}
        '''

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
