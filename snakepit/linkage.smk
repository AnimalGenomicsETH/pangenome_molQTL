from pathlib import PurePath

rule all:
    input:
        expand('LD/{vcf}.tags.csv',vcf=config['vcfs'])

wildcard_constraints:
    chromosome = r'\d+',
    target = r'SV|SNP|random'

localrules: make_targets
rule make_targets:
    input:
        hits = lambda wildcards: config['eQTL'][wildcards.vcf],
        targets = lambda wildcards: [] if wildcards.target == 'SV' else 'LD/{vcf}.SV_targets.txt'
    output:
        'LD/{vcf}.{target,SNP|SV}_targets.txt'
    params:
        size_condition = lambda wildcards: 'length(a[4])>=50||length(a[5])>=50' if wildcards.target == 'SV' else 'length(a[4])<50||length(a[5])<50',
        n_samples = lambda wildcards, input: 9999999 if wildcards.target == 'SV' else f"$(wc -l {input.targets} | awk '{{print $1}}')"
    shell:
        '''
        awk '$22=="1" {{split($8,a,"_"); if({params.size_condition}) {{print $8}} }}' {input} |\
        sort -u |\
        shuf -n {params.n_samples} > {output}
        '''

localrules: random_targets
rule random_targets:
    input:
        vcf = lambda wildcards: config['vcfs'][wildcards.vcf],
        targets = 'LD/{vcf}.SV_targets.txt'
    output: 
        'LD/{vcf}.random_targets.txt' 
    shell:
        '''
        bcftools query -e 'MAF<0.01||abs(ILEN)>=50' -f '%ID\\n' {input.vcf} |\
        sort -u |\
        shuf -n $(wc -l {input.targets} |\
        awk '{{print $1}}') > {output}
        '''

rule plink_make_bed:
    input:
        lambda wildcards: config['vcfs'][wildcards.vcf]
    output:
        'LD/{vcf}.{chromosome}.bed'
    params:
        bfile = lambda wildcards, output: PurePath(output[0]).with_suffix(''),
        memory = lambda wildcards, threads, resources: int(threads*resources.mem_mb)
    threads: 1
    resources:
        mem_mb = 20000,
        walltime = '15'
    shell:
        '''
        plink2 --native --vcf {input} --threads {threads} --memory {params.memory} --maf 0.01:minor --cow --make-bed --chr {wildcards.chromosome} --out {params.bfile}
        '''

rule plink_LD:
    input:
        bfile = rules.plink_make_bed.output[0],
        targets = 'LD/{vcf}.{target}_targets.txt'
    output:
        'LD/{vcf}.{chromosome}.{target}.{r2}.tags.list'
    params:
        bfile = lambda wildcards, input: PurePath(input.bfile).with_suffix(''),
        tags = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        r2 = lambda wildcards: int(wildcards.r2)/100
    threads: 1
    resources:
        mem_mb = lambda wildcards: 40000 if int(wildcards.chromosome) == 1 else 25000,
        walltime = '15'
    shell:
        '''
        plink --bfile {params.bfile} --silent --show-tags {input.targets} --list-all --tag-r2 {params.r2} --tag-kb 1000 --cow --threads {threads} --out {params.tags}
        '''

localrules: gather_LD
rule gather_LD:
    input:
        expand('LD/{{vcf}}.{chromosome}.{{target}}.{{r2}}.tags.list',chromosome=range(1,30))
    output:
        'LD/{vcf}.{target}.{r2}.tags.list'
    shell:
        '''
        awk '{{$1=$1}};1' {input} | grep -v "NTAG" | cut -d' ' -f 1,4 | sort -k1,1V | awk -v C={wildcards.target} -v R={wildcards.r2} '{{print $0,C,R}}' > {output}
        '''

localrules: gather_all
rule gather_all:
    input:
        expand('LD/{vcf}.{target}.{r2}.tags.list',vcf=config['vcfs'],target=('SV','SNP','random'),r2=(70,80,90))
    output:
        'LD/{vcf}.tags.csv'
    shell:
        '''
        echo "ID NTAG SOURCE R2" > {output}
        cat {input} >> {output}
        '''
