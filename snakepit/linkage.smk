from pathlib import PurePath

rule all:
    input:
        expand('LD/{vcf}.{target}.tags.list',vcf=config['vcfs'],target=('SV','random'))

rule make_targets:
    input:
        lambda wildcards: config['eQTL'][wildcards.vcf]
    output:
        'LD/{vcf}.SV_targets.txt'
    shell:
        '''
        awk '$22=="1" {{split($8,a,"_"); if(length(a[4])>=50||length(a[5])>=50) {{print $8}} }}' {input} |\
        sort -u > {output}
        '''

rule plink_make_bed:
    input:
        lambda wildcards: config['vcfs'][wildcards.vcf]
    output:
        'bfiles/{vcf}.{chromosome}.bed'
    params:
        bfile = lambda wildcards, output: PurePath(output[0]).with_suffix(''),
        memory = lambda wildcards, threads, resources: int(threads*resources.mem_mb)
    threads: 2
    resources:
        mem_mb = 10000
    shell:
        '''
        plink2 --native --vcf {input} --threads {threads} --memory {params.memory} --maf 0.05 --cow --make-bed --chr {wildcards.chromosome} --out {params.bfile}
        '''

rule plink_LD:
    input:
        bfile = 'bfiles/{vcf}.{chromosome}.bed',
        targets = 'LD/{vcf}.{target}_targets.txt'
    output:
        'LD/{vcf}.{chromosome}.{target}.tags.list'
    params:
        bfile = lambda wildcards, input: PurePath(input.bfile).with_suffix(''),
        tags = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    threads: 2
    resources:
        mem_mb = 15000
    shell:
        '''
        plink --bfile {params.bfile} --silent --show-tags {input.targets} --list-all --tag-r2 0.7 --tag-kb 1000 --cow --threads {threads} --out {params.tags}
        '''

localrules: gather_LD
rule gather_LD:
    input:
        expand('LD/{{vcf}}.{chromosome}.{{target}}.tags.list',chromosome=range(1,30))
    output:
        'LD/{vcf}.{target}.tags.list'
    shell:
        '''
        echo "ID NTAG SOURCE" > {output}
        awk '{{$1=$1}};1' {input} | grep -v "NTAG" | cut -d' ' -f 1,4 | sort -k1,1V >> {output}
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
        bcftools query -e 'abs(ILEN)>=50' -f '%ID\\n' {input.vcf} | shuf -n $(wc -l {input.targets} | awk '{{print $1}}') > {output}
        '''
