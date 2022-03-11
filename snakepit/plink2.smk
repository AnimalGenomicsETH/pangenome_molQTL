#plink2 --pca --threads 2 --memory 10000 --chr-set 30 --vcf SVg50.vcf.gz --vcf-half-call h --bad-freqs --out SV
#grep -f <(awk 'NR>1{print $1}' full.eigenvec) /cluster/work/pausch/vcf_UCD/2020_09/stats/summary_coverage.txt | awk '{print $1,$3}' | sort -k1,1 | awk '{print $2}'|cat <(echo "breed") - | paste - full.eigenvec > full2.eigenvec


rule bcftools_annotate:
    input:
        config['pangenie_out']
    output:
        'pangenie/samples.annotated.vcf.gz'
    threads: 4
    resources:
        mem_mb = 1000,
        walltime = '24:00'
    shell:
        '''
        bcftools annotate --threads {threads} --set-id +'%CHROM\_%POS\_%TYPE' -o {output} {input}
        '''

rule plink_LD:
    input:
        vcf = 'pangenie/samples.annotated.vcf.gz'
    output:
        'pangenie/plink.ROIs.ld'
    params:
        mem = lambda wildcards, threads, resources: threads*resources.mem_mb,
        ROIs = ' '.join(config['ROIs']),
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 4
    resources:
        mem_mb = 1500,
        walltime = '30'
    envmodules:
        'gcc/8.2.0'
        'plink'
    shell:
        '''
        plink --vcf {input.vcf} --ld-snps {params.ROIs} --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --dprime --with-freqs --threads {threads} --memory {params.mem} --r2 --chr-set 30 --vcf-half-call h --out {params.out}
        '''


def make_region(ROI,window=100000):
    chr, pos = ROI.split('_')[:2]
    return f'{chr}:{pos-window}-{pos+window}'

rule plink_matrix:
    input:
        vcf = 'pangenie/samples.annotated.vcf.gz'
    output:
        'pangenie/plink.{ROI}.ld'
    threads: 2
    resources:
        mem_mb = 2000
    params:
        mem = lambda wildcards, threads, resources: threads*resources.mem_mb,
        out = lambda wildcards, output: PurePath(output[0]).with_suffix(''),
        region = lambda wildcards: make_region(wildcards.ROI)
    shell:
        '''
        plink --vcf <(bcftools --threads {threads} view {input.vcf} {params.region}) --r2 square --chr-set 30 --threads {threads} --memory {params.mem} --out {params.out}
        '''
#plink --vcf test.vcf --ld-snp 17_60436342 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out ex --threads 2 --memory 8000 --r2 --chr-set 30

rule bcftools_get_SV:
    input:
        'pangenie/samples.annotated.vcf.gz'
    output:
        'pangenie/samples.SV.{ILEN}.txt'
    resources:
        mem_mb = 4000
    shell:
        '''
        bcftools query -i 'abs(ILEN)>{wildcards.ILEN}' -f'%ID\n' -H {input} | tail -n +2 > {output}
        '''


rule plink_tag:
    input:
        vcf = 'pangenie/samples.annotated.vcf.gz',
        tags = 'pangenie/samples.SV.{ILEN}.txt'
    output:
        'pangenie/samples.SV.{ILEN}.tags.list'
    threads: 8
    resources:
        mem_mb = 3000
    params:
        mem = lambda wildcards, threads, resources: threads*resources.mem_mb,
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    shell:
        '''
        plink --vcf {input.vcf} --show-tags {input.tags} --tag-kb 1000 --threads {threads} --memory {params.mem} --chr-set 30 --vcf-half-call h --list-all --out {params.out}
        '''

rule match_CHIP_tags:
#zgrep -f <(awk '$3=="1_94218322_INDEL"&&$9>0.8&&$7~/SNP/ {split($7,a,"_"); print a[2]}' plink.ROIs.ld) -H  /cluster/work/pausch/naveen/ALLIMPUTE/hd/CHR1/ref_alt.vcf.gz |  awk 'BEGIN{print "0/0\t0/1\t1/1"}{print gsub(/0\/0/,"")"\t"gsub(/1\/0/,"")+gsub(/0\/1/,"")"\t"gsub(/1\/1/,"")}'
