#plink2 --pca --threads 2 --memory 10000 --chr-set 30 --vcf SVg50.vcf.gz --vcf-half-call h --bad-freqs --out SV
#grep -f <(awk 'NR>1{print $1}' full.eigenvec) /cluster/work/pausch/vcf_UCD/2020_09/stats/summary_coverage.txt | awk '{print $1,$3}' | sort -k1,1 | awk '{print $2}'|cat <(echo "breed") - | paste - full.eigenvec > full2.eigenvec


rule bcftools_annotate:
    input:
        config['panel']
    output:
        'panels/panel.ANNOT.vcf.gz'
    threads: 4
    resources:
        mem_mb = 1000,
        walltime = '30'
    shell:
        '''
        bcftools annotate --threads {threads} --set-id +'%CHROM\_%POS\_%TYPE' -o {output} {input}
        '''

rule plink_LD:
    input:
        vcf = 'panels/panel.ANNOT.vcf.gz'
    output:
        'plink.ROIs.ld'
    params:
        mem = lambda wildcards, threads, resources: threads*resources.mem_mb,
        ROIs = ' '.join(config['ROIs']),
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 6
    resources:
        mem_mb = 5000,
        walltime = '30'
    shell:
        '''
        plink --vcf {input.vcf} --ld-snps {params.ROIs} --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --threads {threads} --memory {params.mem} --r2 --chr-set 30 --vcf-half-call h --out {params.out}
        '''
#plink --vcf test.vcf --ld-snp 17_60436342 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out ex --threads 2 --memory 8000 --r2 --chr-set 30
