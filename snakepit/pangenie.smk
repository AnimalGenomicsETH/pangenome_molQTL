
rule pangenie:
    input:
        reference = config['reference'],
        vcf = config['panel'],
        fastq = lambda wildcards: config['samples'][wildcards.sample]
    output:
        get_dir('VG','{sample}.all.pangenie_genotyping.vcf')
    params:
        lambda wildcards, output: str(PurePath(output[0]).with_suffix('')).replace('_genotyping','')
    threads: 8
    resources:
        mem_mb = 15000
    shell:
        '''
        PanGenie -i {input.fastq} -r {input.reference} -v {input.vcf} -t {threads} -j {threads} -s {wildcards.sample} -p -u -o {params}
        '''
