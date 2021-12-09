rule pbmm2_align:
    input:
        reads = lambda wildcards: config['hifi_samples'][wildcards.sample]
    output:
        temp(get_dir('pbsv','{sample}.ARS.pbmm2.bam'))
    threads: 24
    resources:
        mem_mb = 3000,
        walltime = '4:00'
    shell:
        'pbmm2 align {config[reference]} {input.reads} {output} --sort --preset CCS -j {threads} --sample {wildcards.sample}'

rule pbsv_discover:
    input:
        get_dir('pbsv','{sample}.ARS.pbmm2.bam')
    output:
        get_dir('pbsv','{sample}.ARS.svsig.gz')
    resources:
        mem_mb = 10000
    shell:
        'pbsv discover {input} {output}'

rule pbsv_call:
    input:
        (get_dir('pbsv','{sample}.ARS.svsig.gz',sample=S) for S in config['hifi_samples'])
    output:
        get_dir('pbsv','samples.pbsv.vcf')
    threads: 8
    resources:
        mem_mb = 3000,
        walltime = "24:00"
    shell:
        'pbsv call --ccs -j {threads} {config[reference]} {input} {output}'
