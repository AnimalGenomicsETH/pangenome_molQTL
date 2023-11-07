
include: 'snakepit/pangenie_panel.smk'
include: 'snakepit/pangenie_genotyping.smk'
include: 'snakepit/association_mapping.smk'

rule all:
    input:
        'pangenie_panel/multisample-vcfs/graph-filtered.vcf',
        'pangenie_panel/samples.all.pangenie_genotyping_DV.vcf.gz',
        expand('QTL/{QTL}/Testis_{variants}/conditionals.01.txt.gz',QTL=('eQTL','sQTL'),variants=config['variants'])
