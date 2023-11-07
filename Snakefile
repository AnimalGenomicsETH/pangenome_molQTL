

include: 'snakepit/pangenie_panel.smk'
include: 'snakepit/pangenie_genotyping.smk'
include: 'snakepit/association_mapping.smk'

rule all:
    input:
        'panel',
        'variants',
        'SR_eQTL'
        expand('QTL/{QTL}/Testis_{variants}/conditionals.01.txt.gz',QTL=('eQTL','sQTL'),variants=config['variants'])
