
include: 'snakepit/pangenie_panel.smk'
include: 'snakepit/pangenie_genotyping.smk'
include: 'snakepit/variant_calling.smk'
include: 'snakepit/association_mapping.smk'
include: 'snakepit/LD.smk'

rule all:
    input:
        'pangenie_panel/multisample-vcfs/graph-filtered.vcf',
        'pangenie_panel/samples.all.pangenie_genotyping_DV.vcf.gz',
        'jasmine.vcf',
        'SVs/sizes.gz',
        expand('LD/samples.SV.{r2}.{window}.tags.list',r2=list(range(2,10))+[99],window=(10,100,1000)),
        expand('QTL/{QTL}/Testis_{variants}/conditionals.01.txt.gz',QTL=('eQTL','sQTL'),variants=('PanGenie','SR'))
