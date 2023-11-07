
include: 'snakepit/pangenie_panel.smk'
include: 'snakepit/pangenie_genotyping.smk'
include: 'snakepit/variant_calling.smk'
include: 'snakepit/variant_accuracy.smk'
include: 'snakepit/association_mapping.smk'
include: 'snakepit/LD.smk'

rule all:
    input:
        ## Pangenome panel creation and genotyping
        'pangenie_panel/multisample-vcfs/graph-filtered.vcf',
        'pangenie_panel/samples.all.pangenie_genotyping_DV.vcf.gz',
        ## Variant comparison between pangenome panel, short reads, and HiFi reads
        'jasmine.vcf',
        expand('SVs/{metric}.gz',metric=('sizes','support','GTs')),
        expand('SNPs/{metric}.csv',metric=('F1','isec')),
        ## Linkage disequilibrium analysis of SVs
        expand('LD/samples.SV.{r2}.{window}.tags.list',r2=list(range(2,10))+[99],window=(10,100,1000)),
        ## Molecular QTL mapping
        expand('QTL/{QTL}/Testis_{variants}/conditionals.01.txt.gz',QTL=('eQTL','sQTL'),variants=('PanGenie','SR'))
