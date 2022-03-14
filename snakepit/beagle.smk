grep -f <(bcftools query -l imputed_chip.vcf.gz.vcf.gz) <(bcftools query -l chips_1.vcf.gz) > exclude_samples.txt
zcat impute_region.vcf.gz | perl -pe "s/\s\.:/\t.\/.:/g" | bgzip -c > impute_region2.vcf.gz
java -jar /cluster/work/pausch/alex/software/beagle.08Feb22.fa4.jar gt=impute_region2.vcf.gz out=imputed_chip.vcf.gz ne=100
"java -jar /cluster/work/pausch/alex/software/beagle.08Feb22.fa4.jar ref=imputed_chip.vcf.gz.vcf.gz gt=chips_1.vcf.gz out=imputed_chipX.vcf.gz ne=100 nthreads=4 excludesamples=exclude_samples.txt
