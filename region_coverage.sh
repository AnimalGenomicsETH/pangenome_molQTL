samtools bedcov <(echo -e "14\t1559308\t1605491") $(awk 'NR>205 {print "/cluster/work/pausch/inputs/bam/BTA_eQTL/"$2".bam"}' ../config/large_cohort.yaml)
for i in $(awk 'NR>205 {print $2}' ../config/large_cohort.yaml); do echo "\"${i}\":$(awk '/total/ {print $1*150/2700000000.0}' /cluster/work/pausch/inputs/bam/BTA_eQTL/${i}.stats)"; done |tr '\n' ','
bcftools query -f '{["%SAMPLE":"%GT",]\n' -r 14:1559308 variants.normed.vcf.gz
