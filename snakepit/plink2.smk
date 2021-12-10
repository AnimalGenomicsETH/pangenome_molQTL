plink2 --pca --threads 2 --memory 10000 --chr-set 30 --vcf SVg50.vcf.gz --vcf-half-call h --bad-freqs --out SV
grep -f <(awk 'NR>1{print $1}' full.eigenvec) /cluster/work/pausch/vcf_UCD/2020_09/stats/summary_coverage.txt | awk '{print $1,$3}' | sort -k1,1 | awk '{print $2}'|cat <(echo "breed") - | paste - full.eigenvec > full2.eigenvec
