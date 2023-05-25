

# get LD variants for an ID
#awk -v M=$(LC_ALL=C; grep $1 /cluster/work/pausch/alex/eQTL_GWAS/Variant_accuracy/LD/ANNOTATION_MAP.txt | awk '{print $1}') '$1==M {print $8}' /cluster/work/pausch/alex/eQTL_GWAS/Variant_accuracy/LD/samples.SV.tags.list | sed 's/|/ /g'


for i in $(awk -v M=$(LC_ALL=C; grep $1 /cluster/work/pausch/alex/eQTL_GWAS/Variant_accuracy/LD/ANNOTATION_MAP.txt | awk '{print $1}') '$1==M {print $8}' /cluster/work/pausch/alex/eQTL_GWAS/Variant_accuracy/LD/samples.SV.9.1000.tags.list | sed 's/|/ /g'); do mawk -v M=$i '$1==M {print $2}' /cluster/work/pausch/alex/eQTL_GWAS/Variant_accuracy/LD/ANNOTATION_MAP.txt; done




plink2 --vcf LD/variants.vcf.gz --exclude LD/SV_IDs.txt --allow-extra-chr --cow --threads 4 --pca --out SNPs
plink2 --vcf LD/variants.vcf.gz --extract LD/SV_IDs.txt --allow-extra-chr --cow --threads 4 --pca --out SVs


plink --vcf LD/variants.vcf.gz --allow-extra-chr --cow --threads 4 --ld-snp 19385103 --r2 gz --ld-window 999999999 --ld-window-kb 1000 --ld-window-r2 0 --out STN1


{ echo "sample TPM GT" ; paste <(grep -P "(start|ENSBTAG00000001889)" /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/testis/TPM_matrices/Testis_TPM_normalized_FINAL.tsv | cut -f 7-|  awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{for(i in a)for(j in a[i])printf"%s"(j==NR?RS:FS),a[i][j]}') <(bcftools view -S TPM_sample_order.txt  PanGenie.01_minor.vcf.gz 22:39001060 | bcftools query -f '[%GT\n]') ; } 



bcftools query -i 'abs(ILEN)>=1000' -f '%ID\n' PanGenie.01_minor.vcf.gz | cut -d'_' -f -2 | uniq | sort > big_SV_ID_location.counts
comm -12 SV_ID SNP_ID > SNP_SV_overlaps
