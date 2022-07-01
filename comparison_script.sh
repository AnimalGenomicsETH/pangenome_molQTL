awk 'NR>1&&$21=="1" {print $9"\t"$10"\t"$11"\t"$1"\t"$8}' /cluster/work/pausch/xena/eQTL/cis_all/testis/normal/maf01/conditional/all.txt | sort -k1,1n -k2,2n > eQTL/xena_S.txt
awk -v L=499 '{split($5,a,"_"); if(length(a[4])>L||length(a[5])>L){print $4}}' alex_S.txt >SV_hits.txt
grep -f SV_hits.txt xena_S.txt | cut -f 4 | sort | uniq | wc -l
