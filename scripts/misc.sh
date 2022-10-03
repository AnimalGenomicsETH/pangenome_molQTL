
# find mendelian SV violations
bcftools view  -s OxO,OSire,ODam -U -c 1 graph-filtered.vcf |\
bcftools norm -f /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa -m -any |\
bcftools +setGT -- -t a -n u |\
bcftools +mendelian -t ODam,OSire,OxO -m x |\
bcftools view -i 'abs(ILEN)>=50' > inconsistent_sites.SV.vcf

# find BSW private SVs and AF
bcftools view  -s BSW3,BSW4,BxP -U -c 1 -a graph-filtered.vcf |\
bcftools +setGT -- -t a -n u |\
bcftools +setGT -- -t . -n 0 |\
bcftools norm -m -any |\
bcftools view -a |\
bcftools view -c 1 -m 1 -i 'abs(ILEN)>=50' |\
bcftools +fill-tags - -- -t AF |\
grep -oP "AF=(0\.)?\d+" | sort | uniq -c
