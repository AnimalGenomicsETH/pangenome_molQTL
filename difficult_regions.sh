bcftools isec -n +1 -c some -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 HiFi_DV/all.DeepVariantWGS.vcf.gz SR_DV/all.DeepVariantWGS.vcf.gz > overlap.txt

#grep "centr" ../REF_DATA/ARS.repeats.sorted.bed > genome_regions.Centromeres.bed
#for F in "genome_regions.Centromeres.bed" "/cluster/work/pausch/alex/REF_DATA/genome_regions.Tandem.bed"
for F in "Satellite" "Tandem" "Low" "Normal" "Repetitive"
do
for i in "01" "10" "11"
do
  #awk -v OFS='\t' -v C=$i '$5==C {print $1,$2,$2+1}' overlap.txt | bedtools intersect -a $F -b - | wc -l
  awk -v OFS='\t' -v C=$i '$5==C {print $1,$2,$2+1}' overlap.txt | bedtools intersect -a <(zgrep $F /cluster/work/pausch/alex/REF_DATA/genome_regions.bed.gz) -b - | wc -l
done
done
