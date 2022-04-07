################################ Call variants from haplotype-resolved assemblies ################################
#
# steps:
#
#  1.) align contigs to reference using minimap2
#  2.) determine regions (uniquely) covered contigs haplotypes
#  3.) call variants from assemblies using paftools
#  4.) generate bi-allelic vcf file
#  5.) check mendelian consistency in trios and construct graph (multi-allelic vcf)
#  6.) compute some statistics
#
# output:
#
# variant callset produced from the assemblies (represented as bi-allelic VCF: callset-filtered.vcf)
# pangenome graph produced from the variant calls (represented as multi-allelic VCF: graph-filtered.vcf)
#
##################################################################################################################


configfile: "config.json"
samples = config['assemblies'].keys()
samples_parents = [s for s in samples if not s in config['trios']]
scripts = config['scripts']
outdir = config['outdir']
chromosomes = [config['reference']['prefix'] + str(i) for i in range(1,23)] + [config['reference']['prefix'] + 'X', config['reference']['prefix'] + 'Y']
frac_missing = 0.2 # skip positions with more than this fraction of missing alleles

# paftools skips contig-alignments shorter than this threshold
min_alignment_len = 50000

rule all:
	input:
		outdir + "multisample-vcfs/graph-filtered.vcf",
		outdir + "multisample-vcfs/callset-filtered.vcf",
		outdir + "statistics/vcftools-plots/indel-histogram.pdf",
		outdir + "statistics/vcfstats-stats.txt"



###########################################
#		1) Alignment
###########################################

rule align_assemblies:
    input:
        contigs = lambda wildcards: config['assemblies'][wildcards.sample][int(wildcards.haplotype)],
        reference = config['reference']
    output:
        temp(outdir + "paf/{sample}-hap{haplotype}.sam")
    threads: 8
    resources:
        mem_mb = 4000
    shell:
        '''
        minimap2 -ax asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 -t {threads} {input.reference} {input.contigs} > {output}
        '''

# align assemblies to reference genome
rule align_assemblies_paf:
	input:
	    outdir + "paf/{sample}-hap{haplotype}.sam"
    output:
		temp(outdir + "paf/{sample}-hap{haplotype}.paf")
	resources:
        mem_mb = 5000
    shell:
        '''
        paftools.js sam2paf {input} > {output}
        '''


###########################################
#		2) Callable regions
###########################################

# align assemblies to reference and produce BAM output
rule align_assemblies_bam:
	input:
	    outdir + "paf/{sample}-hap{haplotype}.sam"
    output:
		outdir + "bam/{sample}-hap{haplotype}.bam"
    threads: 4
	resources:
		mem_mb = 8000,
        disk_scratch = 10000
	shell:
		"""
        samtools view --threads {threads} -bS {input} | samtools sort --threads {threads} -T $TMPDIR -o {output} --write-index -
		"""


# compute regions covered by at least one contig
rule compute_covered_regions:
	input:
		outdir + "bam/{sample}-hap{haplotype}.bam"
	output:
		outdir + "bed/{sample}-hap{haplotype}_covered.bed"
	shell:
		"bedtools bamtobed -i {input} | awk '($3-$2) >= {min_alignment_len}' | bedtools merge > {output}"


# compute regions with per-base coverage < 2 
# NOTE: this will NOT remove cases in which there are more than one contig, but all except one contain a deletion)
# 	CCCCCCCCCCCC
# 	C----------C
rule compute_coverage:
	input:
		outdir + "bam/{sample}-hap{haplotype}.bam"
	output:
		outdir + "bed/{sample}-hap{haplotype}_unique.bed"
	resources:
		mem_total_mb=10000,
		runtime_hrs=3,
		runtime_min=1
	shell:
		"bedtools genomecov -bga -ibam {input} | awk '$4 < 2' | bedtools merge > {output}"


rule intersect_beds:
	input:
		covered= outdir + "bed/{sample}-hap{haplotype}_covered.bed",
		unique= outdir + "bed/{sample}-hap{haplotype}_unique.bed"
	output:
		outdir + "bed/{sample}-hap{haplotype}_callable.bed"
	resources:
		mem_total_mb=10000,
		runtime_hrs=3,
		runtime_min=1
	shell:
		"bedtools intersect -a {input.covered} -b {input.unique} > {output}"

rule sort_bed:
	input:
		"{filename}.bed"
	output:
		"{filename}-sorted.bed"
	shell:
		"bedtools sort -i {input} > {output}"


rule callable_regions:
	input:
		expand(outdir + "bed/{sample}-hap{haplotype}_callable-sorted.bed", sample=samples, haplotype=[0,1])
	output:
		outdir + "bed/callable-regions.bed"
	params:
		covered = int((1-frac_missing) * len(samples)*2)
	shell:
		"bedtools multiinter -i {input} | awk '$4 > {params.covered}' | bedtools merge > {output}"



##########################################
#		3) Variant Calling 
##########################################

# call variants from alignments
rule paftools:
	input:
		paf= outdir + "paf/{sample}-hap{haplotype}.paf",
		reference = config['reference']
	output:
		temp(outdir + "calls/{sample}-hap{haplotype, [0,1]}.vcf")
	resources:
		mem_total_mb=50000,
		runtime_hrs=3,
		runtime_min=1
	shell:
		"paftools.js call -L {min_alignment_len} -s {wildcards.sample}_{wildcards.haplotype} -f {input.reference} {input.paf} | sed 's|1/1|1|g' > {output}"


####################################################
#  4) Create bi-allelic VCF and filter it
####################################################


rule compress_vcf:
	input:
		"{filename}.vcf"
	output:
		gz="{filename}.vcf.gz",
		tbi="{filename}.vcf.gz.tbi"
	shell:
		"""
		bgzip -c {input} > {output.gz}
		tabix -p vcf {output.gz}
		"""

# create a multisample VCF containing all haplotypes
rule collect_all_haplotypes:
	input:
		vcfs=expand("{outdir}calls/{sample}-hap{haplotype}.vcf.gz", outdir=outdir, sample=samples, haplotype=[0,1]),
		tbi=expand("{outdir}calls/{sample}-hap{haplotype}.vcf.gz.tbi", outdir=outdir, sample=samples, haplotype=[0,1])
	output:
		outdir + "calls/all-haplotypes.vcf"
	resources:
		mem_total_mb=50000,
		runtime_hrs=3,
		runtime_min=1
	shell:
		"bcftools merge -m none --missing-to-ref {input.vcfs} | python3 {scripts}/assign-variant-ids.py > {output}"


# extract variant ids of callable regions
rule extract_covered_ids:
	input:
		bed = outdir + "bed/{sample}-hap{haplotype}_callable.bed",
		vcf = outdir + "calls/all-haplotypes.vcf"
	output:
		outdir + "bed/{sample}_{haplotype}.txt"
	shell:
		"bedtools intersect -a {input.vcf} -b {input.bed} -wa -f 1.0 | cut -f3 > {output}"


# set alleles outside of callable regions to missing
rule set_to_missing:
	input:
		vcf = outdir + "calls/all-haplotypes.vcf",
		bed = expand("{outdir}bed/{sample}_{haplotype}.txt", outdir=outdir, sample=samples, haplotype=[0,1])
	output:
		outdir + "calls/all-haplotypes-callable.vcf"
	log:
		outdir + "calls/all-haplotypes-callable.log"
	resources:
		mem_total_mb=150000,
		runtime_hrs=10,
		runtime_min=1
	shell:
		"python3 {scripts}/set-to-missing.py -v {input.vcf} -m {frac_missing} -f {input.bed} 2> {log} 1> {output}"


# convert haploid VCF into a diploid one by combining haplotypes of each sample
rule write_input:
	output:
		outdir + "samples.txt"
	run:
		with open(output[0], 'w') as txt_output:
			for sample in samples:
				txt_output.write('\t'.join([sample, sample + '_0', sample + '_1']) + '\n')


rule combine_haplotypes:
	input:
		haps = outdir + "calls/all-haplotypes-callable.vcf",
		samples = outdir + "samples.txt"
	output:
		temp(outdir + "multisample-vcfs/callset.vcf")
	resources:
		mem_total_mb=100000,
		runtime_hrs=10,
		runtime_min=1
	shell:
		'python3 {scripts}/merge_vcfs.py combine_columns -samples {input.samples} -vcf {input.haps} > {output}'


#################################################################
#  5) Check mendelian consistency for trios and construct graph
#################################################################


# generate a file specifying the trio relationships
rule generate_ped_file:
	output:
		"{outdir}trios.ped"
	run:
		with open(output[0], "w") as ped_output:
			for trio in config['trios']:
				father=config['trios'][trio][0]
				mother=config['trios'][trio][1]
				ped_output.write('\t'.join([trio, trio, father, mother]) + '\n')
				
rule generate_samples_file:
	output:
		"{outdir}trio-samples.txt"
	run:
		with open(output[0], "w") as sample_output:
			for trio in config['trios']:
				sample_output.write(trio + '\n')
				for sample in config['trios'][trio]:
					sample_output.write(sample + '\n')


# remove all variants where there is a mendelian conflict in at least one of the trios
# if no trios are given in config, the vcf does not change.
rule check_mendelian_consistency:
	input:
		vcf="{outdir}multisample-vcfs/callset.vcf",
		ped="{outdir}trios.ped",
		samples="{outdir}trio-samples.txt"
	output:
		vcf="{outdir}multisample-vcfs/callset-filtered.vcf",
		tsv="{outdir}multisample-vcfs/mendelian-consistency.tsv"
	log:
		"{outdir}multisample-vcfs/mendelian-consistency.log"
	resources:
		mem_total_mb=100000,
		runtime_hrs=10,
		runtime_min=1
	shell:
		"""
		python3 {scripts}/mendelian-consistency.py filter -vcf {input.vcf} -samples {input.samples} -ped {input.ped} -o {output.tsv} 2> {log} 1> {output.vcf}
		"""

rule merge_haplotypes:
	input:
		vcf = outdir + "multisample-vcfs/callset-filtered.vcf",
		reference = config['reference']
	output:
		tmp = temp(outdir + "multisample-vcfs/graph-filtered-tmp.vcf")
	params:
		chrom = ','.join([c for c in chromosomes])
	log:
		outdir + "multisample-vcfs/assemblies-all-samples-filtered.log"
	resources:
		mem_total_mb=900000,
		runtime_hrs=30,
		runtime_min=1
	shell:
		"""
		python3 {scripts}/merge_vcfs.py merge -vcf {input.vcf} -r {input.reference} -ploidy 2 -chromosomes {params.chrom} 2> {log} 1> {output.tmp}
		"""

#############################################
#  6) Generate statistics and create plots 
#############################################

rule normalize_vcf:
	input:
		vcf = outdir + "multisample-vcfs/graph-filtered-tmp.vcf.gz",
		reference = config['reference']
	output:
		outdir + "multisample-vcfs/graph-filtered.vcf"
	shell:
		"bcftools norm -d all -f {input.reference} {input.vcf} > {output}"


rule vcfstats_statistics:
	input:
		outdir + "multisample-vcfs/graph-filtered.vcf"
	output:
		txt= outdir + "statistics/vcfstats-stats.txt"
	shell:
		"rtg vcfstats {input} > {output}"


rule untypable_ids:
	input:
		outdir + "multisample-vcfs/callset-filtered.vcf.gz"
	output:
		lists=expand(outdir + "statistics/untypable-ids/{sample}-untypable.tsv", sample=samples_parents),
		summary= outdir + "statistics/untypable-ids.tsv"
	params:
		out= outdir + "statistics/untypable-ids"
	shell:
		"zcat {input} | python3 {scripts}/untypable-ids.py {params.out} > {output.summary}"


rule indel_histogram:
	input:
		outdir + "multisample-vcfs/graph-filtered.vcf"
	output:
		histo= outdir + "statistics/vcftools-stats.indel.hist",
		plot= outdir + "statistics/vcftools-plots/indel-histogram.pdf"
	shell:
		"""
		vcftools --vcf {input} --out {outdir}statistics/vcftools-stats --hist-indel-len
		cat {output.histo} | python3 {scripts}/plot-callset-statistics.py length {output.plot} 20000
		"""
