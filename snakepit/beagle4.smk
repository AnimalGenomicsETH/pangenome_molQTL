

reference_file="/cluster/work/pausch/inputs/ref/BTA/UCD1.2/ARS-UCD1.2_Btau5.0.1Y.fa"
OUT_DIR="/cluster/work/pausch/alex/BSW_OBV_SVs/pangenie_imputed/"
IN_DIR="/cluster/work/pausch/alex/BSW_OBV_SVs/pangenie/"


#tools
BEAGLE='/cluster/work/pausch/alex/software/beagle.27Jan18.7e1.jar'


##parameters for split_vcf
#nvar=50_000
#overlap=5_000


#nvar=200_000
#overlap=5_000
nvar=60_000
overlap=5_000



##paramterrs for Beagle
ne=200
nt=6
mem="20g"
window= int(nvar*1.5 + overlap) #beagles default window size is 50k variants [this is the biggest possible window -- by split.vcf]

##chromosomes=list(range (1,30)) + ["X", "Y"]
chromosomes=list (range (1,30))
#chromosomes = [el for el in chromosomes if el !=26]
#chromosomes=["Y"]


rule all:				      
    input:
        expand (OUT_DIR + "CHR{chr}/PHASED/clean_beagle4.1.vcf.gz", chr=chromosomes),
        expand (OUT_DIR + "CHR{chr}/PHASED/clean_beagle4.1.vcf.gz.tbi", chr=chromosomes)

        
checkpoint split_vcf:
    input:
        OUT_DIR + "{chr}.raw.vcf.gz"
    output:
        directory(OUT_DIR + "CHR{chr}/SPLIT/")
    params:
        prefix=OUT_DIR + "CHR{chr}/SPLIT/part",
        nvar=nvar,
        overlap=overlap
    script:
        "split_vcf.py"

rule beagle:
    input:
        vcf=OUT_DIR + "CHR{chr}/SPLIT/part{part}.vcf"
    output:
        vcf=OUT_DIR + "CHR{chr}/PHASED/part{part}.vcf.gz"
    threads: nt
    resources:
        mem_mb = 3000
    params:
        prefix=OUT_DIR + "CHR{chr}/PHASED/part{part}",
	    ne=ne,
        nt=nt,
        mem=mem,
        window=window
    shell:
        "java -Xmx{params.mem} -Xss5m -jar " + BEAGLE +
        " window={params.window} "
        "gl={input.vcf} "
        "ne={params.ne} "
        "gprobs=true "
        "nthreads={threads} "
        "out={params.prefix}"

def list_files (wildcards):
    checkpoint_dir= checkpoints.split_vcf.get (chr=wildcards.chr).output [0]
    return (expand (OUT_DIR + "CHR{chr}/PHASED/part{part}.vcf.gz",
                    chr=wildcards.chr,
                    part=glob_wildcards(os.path.join(checkpoint_dir, "part{mypart}.vcf")).mypart))

rule merge_vcfs:
    input:
        myfiles=list_files
    output:
        OUT_DIR + "CHR{chr}/PHASED/clean_beagle4.1.vcf"
    params:
        prefix=OUT_DIR + "CHR{chr}/PHASED/part",
        nvar=nvar,
        overlap=overlap
    script:
        "merge_vcf.py"

rule zip:
    input:
        OUT_DIR + "CHR{chr}/PHASED/clean_beagle4.1.vcf"
    output:
        OUT_DIR + "CHR{chr}/PHASED/clean_beagle4.1.vcf.gz",
    threads: 6
    resources:
        mem_mb = 3000
    shell:
        'bgzip --threads {threads} {input}'

rule index:
    input:
        OUT_DIR + "CHR{chr}/PHASED/clean_beagle4.1.vcf.gz",
    output:
        OUT_DIR + "CHR{chr}/PHASED/clean_beagle4.1.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"
        


