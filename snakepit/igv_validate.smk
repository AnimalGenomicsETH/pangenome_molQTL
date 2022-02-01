
rule slice_ARS:
    input:
        bam = '/cluster/work/pausch/inputs/bam/BTA_UCD12/{sample}.bam',
        asm = '/cluster/work/pausch/alex/assembly/BSWCHEM120151536851/hifiasmv152_100/hap2.scaffolds.fasta'
    output:
        bam = get_dir('igv','{sample}.TBX3.ARS.bam'),
        asm = get_dir('igv','{sample}.TBX3.BSW.bam')
    threads: 4
    resources:
        mem_mb = 5000,
        disk_scratch = 10,
        walltime = '1:00'
    shell:
        '''
        samtools view -@ {threads} -o {output.bam} {input.bam} 17:60000000-61000000
        samtools fastq -@ {threads} {output.bam} | paste - - - -  | sort -k1,1 -S 3G  | tr '\t' '\n' > $TMPDIR/reads.fq
        minimap2 -t {threads} --cs -axsr {input.asm} <(seqtk seq -1 $TMPDIR/reads.fq) <(seqtk seq -2 $TMPDIR/reads.fq) | samtools sort - -m 1000M -@ {threads} -o {output.asm}
        samtools index -@ {threads} {output.bam}
        samtools index -@ {threads} {output.asm}
        '''
