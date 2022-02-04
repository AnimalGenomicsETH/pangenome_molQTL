
rule slice_ARS:
    input:
        bam = '/cluster/work/pausch/inputs/bam/BTA_UCD12/{sample}.bam',
        asm = '/cluster/work/pausch/alex/assembly/BSWCHEM120151536851/hifiasmv152_100/hap2.scaffolds.fasta'
    output:
        bam = get_dir('igv','{sample}.TBX3.ARS.bam'),
        asm = get_dir('igv','{sample}.TBX3X.BSW.bam')
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

rule slice_BSW:
    input:
        fastq = expand('/cluster/work/pausch/inputs/fastq/BTA/{{sample}}_R{N}.fastq.gz',N=(1,2)),
        asm = '/cluster/work/pausch/alex/assembly/BSWCHEM120151536851/hifiasmv152_100/hap2.scaffolds.fasta'
    output:
        full = multiext(get_dir('igv','{sample}.BSW.bam'),'','.bai'),
        region = multiext(get_dir('igv','{sample}.TBX3.BSW.bam'),'','.bai')
    threads: 24
    resources:
        mem_mb = 5000,
        disk_scratch=100
    shell:
        '''
            minimap2 -t {threads} --cs -axsr {input.asm} {input.fastq} | samtools sort - -m 1000M -@ {threads} -T $TMPDIR -o {output.full[0]}
            samtools index -@ {threads} {output.full[0]}
            samtools view -@ {threads} -o {output.region[0]} {output.full[0]} 17:68200000-68700000
            samtools index -@ {threads} {output.region[0]}
        '''

