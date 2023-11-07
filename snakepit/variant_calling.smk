rule minimap2_align:
    input:
        lambda wildcards: config['HiFi_samples'][wildcards.sample]
    output:
        bam = temp('alignments/{sample}.HiFi.bam'),
        csi = temp('alignments/{sample}.HiFi.bam.csi')
    threads: 12
    resources:
        mem_mb = 6000,
        disk_scratch = 100,
        walltime = '24:00'
    shell:
        '''
        minimap2 -axmap-hifi -t {threads} {config[reference]} {input} |  samtools sort - -m 3000M -@ 4 -T $TMPDIR --write-index -o {output.bam}
        '''

rule sniffles_call:
    input:
        bam = 'alignments/{sample}.HiFi.bam'
    output:
        vcf = temp('variant_calling/{sample}.sniffles.vcf.gz'),
        snf = temp('variant_calling/{sample}.sniffles.snf')
    threads: 4
    resources:
        mem_mb = 2500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.bam} --reference {config[reference]} --sample-id {wildcards.sample} --threads {threads} --vcf {output.vcf} --snf {output.snf}
        '''

rule sniffles_merge:
    input:
        snfs = expand('variant_calling/{sample}.sniffles.snf',sample=config['HiFi_samples'])
    output:
        vcf = 'variant_calling/samples.sniffles.vcf.gz'
    threads: 2
    resources:
        mem_mb = 2500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.snfs} --reference {config[reference]} --threads {threads} --vcf {output.vcf}
        '''

rule bcftools_autosomes:
    input:
        'variant_calling/samples.sniffles.vcf.gz'
    output:
        'variant_calling/samples.sniffles.autosomes.vcf'
    threads: 2
    resources:
        mem_mb = 4000
    params:
        regions = ','.join(map(str,range(1,30)))
    shell:
        '''
        bcftools view --threads {threads} -r {params.regions} -o {output} {input}
        '''

rule bcftools_stuff:
    input:
        rules.bcftools_autosomes.output
    output:
        sizes = 'SVs/sizes.gz',
        support = 'SVs/support.gz',
        GTs = 'SVs/GTs.gz'
    localrule: True
    shell:
        '''
        bcftools query -e 'INFO/SVTYPE=="BND"' -f '%INFO/SVLEN\\n' {input} | pigz -p2 -c > {output.sizes}
        bcftools query -e 'INFO/SVTYPE=="BND"' -f '%INFO/SUPP_VEC\\n' {input} | mawk ' {{ print gsub(1,2,$1) }} ' | pigz -p2 -c > {output.support}
        bcftools query -e 'INFO/SVTYPE=="BND"' -f '[%GT ]\\n' {input} | sed 's"|" "g' | sed 's"/" "g' | sed 's/\./nan/g' | pigz -p2 -c > {output.GTs}
        '''

rule bcftools_split_panel:
    input:
        expand(rules.merge_pangenie.output,pangenie_mode='genotyping',allow_missing=True)
    output:
        SV = 'variant_calling/panel.SV.vcf',
        small = 'variant_calling/panel.small.vcf.gz'
    params:
        SV_size = 50,#config['SV_size'],
        bcf = '$TMPDIR/normed.bcf'
    threads: 2
    resources:
        mem_mb = 4000,
        disk_scratch = 10
    shell:
        '''
        bcftools norm --threads {threads} -f {config[reference]} -m -any -Ou {input} > {params.bcf}
        bcftools view -i 'abs(ILEN)>={params.SV_size}' -o {output.SV} {params.bcf}
        bcftools view -e 'abs(ILEN)>={params.SV_size}' -o {output.small} {params.bcf}
        tabix -fp vcf {output.small}
        '''

rule jasmine_intersect:
    input:
        read = 'variant_calling/samples.sniffles.autosomes.vcf',
        asm = 'variant_calling/panel.SV.vcf'
    output:
        'jasmine.vcf',
        'jasmine_overlaps.txt'
    params:
        _input = lambda wildcards, input: ','.join(input)
    conda:
        'jasmine'
    threads: 2
    resources:
        mem_mb = 3000,
        disk_scratch = 5
    shell:
        '''
        jasmine --comma_filelist file_list={params._input} threads={threads} out_file={output[0]} out_dir=$TMPDIR spec_reads=0 genome_file={config[reference]} min_seq_id=.5 --pre_normalize --ignore_strand --allow_intrasample --normalize_type
        grep -vE "SVTYPE=(INV|TRA)" {output[0]} | grep -oP "(SVLEN=-?\d*|SUPP_VEC=\d{{2}})" | sed 's/[A-Z,=,_]*//g'  | paste -s -d' \n' > {output[1]}
        '''
#jasmine --comma_filelist file_list=smoove_SV/All_filter_type.vcf,eQTL_GWAS/variants/variant_calling/panel.SV.vcf threads=1 out_file=SR_LR.vcf out_dir=$TMPDIR spec_reads=0 genome_file=REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa min_seq_id=0 --pre_normalize --ignore_strand --allow_intrasample max_dist_linear=1 --normalize_type --dup_to_ins
#bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t0\t+\t%POS\t%INFO/END\t%INFO/SUPP_VEC\n' SR_LR.vcf | grep -vE "(INV|TRA)"  | sed s'/10$/50,200,50/g' | sed s'/11$/250,100,150/g' | sed s'/01$/50,20,250/g' >> SR_LR.bed 

## Assembly rules:
rule minimap2_align_asm:
    input:
        lambda wildcards: config['asm'][wildcards.sample]
    output:
        'alignments/{sample}.asm.paf'
    threads: 4
    resources:
        mem_mb = 10000
    shell:
        '''
        minimap2 -cx asm5 -t {threads} --cs {config[reference]} {input} > {output}
        '''

rule paftools_call:
    input:
        'alignments/{sample}.asm.paf'
    output:
        'variant_calling/{sample}.asm.paf'
    threads: 1
    resources:
        mem_mb = 2000
    shell:
        '''
        sort -k6,6 -k8,8n {input} | paftools.js call -f {config[reference]} {input} > {output.vcf}
        '''
