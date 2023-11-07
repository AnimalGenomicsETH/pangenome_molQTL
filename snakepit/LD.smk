from pathlib import PurePath

rule bcftools_annotate:
    input:
        'QTL/PanGenie.vcf.gz'
    output:
        annotation = multiext('LD/annotation.bed.gz','','.tbi'),
        annotated = multiext('LD/variants.vcf.gz','','.tbi'),
        nice = 'LD/nicer_variants.txt'
    threads: 2
    shell:
        '''
        paste <(bcftools query -i 'MAF>0.01' -f '%CHROM\t%POS\t%REF\t%ALT\n' {input}) <(seq 1 $(bcftools index -n {input})) | bgzip -c > {output.annotation[0]}

        tabix -p vcf -b 2 -e 2 {output.annotation[0]}

        zcat {output.annotation[0]} | mawk ' {{ A[$1"_"$2]=$0 }} END {{ for (key in A) {{ print A[key] }} }} ' | cut -f 5 | awk '$1' | sort -n > {output.nice}

        bcftools annotate --threads {threads} -a {output.annotation[0]} -o {output.annotated[0]} -c CHROM,POS,REF,ALT,ID {input}
        tabix -p vcf {output.annotated[0]}
        '''

rule bcftools_query:
    input:
        rules.bcftools_annotate.output['annotated']
    output:
        'LD/SV_IDs.txt'
    shell:
        '''
        bcftools query -i 'abs(ILEN)>=50' -f '%ID\n' {input[0]} > {output}
        '''

rule plink2_LD:
    input:
        vcf = rules.bcftools_annotate.output['annotated'][0],
        nicer_variants = rules.bcftools_annotate.output['nice'],
        tags = rules.bcftools_query.output
    output:
        'LD/samples.SV.{r2}.{window}.tags.list'
    threads: 4
    resources:
        mem_mb = 5000,
        walltime = '60'
    params:
        mem = lambda wildcards, threads, resources: threads*resources.mem_mb,
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    shell:
        '''
        plink --vcf {input.vcf} --show-tags {input.tags} --tag-r2 0.{wildcards.r2} --tag-kb {wildcards.window} --threads {threads} --memory {params.mem} --chr-set 30 --vcf-half-call h --list-all --out {params.out}
        '''
