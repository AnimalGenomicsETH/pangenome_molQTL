from pathlib import PurePath

## Add reference path for happy container
workflow.singularity_args = f'-B $TMPDIR -B {PurePath(config["reference"]).parent}'

rule all:
    input:
        #expand('SVs/{sample}.denovo.sniffles.vcf.gz',sample=config['samples']),
        #'SVs/samples.denovo.sniffles.vcf.gz',
        'SNPs/F1.csv',
        'SNPs/isec.csv'

rule sniffles_call:
    input:
        bam = lambda wildcards: config['SV_samples'][wildcards.sample]#'/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long-read_alignments/PacBio_CCS/eQTL/{sample}.mm2.cram'
    output:
        vcf = temp(multiext('SVs/{sample}.denovo.sniffles.vcf.gz','','.tbi')),
        snf = temp('SVs/{sample}.denovo.sniffles.snf')
    threads: 4
    resources:
        mem_mb = 2500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.bam} --reference {config[reference]} --sample-id {wildcards.sample} --phase --minsvlen 50 --threads {threads} --vcf {output.vcf[0]} --snf {output.snf}
        '''

rule sniffles_merge:
    input:
        snfs = expand(rules.sniffles_call.output['snf'],sample=config['SV_samples'],allow_missing=True)
    output:
        vcf = multiext('SVs/samples.{call}.sniffles.vcf.gz','','.tbi')
    threads: 2
    resources:
        mem_mb = 2500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.snfs} --reference {config[reference]} --threads {threads} --vcf {output.vcf[0]}
        '''

rule sniffles_genotype:
    input:
        bam = '/cluster/work/pausch/alex/CCS_eQTL_cohort/ARS/{sample}.mm2.cram',
        vcf = lambda wildcards: '' if wildcards.call == 'denovo' else '/cluster/work/pausch/alex/eQTL_GWAS/SV_ACCURACY/genotyping.vcf.gz' #'/cluster/work/pausch/alex/eQTL_GWAS/SV_ACCURACY/pangenie.vcf'
    output:
        vcf = temp(multiext('SVs/{sample}.{call}.sniffles.vcf.gz','','.tbi'))
    params:
    threads: 4
    resources:
        mem_mb = 1500
    conda:
        'sniffles'
    shell:
        '''
        sniffles --input {input.bam} --genotype-vcf {input.vcf} --reference {config[reference]} --sample-id {wildcards.sample} --minsvlen 50 --threads {threads} --vcf {output.vcf}
        '''

rule bcftools_stuff:
    input:
        'SVs/samples.denovo.sniffles.vcf.gz'
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


## SNPs
def get_variants(variant,caller):
    input_dict = {
              'SNPs_PG':'/cluster/work/pausch/alex/eQTL_GWAS/pangenie_8_wave/samples.all.pangenie_genotyping.vcf.gz',
              'SNPs_DV':'/cluster/work/pausch/alex/eQTL_GWAS/variants/DV-SR/cohort.autosomes.WGS.imputed.vcf.gz',
              'SVs_PG': '/cluster/work/pausch/alex/eQTL_GWAS/pangenie_panel.vcfwave.vcf.gz',
              'SVs_Sniffles': 'SVs/samples.denovo.sniffles.vcf.gz'
              }
    return input_dict[f'{variant}_{caller}']

rule split_vcf:
    input:
        lambda wildcards: get_variants(wildcards.variant,wildcards.caller)#input_dict[f'{wildcards.variant}_{wildcards.caller}']
        #lambda wildcards: '/cluster/work/pausch/alex/eQTL_GWAS/pangenie_8_wave/samples.all.pangenie_genotyping.vcf.gz' if wildcards.caller == 'PG' else '/cluster/work/pausch/alex/eQTL_GWAS/variants/DV-SR/cohort.autosomes.WGS.imputed.vcf.gz'
    output:
        multiext('{variant,SVs|SNPs}/{sample}.{caller,PG|DV|Sniffles}.vcf.gz','','.tbi')
    params:
        regions = ' '.join(map(str,range(1,30))),
        SVs = lambda wildcards: {'SVs_Sniffles':"-i 'F_MISSING<0.2&&abs(ILEN)>=50&&INFO/SVTYPE!=\"BND\"'",'SVs_PG':"-i 'F_MISSING<0.2&&abs(ILEN)>=50'"}.get(f'{wildcards.variant}_{wildcards.caller}',"-i 'F_MISSING<0.2'"),
        sample = lambda wildcards: f'-s {wildcards.sample}' if wildcards.sample != 'all' else ''
    shell:
        '''
        bcftools view -c 1 -a {params.SVs} {params.sample} -o {output[0]} {input} {params.regions}
        tabix -p vcf {output[0]}
        '''

rule bcftools_isec:
    input:
        vcf_truth = expand(rules.split_vcf.output,caller='DV',variant='SNPs',allow_missing=True),
        vcf_query = expand(rules.split_vcf.output,caller='PG',variant='SNPs',allow_missing=True)
    output:
        'SNPs/{sample}.isec'
    shell:
        '''
        bcftools isec -n +1 {input.vcf_truth[0]} {input.vcf_query[0]} | awk ' {{ A[$5]+=1 }} END {{ print "{wildcards.sample}",A["01"],A["10"],A["11"] }}' > {output}
        '''

rule gather_isec:
    input:
        expand(rules.bcftools_isec.output,sample=config['samples'])
    output:
        'SNPs/isec.csv'
    localrule: True
    shell:
        '''
        {{ echo "sample PG DV Mutual" ; cat {input} ; }} > {output} 
        '''

rule happy:
    input:
        vcf_truth = expand(rules.split_vcf.output,caller='DV',variant='SNPs',allow_missing=True),
        vcf_query = expand(rules.split_vcf.output,caller='PG',variant='SNPs',allow_missing=True),
        reference = config['reference']
    output:
        csv = 'SNPs/{sample}.summary.csv',
        others = temp(multiext('SNPs/{sample}','.bcf','.bcf.csi','.extended.csv','.roc.all.csv.gz','.runinfo.json'))
    params:
        _dir = lambda wildcards, output: PurePath(output.csv).with_suffix('').with_suffix('')
    container: '/cluster/work/pausch/alex/software/images/hap.py_latest.sif'
    threads: 2
    resources:
        mem_mb = 5000,
        scratch = '10G'
    shell:
        '''
        /opt/hap.py/bin/hap.py -r {input.reference} --bcf --usefiltered-truth --no-roc --no-json -L --pass-only --scratch-prefix $TMPDIR -X --threads {threads} -o {params._dir} {input.vcf_truth[0]} {input.vcf_query[0]}
        '''

rule gather_SR_happy:
    input:
        expand(rules.happy.output[0],sample=config['samples'])
    output:
        'SNPs/F1.csv'
    localrule: True
    shell:
        '''
        echo -e "variant truth query recall precision truth_TiTv query_TiTv sample" > {output}
        for i in {input}
        do
          awk -v I=$(basename $i) -F',' '$2=="PASS" {{ split(I,a,"."); print $1,$3,$6,$10,$11,$14,$15,a[1] }}' $i >> {output}
        done
        '''

rule jasmine:
    input:
        vcfs = expand('SVs/{sample}.{caller}.vcf.gz',caller=('PG','Sniffles'),allow_missing=True)
    output:
        'SVs/{sample}.jasmine.txt'
    params:
        _input = lambda wildcards, input: ','.join(input.vcfs).replace('.gz','').replace('SVs','$TMPDIR/SVs')
    conda:
        'jasmine'
    threads: 1
    resources:
        mem_mb= 5000,
        walltime = '60',
        scratch = '5G'
    shell:
        '''
        mkdir -p $TMPDIR/SVs
        pigz -dc {input.vcfs[0]} > $TMPDIR/SVs/{wildcards.sample}.PG.vcf
        pigz -dc {input.vcfs[1]} > $TMPDIR/SVs/{wildcards.sample}.Sniffles.vcf

        java -Xmx6048m -jar /cluster/work/pausch/alex/software/Jasmine/jasmine.jar \
        --comma_filelist file_list={params._input} threads={threads} out_file=/dev/stdout out_dir=$TMPDIR \
        genome_file={config[reference]} --pre_normalize --ignore_strand --allow_intrasample --ignore_type \
        max_dist_linear=1 max_dist=1000 > {output}
        #|\
                #grep -hoP "SUPP_VEC=\K\d+" | awk ' {{ A[$1]+=1 }} END {{ print "{wildcards.sample}",A["01"],A["10"],A["11"] }}' > {output}
        '''

rule gather_jasmine:
    input:
        expand(rules.jasmine.output[0],sample=config['SV_samples'])
    output:
        'SVs/jasmine.csv'
    localrule: True
    shell:
        '''
        {{ echo "sample Sniffles PG Mutual" ; cat {input} ; }} > {output}
        '''
