rule all:
    input:
        'samples.txt'

rule get_samples:
    output:
        'samples.txt'
    run:
        samples = {'OBV':20,'BSW':20,'HOL':2,'FV':2,'SIM':5}
        for breed, quantity in samples.items():
            shell(f"awk '$3==\"{breed}\"&&$4>15&&$4<30 {{{{print $1,$2}}}}' /cluster/work/pausch/vcf_UCD/2020_09/stats/summary_coverage.txt | shuf | head -n {quantity} >> {output}")
        
