rule prep_ref_genome:
    input:
        "ref/human_g1k_v37.fasta",
        "ref/human_g1k_v37.fasta.fai",
        "ref/human_g1k_v37.fasta.bwt",
        "ref/human_g1k_v37.dict"

rule wget_ref:
    output:
        "ref/human_g1k_v37.fasta.gz",
    shell:
        "wget -P ref http://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz"

rule extract_ref:
    input:
        "ref/human_g1k_v37.fasta.gz"
    output:
        "ref/human_g1k_v37.fasta"
    shell: "set +euo pipefail; gunzip --quiet {input}; exit 0"

rule samtools_faidx:
    input:
        "ref/human_g1k_v37.fasta",
    output:
        "ref/human_g1k_v37.fasta.fai"
    conda:
        "../envs/misc.yaml"
    shell:
        "samtools faidx {input}"

rule samtools_dict:
    input:
        "ref/human_g1k_v37.fasta",
    output:
        "ref/human_g1k_v37.dict"
    conda:
        "../envs/misc.yaml"
    shell:
        "samtools dict -o {output} {input}"

rule bwa_index:
    input:
        "ref/human_g1k_v37.fasta",
    output:
        "ref/human_g1k_v37.fasta.bwt"
    conda:
        "../envs/misc.yaml"
    shell:
        "bwa index {input}"

rule downlad_gmap:
    output:
        MAP
    shell:
        " wget https://github.com/odelaneau/GLIMPSE/raw/master/maps/genetic_maps.b37/chr{wildcards.chrom}.b37.gmap.gz -O {output}"
