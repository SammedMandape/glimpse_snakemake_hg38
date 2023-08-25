#SAMPLES=["ERR262997"]
#READS=[1,2]
#COVS=[1,0.1]
SEED="20220503"
CHROMOSOMES=list(range(1,23))
REF_GENOME="/eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
SITEVCF=ancient("Gnomad/1kg.chr{chrom}.sites.vcf.gz")
SITETSV=ancient("Gnomad/1kg.chr{chrom}.sites.tsv.gz")
MAP="Gnomad/chr{chrom}.b38.gmap.gz"
GLIMPSE_VCF="Gnomad/1kg.chr{chrom}.bcf"
SEED="20220503"
REF_GENOME="/eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"

#include: "rules/build_glimpse_panel_sm.smk"
#include: "rules/prep_ref_genome.smk"
include: "chrNameChange.smk"

rule all:
    input:
        expand("chunk/chunks.chr{chrom}.txt",chrom=CHROMOSOMES)

rule make_chunks:
    input:
        chrom_sites="Gnomad/1kg.chr{chrom}.sites.vcf.gz"
    output:
        "chunk/chunks.chr{chrom}.txt"
    singularity: "docker://ghcr.io/signaturescience/glimpse:latest"
    shell:
        "singularity run glimpse.sif chunk --input {input.chrom_sites} --region chr{wildcards.chrom} --window-size 2000000 --buffer-size 200000 --output {output}"
