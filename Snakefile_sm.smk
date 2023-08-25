CHROMOSOMES=list(range(1,23))
SITEVCF=ancient("Gnomad/1kg.chr{chrom}.sites.vcf.gz")
SITETSV=ancient("Gnomad/1kg.chr{chrom}.sites.tsv.gz")
MAP="Gnomad/chr{chrom}.b38.gmap.gz"
GLIMPSE_VCF="Gnomad/1kg.chr{chrom}.bcf"
SEED="20220503"
REF_GENOME="/eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
SAMPLES=["07035-500pg"]
DIR=["MiSeq_run106"]
#include: "rules/build_glimpse_panel_sm.smk"
#include: "rules/prep_ref_genome.smk"

rule all:
    input:
        # bcftools call mpileup
        expand("{dir}/bcfcall/{samplename}_chr{chrom}.la.md.vcf.gz{ext}",dir=DIR,samplename=SAMPLES,chrom=CHROMOSOMES,ext=["", ".tbi"]),
        #expand("bcfcall/{samplename}_chr{chrom}.la.md.vcf.gz", samplename=SAMPLES, chrom=CHROMOSOMES),
        expand("{dir}/glimpse_chunk/{samplename}_chr{chrom}.list",dir=DIR, samplename=SAMPLES, chrom=CHROMOSOMES),
        # GLIMPSE impute outputs for each CHROMOSOMES
        expand("{dir}/glimpse_ligate/{samplename}_chr{chrom}.la.md.vcf.gz",dir=DIR, samplename=SAMPLES, chrom=CHROMOSOMES),
        #expand("glimpse_ligate/{samp}_{chrom}.vcf.gz", samp=SAMPLES, chrom=CHROMOSOMES, suf=['', '.tbi']),
        ancient("Snakefile_sm.smk")

rule bcfcall:
    input:
        bam_in="{dir}/markdup_bams/{samplename}.la.md.bam",
        bai_in="{dir}/markdup_bams/{samplename}.la.md.bam.bai",
        ref=REF_GENOME,
        sites=SITETSV
    output:
        #"bcfcall/{samplename}_chr{chrom}.la.md.vcf.gz"
        vcfout = "{dir}/bcfcall/{samplename}_chr{chrom}.la.md.vcf.gz",
        tbiout = "{dir}/bcfcall/{samplename}_chr{chrom}.la.md.vcf.gz.tbi"
    shell:
        "bcftools mpileup -f {input.ref} -I -E -a 'FORMAT/DP' -T {input.sites} -r chr{wildcards.chrom} {input.bam_in} -Ou "
        " | bcftools call -Aim -C alleles -T {input.sites} -Oz -o {output.vcfout} "
        "&& bcftools index -t {output.vcfout}"

# WIP: GLIMPSE steps here
# File requirements include gmap file and 1000g vcf files
# both are generated from build_glimpse_pane.smk rules
# one glimpse_chunk_phase rule per chunk, per chromosome
# one glimpse_concat rule per chromosome, requires all chunks to be phased


# lambda function to get the positional information from the chunk files.
# This would be super inefficient if chunk files get long, as it reads the file
# 3X per chunk_phase run

checkpoint f:
    input:
        chrom_sites="Gnomad/1kg.chr{chrom}.sites.vcf.gz"
    output:
        "chunk/chunks.chr{chrom}.txt"
    singularity: "docker://ghcr.io/signaturescience/glimpse:latest"
    shell:
        "/GLIMPSE/chunk --input {input.chrom_sites} --region chr{wildcards.chrom} --window-size 2000000 --buffer-size 200000 --output {output}"

def read_chunks(filename):
    return_dict = []
    with open(filename,"r") as fh:
        mylist = fh.read().splitlines()
        for line in mylist:
            short_arry = line.split("\t")
            return_dict.append(short_arry)
    return return_dict

def get_chunk_value(chr,chunk,position):
    chunkfile="chunk/chunks.chr"+str(chr)+".txt"
    with open(chunkfile, "r") as ckf:
        mylist = ckf.read().splitlines()
        for entry in mylist:
            checkvalues = entry.split("\t")
            if checkvalues[0] == chunk:
                return checkvalues[position]

rule glimpse_chunk_phase:
    input:
        chunks="chunk/chunks.chr{chrom}.txt",
        vcf_in="{dir}/bcfcall/{samplename}_chr{chrom}.la.md.vcf.gz",
        index_in="{dir}/bcfcall/{samplename}_chr{chrom}.la.md.vcf.gz.tbi",
        map=MAP,
        ref=GLIMPSE_VCF
    output:
        vcfgzip = "{dir}/glimpse_chunk/{samplename}_{chunk}_chr{chrom}.la.md.vcf.gz",
        indextbi = "{dir}/glimpse_chunk/{samplename}_{chunk}_chr{chrom}.la.md.vcf.gz.tbi"
    singularity:
        "docker://ghcr.io/signaturescience/glimpse:latest"
    params:
        ID= lambda wc: get_chunk_value(wc.chrom,wc.chunk,0),
        IRG=lambda wc: get_chunk_value(wc.chrom,wc.chunk,2),
        ORG=lambda wc: get_chunk_value(wc.chrom,wc.chunk,3)
    log:
        "{dir}/glimpse_chunk/{samplename}_{chunk}_chr{chrom}.log"
    shell:
        "singularity run glimpse.sif phase --input {input.vcf_in} --reference {input.ref} --map {input.map} --input-region {params.IRG} --output-region {params.ORG} --output {output.vcfgzip} "
        "&& bcftools index -t {output.vcfgzip} 2> {log}"

# Function for lambda in concat to generate the number of chunks in a given chromosome
def count_chunks(chr):
    chunkfile="chunk/chunks.chr"+str(chr)+".txt"
    with open(chunkfile, "r") as ckf:
        mylist = ckf.read().splitlines()
    return list(range(0,len(mylist)))

## Intermediate rule to make ligate list files
rule glimpse_ligate_list:
    input:
        gzip=lambda wc: expand("{{dir}}/glimpse_chunk/{{samplename}}_{chunky}_chr{{chrom}}.la.md.vcf.gz", chunky=count_chunks(wc.chrom)),
        index=lambda wc: expand("{{dir}}/glimpse_chunk/{{samplename}}_{chunky}_chr{{chrom}}.la.md.vcf.gz.tbi", chunky=count_chunks(wc.chrom))
    output:
        "{dir}/glimpse_chunk/{samplename}_chr{chrom}.list"
    shell:
        """
        for i in {input.gzip}
        do echo $i >> {output}
        done
        """

# Concats chunks
# inputs of all only ask for these outputs
# temp outputs of phase_chunks are not included in rule all
rule glimpse_concat:
    input:
        list="{dir}/glimpse_chunk/{samplename}_chr{chrom}.list",
        chunks="chunk/chunks.chr{chrom}.txt"
    output:
        "{dir}/glimpse_ligate/{samplename}_chr{chrom}.la.md.vcf.gz"
    container:
        "docker://ghcr.io/signaturescience/glimpse:latest"
    shell:
        """
        singularity run glimpse.sif ligate --input {input.list} --output {output}
        """
