CHROMOSOMES=list(range(1,23))

rule buildchrNameChange:
    input:
        expand("Gnomad/1kg.chr{chrom}.bcf{ext}", chrom=CHROMOSOMES, ext=["", ".csi"]),
        expand("Gnomad/1kg.chr{chrom}.sites.vcf.gz{ext}", chrom=CHROMOSOMES, ext=["", ".tbi"]),
        expand("Gnomad/1kg.chr{chrom}.sites.tsv.gz{ext}",chrom=CHROMOSOMES,ext=["", ".tbi"])

rule create_hg38_bcf:
    input:
        ibcf = "Gnomad/1kg.{chrom}.bcf",
        icsi = "Gnomad/1kg.{chrom}.bcf.csi"
    output:
        obcf = "Gnomad/1kg.chr{chrom}.bcf",
        ocsi = "Gnomad/1kg.chr{chrom}.bcf.csi"
    shell:
        """
         bcftools annotate --rename-chrs chr_name_conv_nonetochr.txt {input.ibcf} -Ob -o {output.obcf}
         bcftools index -f {output.obcf}
        """

rule create_hg38_vcf:
    input:
        ivcf = "Gnomad/1kg.{chrom}.sites.vcf.gz",
        itbi = "Gnomad/1kg.{chrom}.sites.vcf.gz.tbi"
    output:
        ovcf = "Gnomad/1kg.chr{chrom}.sites.vcf.gz",
        otbi = "Gnomad/1kg.chr{chrom}.sites.vcf.gz.tbi"
    shell:
        """
        bcftools annotate --rename-chrs chr_name_conv_nonetochr.txt {input.ivcf} -Oz -o {output.ovcf}
        tabix -f {output.ovcf}
        """

rule create_hg38_tsv:
    input:
        itsv = rules.create_hg38_bcf.output.obcf,
        icsi = rules.create_hg38_bcf.output.ocsi
    output:
        otsv = "Gnomad/1kg.chr{chrom}.sites.tsv.gz",
        ocsi = "Gnomad/1kg.chr{chrom}.sites.tsv.gz.tbi"
    shell:
        """
        bcftools query -f "%CHROM\\t%POS\\t%REF,%ALT\\n" {input.itsv} |  bgzip -c > {output.otsv}
        tabix -f -s1 -b2 -e2 {output.otsv}
        """