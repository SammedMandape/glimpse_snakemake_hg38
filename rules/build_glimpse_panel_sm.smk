CHROMOSOMES=list(range(1,23))

rule build_glimpse_panel:
    input:
        expand("Gnomad/1kg.{chrom}.bcf{ext}", chrom=CHROMOSOMES, ext=["", ".csi"]),
        expand("Gnomad/1kg.{chrom}.sites.vcf.gz{ext}", chrom=CHROMOSOMES, ext=["", ".tbi"]),
        expand("Gnomad/1kg.{chrom}.sites.tsv.gz{ext}", chrom=CHROMOSOMES, ext=["", ".tbi"])

rule create_biallelic_snp_bcf:
    input:
        obcf = "Gnomad/eagle.chr{chrom}.ensembl.bcf",
        ocsi = "Gnomad/eagle.chr{chrom}.ensembl.bcf.csi"
    output:
        bcf="Gnomad/1kg.{chrom}.bcf",
        csi="Gnomad/1kg.{chrom}.bcf.csi"
    # conda:
    #     "../envs/misc.yaml"
    shell:
        """
        bcftools view -r {wildcards.chrom} -m2 -M2 -v snps -Ob -o {output.bcf} {input.obcf}
        bcftools index -f {output.bcf}
        """

rule create_site_vcf:
    input:
        bcf=rules.create_biallelic_snp_bcf.output.bcf,
        csi=rules.create_biallelic_snp_bcf.output.csi
    output:
        vcfgz="Gnomad/1kg.{chrom}.sites.vcf.gz",
        tbi="Gnomad/1kg.{chrom}.sites.vcf.gz.tbi"
    # conda:
    #     "../envs/misc.yaml"
    shell:
        """
        bcftools view -G -m2 -M2 -v snps {input.bcf} -Oz -o {output.vcfgz}
        tabix -f {output.vcfgz}
        """

rule create_site_tsv:
    input:
        bcf=rules.create_biallelic_snp_bcf.output.bcf,
        csi=rules.create_biallelic_snp_bcf.output.csi
    output:
        tsvgz="Gnomad/1kg.{chrom}.sites.tsv.gz",
        tbi="Gnomad/1kg.{chrom}.sites.tsv.gz.tbi"
    # conda:
    #     "../envs/misc.yaml"
    shell:
        """
        bcftools query -f "%CHROM\\t%POS\\t%REF,%ALT\\n" {input.bcf} | bgzip -c > {output.tsvgz}
        tabix -f -s1 -b2 -e2 {output.tsvgz}
        """
