# scripts_vcf
Scripts to handle and filter VCF files (intended for personal use)

#Dependencies

- Platypus callable as *
- javarkit's VcfMultiToOneAllele callable as vcfmulti2oneallele

In order to make these scripts work you need to set the environment variables:
    - SCRIPTSVCF_DIR to indicate the directory where they are located
    - GNOMAD to indicate the tabix VCF file of gnomAD
    - HUMANDB_DIR to indicate where the reference human genome is located for annovar
    - SNPSIFT_DIR to indicate where the file SnpSift.jar from the SnpEff package is located

##Prepare gnomAD data
gunzip -c gnomad.genomes.r2.1.sites.vcf.bgz | head -n 1000 | grep '^#' > gnomad.genomes.r2.1.sites.onlyAFINFO.vcf
gunzip -c gnomad.genomes.r2.1.sites.vcf.bgz | perl -pe 's/\t[^\t]*;(AF=[^;]+).*$/\t$1/' >> gnomad.genomes.r2.1.sites.onlyAFINFO.vcf
