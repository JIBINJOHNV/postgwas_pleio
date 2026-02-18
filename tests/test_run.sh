

docker run --platform linux/amd64 -it jibinjv/postgwas-pleio:1.0 python /opt/mtag/mtag.py --help


docker run --rm --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio --help

docker run --rm -v /Users/JJOHN41/:/Users/JJOHN41/ \
    --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio meta-analysis mtag pipeline \
    --vcfs /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/i42_GRCh37_merged.vcf.gz /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/i65_GRCh37_merged.vcf.gz \
    --run_name kadoorie_biobank \
    --out /Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ 