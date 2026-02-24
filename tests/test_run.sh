

docker run --platform linux/amd64 -it jibinjv/postgwas-pleio:1.0 python /opt/mtag/mtag.py --help


docker run --rm -it --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio --help

docker run --rm -v /Users/JJOHN41/:/Users/JJOHN41/ -it --platform linux/amd64 jibinjv/postgwas-pleio:1.1 bash 


docker run --rm -v /Users/JJOHN41/:/Users/JJOHN41/ \
    --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio meta-analysis mtag pipeline \
    --inputfile /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/pleio_input.txt \
    --run_name mtag_test \
    --cores 8 \
    --out /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/mtag_test \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ 


docker run --rm -v /Users/JJOHN41/:/Users/JJOHN41/ \
    --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio meta-analysis pleio pipeline \
    --inputfile /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/pleio_input.txt \
    --out /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/pleio_test/ \
    --run_name pleio_test \
    --cores 8 \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --nis 100000 \
    --flattening_p_values

docker run --rm -it -v /Users/JJOHN41/:/Users/JJOHN41/ \
    --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio meta-analysis asset pipeline  \
    --inputfile /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/pleio_input.txt \
    --run_name fastasset_testrun \
    --out /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/asset_test/ \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --cores 8 \
    --hm3 /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
    --info_filter 0.9 \
    --maf_filter 0.01 \
    --chunk_size 100000 \
    --scr_pthr 0.05 \
    --meth_pval DLM 

docker run --platform linux/amd64 --rm -it  -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2  \
    postgwas-pleio meta-analysis metal pipeline \
    --inputfile /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/pleio_input.txt \
    --run_name metal_testrun\
    --out /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/metal_test/ \
    --scheme SAMPLESIZE \
    --heterogeneity \
    --genomic-control \
    --overlap-correction \
    --track-freq  

docker run --platform linux/amd64 --rm -it  -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2  \
    postgwas-pleio meta-analysis metal pipeline \
    --inputfile /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/pleio_input.txt \
    --run_name metal_testrun_stderr \
    --out /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/metal_test/ \
    --scheme STDERR \
    --heterogeneity \
    --genomic-control \
    --track-freq 


docker run --platform linux/amd64 --rm -it  \
    -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis placo  pipeline \
    --inputfile /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/placo_input.txt \
    --run_name placo_testrun \
    --cores 8 \
    --out /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/placo_test/ \
    --method placo.plus 







docker run --rm -v /Users/JJOHN41/:/Users/JJOHN41/ \
    --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio meta-analysis pleio direct \
    --out /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/pleio_output_direct/ \
    --metain /Users/JJOHN41/Downloads/kadoorie_biobank/pleio_output/pleio_ldsc/metain.txt.gz \
    --ce /Users/JJOHN41/Downloads/kadoorie_biobank/pleio_output/pleio_ldsc/ce.txt.gz \
    --sg /Users/JJOHN41/Downloads/kadoorie_biobank/pleio_output/pleio_ldsc/sg.txt.gz \
    --blup \
    --parallel 

docker run --rm -it -v /Users/JJOHN41/:/Users/JJOHN41/ \
    --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio meta-analysis asset pipeline  \
    --inputfile /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/pleio_input.txt \
    --run_name fastasset_testrun \
    --out /Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output/fastasset_pipeline/ \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --ncores 8 \
    --hm3 /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
    --info_filter 0.9 \
    --maf_filter 0.01 \
    --chunk_size 100000 \
    --scr_pthr 0.05 \
    --meth_pval DLM 




docker run --rm -it -v /Users/JJOHN41/:/Users/JJOHN41/ \
    --platform linux/amd64 jibinjv/postgwas-pleio:1.1 Rscript src/postgwas_pleio/meta_analysis/fastasset/asset_cli.R \
    --input_manifest /Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output/fastasset/1_fastasset_inputs/fastasset_testrun_fastasset_ldsc_input.txt \
    --fastasset_input /Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output/fastasset/1_fastasset_inputs/fastasset_testrun_fastasset_input.tsv \
    --output_dir /Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output/fastasset/ \
    --hm3 /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
    --ld_ref /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/




input_manifest_df= pd.read_csv("/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/pleio_input.txt",sep="\t")

tsv_file="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/kadoorie_biobank_test.tsv"
joined_out_tsv="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/kadoorie_biobank_test_joined.tsv"
sample_specific_directory="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/"
joined_remove_duplicate_ids=True,
joined_create_unique_id=False,
joined_sample_retained_cols=["ID","CHROM","POS","REF","ALT","AF", "ES","SI","SE"]
joined_header_map={"CHROM":"CHR","POS":"Start"}
joined_suffix_replace={":ES":"_Beta","SI":"_INFO"}
sample_remove_duplicate_ids=True
sample_create_unique_id=False
sample_specific_retained_cols=["ID","CHROM","POS","REF","ALT","AF", "ES","SI","SE"]
sample_header_map={"CHROM":"CHR","POS":"Start"}
sample_suffix_replace={":ES":"_Beta","SI":"_INFO"}
log_file="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/kadoorie_biobank_test.tsv.log"


run_postprocess_pipeline(
        mode={"joined_samples","sample_specific"},
        tsv_file=tsv_file,
        input_manifest_df=input_manifest_df,
        joined_out_tsv=joined_out_tsv,
        sample_specific_directory=sample_specific_directory,
        joined_remove_duplicate_ids=joined_remove_duplicate_ids,
        joined_create_unique_id=joined_create_unique_id,
        joined_sample_retained_cols=joined_sample_retained_cols,
        joined_header_map=joined_header_map,
        joined_suffix_replace=joined_suffix_replace,
        sample_remove_duplicate_ids=sample_remove_duplicate_ids,
        sample_create_unique_id=sample_create_unique_id,
        sample_specific_retained_cols=sample_specific_retained_cols,
        sample_header_map=sample_header_map,
        sample_suffix_replace=sample_suffix_replace,
        log_file=log_file )
