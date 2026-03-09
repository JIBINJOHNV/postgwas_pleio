

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





docker run --platform linux/amd64 --rm -it \
    -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2 \
        postgwas-pleio meta-analysis  genomicpca pipeline \
        --inputfile /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/pleio_input.txt \
        --run_name genomicpca_test \
        --cores 8 \
        --hm3 /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
        --ld_ref /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --out /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/genomicpca_test \
        --info_filter 0.7 \
        --maf_filter 0.01 \
        --approach both 













docker run --platform linux/amd64 --rm -it  -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2  bash 

out_path="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/gpca_test/"
out_path = Path(out_path)
val_input = validate_manifest("/Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/placo_input.txt", str(out_path))
run_name="gpca_test"
cores=4
pre_process_dir = out_path / f"0_{run_name}_pre_process"
pre_process_dir.mkdir(parents=True, exist_ok=True)
out_tsv = pre_process_dir / f"{run_name}_master_sumstats.tsv"
preprocess_log_file = pre_process_dir / f"{run_name}_master_sumstats.log"


master_vcf_tsv, success_2 = run_merge_vcf_to_tsv_pipeline(
        input_df=val_input["manifest_df"],
        out_tsv=str(out_tsv),
        fixed_fields=["ID", "CHROM", "POS", "REF", "ALT"],
        format_fields=["AF", "ES", "SE", "NEF", "LP", "SI", "EZ"],
        log_file=str(preprocess_log_file),
        n_threads=int(cores),
        qc_n_rows=5_000_000,
        bcftools_exec="bcftools",
        fail_fast_on_merge=True,
    )

asset_input_dir = out_path / f"1_{run_name}_asset_input"
ldsc_input_dir = out_path / f"2_{run_name}_ldsc_input"
asset_input_dir.mkdir(parents=True, exist_ok=True)
ldsc_input_dir.mkdir(parents=True, exist_ok=True)


asset_inputs = run_mastertsv_to_toolformat_pipeline(
        mode={"joined_samples", "sample_specific","sample_specific_ldsc"},
        tsv_file=master_vcf_tsv,
        input_manifest_df=val_input["manifest_df"],
        joined_out_tsv=str(asset_input_dir / f"{run_name}_asset_joined.tsv"),
        sample_specific_directory=str(asset_input_dir),
        sample_specific_ldsc_directory=str(ldsc_input_dir),
        joined_remove_duplicate_ids=True,
        joined_create_unique_id=True,
        joined_sample_retained_cols=["ID", "NEF", "ES", "SE"],
        joined_header_map=None,
        joined_suffix_replace={":ES": ".Beta", ":SE": ".SE", ":NEF": ".N"},
        joined_convert_to_raw_p_value=False,
        sample_remove_duplicate_ids=True,
        sample_create_unique_id=True,
        sample_specific_retained_cols=["ID", "CHROM", "POS", "ALT","REF", "AF", "NEF","EZ","LP"],
        sample_header_map={"ID": "SNPID", "CHROM": "CHR", "POS": "BP","ALT": "EA", "REF": "OA"},
        sample_suffix_replace={"AF":"EAF","NEF":"N","EZ":"Z", "LP":"P"},
        sample_convert_to_raw_p_value=True,
        sample_ldsc_remove_duplicate_ids=False,
        sample_ldsc_create_unique_id=False,
        sample_ldsc_specific_retained_cols=["ID","REF", "ALT", "AF", "NEF", "ES", "LP", "SI"],
        sample_ldsc_header_map={"ID": "SNPID", "REF": "A2", "ALT": "A1"},
        sample_ldsc_suffix_replace={"NEF":"N","AF":"MAF","ES":"effect", "LP":"P","SI":"INFO" },
        sample__ldscconvert_to_raw_p_value = True,
        log_file=str(preprocess_log_file),  
    )


# Build ASSET input manifest (NAME + FILE) from sample_specific table
asset_input_df = asset_inputs["sample_specific"].rename({"matrix_tsv": "FILE", "sample_id": "NAME"})
asset_ldsc_input_df = asset_inputs["sample_specific_ldsc"].rename({"ldsc_matrix_tsv": "ldsc_FILE", "sample_id": "NAME"}).select(["sumstat_vcf", "NAME", "ldsc_FILE"])
asset_input_df = asset_input_df.join( asset_ldsc_input_df, on=["sumstat_vcf", "NAME"], how="left")

# IMPORTANT: Polars does not support pandas-style assignment.
# If your R script requires PPREV/SPREV columns, they must already exist in asset_input_df,
# or you should add them using Polars expressions. Here we only fill if columns exist.
if "PPREV" in asset_input_df.columns:
    asset_input_df = asset_input_df.with_columns(
        asset_input_df["PPREV"].fill_null("nan").alias("PPREV")
    )
if "SPREV" in asset_input_df.columns:
    asset_input_df = asset_input_df.with_columns(
        asset_input_df["SPREV"].fill_null("nan").alias("SPREV")
    )


asset_input_file = asset_input_dir / f"{run_name}_asset_input_manifest.tsv"
asset_input_df.write_csv(asset_input_file, separator="\t")


munge_df<-run_munge_stage(input_manifest='/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/gpca_test/1_gpca_test_asset_input/gpca_test_asset_input_manifest.tsv'
                            hm3="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist",
                            output_dir="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/gpca_test/3_ldsc_analysis",
                            info_filter =0.7,
                            maf_filter =0.3,
                            run_name = "munge")



ldsc=run_ldsc_stage(df=munge_df, 
        ld_ref="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/", 
        output_dir="/Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/gpca_test/",
        run_name="ldsc")









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
