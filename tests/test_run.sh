

docker run --platform linux/amd64 -it jibinjv/postgwas-pleio:1.0 python /opt/mtag/mtag.py --help


docker run --rm -it --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio --help

docker run --rm -v /Users/JJOHN41/:/Users/JJOHN41/ -it --platform linux/amd64 jibinjv/postgwas-pleio:1.1 bash 


docker run --rm -v /Users/JJOHN41/:/Users/JJOHN41/ \
    --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio meta-analysis mtag pipeline \
    --vcfs /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/i42_GRCh37_merged.vcf.gz /Users/JJOHN41/Downloads/kadoorie_biobank/vcf_files/i65_GRCh37_merged.vcf.gz \
    --run_name kadoorie_biobank \
    --cores 8 \
    --out /Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ 




docker run --rm -v /Users/JJOHN41/:/Users/JJOHN41/ \
    --platform linux/amd64 jibinjv/postgwas-pleio:1.1 postgwas-pleio meta-analysis pleio pipeline \
    --inputfile /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/pleio_input.txt \
    --out /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/pleio_output/ \
    --run_name osteoarthritis_erectile_dysfunction \
    --cores 8 \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --nis 100000 \
    --flattening_p_values


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




'micromamba', 'run', '-n', 'postgwas', 'Rscript', 
'/opt/conda/envs/postgwas/lib/python3.11/site-packages/postgwas_pleio/meta_analysis/asset/R/asset_cli.R', '--input_manifest', 
'/Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output/fastasset_pipeline/1_fastasset_inputs/fastasset_testrun_fastasset_ldsc_input.txt', '--fastasset_input', 
'/Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output/fastasset_pipeline/1_fastasset_inputs/fastasset_testrun_fastasset_input.tsv', '--output_dir', 
'/Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output/fastasset_pipeline', '--hm3', 
'/Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist', '--ld_ref', 
'/Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/', '--info_filter', '0.9', '--maf_filter', '0.01', '--chunk_size', 
'100000', '--scr_pthr', '0.05', '--meth_pval', 'DLM', '--ncores', '8', '--run_name', 'fastasset_testrun'

 
inputfile='/Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/pleio_input.txt'
output_folder='/Users/JJOHN41/Downloads/kadoorie_biobank/mtag_output/fastasset/'
run_name='kadoorie_biobank'

formatter_out=fastasset_formatter(inputfile=inputfile,
         output_folder=output_folder, 
         run_name="fastasset_testrun", 
         n_threads=8)  