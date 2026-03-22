### Downloaded git clone https://github.com/JIBINJOHNV/postgwas.git 
    cd postgwas/tests 
    edit the file path mentioned in harmonisation.yaml



##  pull postgwas docker iamges 
docker run --rm --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -it jibinjv/postgwas:1.3 postgwas --help


## create config file 
docker run --rm --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/:/Users/JJOHN41/ \
  -it jibinjv/postgwas:1.3 python /opt/postgwas/src/postgwas/scripts/create_sumstat_map_pl.py \
  --input /Users/JJOHN41/Downloads/meta_analysis_testing/raw_data/ \
  --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
  --output-path /Users/JJOHN41/Downloads/meta_analysis_testing/harmonised/harmonisatio_example_input_file.csv \
  --harmonisation-output-path /Users/JJOHN41/Downloads/meta_analysis_testing/harmonised/


## open the newly created config file and check  : /Users/JJOHN41/Downloads/meta_analysis_testing/harmonised/harmonisatio_example_input_file.csv

## Perform harmormonisation
docker run --rm  --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    -it jibinjv/postgwas:1.3 postgwas harmonisation \
    --nthreads 10 \
    --max-mem 50G \
    --config /Users/JJOHN41/Downloads/meta_analysis_testing/harmonised/harmonisatio_example_input_file.csv \
    --defaults /Users/JJOHN41/Downloads/meta_analysis_testing/postgwas/tests/harmonisation.yaml

docker run --rm  --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    -it jibinjv/postgwas:1.3 postgwas harmonisation \
    --nthreads 10 \
    --max-mem 50G \
    --config /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/pleio/4_harmonisation_input_files/pleio_sczbip_harmonisation_manifest.tsv \
    --defaults /Users/JJOHN41/Downloads/meta_analysis_testing/postgwas/tests/harmonisation.yaml



/Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/pleio/4_harmonisation_input_files/pleio_sczbip_harmonisation_manifest.tsv
####-------------------------------------------Meta Analysis----------------------------------------
## Prepare your metaanalysis_input_file.tsv it should contain four columns namely 
# sumstat_vcf	TYPE	SPREV	PPREV	sample_id 
# Make sure sample_id should match with sample name inside vcf file , otherwise it will through error 



## Metal analysis 

docker run --platform linux/amd64 --rm -it \
    -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis metal pipeline --help 

docker run --platform linux/amd64 \
    --rm -it \
    --user root \
    -e HOME=/Users/JJOHN41 \
    -e DOCKER_HOST=unix:///var/run/docker.sock \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis metal pipeline \
    --inputfile /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis_input_file.tsv \
    --run_name metal_sczbip \
    --out /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/metal/ \
    --scheme SAMPLESIZE \
    --heterogeneity \
    --overlap-correction \
    --track-freq \
    --harmonise \
    --defaults_config /Users/JJOHN41/Downloads/meta_analysis_testing/postgwas/tests/harmonisation.yaml \
    --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
    --sample_size_approach weight \
    --info_method mean \
    --nthreads 12  

  #     --genomic-control \


docker run --platform linux/amd64 \
    --rm -it \
    --user root \
    -e HOME=/Users/JJOHN41 \
    -e DOCKER_HOST=unix:///var/run/docker.sock \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis metal pipeline \
    --inputfile /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis_input_file.tsv \
    --run_name metal_sczbip_stderror \
    --out /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/metal_stderr/ \
    --scheme STDERR \
    --heterogeneity \
    --track-freq \
    --harmonise \
    --defaults_config /Users/JJOHN41/Downloads/meta_analysis_testing/postgwas/tests/harmonisation.yaml \
    --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
    --sample_size_approach totalnef \
    --info_method mean 




## mtag analysis
docker run --platform linux/amd64 \
    --rm -it \
    --user root \
    -e HOME=/Users/JJOHN41 \
    -e DOCKER_HOST=unix:///var/run/docker.sock \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis mtag pipeline \
    --inputfile /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis_input_file.tsv \
    --run_name mtag_sczbip \
    --nthreads 8 \
    --out /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/mtag/ \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --harmonise \
    --defaults_config /Users/JJOHN41/Downloads/meta_analysis_testing/postgwas/tests/harmonisation.yaml \
    --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
    --info_method mean  


docker run --platform linux/amd64 \
    --rm -it \
    --user root \
    -e HOME=/Users/JJOHN41 \
    -e DOCKER_HOST=unix:///var/run/docker.sock \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis mtag pipeline \
    --inputfile /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis_input_file.tsv \
    --run_name mtag_sczbip_perfectgencov_equalh \
    --nthreads 8 \
    --perfect_gencov \
    --equal_h2 \
    --out /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/mtag_sczbip_perfectgencov_equalh/ \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --harmonise \
    --defaults_config /Users/JJOHN41/Downloads/meta_analysis_testing/postgwas/tests/harmonisation.yaml \
    --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
    --info_method mean  

## Pleio test 

docker run --platform linux/amd64 \
    --rm -it \
    --user root \
    -e HOME=/Users/JJOHN41 \
    -e DOCKER_HOST=unix:///var/run/docker.sock \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis pleio pipeline \
    --inputfile /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis_input_file.tsv \
    --out /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/pleio/  \
    --run_name pleio_sczbip \
    --nthreads 8 \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --nis 100000 \
    --harmonise \
    --defaults_config /Users/JJOHN41/Downloads/meta_analysis_testing/postgwas/tests/harmonisation.yaml \
    --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
    --info_method mean  

        # --flattening_p_values


docker run --platform linux/amd64 \
    --rm -it \
    --user root \
    -e HOME=/Users/JJOHN41 \
    -e DOCKER_HOST=unix:///var/run/docker.sock \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis pleio pipeline \
    --inputfile /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis_input_file.tsv \
    --out /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/pleio_fltap/  \
    --run_name pleio_sczbip \
    --nthreads 8 \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --nis 100000 \
    --harmonise \
    --defaults_config /Users/JJOHN41/Downloads/meta_analysis_testing/postgwas/tests/harmonisation.yaml \
    --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
    --info_method mean \
    --flattening_p_values


## genomicpca
docker run --platform linux/amd64 --rm -it \
    -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2 \
        postgwas-pleio meta-analysis  genomicpca pipeline \
        --inputfile /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis_input_file.tsv \
        --run_name genomicpca_sczbip  \
        --cores 8 \
        --hm3 /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
        --ld_ref /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --out /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/genomicpca \
        --info_filter 0.7 \
        --maf_filter 0.01 \
        --approach both 


docker run --platform linux/amd64 -it \
    -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2 bash 

## placo
docker run --platform linux/amd64 --rm -it  \
    -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis placo  pipeline \
    --inputfile /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis_input_file.tsv \
    --run_name placo_sczbip \
    --cores 8 \
    --out /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/placo \
    --method placo.plus 



## asset 
docker run --platform linux/amd64 --rm -it  \
    -v /Users/JJOHN41/:/Users/JJOHN41/ jibinjv/postgwas-pleio:1.2 \
    postgwas-pleio meta-analysis asset pipeline  \
    --inputfile /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis_input_file.tsv \
    --run_name fastasset_testrun \
    --out /Users/JJOHN41/Downloads/meta_analysis_testing/metaanalysis/asset \
    --ld_ref_panel /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --cores 8 \
    --hm3 /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
    --info_filter 0.9 \
    --maf_filter 0.01 \
    --chunk_size 100000 \
    --scr_pthr 0.05 \
    --meth_pval DLM 