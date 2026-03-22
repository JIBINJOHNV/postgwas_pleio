

def parse_metal_pipeline_args(subparser):
    """
    Defines command-line arguments for the METAL pipeline mode.
    """
    description = (
        "METAL Pipeline: Automated Meta-Analysis Workflow\n\n"
        "1. Reads validated GWAS summary statistic VCFs from manifest.\n"
        "2. Converts and harmonizes them to METAL-ready format.\n"
        "3. Generates a METAL script with consistent column mapping.\n"
        "4. Performs meta-analysis using METAL.\n"
    )

    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis metal pipeline \\\n"
        "    --inputfile input_manifest.tsv \\\n"
        "    --run_name psych_meta \\\n"
        "    --out ./meta_results \\\n"
        "    --scheme STDERR \\\n"
        "    --heterogeneity\n"
        "For detailed explanation of METAL options, visit:\n"
         "https://genome.sph.umich.edu/wiki/METAL_Documentation\n"
    )

    pipe = subparser.add_parser(
        "pipeline",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter
    )

    # =====================================
    # Input Options
    # =====================================

    in_opts = pipe.add_argument_group(title="Input Options")

    in_opts.add_argument(
        "--inputfile",
        required=True,
        metavar="FILE",
        help=(
            "TSV file with columns:\n"
            "  sumstat_vcf, TYPE, SPREV, PPREV, sample_id\n\n"
            "sumstat_vcf : Path to the GWAS summary statistics VCF file.\n"
            "TYPE        : Trait type (binary | quantitative).\n"
            "SPREV       : Sample prevalence (cases / total samples). "
                          "Since this pipeline uses effective sample size for binary traits, "
                          "SPREV is typically set to 0.5. For quantitative traits, this field can be empty.\n"
            "PPREV       : Population prevalence (binary traits only). "
                          "For quantitative traits, this field can be empty. "
                          "If liability-scale conversion is not required, this field may also be left empty.\n"
            "sample_id   : Unique trait identifier. Make sure the sample_id matches "
                           "the sample name inside the VCF file.\n"
        )
    )

    in_opts.add_argument(
        "--run_name",
        required=True,
        metavar="NAME",
        help="Prefix for naming intermediate files, output files, and logs."
    )

    in_opts.add_argument(
        "--nthreads",
        type=int,
        default=max(1, mp.cpu_count() // 2),
        metavar="INT",
        help="Number of parallel threads to use. [default: %(default)s]"
    )

    # =====================================
    # Output Options
    # =====================================

    out_opts = pipe.add_argument_group(title="Output Options")

    out_opts.add_argument(
        "--out",
        default="./metal_results",
        metavar="DIR",
        help="Output directory. [default: %(default)s]"
    )

    # =====================================
    # Meta-Analysis Scheme
    # =====================================

    scheme_opts = pipe.add_argument_group(title="Meta-Analysis Scheme")

    scheme_opts.add_argument(
        "--scheme",
        default="STDERR",
        choices=["STDERR", "SAMPLESIZE"],
        metavar="TYPE",
        help=(
            "Meta-analysis scheme.\n"
            "  STDERR      → Inverse-variance meta-analysis using BETA and SE.\n"
            "                Use when effect sizes are directly comparable across studies\n"
            "                (same phenotype, same units, same transformations).\n\n"
            "  SAMPLESIZE  → Z-score based meta-analysis using Z and sample size.\n"
            "                Use when effect sizes are not directly comparable across studies\n"
            "                (different measurement scales or transformations).\n\n"
            "[default: %(default)s]"
        )
    )

    scheme_opts.add_argument(
        "--heterogeneity",
        action="store_true",
        default=False,
        help="Enable heterogeneity test (METAL HETEROGENEITY ON). [default: %(default)s]"
    )

    scheme_opts.add_argument(
        "--genomic-control",
        action="store_true",
        default=False,
        help="Enable genomic control correction (GENOMICCONTROL ON). [default: %(default)s]"
    )

    scheme_opts.add_argument(
        "--overlap-correction",
        dest="overlap_correction",
        action="store_true",
        default=False,
        help=(
            "Enable sample overlap correction (OVERLAP ON).\n"
            "Valid only when --scheme SAMPLESIZE is used.\n"
            "[default: %(default)s]"
        )
    )

    scheme_opts.add_argument(
        "--zcutoff",
        type=float,
        default=None,
        metavar="FLOAT",
        help="ZCUTOFF threshold used for overlap estimation. [default: %(default)s]"
    )

    # =====================================
    # Quality & Debug Options
    # =====================================

    qc_opts = pipe.add_argument_group(title="Quality & Debug Options")

    qc_opts.add_argument(
        "--column-counting",
        default="STRICT",
        choices=["STRICT", "LENIENT"],
        metavar="MODE",
        help=(
            "COLUMNCOUNTING mode.\n"
            "  STRICT   → exact column match (recommended).\n"
            "  LENIENT  → allow variable column counts.\n"
            "[default: %(default)s]"
        )
    )

    qc_opts.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Enable VERBOSE ON (large output; useful for debugging alignment issues). "
             "[default: %(default)s]"
    )

    qc_opts.add_argument(
        "--track-freq",
        action="store_true",
        default=False,
        help=(
            "Track allele frequency across studies "
            "(AVERAGEFREQ + MINMAXFREQ in METAL).\n"
            "[default: %(default)s]"
        )
    )

    return pipe





def parse_mtag_pipeline_args(subparser):
    """
    Defines command-line arguments for the MTAG Pipeline mode using standard argparse.
    """
    description = (
        "mtag-pipeline: Automated End-to-End Analysis\n\n"
        "1. Reads GWAS summary statistics VCFs listed in a manifest.\n"
        "2. Filters to common variants (inner join) and biallelic SNPs.\n"
        "3. Converts and formats data into MTAG-ready TSVs.\n"
        "4. Runs MTAG.\n"
    )

    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis mtag pipeline \\\n"
        "    --inputfile input_manifest.tsv \\\n"
        "    --run_name heart_disease_study \\\n"
        "    --ld_ref_panel ./eur_ld_ref/ \\\n"
        "    --cores 4\n\n"
        "For detailed explanation of MTAG options, visit:\n"
        "  https://github.com/JonJala/mtag/wiki/Tutorial-1:-The-Basics\n"
    )

    pipe = subparser.add_parser(
        "pipeline",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # =========================
    # Input Options
    # =========================
    in_opts = pipe.add_argument_group(title="Input Options")

    in_opts.add_argument(
        "--inputfile",
        required=True,
        metavar="FILE",
        help=(
            "TSV manifest with columns:\n"
            "  sumstat_vcf, TYPE, SPREV, PPREV, sample_id\n\n"
            "sumstat_vcf : Path to the GWAS summary statistics VCF file.\n"
            "TYPE        : Trait type (binary | quantitative).\n"
            "SPREV       : Sample prevalence (cases / total). For binary traits, if you are using\n"
            "              effective sample size, SPREV is typically set to 0.5. For quantitative\n"
            "              traits, this field can be empty.\n"
            "PPREV       : Population prevalence (binary traits only). For quantitative traits, this\n"
            "              field can be empty. If liability-scale conversion is not required, this\n"
            "              field may also be left empty.\n"
            "sample_id   : Unique trait identifier. Make sure sample_id matches the sample name in\n"
            "              the VCF header.\n"
        ),
    )

    in_opts.add_argument(
        "--run_name",
        type=str,
        required=True,
        metavar="NAME",
        help="Prefix for naming intermediate files, output files, and logs.",
    )

    in_opts.add_argument(
        "--out",
        metavar="DIR",
        default="./pipeline_results",
        help="Output directory. [default: %(default)s]",
    )

    # =========================
    # Genetic / MTAG Parameters
    # =========================
    gen_opts = pipe.add_argument_group(title="Genetic Parameters")

    gen_opts.add_argument(
        "--ld_ref_panel",
        metavar="PATH",
        required=True,
        type=str,
        help="Path to the LD reference panel directory/files required by MTAG.",
    )

    gen_opts.add_argument(
        "--no_overlap",
        action="store_true",
        help=(
            "Assume there is no sample overlap between any pair of GWAS studies. "
            "This sets the off-diagonal elements of the MTAG residual covariance matrix (Sigma) to 0."
        ),
    )

    gen_opts.add_argument(
        "--perfect_gencov",
        action="store_true",
        help=(
            "Assume the GWAS are different measures of the same underlying trait (genetic correlation = 1). "
            "Useful when traits are essentially the same phenotype measured differently (possibly with different "
            "measurement error)."
        ),
    )

    gen_opts.add_argument(
        "--equal_h2",
        action="store_true",
        help=(
            "Assume all GWAS have the same SNP heritability (equal h2). "
            "Typically appropriate only when the GWAS are on the same underlying trait."
        ),
    )

    # =========================
    # Execution Flags
    # =========================
    flag_opts = pipe.add_argument_group(title="MTAG Execution Flags")

    flag_opts.add_argument(
        "--fdr",
        action="store_true",
        help="Compute the approximate upper bound on the FDR (maxFDR) for MTAG results.",
    )

    flag_opts.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging (more detailed MTAG/pipeline output).",
    )

    flag_opts.add_argument(
        "--nthreads",
        type=int,
        default=max(1, mp.cpu_count() // 2),
        metavar="INT",
        help="Number of parallel threads to use. [default: %(default)s]",
    )

    flag_opts.add_argument(
        "--chunksize",
        type=int,
        default=10000,
        metavar="INT",
        help="Number of variants per MTAG chunk. [default: %(default)s]",
    )

    flag_opts.add_argument(
        "--std_betas",
        action="store_true",
        default=False,
        help=(
            "Output standardized effect sizes (beta) from MTAG. "
            "If set, MTAG will not rescale betas by 1/sqrt(2*MAF*(1-MAF))."
        ),
    )

    return pipe
    





def parse_pleio_pipeline_args(subparser):
    """
    Defines command-line arguments for the PLEIO Pipeline mode.
    """

    description = (
        "pleio pipeline: Automated End-to-End Analysis\n\n"
        "Steps performed automatically:\n"
        "1. Merge multiple GWAS summary statistics VCF files.\n"
        "2. Filter for common variants (inner join) and biallelic SNPs.\n"
        "3. Format merged data into PLEIO-ready TSV files.\n"
        "4. Execute the PLEIO multi-trait meta-analysis.\n"
    )

    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis pleio pipeline \\\n"
        "      --inputfile inputfile.tsv \\\n"
        "      --run_name psych_study \\\n"
        "      --ld_ref_panel ./eur_w_ld_chr/ \\\n"
        "      --cores 4\n"
        "For detailed explanation of pleio options, visit:\n"
        "  https://github.com/cuelee/pleio/wiki\n"

    )

    pipe = subparser.add_parser(
        "pipeline",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter
    )

    # -----------------------------
    # Input Options
    # -----------------------------
    in_opts = pipe.add_argument_group(title='Input Options')

    in_opts.add_argument(
        "--inputfile",
        required=True,
        metavar="FILE",
        help=(
            "TSV manifest containing the following columns:\n\n"
            "  sumstat_vcf : Path to the GWAS summary statistics VCF file\n"
            "  TYPE        : Trait type (binary | quantitative)\n"
            "  SPREV       : Sample prevalence (cases / total).\n"
            "                For binary traits using effective sample size,\n"
            "                SPREV is typically set to 0.5.\n"
            "                For quantitative traits this field can be empty.\n"
            "  PPREV       : Population prevalence (binary traits only).\n"
            "                Leave empty for quantitative traits or when\n"
            "                liability-scale conversion is not required.\n"
            "  sample_id   : Unique trait identifier. Must match the sample\n"
            "                name in the VCF header.\n"
        )
    )

    in_opts.add_argument(
        "--run_name",
        required=True,
        type=str,
        metavar="NAME",
        help="Prefix used for naming intermediate files and logs."
    )

    # -----------------------------
    # Workflow Automation
    # -----------------------------
    flow_opts = pipe.add_argument_group(title='Workflow Automation')

    flow_opts.add_argument(
        "--out",
        metavar="DIR",
        default="./pipeline_results",
        help="Output directory for pipeline results (default: %(default)s)."
    )

    flow_opts.add_argument(
        "--cores",
        type=int,
        default=max(1, mp.cpu_count() // 2),
        metavar="INT",
        help="Number of CPU threads used for parallel processing (default: %(default)s)."
    )

    # -----------------------------
    # Genetic Parameters
    # -----------------------------
    gen_opts = pipe.add_argument_group(title='Genetic Parameters')

    gen_opts.add_argument(
        "--ld_ref_panel",
        metavar="PATH",
        required=True,
        type=str,
        help="Path to LD reference panel directory (e.g. eur_w_ld_chr/)."
    )

    gen_opts.add_argument(
        "--nis",
        type=int,
        default=100000,
        metavar="INT",
        help="Number of samples used in importance sampling for estimating "
             "the null distribution (default: %(default)s)."
    )

    gen_opts.add_argument(
        "--flattening_p_values",
        default=False,
        action="store_true",
        help=(
            "Apply p-value flattening to recalibrate the PLEIO p-value "
            "distribution so that null p-values follow a uniform distribution. "
            "Useful for evaluating genomic inflation (λGC) or producing QQ plots. "
            "This option modifies the reported p-values and is generally "
            "not recommended for discovery analyses."
        )
    )

    return pipe





def parse_genomicpca_pipeline_args(subparser):
    """
    Defines command-line arguments for the GenomicPCA pipeline.
    """

    description = (
        "GenomicPCA Pipeline: Multivariate GWAS Meta-Analysis\n\n"
        "Pipeline Steps:\n"
        "1. Munge GWAS summary statistics using GenomicSEM.\n"
        "2. Estimate genetic covariance/correlation matrices using LDSC.\n"
        "3. Perform GenomicPCA decomposition of the LDSC matrix.\n"
        "4. Run multivariate GWAMA using PCA-derived loadings.\n"
    )

    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis genomicpca pipeline \\\n"
        "    --inputfile gwas_manifest.tsv \\\n"
        "    --run_name psych_genomicpca \\\n"
        "    --out ./genomicpca_results \\\n"
        "    --hm3 w_hm3.snplist \\\n"
        "    --ld_ref eur_w_ld_chr/ \\\n"
        "    --cores 8 \\\n"
        "    --approach correlation\n"
        "For detailed explanation of pleio options, visit:\n"
        "  https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html\n"
    )

    pipe = subparser.add_parser(
        "pipeline",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter
    )

    # =========================================================
    # Input Options
    # =========================================================
    in_opts = pipe.add_argument_group(title="Input Options")

    in_opts.add_argument(
        "--inputfile",
        required=True,
        metavar="FILE",
        help=(
            "TSV manifest containing the following columns:\n\n"
            "  sumstat_vcf : Path to the GWAS summary statistics VCF file\n"
            "  TYPE        : Trait type (binary | quantitative)\n"
            "  SPREV       : Sample prevalence (cases / total).\n"
            "                For binary traits using effective sample size,\n"
            "                SPREV is typically set to 0.5.\n"
            "                For quantitative traits this field can be empty.\n"
            "  PPREV       : Population prevalence (binary traits only).\n"
            "                Leave empty for quantitative traits or when\n"
            "                liability-scale conversion is not required.\n"
            "  sample_id   : Unique trait identifier. Must match the sample\n"
            "                name in the VCF header.\n"
        )
    )

    in_opts.add_argument(
        "--run_name",
        required=True,
        metavar="NAME",
        help="Prefix used for naming intermediate files, output files, and logs."
    )

    in_opts.add_argument(
        "--cores",
        type=int,
        default=max(1, mp.cpu_count() // 2),
        metavar="INT",
        help="Number of parallel CPU cores to use (default: %(default)s)."
    )

    in_opts.add_argument(
        "--hm3",
        required=True,
        metavar="FILE",
        help=(
            "Path to HapMap3 SNP list used during the GenomicSEM munge stage.\n"
            "Example: w_hm3.snplist"
        )
    )

    in_opts.add_argument(
        "--ld_ref",
        required=True,
        metavar="DIR",
        help=(
            "Directory containing the LD Score Regression reference panel.\n\n"
            "Example structure:\n"
            "  eur_w_ld_chr/\n"
            "      chr1.l2.ldscore.gz\n"
            "      chr2.l2.ldscore.gz\n"
        )
    )

    # =========================================================
    # Output Options
    # =========================================================
    out_opts = pipe.add_argument_group(title="Output Options")

    out_opts.add_argument(
        "--out",
        default="./genomicpca_results",
        metavar="DIR",
        help="Output directory for pipeline results (default: %(default)s)."
    )

    # =========================================================
    # Quality Control Options
    # =========================================================
    qc_opts = pipe.add_argument_group(title="Quality Control Options")

    qc_opts.add_argument(
        "--info_filter",
        type=float,
        default=0.7,
        metavar="FLOAT",
        help=(
            "INFO score threshold for SNP filtering applied during preprocessing "
            "(default: %(default)s)."
        )
    )

    qc_opts.add_argument(
        "--maf_filter",
        type=float,
        default=0.01,
        metavar="FLOAT",
        help=(
            "Minor allele frequency threshold applied during preprocessing "
            "(default: %(default)s)."
        )
    )

    # =========================================================
    # GenomicPCA Options
    # =========================================================
    gpca_opts = pipe.add_argument_group(title="GenomicPCA Options")

    gpca_opts.add_argument(
        "--approach",
        default="both",
        choices=["correlation", "covariance", "both"],
        metavar="TYPE",
        help=(
            "GenomicPCA method:\n"
            "  correlation → PCA on LDSC genetic correlation matrix (recommended)\n"
            "  covariance  → PCA on LDSC genetic covariance matrix\n"
            "  both        → run both analyses\n"
            "(default: %(default)s)"
        )
    )

    return pipe





def parse_placo_pipeline_args(subparser):
    """
    Defines command-line arguments for the PLACO pipeline mode.
    """

    description = (
        "PLACO Pipeline: Variant-level Pleiotropy Analysis\n\n"
        "Pipeline Steps:\n"
        "1. Reads GWAS summary statistic VCF files from a manifest.\n"
        "2. Harmonizes variants across the two traits.\n"
        "3. Computes Z-scores and p-values for each variant.\n"
        "4. Estimates variance parameters (VarZ) and trait correlation (CorZ).\n"
        "5. Performs pleiotropic association testing using PLACO or PLACO+.\n\n"
        "PLACO tests whether a genetic variant is associated with BOTH traits\n"
        "under a composite null hypothesis that at most one trait is associated.\n"
    )

    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis placo pipeline \\\n"
        "    --inputfile input_manifest.tsv \\\n"
        "    --run_name t2d_prostate \\\n"
        "    --out ./placo_results \\\n"
        "    --method placo.plus\n\n"
        "For detailed explanation of PLACO options, visit:\n"
        "  https://github.com/RayDebashree/PLACO\n"
    )

    pipe = subparser.add_parser(
        "pipeline",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter
    )

    # =====================================
    # Input Options
    # =====================================
    in_opts = pipe.add_argument_group(title="Input Options")

    in_opts.add_argument(
        "--inputfile",
        required=True,
        metavar="FILE",
        help=(
            "TSV manifest containing the following columns:\n\n"
            "  sumstat_vcf : Path to the GWAS summary statistics VCF file\n"
            "  TYPE        : Trait type (binary | quantitative)\n"
            "  SPREV       : Sample prevalence (cases / total).\n"
            "                For binary traits using effective sample size,\n"
            "                SPREV is typically set to 0.5.\n"
            "                For quantitative traits this field can be empty.\n"
            "  PPREV       : Population prevalence (binary traits only).\n"
            "                Leave empty for quantitative traits or when\n"
            "                liability-scale conversion is not required.\n"
            "  sample_id   : Unique trait identifier. Must match the sample\n"
            "                name in the VCF header.\n"
        )
    )

    in_opts.add_argument(
        "--run_name",
        required=True,
        metavar="NAME",
        help="Prefix used for naming intermediate files, output files, and logs."
    )

    in_opts.add_argument(
        "--cores",
        type=int,
        default=max(1, mp.cpu_count() // 2),
        metavar="INT",
        help="Number of parallel CPU cores to use. [default: %(default)s]"
    )

    # =====================================
    # Output Options
    # =====================================
    out_opts = pipe.add_argument_group(title="Output Options")

    out_opts.add_argument(
        "--out",
        default="./placo_results",
        metavar="DIR",
        help="Output directory for PLACO results. [default: %(default)s]"
    )

    # =====================================
    # PLACO Analysis Options
    # =====================================
    scheme_opts = pipe.add_argument_group(title="PLACO Analysis Options")

    scheme_opts.add_argument(
        "--method",
        default="placo.plus",
        choices=["placo.plus", "placo"],
        metavar="TYPE",
        help=(
            "PLACO testing method.\n"
            "  placo.plus  → General method for ANY two traits (recommended).\n"
            "                 Accounts for correlation between traits and\n"
            "                 potential sample overlap.\n\n"
            "  placo       → Original PLACO method assuming the two traits\n"
            "                 are independent (no correlation).\n\n"
            "[default: %(default)s]"
        )
    )

    pipe.add_argument(
        "--pthreshold",
        type=float,
        default=1e-4,
        metavar="FLOAT",
        help=(
            "P-value threshold used when estimating null parameters "
            "(variance of Z-scores and correlation between traits).\n"
            "Variants with p-values below this threshold are excluded\n"
            "from the null estimation step.\n"
            "[default: %(default)s]"
        )
    )

    return pipe
