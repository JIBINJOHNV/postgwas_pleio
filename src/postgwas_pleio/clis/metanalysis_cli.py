import argparse
import sys
from rich_argparse import RichHelpFormatter
import argparse
import textwrap
from rich_argparse import RichHelpFormatter
import argparse
import multiprocessing as mp
import psutil 

AUTO_THREADS = mp.cpu_count()
AUTO_MAX_MEM = f"{int(psutil.virtual_memory().total / 1024**3)}G"


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
        "For detailed explanation of genomicpca options, visit:\n"
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


def parse_placo_direct_args(subparser):
    """
    Defines command-line arguments for the METAL pipeline mode.
    """
    description = (
        "PLCAO Direct: Automated Meta-Analysis Workflow\n\n"
        "1. Reads validated GWAS summary statistic VCFs from manifest.\n"
        "2. Converts and harmonizes to METAL-ready format.\n"
        "3. Generates METAL script with consistent column mapping.\n"
        "4. Performs meta-analysis using METAL.\n"
    )
    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis metal pipeline \\\n"
        "    --inputfile input_manifest.tsv \\\n"
        "    --run_name psych_meta \\\n"
        "    --out ./meta_results \\\n"
        "    --scheme STDERR \\\n"
    )
    pipe = subparser.add_parser(
        "direct",
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
            "sumstat_vcf : Path to GWAS sumstat VCF file\n"
            "TYPE        : binary | quantitative \n"
            "SPREV       : Sample prevalence (0.5 for binary effective N);  This can be empty\n"
            "PPREV       : Population prevalence (binary traits only) ;  This can be empty \n"
            "sample_id   : Unique trait identifier:; Make sure the sample_id matches the sample name inside the VCF file.”\n" 
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
                          "SPREV can set empty, in that case it will try to take mean of prevalence using cases and controls. "
                          "For quantitative traits, this field can be empty.\n"
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


    harmonisation_group = pipe.add_argument_group(
        title="Harmonisation Options"
    )

    harmonisation_group.add_argument(
        "--harmonise",
        action="store_true",
        default=False,
        help=(
            "Perform harmonisation of the METAL output and generate VCF files. "
            "If enabled, you must provide the PostGWAS resource folder and a "
            "defaults configuration file using '--defaults_config'. "
            "[default: %(default)s]"
        )
    )

    harmonisation_group.add_argument(
        "--defaults_config",
        metavar="FILE",
        help=(
            "Full path to the harmonisation defaults YAML configuration file. "
            "The file should contain paths to reference resources and population "
            "settings used by the PostGWAS harmonisation step. "
            "Example configuration:\n"
            "https://github.com/JIBINJOHNV/postgwas/blob/main/tests/harmonisation.yaml"
        )
    )

    harmonisation_group.add_argument(
        "--resource-folder",
        metavar="FILE",
        help=("Postgwas harmonisation resource/reference folder, absolute path require")
    )

    harmonisation_group.add_argument(
        "--sample_size_approach",
        choices=["weight", "sample_overlap_corrected", "totalnef"], 
        default="totalnef",
        help=(
            "Which sample size estimate to report in the harmonised METAL output.\n"
            "weight : METAL internal weight (sum of study weights).\n"
            "sample_overlap_corrected : effective sample size after overlap correction (METAL 'N').\n"
            "totalnef : total effective sample size summed across studies (recommended).\n"
            "If --scheme STDERR is used, only totalnef is used.\n"
            "[default: %(default)s]"
        )
    )

    harmonisation_group.add_argument(
        "--info_method",
        choices=["mean", "median", "min", "max"],
        default="mean",
        help=(
            "Method used to aggregate INFO scores across studies when creating the "
            "harmonised METAL output.\n"
            "mean   : average INFO score across studies.\n"
            "median : median INFO score across studies.\n"
            "min    : minimum INFO score across studies (most conservative).\n"
            "max    : maximum INFO score across studies.\n"
            "[default: %(default)s]"
        )
    )


    # Create a dedicated group to avoid showing under "Optional Arguments"
    resource_group = pipe.add_argument_group("SYSTEM RESOURCES")

    # -------- THREADS ----------
    resource_group.add_argument(
        "--nthreads",
        type=int,
        metavar='',
        default=AUTO_THREADS,
        help="Number of parallel threads to use. [default: %(default)s]"
    )

    # -------- MEMORY ----------
    resource_group.add_argument(
        "--max-mem",
        metavar='',
        default=AUTO_MAX_MEM,
        help=(
            f"Maximum memory allowed. "
            "Formats accepted: 4G, 800M, 1200M."
            f"[bold green]Default:[/bold green] [cyan]{AUTO_MAX_MEM}[/cyan]"
        )
    )

    return pipe




def parse_asset_direct_args(subparsers):
    parser = subparsers.add_parser(
        "direct",
        prog="postgwas-pleio meta-analysis asset direct",
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent("""
        ASSET MULTI-TRAIT PLEIOTROPY PIPELINE
        =====================================

        Pipeline stages:
          1) Convert summary statistics into LDSC format using GenomicSEM::munge
          2) LDSC analysis (genetic correlation and intercept estimation)
          3) fastASSET analysis (fastASSET::fast_asset)
          4) Full ASSET analysis (ASSET::h.traits)

        ────────────────────────────────────────
        INPUT FILE REQUIREMENTS
        ────────────────────────────────────────

        1) LDSC MANIFEST (--input_manifest)

          Each row represents ONE trait.

          Required columns:

            FILE   → Path to GWAS summary statistics
            NAME   → Trait name (must be unique)
            SPREV  → Sample prevalence (case fraction)
            PPREV  → Population prevalence

          Required summary statistic columns:
            • SNP
            • A1
            • A2
            • Z
            • N

          Optional columns:
            INFO, P, BETA, SE

        2) fastASSET SNP INPUT (--fastasset_input)

          Required:
            • ID column
            • Beta column per trait
            • SE column per trait
            • Sample size column per trait

          Naming rule:
            TraitName.Beta
            TraitName.SE
            TraitName.N

        ────────────────────────────────────────
        OUTPUT STRUCTURE
        ────────────────────────────────────────

          3_munge_output/
          4_ldsc_output/
          5_fastasset_output/
          6_htraits_output/
          environment metadata
        """)
    )

    # ===============================
    # REQUIRED
    # ===============================
    parser.add_argument(
        "--input_manifest",
        required=True, metavar="MANIFEST",
        help="Trait manifest describing GWAS inputs"
    )

    parser.add_argument(
        "--fastasset_input",
        required=True, metavar="TABLE",
        help="SNP effect-size matrix"
    )

    parser.add_argument(
        "--output_dir",
        required=True, metavar="DIR",
        help="Pipeline output directory"
    )

    parser.add_argument(
        "--run_name",
        required=True, metavar="NAME",
        help="Prefix for naming intermediate files and logs."
    )

    parser.add_argument(
        "--hm3",
        required=True, metavar="SNPLIST",
        help="HapMap3 SNP list"
    )

    parser.add_argument(
        "--ld_ref_panel",
        required=True, metavar="DIR",
        help="LDSC LD reference directory"
    )

    # ===============================
    # OPTIONAL
    # ===============================
    parser.add_argument(
        "--info_filter",
        type=float, default=0.9, metavar="FLOAT",
        help="Munge INFO filter threshold [default: %(default)s]"
    )

    parser.add_argument(
        "--maf_filter",
        type=float, default=0.01, metavar="FLOAT",
        help="Munge MAF filter threshold [default: %(default)s]"
    )

    parser.add_argument(
        "--chunk_size",
        type=int, default=100000, metavar="INT",
        help="Chunk size for SNP processing [default: %(default)s]"
    )

    parser.add_argument(
        "--scr_pthr",
        type=float, default=0.999999, metavar="FLOAT",
        help="fastASSET screening threshold [default: %(default)s]"
    )

    parser.add_argument(
        "--meth_pval",
        default="DLM", choices=["DLM", "IS", "B"],
        metavar="METHOD",
        help="ASSET p-value method [choices: %(choices)s] (default: %(default)s)"
    )

    parser.add_argument(
        "--cores",
        type=int,
        default=max(1, mp.cpu_count() // 2),
        metavar="INT",
        help="Parallel cores [default: %(default)s]"
    )

    return parser


def parse_asset_pipeline_args(subparser):
    """
    Defines command-line arguments for the asset Pipeline mode.
    """
    description = (
        "asset pipeline: Automated End-to-End Analysis\n\n"
        "1. Merges multiple VCF summary statistics.\n"
        "2. Filters for common variants (inner join) and biallelic SNPs.\n"
        "3. Automatically splits and formats data into asset-ready TSV files.\n"
        "4. LDSC analysis (genetic correlation and intercept estimation)\n"
        "5. fastASSET analysis (fastASSET::fast_asset)\n"
        "6. Full ASSET analysis (ASSET::h.traits)\n"
    )
    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis asset pipeline \\\n"
        "    --inputfile inputfile.tsv \\\n"
        "    --run_name psych_study \\\n"
        "    --ld_ref_panel ./eur_w_ld_chr/ \\\n"
        "    --cores 4"
    )

    pipe = subparser.add_parser(
        "pipeline",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter
    )

    # --- Input Options ---
    in_opts = pipe.add_argument_group(title='Input Options')
    in_opts.add_argument("--inputfile", required=True, metavar="FILE",
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
        ))
    in_opts.add_argument("--run_name", type=str, required=True, metavar="NAME",
                        help='Prefix for naming intermediate files and logs.')

    # --- Workflow Automation ---
    flow_opts = pipe.add_argument_group(title='Workflow Automation')
    flow_opts.add_argument("--out", metavar='DIR', default='./pipeline_results', 
                          help='Output directory (default: %(default)s).')

    # --- Genetic Context ---
    gen_opts = pipe.add_argument_group(title='Genetic Parameters')
    gen_opts.add_argument('--ld_ref_panel', metavar="PATH", required=True, type=str, 
                         help='Path to LD reference panel (e.g., eur_w_ld_chr/).')

    gen_opts.add_argument("--cores", type=int, default=max(1, mp.cpu_count() // 2), metavar="INT",
                          help='Number of CPU threads to use (default: %(default)s).')

    gen_opts.add_argument("--hm3", required=True, metavar="SNPLIST",
                          help="HapMap3 SNP list")

    # ===============================
    # OPTIONAL
    # ===============================

    gen_opts.add_argument("--info_filter", type=float, default=0.9, metavar="FLOAT",
                          help="Munge INFO filter threshold [default: %(default)s]")

    gen_opts.add_argument("--maf_filter", type=float, default=0.01, metavar="FLOAT",
                          help="Munge MAF filter threshold [default: %(default)s]")

    gen_opts.add_argument("--chunk_size", type=int, default=100000, metavar="INT",
                          help="Chunk size for SNP processing [default: %(default)s]")

    gen_opts.add_argument("--scr_pthr", type=float, default=0.999999, metavar="FLOAT",
                          help="fastASSET screening threshold [default: %(default)s]")

    gen_opts.add_argument("--meth_pval", default="DLM", metavar="METHOD",
                          choices=["DLM", "IS", "B"],
                          help="ASSET p-value method (DLM|IS|B) [default: %(default)s]; choices : DLM, IS, B : https://rdrr.io/bioc/ASSET/man/h.traits.html")

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
        "      --nthreads 4\n"
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
        "--nthreads",
        type=int,
        default=max(1, mp.cpu_count() // 2),
        metavar="INT",
        help="Number of CPU threads used for parallel processing (default: %(default)s)."
    )

    flow_opts.add_argument(
        "--max-mem",
        metavar='',
        default=AUTO_MAX_MEM,
        help=(
            f"Maximum memory allowed. "
            "Formats accepted: 4G, 800M, 1200M."
            f"[bold green]Default:[/bold green] [cyan]{AUTO_MAX_MEM}[/cyan]"
        )
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


    harmonisation_group = pipe.add_argument_group(title="Harmonisation Options")

    harmonisation_group.add_argument(
        "--harmonise",
        action="store_true",
        default=False,
        help=(
            "Perform harmonisation of the METAL output and generate VCF files. "
            "If enabled, you must provide the PostGWAS resource folder and a "
            "defaults configuration file using '--defaults_config'. "
            "[default: %(default)s]"
        )
    )

    harmonisation_group.add_argument(
        "--defaults_config",
        metavar="FILE",
        help=(
            "Full path to the harmonisation defaults YAML configuration file. "
            "The file should contain paths to reference resources and population "
            "settings used by the PostGWAS harmonisation step. "
            "Example configuration:\n"
            "https://github.com/JIBINJOHNV/postgwas/blob/main/tests/harmonisation.yaml"
        )
    )

    harmonisation_group.add_argument(
        "--resource-folder",
        metavar="FILE",
        help=("Postgwas harmonisation resource/reference folder, absolute path require")
    )

    harmonisation_group.add_argument(
        "--info_method",
        choices=["mean", "median", "min", "max"],
        default="mean",
        help=(
            "Method used to aggregate INFO scores across studies when creating the "
            "harmonised METAL output.\n"
            "mean   : average INFO score across studies.\n"
            "median : median INFO score across studies.\n"
            "min    : minimum INFO score across studies (most conservative).\n"
            "max    : maximum INFO score across studies.\n"
            "[default: %(default)s]"
        )
    )

    harmonisation_group.add_argument(
        "--n_method",
        choices=["mean", "median", "min", "max", "sum"],
        default="sum",
        help=(
            "Aggregation method for effective sample size (NEF) across studies when generating "
            "harmonised MTAG inputs. This option is used only when all of the following flags are enabled: "
            "--harmonise, --equal_h2, and --perfect_gencov.\n\n"
            "Available methods:\n"
            "  mean   : Average effective sample size across studies.\n"
            "  median : Median effective sample size across studies.\n"
            "  min    : Minimum effective sample size (most conservative).\n"
            "  max    : Maximum effective sample size.\n"
            "  sum    : Sum of effective sample sizes across studies.\n\n"
            "[default: %(default)s]"
        )
    )

    return pipe


def parse_pleio_direct_args(subparser):
    """
    Defines command-line arguments for the PLEIO Pipeline mode.
    """
    description = (
        "pleio direct: Direct PLEIO Analysis. It need following input files. \n\n"
        "1. meta input data created by ldsc_preprocess.py \n"
        "2. genetic covariance matrix created by ldsc_preprocess.py \n"
        "3. non-genetic correlation matrix created by ldsc_preprocess.py \n"
        "4. importance sampling files created by ldsc_preprocess.py \n"
    )
    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis pleio direct \\\n"
        "    --metain inputfile.tsv \\\n"
        "    --run_name psych_study \\\n"
        "    --ld_ref_panel ./eur_w_ld_chr/ \\\n"
        "    --cores 4"
    )

    direct = subparser.add_parser(
        "direct",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    direct.add_argument('--out', default='pleio',type=str,
        help='Path to the output. If --out is not set, PLEIO will use pleio as the default output directory.')
    direct.add_argument('--metain', default=None, type=str,
        help='input file: file prefix of the meta input data.')
    direct.add_argument('--sg', default=None, type=str,
        help='input file: file prefix of the genetic covariance matrix.')
    direct.add_argument('--ce', default=None, type=str,
        help='input file: file prefix of the non-genetic correlation matrix.')
    direct.add_argument('--isf', default=None, type=str,
        help='Filename Prefix for Estimated null-distribution cumulative densities. ')
    direct.add_argument('--create', default = False, action='store_true',
        help='If this flag is set, PLEIO will run importance sampling and create new isf. ')
    direct.add_argument('--nis', default = int(100000), type = int,
        help='Number of samples for importance sampling. ')
    direct.add_argument('--blup', default = False, action='store_true',
        help='If this flag is set, PLEIO will estimate Best Linear Unbiased Prediction (BLUP)'
        'and write [output].blup.gz')
    direct.add_argument('--flattening_p_values', default = False, action= 'store_true',
        help='Flattening p-value distribution. This is necessary if you want to check lambda GC')
    direct.add_argument('--parallel', default = False, action='store_true',
        help='If this flag is set, PLEIO will run parallel computing ')
    direct.add_argument('--cores', default = mp.cpu_count() - 2, type = int, 
        help='Number of cpu cores for parallel computing. If --cores is not set, PLEIO will use n_cores - 1 as the default.')
    direct.add_argument('--snp', default='SNP', type=str,
        help='Name of SNP column (if not a name that ldsc understands). NB: case insensitive.')



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
        "    --nthreads 4\n\n"
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
            "TSV file with columns:\n"
            "  sumstat_vcf, TYPE, SPREV, PPREV, sample_id\n\n"
            "sumstat_vcf : Path to the GWAS summary statistics VCF file.\n"
            "TYPE        : Trait type (binary | quantitative).\n"
            "SPREV       : Sample prevalence (cases / total samples). "
                          "SPREV can set empty, in that case it will try to take mean of prevalence using cases and controls. "
                          "For quantitative traits, this field can be empty.\n"
            "PPREV       : Population prevalence (binary traits only). "
                          "For quantitative traits, this field can be empty. "
                          "If liability-scale conversion is not required, this field may also be left empty.\n"
            "sample_id   : Unique trait identifier. Make sure the sample_id matches "
                           "the sample name inside the VCF file.\n"
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
            "--incld_ambig_snps", 
            default=False, 
            action="store_true", 
            help=("Include strand-ambiguous SNPs (A/T or C/G) in the MTAG analysis. "
                  "By default these SNPs are removed. Even when included, they are "
                  "excluded from the estimation of Omega and Sigma and are only used "
                  "in the final MTAG results."
                ),
            )

    flag_opts.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging (more detailed MTAG/pipeline output).",
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


    harmonisation_group = pipe.add_argument_group(title="Harmonisation Options")

    harmonisation_group.add_argument(
        "--harmonise",
        action="store_true",
        default=False,
        help=(
            "Perform harmonisation of the METAL output and generate VCF files. "
            "If enabled, you must provide the PostGWAS resource folder and a "
            "defaults configuration file using '--defaults_config'. "
            "[default: %(default)s]"
        )
    )

    harmonisation_group.add_argument(
        "--defaults_config",
        metavar="FILE",
        help=(
            "Full path to the harmonisation defaults YAML configuration file. "
            "The file should contain paths to reference resources and population "
            "settings used by the PostGWAS harmonisation step. "
            "Example configuration:\n"
            "https://github.com/JIBINJOHNV/postgwas/blob/main/tests/harmonisation.yaml"
        )
    )

    harmonisation_group.add_argument(
        "--resource-folder",
        metavar="FILE",
        help=("Postgwas harmonisation resource/reference folder, absolute path require")
    )

    harmonisation_group.add_argument(
        "--info_method",
        choices=["mean", "median", "min", "max"],
        default="mean",
        help=(
            "Method used to aggregate INFO scores across studies when creating the "
            "harmonised METAL output.\n"
            "mean   : average INFO score across studies.\n"
            "median : median INFO score across studies.\n"
            "min    : minimum INFO score across studies (most conservative).\n"
            "max    : maximum INFO score across studies.\n"
            "[default: %(default)s]"
        )
    )


    harmonisation_group.add_argument(
        "--n_method",
        choices=["mean", "median", "min", "max", "sum"],
        default="sum",
        help=(
            "Aggregation method for effective sample size (NEF) across studies when generating "
            "harmonised MTAG inputs. This option is used only when all of the following flags are enabled: "
            "--harmonise, --equal_h2, and --perfect_gencov.\n\n"
            "Available methods:\n"
            "  mean   : Average effective sample size across studies.\n"
            "  median : Median effective sample size across studies.\n"
            "  min    : Minimum effective sample size (most conservative).\n"
            "  max    : Maximum effective sample size.\n"
            "  sum    : Sum of effective sample sizes across studies.\n\n"
            "[default: %(default)s]"
        )
    )


    # Create a dedicated group to avoid showing under "Optional Arguments"
    resource_group = pipe.add_argument_group("SYSTEM RESOURCES")

    # -------- THREADS ----------
    resource_group.add_argument(
        "--nthreads",
        type=int,
        metavar='',
        default=AUTO_THREADS,
        help=f"Threads to use [bold green]Default:[/bold green] [cyan]{AUTO_THREADS}[/cyan]"
    )

    # -------- MEMORY ----------
    resource_group.add_argument(
        "--max-mem",
        metavar='',
        default=AUTO_MAX_MEM,
        help=(
            f"Maximum memory allowed. "
            "Formats accepted: 4G, 800M, 1200M."
            f"[bold green]Default:[/bold green] [cyan]{AUTO_MAX_MEM}[/cyan]"
        )
    )
    return pipe


def parse_mtag_direct_args(subparser):
    """
    Defines command-line arguments for MTAG Direct mode.
    Note: We now take 'subparser' as an argument to attach to the main chain.
    """
    # Create the 'direct' sub-command
    # RawDescriptionHelpFormatter preserves your manual newlines
    parser = subparser.add_parser(
        "direct",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Direct implementation of MTAG method",
        description=(
            "mtag: Multitrait Analysis of GWAS\n"
            "This program is the implementation of MTAG method described by Turley et. al.\n"
            "Requires the input of a comma-separated list of GWAS summary statistics with identical columns.\n"
            "It is recommended to pass the column names manually to the program using the options below.\n"
            "The implementation of MTAG makes use of the LD Score Regression (ldsc) for cleaning the data\n"
            "and estimating residual variance-covariance matrix, so the input must also be compatible\n"
            "./munge_sumstats.py command in the ldsc distribution included with mtag.\n"
            "The default estimation method for the genetic covariance matrix Omega is GMM\n"
            "Note below: any list of passed to the options below must be comma-separated without whitespace."
        )
    )

    # --- Input Files Group ---
    in_opts = parser.add_argument_group(title='Input Files', description="Input files to be used by MTAG.")
    in_opts.add_argument("--sumstats", metavar="{File1},{File2}...", type=str, required=False, 
                        help='Specify the list of summary statistics files. Separate by ",".')
    in_opts.add_argument("--gencov_path", metavar="FILE_PATH", default=None, action="store", 
                        help="Read in the genetic covariance matrix and skip estimation.")
    in_opts.add_argument("--residcov_path", metavar="FILE_PATH", default=None, action="store", 
                        help="Read in the residual covariance matrix and skip estimation.")

    # --- Output Formatting Group ---
    out_opts = parser.add_argument_group(title="Output formatting", description="Set the output directory and prefix.")
    out_opts.add_argument("--out", metavar='DIR/PREFIX', default='./mtag_results', type=str, 
                         help='Specify directory and name prefix. Default is ./mtag_results')
    out_opts.add_argument("--make_full_path", default=False, action="store_true", 
                         help="Make output path if it does not exist.")
    out_opts.add_argument("--meta_format", default=False, action="store_true", 
                         help="Creates a file of the union of SNPs across traits.")

    # --- Column Names Group ---
    input_formatting = parser.add_argument_group(title="Column names of input files", 
                                                description="Manually pass names of relevant summary statistics columns.")
    input_formatting.add_argument("--snp_name", default="snpid", type=str, help="Unique identifier for SNPs. Default 'snpid'.")
    input_formatting.add_argument("--z_name", default="z", help="Common name of Z scores column. Default 'z'.")
    input_formatting.add_argument("--beta_name", default="beta", help="Common name of beta coefficients column.")
    input_formatting.add_argument("--se_name", default="se", help="Common name of standard errors column. Default 'se'.")
    input_formatting.add_argument("--n_name", default="n", help="Common name of sample sizes column. Default 'n'.")
    input_formatting.add_argument("--n_value", default=None, metavar="N1, N2,...", type=str, 
                                 help="Comma separated sample size values for GWAS without N column.")
    input_formatting.add_argument('--eaf_name', default="freq", help="Common name of effect allele frequencies. Default 'freq'.")
    input_formatting.add_argument('--no_chr_data', default=False, action='store_true', help="Do not use chr/pos data.")
    input_formatting.add_argument('--chr_name', default='chr', type=str, help="Column name for chromosome. Default 'chr'.")
    input_formatting.add_argument('--bpos_name', default='bpos', type=str, help="Column name for base pair. Default 'bpos'.")
    input_formatting.add_argument('--a1_name', default='a1', type=str, help="Column name for effect allele. Default 'a1'.")
    input_formatting.add_argument('--a2_name', default='a2', type=str, help="Column name for non-effect allele. Default 'a2'.")
    input_formatting.add_argument('--p_name', default='p', type=str, help="Column name for p-value. Default 'p'.")

    # --- Filter Options Group ---
    filter_opts = parser.add_argument_group(title="Filter Options", description="Filtering for input summary statistics.")
    filter_opts.add_argument("--include", default=None, metavar="SNPLIST1,SNPLIST2..", type=str, help="Restrict to union of SNPs in lists.")
    filter_opts.add_argument("--exclude", "--excludeSNPs", default=None, metavar="SNPLIST1,SNPLIST2..", type=str, help="Exclude union of SNPs in lists.")
    filter_opts.add_argument('--only_chr', metavar="CHR_A,CHR_B..", default=None, type=str, help="Restrict MTAG to listed chromosomes.")
    filter_opts.add_argument("--homogNs_frac", default=None, type=str, metavar="FRAC", help="Filter by (N-Mode)/Mode < FRAC.")
    filter_opts.add_argument("--homogNs_dist", default=None, type=str, metavar="D", help="Filter by sample size distance from mode.")
    filter_opts.add_argument('--maf_min', default='0.01', type=str, help="Minimum MAF threshold. Default 0.01.")
    filter_opts.add_argument('--n_min', default=None, type=str, help="Minimum threshold for SNP sample size.")
    filter_opts.add_argument('--n_max', default=None, type=str, help="Maximum threshold for SNP sample size.")
    filter_opts.add_argument("--info_min", default=None, type=str, help="Minimum info score for filtering.")
    filter_opts.add_argument("--incld_ambig_snps", default=False, action="store_true", help="Include strand ambiguous SNPs.")
    filter_opts.add_argument("--no_allele_flipping", default=False, action="store_true", help="Prevents flipping effect sizes.")

    # --- Special Cases Group ---
    special_cases = parser.add_argument_group(title="Special Cases", description="Improve runtime based on specific assumptions.")
    special_cases.add_argument('--use_beta_se', default=False, action='store_true', help='Use BETA and SE instead of Z-statistic.')
    special_cases.add_argument('--no_overlap', default=False, action='store_true', help='Assume no sample overlap.')
    special_cases.add_argument('--perfect_gencov', default=False, action='store_true', help='Assume traits are perfectly genetically correlated.')
    special_cases.add_argument('--equal_h2', default=False, action='store_true', help='Assume all traits have equal heritability.')
    special_cases.add_argument('--force', default=False, action='store_true', help='Force estimation despite small mean chi2.')

    # --- FDR Options Group ---
    fdr_opts = parser.add_argument_group(title='Max FDR calculation', description="Calculation of upper bound on false discovery.")
    fdr_opts.add_argument('--fdr', default=False, action='store_true', help='Perform max FDR calculations.')
    fdr_opts.add_argument('--skip_mtag', default=False, action='store_true', help='Skip MTAG and perform FDR only.')
    fdr_opts.add_argument('--grid_file', default=None, action='store', help='Pre-set list of grid points.')
    fdr_opts.add_argument('--fit_ss', default=False, action='store_true', help='Restrict grid search using prior null estimates.')
    fdr_opts.add_argument('--intervals', default=10, type=int, help='Number of intervals to partition [0,1].')
    fdr_opts.add_argument('--cores', default=1, type=int, help='Number of threads/cores for FDR computation.')
    fdr_opts.add_argument('--p_sig', default=5.0e-8, type=float, help='Significance threshold. Default 5.0e-8.')
    fdr_opts.add_argument('--n_approx', default=True, action='store_true', help='Use mean sample size across SNPs for FDR speedup.')

    # --- Miscellaneous Group ---
    misc = parser.add_argument_group(title="Miscellaneous")
    misc.add_argument('--ld_ref_panel', default=None, action='store', metavar="FOLDER_PATH", type=str, help='Specify folder of LD reference panel.')
    misc.add_argument('--time_limit', default=100., type=float, action="store", help="Set time limit (hours) on variance covariance estimation.")
    misc.add_argument('--std_betas', default=False, action='store_true', help="Results files will have standardized effect sizes.")
    misc.add_argument("--tol", default=1e-6, type=float, help="Relative (x) tolerance for genetic VCV estimation.")
    misc.add_argument('--numerical_omega', default=False, action='store_true', help='Use MLE estimator of genetic VCV matrix.')
    misc.add_argument('--verbose', default=False, action='store_true', help='Include verbose output from ldsc and optimization.')
    misc.add_argument('--chunksize', default=int(1e7), type=int, help='Chunksize for reading in data.')
    misc.add_argument('--stream_stdout', default=False, action='store_true', help='Stream processing to console and log.')
    misc.add_argument('--median_z_cutoff', default=0.0, type=float, help='Maximum allowed median Z-score for sumstats QC.')

    return parser 