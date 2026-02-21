import argparse
import sys
from rich_argparse import RichHelpFormatter
import argparse
import textwrap
from rich_argparse import RichHelpFormatter
import argparse
import multiprocessing as mp



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
        "--ncores",
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
            "A TSV file with five columns: sumstat_vcf, TYPE, SPREV, PPREV, sample_id.\n\n"
            "sumstat_vcf: Path to the input GWAS sumstat VCF file.\n"
            "TYPE:        'binary' or 'quantitative'.\n"
            "SPREV:       Sample prevalence. Should be 0.5 for binary traits\n"
            "             because we are using the effective sample size.\n"
            "PPREV:       Population prevalence. Required for binary traits.\n"
            "sample_id:   Unique trait identifier matching VCF samples.\n"
            "\n"
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

    gen_opts.add_argument("--ncores", type=int, default=max(1, mp.cpu_count() // 2), metavar="INT",
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
        "1. Merges multiple VCF summary statistics.\n"
        "2. Filters for common variants (inner join) and biallelic SNPs.\n"
        "3. Automatically splits and formats data into pleio-ready TSV files.\n"
        "4. Executes the PLEIO meta-analysis.\n"
    )
    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis pleio pipeline \\\n"
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
            "A TSV file with five columns: sumstat_vcf, TYPE, SPREV, PPREV, sample_id.\n\n"
            "sumstat_vcf: Path to the input GWAS sumstat VCF file.\n"
            "TYPE:        'binary' or 'quantitative'.\n"
            "SPREV:       Sample prevalence. Should be 0.5 for binary traits\n"
            "             because we are using the effective sample size.\n"
            "PPREV:       Population prevalence. Required for binary traits.\n"
            "sample_id:   Unique trait identifier matching VCF samples.\n"
            "\n"
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
    # Corrected from flag_opts to gen_opts (or flow_opts) 
    gen_opts.add_argument('--cores', default=1, type=int, 
                          help='Number of CPU threads to use (default: %(default)s).')
    gen_opts.add_argument('--nis', default = int(100000), type = int,
                        help='Number of samples for importance sampling. (default: %(default)s).')
    gen_opts.add_argument('--flattening_p_values', default = False, action= 'store_true',
        help='Flattening p-value distribution. This is necessary if you want to check lambda GC')
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
    direct.add_argument('--ncores', default = mp.cpu_count() - 2, type = int, 
        help='Number of cpu cores for parallel computing. If --ncores is not set, PLEIO will use n_cores - 1 as the default.')
    direct.add_argument('--snp', default='SNP', type=str,
        help='Name of SNP column (if not a name that ldsc understands). NB: case insensitive.')



# def parse_pleio_direct_args(subparser):
#     """
#     Defines command-line arguments for the PLEIO Pipeline mode.
#     """
#     description = (
#         "pleio-direct: Direct PLEIO Analysis\n\n"
#         "1. Merges multiple VCF summary statistics.\n"
#         "2. Filters for common variants (inner join) and biallelic SNPs.\n"
#         "3. Automatically splits and formats data into pleio-ready TSV files.\n"
#         "4. Executes the PLEIO meta-analysis.\n"
#     )
#     epilog = (
#         "Example Usage:\n"
#         "  postgwas-pleio meta-analysis pleio direct \\\n"
#         "    --inputfile inputfile.tsv \\\n"
#         "    --run_name psych_study \\\n"
#         "    --ld_ref_panel ./eur_w_ld_chr/ \\\n"
#         "    --cores 4"
#     )

#     pipe = subparser.add_parser(
#         "direct",
#         description=description,
#         epilog=epilog,
#         formatter_class=argparse.RawDescriptionHelpFormatter
#     )
#     # --- Input Options ---
#     in_opts = pipe.add_argument_group(title='Input Options')
#     in_opts.add_argument("--inputfile", required=True, metavar="FILE",
#         help=(
#             "A TSV file with five columns: sumstat_file, TYPE, SPREV, PPREV, sample_id.\n"
#             "sumstat_file: Path to the input GWAS sumstat file (post-harmonization).\n"
#             "TYPE: 'binary' or 'quantitative'.\n"
#             "SPREV: Sample prevalence (cases/total). Can be blank for quantitative.\n"
#             "       Note: PLEIO requires Effective Sample Size (Neff). If TYPE is 'binary'\n"
#             "       and SPREV is blank, it defaults to 0.5 to match Neff assumptions.\n"
#             "PPREV: Population prevalence. Required for binary traits; if empty, \n"
#             "       the pipeline will throw a warning and may fail ldsc_preprocess.\n"
#             "sample_id: Unique trait identifier for output files."
#         ))
#     in_opts.add_argument("--run_name", type=str, required=True, metavar="NAME",
#                         help='Prefix for naming intermediate files and logs.')
#     # --- Workflow Automation ---
#     flow_opts = pipe.add_argument_group(title='Workflow Automation')
#     flow_opts.add_argument("--out", metavar='DIR', default='./pipeline_results', 
#                           help='Output directory (default: %(default)s).')
#     # --- Genetic Context --- 
#     gen_opts = pipe.add_argument_group(title='Genetic Parameters')
#     gen_opts.add_argument('--ld_ref_panel', metavar="PATH", required=True, type=str, 
#                          help='Path to LD reference panel (e.g., eur_w_ld_chr/).') 
#     # Corrected from flag_opts to gen_opts (or flow_opts)
#     gen_opts.add_argument('--cores', default=1, type=int, 
#                           help='Number of CPU threads to use (default: %(default)s).')
#     return pipe




def parse_mtag_pipeline_args(subparser):
    """
    Defines command-line arguments for the MTAG Pipeline mode using standard argparse.
    """
    
    # Standard text description (Rich tags removed for standard argparse compatibility)
    description = (
        "mtag-pipeline: Automated End-to-End Analysis\n\n"
        "1. Merges multiple VCF summary statistics\n"
        "2. Filters for common variants (inner join) and biallelic SNPs.\n"
        "3. Automatically splits and formats data into MTAG-ready TSVs.\n"
        "4. Executes the MTAG meta analysis\n\n"
    )
    epilog = (
        "Example Usage:\n"
        "  postgwas-pleio meta-analysis mtag pipeline \n"
        "    --vcfs trait1.vcf.gz trait2.vcf.gz \n"
        "    --run_name heart_disease_study \n"
        "    --ld_ref_panel ./eur_ld_ref/ \n"
        "    --cores 4"
    )
    # Use the built-in RawDescriptionHelpFormatter to preserve the newlines in description/epilog
    pipe = subparser.add_parser(
        "pipeline",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # --- Input VCF Options ---
    in_opts = pipe.add_argument_group(title='VCF Input Options')
    in_opts.add_argument("--vcfs", nargs="+", required=True, metavar="FILE.vcf.gz",
                        help='List of raw VCF.gz files to process.')
    in_opts.add_argument("--run_name", type=str, required=True, metavar="NAME",
                        help='Prefix for naming intermediate files and logs.')
    # --- Workflow Automation ---
    flow_opts = pipe.add_argument_group(title='Workflow Automation')
    flow_opts.add_argument("--out", metavar='DIR', default='./pipeline_results', 
                          help='Output directory. (default: %(default)s)')
    # --- Genetic Context --- 
    gen_opts = pipe.add_argument_group(title='Genetic Parameters')
    gen_opts.add_argument('--no_overlap', action='store_true', help='Assume no sample overlap.')
    gen_opts.add_argument('--perfect_gencov', action='store_true', help='Assume genetic correlations = 1.0.')
    gen_opts.add_argument('--equal_h2', action='store_true', help='Assume equal heritability.')
    gen_opts.add_argument('--ld_ref_panel', metavar="PATH", required=True, type=str, 
                         help='Path to LD reference panel.') 
    # --- Execution Flags ---
    flag_opts = pipe.add_argument_group(title='MTAG Execution Flags')
    flag_opts.add_argument('--fdr', action='store_true', help='Perform max FDR calculations.')
    flag_opts.add_argument('--verbose', action='store_true', help='Show detailed logs.')
    flag_opts.add_argument('--cores', default=1, type=int, 
                          help='CPU threads. (default: %(default)s)')
    flag_opts.add_argument('--chunksize', default=10000, type=int, 
                          help='Variants per MTAG chunk. (default: %(default)s)')
    flag_opts.add_argument('--std_betas', default=False, action='store_true', help="Results files will have standardized effect sizes.")
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