from pathlib import Path
from typing import List, Union, Optional

PathLike = Union[str, Path]


def generate_metal_script(
    input_files: List[PathLike],
    output_prefix: str,
    script_path: PathLike,
    *,
    scheme: str = "STDERR",
    genomic_control: bool = False,
    overlap_correction: bool = False,
    zcutoff: Optional[float] = None,
    verbose: bool = False,
    column_counting: str = "STRICT",
    track_freq: bool = True,
    heterogeneity: bool = True,
):
    """
    Generate a METAL script for fixed-format GWAS summary statistics.

    Expected columns:
        ID, CHROM, POS, REF, ALT, AF, ES, SE, NEF, LP, SI, EZ
    """

    scheme = scheme.upper()
    column_counting = column_counting.upper()

    if overlap_correction and scheme != "SAMPLESIZE":
        raise ValueError("OVERLAP correction only valid with SCHEME SAMPLESIZE")

    if column_counting not in {"STRICT", "LENIENT"}:
        raise ValueError("column_counting must be STRICT or LENIENT")

    input_files = [str(f) for f in input_files]
    output_prefix = str(output_prefix)
    script_path = Path(script_path)

    lines = []
    lines.append("#############################################")
    lines.append("# AUTO-GENERATED METAL SCRIPT")
    lines.append("#############################################\n")

    # -------------------------
    # Core settings
    # -------------------------

    lines.append(f"SCHEME {scheme}")

    if genomic_control:
        lines.append("GENOMICCONTROL ON")

    if verbose:
        lines.append("VERBOSE ON")

    lines.append(f"COLUMNCOUNTING {column_counting}")

    # -------------------------
    # Column definitions
    # -------------------------

    lines.append("\nMARKER ID")
    lines.append("ALLELE ALT REF")

    lines.append("EFFECT ES")
    lines.append("STDERR SE")
    lines.append("WEIGHT NEF") 
    lines.append("PVALUE LP") 

    if scheme == "SAMPLESIZE":
        if overlap_correction:
            lines.append("OVERLAP ON")
            if zcutoff is not None:
                lines.append(f"ZCUTOFF {zcutoff}")

    # -------------------------
    # Allele frequency tracking
    # -------------------------

    if track_freq:
        lines.append("FREQLABEL AF")
        lines.append("AVERAGEFREQ ON")
        lines.append("MINMAXFREQ ON")

    # -------------------------
    # Custom variable: Total NEF
    # -------------------------

    lines.append("\nCUSTOMVARIABLE TotalNEF")
    lines.append("LABEL TotalNEF as NEF")

    lines.append("")

    # -------------------------
    # Process input files
    # -------------------------

    for f in input_files:
        lines.append(f"PROCESS {f}")

    # -------------------------
    # Final analysis
    # -------------------------

    lines.append("")
    lines.append(f"OUTFILE {output_prefix}_ .tbl")
    lines.append("")
    
    if heterogeneity:
        lines.append("ANALYZE HETEROGENEITY")
        lines.append("")

    lines.append("QUIT")

    script_path.write_text("\n".join(lines))

    return script_path