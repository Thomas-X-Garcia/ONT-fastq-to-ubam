#!/usr/bin/env python3
"""
ONT-fastq-to-ubam: Convert ONT FASTQ to unmapped BAM while preserving metadata.

Author: Thomas-X-Garcia, PhD, HCLD
Repository: https://github.com/Thomas-X-Garcia/ONT-fastq-to-ubam
License: MIT
"""

import argparse
import json
import logging
import shlex
import sys
from pathlib import Path
from typing import Dict, Any, Tuple, Optional

try:
    import pysam
except ImportError:
    print("Error: pysam is not installed. Please run: pip install pysam", file=sys.stderr)
    sys.exit(1)

__version__ = "1.0.0"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Map ONT header keys -> custom two-letter SAM tags and types.
# Avoids collisions with standard tags (NM, AS, XS, SA, RG, etc.).
ONT_TO_SAM: Dict[str, Tuple[str, str]] = {
    "runid": ("XR", "Z"),                         # run id
    "ch": ("XH", "i"),                            # channel (int)
    "start_time": ("XT", "Z"),                    # ISO datetime
    "flow_cell_id": ("XF", "Z"),                  # flow cell id
    "basecall_gpu": ("XG", "Z"),                  # GPU model
    "protocol_group_id": ("XP", "Z"),             # protocol group
    "sample_id": ("XQ", "Z"),                     # original sample label
    "parent_read_id": ("XI", "Z"),                # parent read id
    "basecall_model_version_id": ("XV", "Z"),     # basecaller model
}

def cast_value(v: str) -> Any:
    """Convert string to appropriate numeric type if possible."""
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v

def parse_kv_blob(comment: str) -> Tuple[str, Dict[str, Any]]:
    """
    Parse ONT-style 'k=v' tokens from the FASTQ comment using shlex for robustness.
    Returns (raw_comment, dict_of_parsed).
    Non k=v tokens are ignored.
    """
    kv = {}
    if not comment:
        return comment, kv
    
    try:
        tokens = shlex.split(comment, posix=True)
    except ValueError as e:
        logger.warning(f"Failed to parse comment with shlex: {e}. Using simple split.")
        tokens = comment.split()
    
    for tok in tokens:
        if "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        kv[k] = cast_value(v)
    return comment, kv

def build_header(sample_name: str, description: Optional[str] = None) -> Dict[str, Any]:
    """Build BAM header with proper @HD and @RG lines."""
    header = {"HD": {"VN": "1.6", "SO": "unknown"}}
    
    rg_dict = {"ID": sample_name, "SM": sample_name, "PL": "ONT"}
    if description:
        rg_dict["DS"] = description
    
    header["RG"] = [rg_dict]
    header["PG"] = [{
        "ID": "ONT-fastq-to-ubam",
        "PN": "ONT-fastq-to-ubam",
        "VN": __version__,
        "CL": " ".join(sys.argv)
    }]
    
    return header

def validate_inputs(input_fastq: Path, output_bam: Path, sample_name: str) -> None:
    """Validate input parameters before processing."""
    if not input_fastq.exists():
        raise FileNotFoundError(f"Input FASTQ file does not exist: {input_fastq}")
    
    if not input_fastq.is_file():
        raise ValueError(f"Input path is not a file: {input_fastq}")
    
    valid_extensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
    if not str(input_fastq).lower().endswith(valid_extensions):
        logger.warning(f"Input file may not be a FASTQ file: {input_fastq}")
    
    if output_bam.exists():
        logger.warning(f"Output file already exists and will be overwritten: {output_bam}")
    
    output_dir = output_bam.parent
    if not output_dir.exists():
        raise FileNotFoundError(f"Output directory does not exist: {output_dir}")
    
    if not sample_name or not sample_name.strip():
        raise ValueError("Sample name cannot be empty")
    
    if not sample_name.replace('_', '').replace('-', '').replace('.', '').isalnum():
        logger.warning(f"Sample name contains special characters: {sample_name}")

def write_ubam(
    input_fastq: Path, 
    output_bam: Path, 
    sample_name: str, 
    emit_individual_tags: bool = True,
    description: Optional[str] = None,
    validate: bool = True
) -> int:
    """Convert FASTQ to unmapped BAM.
    
    Returns the number of reads processed.
    """
    if validate:
        validate_inputs(input_fastq, output_bam, sample_name)
    
    read_count = 0
    logger.info(f"Starting conversion: {input_fastq} -> {output_bam}")
    logger.info(f"Sample name: {sample_name}")
    
    try:
        with pysam.FastxFile(str(input_fastq)) as infile:
            with pysam.AlignmentFile(
                str(output_bam), 
                "wb", 
                header=build_header(sample_name, description)
            ) as outbam:
                for rec in infile:
                    a = pysam.AlignedSegment()
                    a.query_name = rec.name                               # up to first whitespace
                    a.query_sequence = rec.sequence
                    a.query_qualities = (
                        pysam.qualitystring_to_array(rec.quality) if rec.quality else None
                    )

                    # Unmapped template scaffold
                    a.flag = 4
                    a.reference_id = -1
                    a.reference_start = -1
                    a.mapping_quality = 0
                    a.next_reference_id = -1
                    a.next_reference_start = -1
                    a.template_length = 0

                    # Read group
                    a.set_tag("RG", sample_name, "Z")

                    raw, kv = parse_kv_blob(rec.comment or "")

                    # Store exact comment + parsed kv as JSON
                    a.set_tag("OX", json.dumps({"raw": raw, "kv": kv}, ensure_ascii=False), "Z")

                    # Also project selected keys into stable per-read tags
                    if emit_individual_tags and kv:
                        for k, (tag, typ) in ONT_TO_SAM.items():
                            if k in kv:
                                val = kv[k]
                                # Coerce type to SAM tag type
                                if typ == "i":
                                    try:
                                        val = int(val)
                                    except Exception:
                                        # fall back to string if not integer-like
                                        typ = "Z"
                                a.set_tag(tag, val, typ)

                    outbam.write(a)
                    read_count += 1
                    
                    if read_count % 10000 == 0:
                        logger.info(f"Processed {read_count:,} reads...")
    
    except Exception as e:
        logger.error(f"Error during conversion: {e}")
        if output_bam.exists():
            output_bam.unlink()
        raise
    
    logger.info(f"Successfully converted {read_count:,} reads")
    return read_count

def main() -> int:
    """Main entry point."""
    ap = argparse.ArgumentParser(
        description="Convert ONT FASTQ to unmapped BAM while preserving metadata.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.fastq output.bam --sample-name SAMPLE001
  %(prog)s reads.fq.gz reads.ubam --sample-name HG002 --description "ONT PromethION run"
  
For use with epi2me-labs/wf-human-variation:
  %(prog)s ont_reads.fastq.gz sample.unaligned.bam --sample-name SAMPLE_ID
  
The --sample-name value will be used by the workflow to name all output files.
"""
    )
    ap.add_argument(
        "fastq", 
        type=Path, 
        help="Input FASTQ file (.fastq, .fq, .fastq.gz, or .fq.gz)"
    )
    ap.add_argument(
        "out_bam", 
        type=Path, 
        help="Output unmapped BAM file path"
    )
    ap.add_argument(
        "--sample-name",
        required=True,
        help="Sample name for @RG:ID and @RG:SM (REQUIRED). This will be used by downstream workflows.",
    )
    ap.add_argument(
        "--description",
        help="Optional description for @RG:DS field",
    )
    ap.add_argument(
        "--no-individual-tags",
        action="store_true",
        help="Do not emit per-read custom tags for ONT fields; keep only OX JSON.",
    )
    ap.add_argument(
        "--no-validation",
        action="store_true",
        help="Skip input validation checks",
    )
    ap.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress messages",
    )
    ap.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}"
    )
    
    args = ap.parse_args()
    
    if args.quiet:
        logger.setLevel(logging.WARNING)

    try:
        read_count = write_ubam(
            args.fastq,
            args.out_bam,
            args.sample_name,
            emit_individual_tags=not args.no_individual_tags,
            description=args.description,
            validate=not args.no_validation
        )
        
        if not args.quiet:
            print(f"\nConversion complete:")
            print(f"  Input:  {args.fastq}")
            print(f"  Output: {args.out_bam}")
            print(f"  Reads:  {read_count:,}")
            print(f"  Sample: {args.sample_name}")
            print(f"\nVerify with: samtools view -H {args.out_bam} | head -20")
        
        return 0
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        return 1
    except ValueError as e:
        logger.error(f"Invalid input: {e}")
        return 1
    except BrokenPipeError:
        try:
            sys.stderr.close()
        except Exception:
            pass
        return 0
    except KeyboardInterrupt:
        logger.info("\nInterrupted by user")
        return 130
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())