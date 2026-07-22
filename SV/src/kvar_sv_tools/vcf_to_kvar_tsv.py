#!/usr/bin/env python3
"""Public CLI for converting an SV VCF submission to Variant_Call TSV."""

import argparse
import os
import sys

try:
    from .KVar2TSV import KVarTSVConverter
except ImportError:
    from KVar2TSV import KVarTSVConverter


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate an SV VCF file and convert it to a KVar Variant_Call TSV file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input SV VCF path (.gz supported)")
    parser.add_argument(
        "-f",
        "-r",
        "--reference-fasta",
        dest="reference_fasta",
        required=True,
        help="Reference FASTA path used for CHROM/POS/REF validation; .fai index is required",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        required=True,
        help="Metadata text file containing ##SampleSet_id, ##Experiment_id, ##reference, and optional ##organism_taxid lines",
    )
    parser.add_argument("-t", "--call-tsv", required=True, help="Output Variant_Call TSV path")
    parser.add_argument("-e", "--error-report", required=True, help="Validation report path")
    parser.add_argument(
        "-c",
        "--call-accession-start",
        type=_positive_int,
        required=True,
        metavar="N",
        help="Starting number for Variant Call accessions, written as kssvN",
    )
    parser.add_argument(
        "--sanitize-error-report",
        action="store_true",
        help="Redact absolute paths and raw row content from the validation report",
    )
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    for label, path in (
        ("VCF", args.vcf),
        ("reference FASTA", args.reference_fasta),
        ("metadata", args.metadata),
    ):
        if not os.path.exists(path):
            parser.error(f"{label} file not found: {path}")

    converter = KVarTSVConverter(reference_fasta_path=args.reference_fasta)
    converter.convert_vcf_to_tsv(
        args.vcf,
        args.call_tsv,
        error_report_path=args.error_report,
        metadata_file_path=args.metadata,
        call_accession_start=args.call_accession_start,
        sanitize_error_report=args.sanitize_error_report,
    )
    print(f"Variant_Call TSV: {args.call_tsv}")
    print(f"Validation report: {args.error_report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
