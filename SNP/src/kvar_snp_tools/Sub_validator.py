#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Public CLI for SNP VCF validation and dbSNP VCF creation."""

import argparse
import os
import sys
from typing import Optional

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from kvar_snp_tools.VCF2dbSNP import VCF2dbSNPConverter
    from kvar_snp_tools.dbsnp_vcf_cleaner import DbSNPVCFCleaner
else:
    from .VCF2dbSNP import VCF2dbSNPConverter
    from .dbsnp_vcf_cleaner import DbSNPVCFCleaner


def _run_reference_check(vcf_path: str, reference_path: Optional[str], report_path: Optional[str]) -> None:
    """Run reference validation when a reference FASTA path is provided."""
    if not reference_path:
        return

    if __package__ in (None, ""):
        from kvar_snp_tools.VCF_ref_check import VCFRefChecker
    else:
        from .VCF_ref_check import VCFRefChecker

    checker = VCFRefChecker()
    checker.check_vcf_against_fasta(
        vcf_path,
        reference_path,
        output_report_path=report_path,
    )


def convert_generic_to_dbsnp(args: argparse.Namespace) -> None:
    """Convert a generic VCF into a dbSNP-formatted VCF."""
    _run_reference_check(args.vcf, args.reference, args.reference_report)
    converter = VCF2dbSNPConverter()
    converter.convert_vcf_to_dbsnp(
        vcf_file_path=args.vcf,
        output_file_path=args.output,
        metadata_file_path=args.metadata,
        error_report_path=args.error_report,
    )


def validate_dbsnp(args: argparse.Namespace) -> None:
    """Validate and rewrite a cleaned dbSNP VCF."""
    _run_reference_check(args.vcf, args.reference, args.reference_report)
    cleaner = DbSNPVCFCleaner()
    cleaner.clean(
        vcf_file_path=args.vcf,
        output_file_path=args.output,
        metadata_file_path=args.metadata,
        error_report_path=args.error_report,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Validate SNP VCF input and create cleaned dbSNP VCF output"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    generic_parser = subparsers.add_parser(
        "generic-to-dbsnp",
        help="Convert a generic VCF into a dbSNP-formatted VCF",
    )
    _add_common_arguments(generic_parser)
    generic_parser.set_defaults(func=convert_generic_to_dbsnp)

    strict_parser = subparsers.add_parser(
        "validate-dbsnp",
        help="Validate and rewrite an input dbSNP VCF as cleaned dbSNP VCF",
    )
    _add_common_arguments(strict_parser)
    strict_parser.set_defaults(func=validate_dbsnp)

    return parser


def _add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared by both modes."""
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF path")
    parser.add_argument("-o", "--output", required=True, help="Output dbSNP VCF path")
    parser.add_argument("-m", "--metadata", help="Optional metadata file path")
    parser.add_argument("-e", "--error-report", help="Optional validation report path")
    parser.add_argument("-r", "--reference", help="Optional reference FASTA for REF allele validation")
    parser.add_argument(
        "-rr",
        "--reference-report",
        help="Optional reference validation report path when --reference is used",
    )


def main() -> None:
    """Run the public submission validation CLI."""
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
