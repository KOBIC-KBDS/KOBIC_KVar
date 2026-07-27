#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Public CLI for SNP VCF validation and dbSNP VCF creation."""

import argparse
import os
import sys
from typing import Optional

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from kvar_snp_tools.VCF2dbSNP import VCF2dbSNPConverter, default_error_report_path
    from kvar_snp_tools.dbsnp_vcf_cleaner import DbSNPVCFCleaner
    from kvar_snp_tools.error_handler import ErrorHandler
else:
    from .VCF2dbSNP import VCF2dbSNPConverter, default_error_report_path
    from .dbsnp_vcf_cleaner import DbSNPVCFCleaner
    from .error_handler import ErrorHandler


def _run_reference_check(
    vcf_path: str,
    reference_path: Optional[str],
    error_handler: ErrorHandler,
) -> None:
    """Run reference validation when a reference FASTA path is provided."""
    if not reference_path:
        return

    if __package__ in (None, ""):
        from kvar_snp_tools.VCF_ref_check import VCFRefChecker
    else:
        from .VCF_ref_check import VCFRefChecker

    checker = VCFRefChecker(error_handler)
    checker.check_vcf_against_fasta(
        vcf_path,
        reference_path,
    )


def _effective_report_path(args: argparse.Namespace) -> str:
    """Return the explicit or command-specific default report path."""
    if args.error_report:
        return args.error_report
    if args.command == "generic-to-dbsnp":
        return default_error_report_path(args.output)
    return DbSNPVCFCleaner._default_report_path(args.output)


def _default_output_path(command: str, vcf_path: str) -> str:
    """Derive a command-specific dbSNP VCF output path from the input path."""
    for suffix in (".vcf.gz", ".vcf"):
        if vcf_path.endswith(suffix):
            stem = vcf_path[: -len(suffix)]
            output_suffix = suffix
            break
    else:
        stem = vcf_path
        output_suffix = ".vcf"

    qualifier = ".dbsnp" if command == "generic-to-dbsnp" else ".cleaned"
    return f"{stem}{qualifier}{output_suffix}"


def _write_failure_report(args: argparse.Namespace, error_handler: ErrorHandler) -> None:
    """Write collected issues when validation exits before normal reporting."""
    error_handler.generate_report(
        _effective_report_path(args),
        vcf_file_path=args.vcf,
        output_tsv_path=args.output,
    )


def convert_generic_to_dbsnp(args: argparse.Namespace) -> None:
    """Convert a generic VCF into a dbSNP-formatted VCF."""
    error_handler = ErrorHandler()
    try:
        _run_reference_check(args.vcf, args.reference, error_handler)
        converter = VCF2dbSNPConverter(error_handler)
        converter.convert_vcf_to_dbsnp(
            vcf_file_path=args.vcf,
            output_file_path=args.output,
            metadata_file_path=args.metadata,
            error_report_path=args.error_report,
        )
    except Exception:
        _write_failure_report(args, error_handler)
        raise


def validate_dbsnp(args: argparse.Namespace) -> None:
    """Validate and rewrite a cleaned dbSNP VCF."""
    error_handler = ErrorHandler()
    try:
        _run_reference_check(args.vcf, args.reference, error_handler)
        cleaner = DbSNPVCFCleaner(error_handler)
        cleaner.clean(
            vcf_file_path=args.vcf,
            output_file_path=args.output,
            metadata_file_path=args.metadata,
            error_report_path=args.error_report,
        )
    except Exception:
        _write_failure_report(args, error_handler)
        raise


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
    parser.add_argument(
        "-o",
        "--output",
        help="Optional output dbSNP VCF path (.vcf or .vcf.gz); derived from the input when omitted",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        required=True,
        help="Metadata file path containing Experiment_id and exactly one SampleSet_id",
    )
    parser.add_argument(
        "-e",
        "--error-report",
        help="Optional integrated validation report path, including reference results when used",
    )
    parser.add_argument("-r", "--reference", help="Optional reference FASTA for REF allele validation")


def _paths_refer_to_same_file(first_path: str, second_path: str) -> bool:
    """Compare paths by resolved name and, when possible, filesystem identity."""
    if os.path.realpath(first_path) == os.path.realpath(second_path):
        return True
    try:
        return os.path.samefile(first_path, second_path)
    except OSError:
        return False


def _validate_generated_paths(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Prevent generated files from overwriting each other or an input file."""
    report_path = _effective_report_path(args)
    if _paths_refer_to_same_file(args.output, report_path):
        parser.error("--output and --error-report must use different paths")

    input_paths = (
        ("--vcf", args.vcf),
        ("--metadata", args.metadata),
        ("--reference", args.reference),
    )
    generated_paths = (
        ("--output", args.output),
        ("--error-report", report_path),
    )
    for generated_option, generated_path in generated_paths:
        for input_option, input_path in input_paths:
            if input_path and _paths_refer_to_same_file(generated_path, input_path):
                parser.error(
                    f"{generated_option} must use a different path from {input_option}"
                )


def _validate_output_suffix(parser: argparse.ArgumentParser, output_path: str) -> None:
    """Require an explicit VCF filename for generated SNP output."""
    if not output_path.endswith((".vcf", ".vcf.gz")):
        parser.error("--output must end with .vcf or .vcf.gz")


def main() -> None:
    """Run the public submission validation CLI."""
    parser = build_parser()
    args = parser.parse_args()
    if not args.output:
        args.output = _default_output_path(args.command, args.vcf)
    _validate_generated_paths(parser, args)
    _validate_output_suffix(parser, args.output)
    args.func(args)


if __name__ == "__main__":
    main()
