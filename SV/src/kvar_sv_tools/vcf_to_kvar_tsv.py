#!/usr/bin/env python3
"""Public CLI for converting an SV VCF submission to a Call TSV."""

import argparse
import os
import sys

try:
    from .KVar2TSV import KVarTSVConverter
    from .error_handler import ErrorCode, ErrorHandler
except ImportError:
    from KVar2TSV import KVarTSVConverter
    from error_handler import ErrorCode, ErrorHandler


def default_call_tsv_path(vcf_path: str) -> str:
    """Return the default Call TSV path beside the input VCF."""
    for suffix in (".vcf.gz", ".vcf"):
        if vcf_path.endswith(suffix):
            return vcf_path[: -len(suffix)] + ".call.tsv"
    return vcf_path + ".call.tsv"


def default_error_report_path(call_tsv_path: str) -> str:
    """Return the default validation report path for a Call TSV output."""
    for suffix in (".tsv.gz", ".tsv"):
        if call_tsv_path.endswith(suffix):
            return call_tsv_path[: -len(suffix)] + ".errors.txt"
    return call_tsv_path + ".errors.txt"


def _paths_refer_to_same_file(first_path: str, second_path: str) -> bool:
    """Compare full paths and, when possible, existing filesystem identity."""
    if os.path.realpath(first_path) == os.path.realpath(second_path):
        return True
    try:
        return os.path.samefile(first_path, second_path)
    except OSError:
        return False


def _validate_paths(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Prevent generated files from overwriting each other or an input file."""
    if _paths_refer_to_same_file(args.call_tsv, args.error_report):
        parser.error("--call-tsv and --error-report must use different file paths")

    input_paths = [("--vcf", args.vcf)]
    if args.reference:
        input_paths.extend(
            [
                ("--reference", args.reference),
                ("reference FASTA index", f"{args.reference}.fai"),
            ]
        )

    generated_paths = (
        ("--call-tsv", args.call_tsv),
        ("--error-report", args.error_report),
    )
    for generated_option, generated_path in generated_paths:
        for input_option, input_path in input_paths:
            if _paths_refer_to_same_file(generated_path, input_path):
                parser.error(
                    f"{generated_option} must use a different file path from {input_option}"
                )


def _validate_call_tsv_suffix(
    parser: argparse.ArgumentParser,
    call_tsv_path: str,
) -> None:
    """Require a plain TSV filename for generated Call output."""
    if not call_tsv_path.endswith(".tsv"):
        parser.error("--call-tsv must end with .tsv")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate an SV VCF file and convert it to a KVar submission Call TSV file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input SV VCF path (.gz supported)")
    parser.add_argument(
        "-r",
        "--reference",
        help="Optional reference FASTA for chromosome, coordinate, and REF validation; its .fai index is created automatically",
    )
    parser.add_argument(
        "-t",
        "--call-tsv",
        help="Optional output Call TSV path; derived from the input when omitted",
    )
    parser.add_argument(
        "-e",
        "--error-report",
        help="Optional validation report path; defaults to <call-tsv-base>.errors.txt",
    )
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.call_tsv:
        args.call_tsv = default_call_tsv_path(args.vcf)
    if not args.error_report:
        args.error_report = default_error_report_path(args.call_tsv)
    _validate_paths(parser, args)
    _validate_call_tsv_suffix(parser, args.call_tsv)

    if not os.path.exists(args.vcf):
        parser.error(f"VCF file not found: {args.vcf}")

    error_handler = ErrorHandler()
    try:
        converter = KVarTSVConverter(
            error_handler=error_handler,
            reference_fasta_path=args.reference,
        )
    except Exception as exc:
        if not error_handler.has_errors():
            error_handler.create_error(
                ErrorCode.FASTA_READ_ERROR,
                additional_info={
                    "file_path": args.reference,
                    "error": str(exc),
                },
            )
        try:
            error_handler.generate_report(
                args.error_report,
                vcf_file_path=args.vcf,
                output_tsv_path=args.call_tsv,
            )
        except Exception as report_exc:
            print(
                f"Reference initialization failed and the validation report "
                f"could not be written: {report_exc}",
                file=sys.stderr,
            )
        print(f"Reference initialization failed: {exc}", file=sys.stderr)
        return 1

    try:
        converter.convert_vcf_to_tsv(
            args.vcf,
            args.call_tsv,
            error_report_path=args.error_report,
        )
    except Exception as exc:
        print(f"Conversion failed: {exc}", file=sys.stderr)
        return 1

    print(f"Call TSV: {args.call_tsv}")
    print(f"Validation report: {args.error_report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
