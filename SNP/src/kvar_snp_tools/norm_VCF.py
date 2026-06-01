#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF normalization using bcftools.
Runs: norm -f REF -m -both (normalize + biallelic), then norm -d exact (dedup),
then assigns deterministic local IDs, indexes, and runs bcftools stats.
"""

import sys
import os
import argparse
import gzip
import subprocess
import tempfile
from typing import Optional


def _vcf_basename(path: str) -> str:
    """Return path with .vcf.gz or .vcf stripped (for .stats.txt and .norm.log)."""
    if path.endswith(".vcf.gz"):
        return path[:-7]  # strip .vcf.gz
    if path.endswith(".vcf"):
        return path[:-4]  # strip .vcf
    return path


def _stats_path(vcf_path: str) -> str:
    """Path for bcftools stats output: same dir as VCF, basename.stats.txt."""
    base = _vcf_basename(vcf_path)
    return base + ".stats.txt"


def _log_path(output_vcf_path: str) -> str:
    """Log path: same dir as output, basename.norm.log (avoid .norm.norm.log when output is already *.norm.vcf.gz)."""
    base = _vcf_basename(output_vcf_path)
    if base.endswith(".norm"):
        return base + ".log"
    return base + ".norm.log"


def _run_cmd(cmd: list, log_file, section_title: str) -> int:
    """Run command, write section header and stdout/stderr to log_file. Return returncode."""
    print(section_title, file=log_file)
    print("-" * 60, file=log_file)
    log_file.flush()
    try:
        result = subprocess.run(
            cmd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
        )
        print("", file=log_file)
        log_file.flush()
        return result.returncode
    except Exception as e:
        print(f"Error running command: {e}", file=log_file)
        log_file.flush()
        return -1


def _vcf_has_contig_header(vcf_path: str) -> bool:
    """Return True when the VCF header contains at least one contig line."""
    opener = gzip.open if vcf_path.endswith(".gz") else open
    mode = "rt"
    with opener(vcf_path, mode, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##contig="):
                return True
            if line.startswith("#CHROM"):
                return False
    return False


def _write_vcf_with_sequential_ids(
    bcftools: str,
    input_vcf: str,
    output_vcf: str,
    log_file,
) -> int:
    """Rewrite VCF IDs as 1-based sequential local IDs while preserving headers."""
    out_dir = os.path.dirname(output_vcf) or "."
    tmp_plain = None
    print("=== Step 3: assign sequential local IDs ===", file=log_file)
    print("-" * 60, file=log_file)
    log_file.flush()

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            suffix=".ids.vcf",
            dir=out_dir,
            delete=False,
        ) as tmp_file:
            tmp_plain = tmp_file.name

            process = subprocess.Popen(
                [bcftools, "view", input_vcf],
                stdout=subprocess.PIPE,
                stderr=log_file,
                text=True,
            )

            if process.stdout is None:
                print("Error: Could not read bcftools view output.", file=log_file)
                return -1

            record_index = 0
            for line in process.stdout:
                if line.startswith("#"):
                    tmp_file.write(line)
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    tmp_file.write(line)
                    continue

                record_index += 1
                fields[2] = str(record_index)
                tmp_file.write("\t".join(fields) + "\n")

            view_rc = process.wait()
            if view_rc != 0:
                print(f"bcftools view failed with exit code {view_rc}.", file=log_file)
                return view_rc

        print("", file=log_file)
        log_file.flush()
        return _run_cmd(
            [bcftools, "view", "-Oz", "-o", output_vcf, tmp_plain],
            log_file,
            "=== Step 4: compress ID-normalized VCF ===",
        )
    except Exception as e:
        print(f"Error assigning local IDs: {e}", file=log_file)
        log_file.flush()
        return -1
    finally:
        if tmp_plain and os.path.isfile(tmp_plain):
            try:
                os.remove(tmp_plain)
            except OSError as e:
                print(f"Warning: Could not remove temporary file {tmp_plain}: {e}", file=log_file)
                log_file.flush()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalize VCF with bcftools, split to biallelic records, assign local IDs, index, and run stats.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-b",
        dest="bcftools",
        default="bcftools",
        help="Path to bcftools executable (default: bcftools from PATH)",
    )
    parser.add_argument(
        "-r",
        dest="reference",
        required=True,
        help="Reference FASTA path (required)",
    )
    parser.add_argument(
        "-v",
        dest="vcf",
        required=True,
        help="Input VCF path (.vcf or .vcf.gz)",
    )
    parser.add_argument(
        "-o",
        dest="output",
        default=None,
        help="Output VCF path (default: same dir as input, <input_basename>.norm.vcf.gz)",
    )
    args = parser.parse_args()

    ref = os.path.abspath(args.reference)
    vcf_in = os.path.abspath(args.vcf)

    # Require reference to exist
    if not os.path.isfile(ref):
        print(f"Error: Reference file not found: {ref}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(vcf_in):
        print(f"Error: Input VCF not found: {vcf_in}", file=sys.stderr)
        sys.exit(1)

    if not (vcf_in.endswith(".vcf.gz") or vcf_in.endswith(".vcf")):
        print("Error: Input VCF must have extension .vcf or .vcf.gz", file=sys.stderr)
        sys.exit(1)

    # Output path: always .vcf.gz
    if args.output:
        out_vcf = os.path.abspath(args.output)
        if not out_vcf.endswith(".vcf.gz"):
            out_vcf = out_vcf.rstrip("/")
            if out_vcf.endswith(".vcf"):
                out_vcf = out_vcf[:-4] + ".vcf.gz"
            else:
                out_vcf = out_vcf + ".vcf.gz"
    else:
        base = _vcf_basename(vcf_in)
        out_vcf = base + ".norm.vcf.gz"

    out_dir = os.path.dirname(out_vcf)
    if out_dir and not os.path.isdir(out_dir):
        print(f"Error: Output directory does not exist: {out_dir}", file=sys.stderr)
        sys.exit(1)

    log_path = _log_path(out_vcf)
    out_base = _vcf_basename(out_vcf)
    reheader_source_intermediate = os.path.join(out_dir, os.path.basename(out_base) + ".tmp.reheader_source.vcf.gz")
    reheader_intermediate = os.path.join(out_dir, os.path.basename(out_base) + ".tmp.reheader.vcf.gz")
    normalized_intermediate = os.path.join(out_dir, os.path.basename(out_base) + ".tmp.norm.vcf.gz")
    dedup_intermediate = os.path.join(out_dir, os.path.basename(out_base) + ".tmp.dedup.vcf.gz")

    bcftools = args.bcftools

    with open(log_path, "w", encoding="utf-8") as log_file:
        norm_input_vcf = vcf_in
        if not _vcf_has_contig_header(vcf_in):
            fai_path = ref + ".fai"
            if not os.path.isfile(fai_path):
                print(f"Error: Reference FASTA index not found: {fai_path}", file=log_file)
                print(f"Contig-less VCF input requires a FASTA index for reheadering.", file=log_file)
                print(f"Step 0 failed (see {log_path}).", file=sys.stderr)
                sys.exit(1)

            cmd_bgzip = [bcftools, "view", "-Oz", "-o", reheader_source_intermediate, vcf_in]
            rc_bgzip = _run_cmd(
                cmd_bgzip,
                log_file,
                "=== Step 0a: bcftools view -Oz (prepare for reheader) ===",
            )
            if rc_bgzip != 0:
                print(f"Step 0a failed (see {log_path}).", file=sys.stderr)
                sys.exit(1)

            cmd_reheader = [bcftools, "reheader", "-f", fai_path, "-o", reheader_intermediate, reheader_source_intermediate]
            rc_reheader = _run_cmd(
                cmd_reheader,
                log_file,
                "=== Step 0b: bcftools reheader -f REF.fai (add contig headers) ===",
            )
            if rc_reheader != 0:
                print(f"Step 0b failed (see {log_path}).", file=sys.stderr)
                sys.exit(1)
            norm_input_vcf = reheader_intermediate

        # Step 1: norm -f REF -m -both
        cmd1 = [bcftools, "norm", "-f", ref, "-m", "-both", norm_input_vcf, "-Oz", "-o", normalized_intermediate, "--force"]
        rc1 = _run_cmd(cmd1, log_file, "=== Step 1: bcftools norm -f REF -m -both (normalize + biallelic) ===")
        if rc1 != 0:
            print(f"Step 1 failed (see {log_path}).", file=sys.stderr)
            sys.exit(1)

        # Step 2: norm -d exact
        cmd2 = [bcftools, "norm", "-d", "exact", normalized_intermediate, "-Oz", "-o", dedup_intermediate, "--force"]
        rc2 = _run_cmd(cmd2, log_file, "=== Step 2: bcftools norm -d exact (deduplicate) ===")
        if rc2 != 0:
            if os.path.isfile(normalized_intermediate):
                try:
                    os.remove(normalized_intermediate)
                except OSError:
                    pass
            print(f"Step 2 failed (see {log_path}).", file=sys.stderr)
            sys.exit(1)

        # Step 3: assign deterministic local IDs without keeping original IDs in the output.
        rc3 = _write_vcf_with_sequential_ids(bcftools, dedup_intermediate, out_vcf, log_file)
        if rc3 != 0:
            print(f"ID reassignment step failed (see {log_path}).", file=sys.stderr)
            sys.exit(1)

        # Remove intermediates
        for intermediate in (reheader_source_intermediate, reheader_intermediate, normalized_intermediate, dedup_intermediate):
            if not os.path.isfile(intermediate):
                continue
            try:
                os.remove(intermediate)
            except OSError as e:
                print(f"Warning: Could not remove intermediate file {intermediate}: {e}", file=log_file)
                log_file.flush()

        # Index final output
        cmd_index = [bcftools, "index", "-t", out_vcf]
        rc_index = _run_cmd(cmd_index, log_file, "=== bcftools index -t (final VCF) ===")
        if rc_index != 0:
            print(f"Index step failed (see {log_path}).", file=sys.stderr)
            sys.exit(1)

        # bcftools stats on input VCF
        input_stats_path = _stats_path(vcf_in)
        cmd_stats_in = [bcftools, "stats", vcf_in]
        print(f"=== bcftools stats (input VCF) -> {input_stats_path} ===", file=log_file)
        print("-" * 60, file=log_file)
        log_file.flush()
        try:
            with open(input_stats_path, "w", encoding="utf-8") as f_stats:
                subprocess.run(cmd_stats_in, stdout=f_stats, stderr=log_file, text=True, check=False)
        except Exception as e:
            print(f"Error writing input stats: {e}", file=log_file)
            log_file.flush()
        print("", file=log_file)
        log_file.flush()

        # bcftools stats on output VCF
        output_stats_path = _stats_path(out_vcf)
        cmd_stats_out = [bcftools, "stats", out_vcf]
        print(f"=== bcftools stats (output VCF) -> {output_stats_path} ===", file=log_file)
        print("-" * 60, file=log_file)
        log_file.flush()
        try:
            with open(output_stats_path, "w", encoding="utf-8") as f_stats:
                subprocess.run(cmd_stats_out, stdout=f_stats, stderr=log_file, text=True, check=False)
        except Exception as e:
            print(f"Error writing output stats: {e}", file=log_file)
            log_file.flush()

    print(f"Done. Output: {out_vcf}")
    print(f"Log: {log_path}")
    print(f"Input stats: {input_stats_path}")
    print(f"Output stats: {output_stats_path}")


if __name__ == "__main__":
    main()
