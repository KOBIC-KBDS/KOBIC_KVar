#!/usr/bin/env python3
"""Smoke tests for the public SNP VCF validation tools."""

import subprocess
import shutil
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "src" / "kvar_snp_tools" / "Sub_validator.py"
EXAMPLES = ROOT / "examples"


def run_command(args):
    """Run a command and fail with useful output."""
    result = subprocess.run(args, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(
            f"Command failed: {' '.join(str(arg) for arg in args)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result


def test_generic_to_dbsnp_writes_population_id_only(tmp_dir):
    """Generic VCF conversion writes population_id and no SampleSet_id header."""
    output_vcf = tmp_dir / "generic.cleaned.dbsnp.vcf"
    error_report = tmp_dir / "generic.errors.txt"
    run_command(
        [
            sys.executable,
            str(CLI),
            "generic-to-dbsnp",
            "-v",
            str(EXAMPLES / "toy.generic.vcf"),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    output_text = output_vcf.read_text(encoding="utf-8")
    assert "##population_id=POP1" in output_text
    assert "##contig=<ID=chr1,length=20>" in output_text
    assert "##SampleSet_id=" not in output_text
    assert "VRT=1" in output_text
    assert "AF=0.2" in output_text
    assert error_report.exists()


def test_validate_dbsnp_writes_cleaned_vcf(tmp_dir):
    """dbSNP VCF validation rewrites a cleaned dbSNP VCF."""
    input_vcf = tmp_dir / "toy.annotated.dbsnp.vcf"
    output_vcf = tmp_dir / "dbsnp.cleaned.vcf"
    error_report = tmp_dir / "dbsnp.errors.txt"
    input_text = (EXAMPLES / "toy.dbsnp.vcf").read_text(encoding="utf-8")
    input_text = input_text.replace(
        "##FORMAT=<ID=NA",
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. Format: Allele|Consequence\">\n"
        "##INFO=<ID=61KJPN_AC,Number=A,Type=Integer,Description=\"61KJPN allele count\">\n"
        "##FORMAT=<ID=NA",
    )
    input_text = input_text.replace(
        "VRT=1;AC=2;AN=10;AF=0.2",
        "VRT=1;AC=2;AN=10;AF=0.2;CSQ=T|missense_variant;61KJPN_AC=7",
    )
    input_vcf.write_text(input_text, encoding="utf-8")
    run_command(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(input_vcf),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    output_text = output_vcf.read_text(encoding="utf-8")
    assert "##population_id=POP1" in output_text
    assert "##contig=<ID=chr1,length=20>" in output_text
    assert "##SampleSet_id=" not in output_text
    assert output_text.count("##INFO=<ID=VRT") == 1
    assert "\tNA:FRQ\t10:0.2" in output_text
    assert "CSQ" not in output_text
    assert "61KJPN_AC" not in output_text
    assert "JPN61K_AC" not in output_text
    assert "UNSUPPORTED_DBSNP_INFO_TAG" in error_report.read_text(encoding="utf-8")
    assert error_report.exists()


def test_validate_dbsnp_with_reference_check(tmp_dir):
    """Reference checking can run before dbSNP VCF cleaning."""
    input_vcf = tmp_dir / "toy.dbsnp.vcf"
    reference_fasta = tmp_dir / "toy.reference.fa"
    output_vcf = tmp_dir / "dbsnp.refchecked.cleaned.vcf"
    error_report = tmp_dir / "dbsnp.refchecked.errors.txt"
    reference_report = tmp_dir / "dbsnp.reference.txt"
    shutil.copyfile(EXAMPLES / "toy.dbsnp.vcf", input_vcf)
    shutil.copyfile(EXAMPLES / "toy.reference.fa", reference_fasta)

    run_command(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(input_vcf),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-r",
            str(reference_fasta),
            "-rr",
            str(reference_report),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    assert output_vcf.exists()
    assert error_report.exists()
    assert reference_report.exists()


def main():
    """Run smoke tests without requiring pytest."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        test_generic_to_dbsnp_writes_population_id_only(tmp_dir)
        test_validate_dbsnp_writes_cleaned_vcf(tmp_dir)
        test_validate_dbsnp_with_reference_check(tmp_dir)
    print("Public CLI smoke tests passed.")


if __name__ == "__main__":
    main()
