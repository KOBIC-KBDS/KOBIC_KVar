#!/usr/bin/env python3
"""Focused tests for streaming dbSNP VCF cleaning and atomic output publishing."""

import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from kvar_snp_tools.dbsnp_vcf_cleaner import DbSNPVCFCleaner  # noqa: E402


def _write_vcf(path: Path, row: str) -> None:
    """Write a minimal dbSNP VCF containing one caller-provided data row."""
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##fileDate=20260601",
                "##handle=KVar",
                "##batch=EXP001",
                "##bioproject_id=PRJNA000001",
                "##biosample_id=SAMD000001",
                "##reference=toy_ref",
                "##contig=<ID=chr1,length=20>",
                '##INFO=<ID=VRT,Number=1,Type=Integer,Description="Variation type, 1 - SNV: single nucleotide variation, 2 - DIV: deletion/insertion variation, 3 - HETEROZYGOUS: variable, but undefined at nucleotide level, 4 - STR: short tandem repeat (microsatellite) variation, 5 - NAMED: insertion/deletion variation of named repetitive element, 6 - NO VARIATION: sequence scanned for variation, but none observed, 7 - MIXED: cluster contains submissions from 2 or more allelic classes, 8 - MNV: multiple nucleotide variation with alleles of common length greater than 1">',
                '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
                '##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">',
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Unsupported annotation">',
                '##FORMAT=<ID=NA,Number=1,Type=Integer,Description="Number of alleles">',
                "##population_id=POP1",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPOP1",
                row,
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_cleaner_streams_rows_and_preserves_output_semantics(tmp_path: Path) -> None:
    """Validation does not retain rows while output still applies corrections."""
    input_vcf = tmp_path / "streaming.dbsnp.vcf"
    output_vcf = tmp_path / "streaming.cleaned.vcf"
    error_report = tmp_path / "streaming.errors.txt"
    _write_vcf(
        input_vcf,
        "chr1\t2\tvar1\tC\tCT\t.\tPASS\tVRT=1;AC=2;AN=10;AF=0.2;CSQ=C|insertion\tNA\t10",
    )

    cleaner = DbSNPVCFCleaner()
    cleaner.clean(str(input_vcf), str(output_vcf), error_report_path=str(error_report))

    output_text = output_vcf.read_text(encoding="utf-8")
    report_text = error_report.read_text(encoding="utf-8")
    assert cleaner.parser.data_rows == []
    assert "##contig=<ID=chr1,length=20>" in output_text
    assert "VRT=2;AC=2;AN=10;AF=0.2" in output_text
    assert "CSQ" not in output_text
    assert "VRT_REF_ALT_MISMATCH" in report_text
    assert "UNSUPPORTED_DBSNP_INFO_TAG" in report_text


def test_cleaner_does_not_publish_output_when_validation_blocks(tmp_path: Path) -> None:
    """Blocking validation errors leave no final or temporary cleaned VCF."""
    input_vcf = tmp_path / "blocked.dbsnp.vcf"
    output_vcf = tmp_path / "blocked.cleaned.vcf"
    error_report = tmp_path / "blocked.errors.txt"
    _write_vcf(
        input_vcf,
        "chr1\t2\tvar1\tC\tT\t.\tPASS\tVRT=1;AC=2;AN=10;AF=2.0\tNA\t10",
    )

    try:
        DbSNPVCFCleaner().clean(str(input_vcf), str(output_vcf), error_report_path=str(error_report))
    except RuntimeError as exc:
        assert "blocked by" in str(exc)
    else:
        raise AssertionError("Expected blocking validation to raise RuntimeError")

    assert not output_vcf.exists()
    assert not list(tmp_path.glob(f".{output_vcf.name}.*.tmp"))
    assert "INVALID_ALLELE_FREQUENCY" in error_report.read_text(encoding="utf-8")


def main() -> None:
    """Run tests without requiring pytest."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        test_cleaner_streams_rows_and_preserves_output_semantics(tmp_dir)
        test_cleaner_does_not_publish_output_when_validation_blocks(tmp_dir)
    print("Public dbSNP cleaner streaming tests passed.")


if __name__ == "__main__":
    main()
