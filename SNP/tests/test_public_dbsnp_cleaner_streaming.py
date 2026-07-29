#!/usr/bin/env python3
"""Focused tests for streaming dbSNP VCF cleaning and atomic output publishing."""

import gzip
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
METADATA = ROOT / "examples" / "toy.human.snp.metadata.txt"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from kvar_snp_tools.dbsnp_vcf_cleaner import DbSNPVCFCleaner  # noqa: E402
from kvar_snp_tools.VCF2dbSNP import VCF2dbSNPConverter  # noqa: E402


def _write_vcf(path: Path, row: str) -> None:
    """Write a minimal dbSNP VCF containing one caller-provided data row."""
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##fileDate=20260601",
                "##handle=KVar",
                "##batch=EXP001",
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
    cleaner.clean(
        str(input_vcf),
        str(output_vcf),
        metadata_file_path=str(METADATA),
        error_report_path=str(error_report),
    )

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
        DbSNPVCFCleaner().clean(
            str(input_vcf),
            str(output_vcf),
            metadata_file_path=str(METADATA),
            error_report_path=str(error_report),
        )
    except RuntimeError as exc:
        assert "blocked by" in str(exc)
    else:
        raise AssertionError("Expected blocking validation to raise RuntimeError")

    assert not output_vcf.exists()
    assert not list(tmp_path.glob(f".{output_vcf.name}.*.tmp"))
    assert "INVALID_ALLELE_FREQUENCY" in error_report.read_text(encoding="utf-8")


def test_cleaner_never_relabels_input_population_from_metadata(tmp_path: Path) -> None:
    """A metadata mismatch blocks output without changing the parsed VCF ID."""
    input_vcf = tmp_path / "population-mismatch.dbsnp.vcf"
    output_vcf = tmp_path / "population-mismatch.cleaned.vcf"
    error_report = tmp_path / "population-mismatch.errors.txt"
    metadata = tmp_path / "population-mismatch.metadata.txt"
    _write_vcf(
        input_vcf,
        "chr1\t2\tvar1\tC\tT\t.\tPASS\tVRT=1;AC=2;AN=10;AF=0.2\tNA\t10",
    )
    metadata.write_text(
        "##Experiment_id=EXP001\n"
        "##reference=toy_ref\n"
        "##SampleSet_id=OTHER_POPULATION\n",
        encoding="utf-8",
    )
    cleaner = DbSNPVCFCleaner()

    try:
        cleaner.clean(
            str(input_vcf),
            str(output_vcf),
            metadata_file_path=str(metadata),
            error_report_path=str(error_report),
        )
    except RuntimeError as exc:
        assert "blocked by" in str(exc)
    else:
        raise AssertionError("Expected population mismatch to block output")

    assert cleaner.parser.header.population_ids == ["POP1"]
    assert not output_vcf.exists()
    assert "METADATA_POPULATION_MISMATCH" in error_report.read_text(encoding="utf-8")


def test_generic_converter_publishes_output_atomically(tmp_path: Path) -> None:
    """A write failure preserves an existing output and removes staged data."""
    output_vcf = tmp_path / "generic.atomic.vcf"
    error_report = tmp_path / "generic.atomic.errors.txt"
    output_vcf.write_text("sentinel\n", encoding="utf-8")

    class FailingConverter(VCF2dbSNPConverter):
        def _write_metadata(self, handle) -> None:
            handle.write("partial\n")
            raise OSError("forced output failure")

    try:
        FailingConverter().convert_vcf_to_dbsnp(
            str(ROOT / "examples" / "toy.human.snp.generic.vcf"),
            str(output_vcf),
            metadata_file_path=str(METADATA),
            error_report_path=str(error_report),
        )
    except OSError as exc:
        assert "forced output failure" in str(exc)
    else:
        raise AssertionError("Expected the forced output write to fail")

    assert output_vcf.read_text(encoding="utf-8") == "sentinel\n"
    assert not list(tmp_path.glob(f".{output_vcf.name}.*.tmp"))


def test_generic_converter_preserves_gzip_output(tmp_path: Path) -> None:
    """Atomic staging still compresses output selected by the final suffix."""
    output_vcf = tmp_path / "generic.atomic.vcf.gz"
    error_report = tmp_path / "generic.atomic-gzip.errors.txt"

    VCF2dbSNPConverter().convert_vcf_to_dbsnp(
        str(ROOT / "examples" / "toy.human.snp.generic.vcf"),
        str(output_vcf),
        metadata_file_path=str(METADATA),
        error_report_path=str(error_report),
    )

    with gzip.open(output_vcf, "rt", encoding="utf-8") as handle:
        output_text = handle.read()
    assert output_text.startswith("##fileformat=VCFv4.1\n")
    assert "##population_id=POP1" in output_text
    assert not list(tmp_path.glob(f".{output_vcf.name}.*.tmp"))


def main() -> None:
    """Run tests without requiring pytest."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        test_cleaner_streams_rows_and_preserves_output_semantics(tmp_dir)
        test_cleaner_does_not_publish_output_when_validation_blocks(tmp_dir)
        test_cleaner_never_relabels_input_population_from_metadata(tmp_dir)
        test_generic_converter_publishes_output_atomically(tmp_dir)
        test_generic_converter_preserves_gzip_output(tmp_dir)
    print("Public dbSNP cleaner streaming tests passed.")


if __name__ == "__main__":
    main()
