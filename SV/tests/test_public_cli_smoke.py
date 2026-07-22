#!/usr/bin/env python3
"""Smoke test for the public KVar SV VCF-to-TSV CLI."""

import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "src" / "kvar_sv_tools" / "vcf_to_kvar_tsv.py"


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def main() -> None:
    help_result = subprocess.run(
        [sys.executable, str(CLI), "--help"],
        cwd=str(ROOT),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert help_result.returncode == 0
    for option in ("--call-tsv", "--error-report", "--call-accession-start"):
        assert option in help_result.stdout
    for option in ("--region-tsv", "--region-accession-start", "--id-mapping"):
        assert option not in help_result.stdout

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        reference = tmp / "toy.fa"
        reference_index = tmp / "toy.fa.fai"
        metadata = tmp / "metadata.txt"
        vcf = tmp / "input.vcf"
        call_tsv = tmp / "Variant_Call.tsv"
        errors = tmp / "validation_report.txt"
        invalid_vcf = tmp / "invalid_ref.vcf"
        invalid_call_tsv = tmp / "invalid.Variant_Call.tsv"
        invalid_errors = tmp / "invalid.validation_report.txt"
        report_failure_call_tsv = tmp / "report_failure.Variant_Call.tsv"
        missing_report = tmp / "missing" / "validation_report.txt"

        sequence = "A" * 200
        reference_header = ">chr1 AC:CM000663.2 AS:GRCh38\n"
        write_text(reference, f"{reference_header}{sequence}\n")
        write_text(reference_index, f"chr1\t200\t{len(reference_header)}\t200\t201\n")
        write_text(
            metadata,
            "##SampleSet_id=toy_samples\n"
            "##Experiment_id=toy_experiment\n"
            "##reference=toy_reference\n",
        )
        write_text(
            vcf,
            "##fileformat=VCFv4.2\n"
            "##reference=toy_reference\n"
            "##SampleSet_id=toy_samples\n"
            "##Experiment_id=toy_experiment\n"
            "##contig=<ID=chr1,length=200>\n"
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate">\n'
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n'
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n'
            '##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t20\tsv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30;SVLEN=-10;AC=1;AN=2\n",
        )
        write_text(
            invalid_vcf,
            vcf.read_text(encoding="utf-8").replace("\tsv1\tA\t", "\tsv_bad\tC\t"),
        )

        result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "--reference-fasta",
                str(reference),
                "--metadata",
                str(metadata),
                "--call-tsv",
                str(call_tsv),
                "--error-report",
                str(errors),
                "--call-accession-start",
                "7",
                "--sanitize-error-report",
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if result.returncode != 0:
            raise AssertionError(
                "CLI failed\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}\n"
            )

        call_text = call_tsv.read_text(encoding="utf-8")
        report_text = errors.read_text(encoding="utf-8")

        call_header = next(
            line.lstrip("#").split("\t")
            for line in call_text.splitlines()
            if line.startswith("#") and not line.startswith("##")
        )

        assert "Variant_Call_ID" in call_header
        assert "Variant_Call_Type" in call_header
        assert "Outer_Start" in call_header
        assert "HGVSG" in call_header
        assert "Variant Call ID" not in call_header
        assert "kssv7" in call_text
        assert "NC_000001.11:g.20_30del" in call_text
        assert "No errors." in report_text

        invalid_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(invalid_vcf),
                "--reference-fasta",
                str(reference),
                "--metadata",
                str(metadata),
                "--call-tsv",
                str(invalid_call_tsv),
                "--error-report",
                str(invalid_errors),
                "--call-accession-start",
                "1",
                "--sanitize-error-report",
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        assert invalid_result.returncode != 0
        assert not invalid_call_tsv.exists()
        assert invalid_errors.exists()
        assert "REF_MISMATCH" in invalid_errors.read_text(encoding="utf-8")

        report_failure_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "--reference-fasta",
                str(reference),
                "--metadata",
                str(metadata),
                "--call-tsv",
                str(report_failure_call_tsv),
                "--error-report",
                str(missing_report),
                "--sanitize-error-report",
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        assert report_failure_result.returncode != 0
        assert not report_failure_call_tsv.exists()
        assert not missing_report.exists()

    print("public CLI smoke test passed")


if __name__ == "__main__":
    main()
