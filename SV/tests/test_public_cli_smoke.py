#!/usr/bin/env python3
"""Smoke test for the public KVar SV VCF-to-TSV CLI."""

import inspect
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "src" / "kvar_sv_tools" / "vcf_to_kvar_tsv.py"
EXAMPLE_VCF = ROOT / "examples" / "toy.human.sv.vcf"
EXAMPLE_REFERENCE = ROOT / "examples" / "toy.human.sv.reference.fasta"
sys.path.insert(0, str(ROOT / "src"))

EXPECTED_SUBMISSION_CALL_HEADER = [
    "Variant_Call_ID",
    "Variant_Call_Type",
    "Chr",
    "Outer_Start",
    "Start",
    "Inner_Start",
    "Inner_Stop",
    "Stop",
    "Outer_Stop",
    "Insertion_Length",
    "Allele_Count",
    "Allele_Frequency",
    "Allele_Number",
    "Copy_Number",
    "Description",
    "Validation",
    "Zygosity",
    "Origin",
    "Phenotype",
    "External_Links",
    "Evidence",
    "Sequence",
    "From_Chr",
    "From_Coord",
    "From_Strand",
    "To_Chr",
    "To_Coord",
    "To_Strand",
    "Mutation_ID",
    "Mutation_Order",
    "Mutation_Molecule",
    "BND_Source_VCF_IDs",
]


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def assert_submission_call_schema(text: str):
    lines = text.splitlines()
    header = next(
        line.lstrip("#").split("\t")
        for line in lines
        if line.startswith("#") and not line.startswith("##")
    )
    rows = [line.split("\t") for line in lines if line and not line.startswith("#")]

    assert header == EXPECTED_SUBMISSION_CALL_HEADER
    assert len(header) == 32
    assert all(len(row) == len(header) for row in rows)
    return header, rows


def main() -> None:
    from kvar_sv_tools.KVar2TSV import KVarTSVConverter

    help_result = subprocess.run(
        [sys.executable, str(CLI), "--help"],
        cwd=str(ROOT),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert help_result.returncode == 0
    for option in ("--reference", "--call-tsv", "--error-report"):
        assert option in help_result.stdout
    for option in (
        "--reference-fasta",
        "--call-accession-start",
        "--metadata",
        "--region-tsv",
        "--region-accession-start",
        "--id-mapping",
        "--sanitize-error-report",
    ):
        assert option not in help_result.stdout
    assert "call_accession_start" not in inspect.signature(
        KVarTSVConverter.convert_vcf_to_tsv
    ).parameters
    assert "metadata_file_path" not in inspect.signature(
        KVarTSVConverter.convert_vcf_to_tsv
    ).parameters

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        reference = tmp / "toy.fasta"
        reference_index = tmp / "toy.fasta.fai"
        vcf = tmp / "input.vcf"
        call_tsv = tmp / "submission.call.tsv"
        errors = tmp / "validation_report.txt"
        automatic_vcf = tmp / "automatic.vcf"
        automatic_call_tsv = tmp / "automatic.call.tsv"
        automatic_errors = tmp / "automatic.call.errors.txt"
        invalid_vcf = tmp / "invalid_ref.vcf"
        invalid_call_tsv = tmp / "invalid.call.tsv"
        invalid_errors = tmp / "invalid.validation_report.txt"
        unchecked_call_tsv = tmp / "unchecked.call.tsv"
        unchecked_errors = tmp / "unchecked.validation_report.txt"
        ignored_context_vcf = tmp / "ignored_context.vcf"
        ignored_context_call_tsv = tmp / "ignored_context.call.tsv"
        ignored_context_errors = tmp / "ignored_context.validation_report.txt"
        report_failure_call_tsv = tmp / "report_failure.call.tsv"
        missing_report = tmp / "missing" / "validation_report.txt"
        insertion_without_end_vcf = tmp / "insertion_without_end.vcf"
        insertion_without_end_tsv = tmp / "insertion_without_end.call.tsv"
        insertion_without_end_errors = tmp / "insertion_without_end.errors.txt"

        sequence = "A" * 200
        reference_header = ">chr1 AC:CM000663.2 AS:GRCh38\n"
        write_text(reference, f"{reference_header}{sequence}\n")
        write_text(reference_index, f"chr1\t200\t{len(reference_header)}\t200\t201\n")
        write_text(
            vcf,
            "##fileformat=VCFv4.2\n"
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
        write_text(
            insertion_without_end_vcf,
            vcf.read_text(encoding="utf-8").replace(
                "chr1\t20\tsv1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30;SVLEN=-10;AC=1;AN=2\n",
                "chr1\t20\tins_without_end\tA\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=12\n",
            ),
        )
        shutil.copyfile(vcf, automatic_vcf)
        write_text(
            ignored_context_vcf,
            vcf.read_text(encoding="utf-8")
            .replace(
                "##fileformat=VCFv4.2\n",
                "##fileformat=VCFv4.2\n"
                "##reference=unused_reference_a\n"
                "##reference=unused_reference_b\n"
                "##SampleSet_id=header_set_a\n"
                "##population_id=header_set_b\n"
                "##Experiment_id=header_experiment_a\n"
                "##batch=header_experiment_b\n"
                '##INFO=<ID=SAMPLESET,Number=1,Type=String,Description="Ignored context">\n'
                '##INFO=<ID=EXPERIMENT,Number=1,Type=String,Description="Ignored context">\n',
            )
            .replace(
                "SVTYPE=DEL;END=30;SVLEN=-10;AC=1;AN=2\n",
                "SVTYPE=DEL;END=30;SVLEN=-10;AC=1;AN=2;"
                "SAMPLESET=row_set;EXPERIMENT=row_experiment\n",
            ),
        )

        for legacy_option in ("-f", "--reference-fasta"):
            legacy_reference_result = subprocess.run(
                [
                    sys.executable,
                    str(CLI),
                    "--vcf",
                    str(vcf),
                    legacy_option,
                    str(reference),
                ],
                cwd=str(ROOT),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            assert legacy_reference_result.returncode != 0
            assert "unrecognized arguments" in legacy_reference_result.stderr

        result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "-r",
                str(reference),
                "--call-tsv",
                str(call_tsv),
                "--error-report",
                str(errors),
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

        call_header, call_rows = assert_submission_call_schema(call_text)

        assert "Variant_Call_ID" in call_header
        assert "Variant_Call_Type" in call_header
        assert "Outer_Start" in call_header
        assert "BND_Source_VCF_IDs" in call_header
        assert "Submitted_Variant_Call_IDs" not in call_header
        assert "Variant Call ID" not in call_header
        assert call_rows
        call_row = dict(zip(call_header, call_rows[0]))
        assert call_row["BND_Source_VCF_IDs"] == "."
        assert "sv1" in call_text
        assert "kssv" not in call_text
        assert "KVar SV Validation Report" in report_text
        assert "=== Validation Summary ===" in report_text
        assert "Total issues: 0" in report_text
        assert "No issues." in report_text

        KVarTSVConverter().convert_vcf_to_tsv(
            str(insertion_without_end_vcf),
            str(insertion_without_end_tsv),
            str(insertion_without_end_errors),
        )
        insertion_header, insertion_rows = assert_submission_call_schema(
            insertion_without_end_tsv.read_text(encoding="utf-8")
        )
        insertion_row = dict(zip(insertion_header, insertion_rows[0]))
        assert insertion_row["Variant_Call_Type"] == "insertion"
        assert insertion_row["Start"] == "20"
        assert insertion_row["Stop"] == "20"

        automatic_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(automatic_vcf),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        assert automatic_result.returncode == 0
        assert automatic_call_tsv.exists()
        assert automatic_errors.exists()
        automatic_report_text = automatic_errors.read_text(encoding="utf-8")
        assert "KVar SV Validation Report" in automatic_report_text
        assert "=== Validation Summary ===" in automatic_report_text
        assert "Total issues: 0" in automatic_report_text
        assert "No issues." in automatic_report_text
        assert str(automatic_errors) in automatic_result.stdout

        original_vcf = vcf.read_text(encoding="utf-8")
        input_collision_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "--call-tsv",
                str(vcf),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert input_collision_result.returncode != 0
        assert "--call-tsv must use a different file path from --vcf" in input_collision_result.stderr
        assert vcf.read_text(encoding="utf-8") == original_vcf

        shared_output = tmp / "shared-output.txt"
        output_collision_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "--call-tsv",
                str(shared_output),
                "--error-report",
                str(shared_output),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert output_collision_result.returncode != 0
        assert (
            "--call-tsv and --error-report must use different file paths"
            in output_collision_result.stderr
        )
        assert not shared_output.exists()

        invalid_suffix_output = tmp / "invalid-call-output.txt"
        invalid_suffix_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "--call-tsv",
                str(invalid_suffix_output),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert invalid_suffix_result.returncode != 0
        assert "--call-tsv must end with .tsv" in invalid_suffix_result.stderr
        assert not invalid_suffix_output.exists()

        original_reference = reference.read_text(encoding="utf-8")
        reference_collision_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "--reference",
                str(reference),
                "--call-tsv",
                str(reference),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert reference_collision_result.returncode != 0
        assert (
            "--call-tsv must use a different file path from --reference"
            in reference_collision_result.stderr
        )
        assert reference.read_text(encoding="utf-8") == original_reference

        invalid_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(invalid_vcf),
                "--reference",
                str(reference),
                "--call-tsv",
                str(invalid_call_tsv),
                "--error-report",
                str(invalid_errors),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        assert invalid_result.returncode != 0
        assert not invalid_call_tsv.exists()
        assert invalid_errors.exists()
        invalid_report_text = invalid_errors.read_text(encoding="utf-8")
        assert "KVar SV Validation Report" in invalid_report_text
        assert "=== Validation Summary ===" in invalid_report_text
        assert "REF_MISMATCH" in invalid_report_text
        assert "Total issues:" in invalid_report_text
        assert "Total errors:" not in invalid_report_text

        unchecked_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(invalid_vcf),
                "--call-tsv",
                str(unchecked_call_tsv),
                "--error-report",
                str(unchecked_errors),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        assert unchecked_result.returncode == 0
        assert unchecked_call_tsv.exists()
        unchecked_report_text = unchecked_errors.read_text(encoding="utf-8")
        assert "Total issues: 0" in unchecked_report_text
        assert "No issues." in unchecked_report_text
        assert "REF_MISMATCH" not in unchecked_report_text

        ignored_context_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(ignored_context_vcf),
                "--reference",
                str(reference),
                "--call-tsv",
                str(ignored_context_call_tsv),
                "--error-report",
                str(ignored_context_errors),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        assert ignored_context_result.returncode == 0
        assert ignored_context_call_tsv.exists()
        ignored_context_report = ignored_context_errors.read_text(encoding="utf-8")
        assert "Total issues: 0" in ignored_context_report
        assert "No issues." in ignored_context_report

        accession_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "--reference",
                str(reference),
                "--call-tsv",
                str(tmp / "accession.call.tsv"),
                "--error-report",
                str(tmp / "accession.validation_report.txt"),
                "--call-accession-start",
                "1",
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        assert accession_result.returncode != 0
        assert "unrecognized arguments: --call-accession-start 1" in accession_result.stderr

        metadata_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "--reference",
                str(reference),
                "--call-tsv",
                str(tmp / "metadata.call.tsv"),
                "--error-report",
                str(tmp / "metadata.validation_report.txt"),
                "--metadata",
                str(tmp / "metadata.txt"),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        assert metadata_result.returncode != 0
        assert "unrecognized arguments: --metadata" in metadata_result.stderr

        example_call_tsv = tmp / "example.call.tsv"
        example_errors = tmp / "example.validation_report.txt"
        example_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(EXAMPLE_VCF),
                "--reference",
                str(EXAMPLE_REFERENCE),
                "--call-tsv",
                str(example_call_tsv),
                "--error-report",
                str(example_errors),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if example_result.returncode != 0:
            raise AssertionError(
                "Synthetic human example failed\n"
                f"STDOUT:\n{example_result.stdout}\n"
                f"STDERR:\n{example_result.stderr}\n"
            )

        example_text = example_call_tsv.read_text(encoding="utf-8")
        example_header, example_rows = assert_submission_call_schema(example_text)
        example_calls = {
            row["Variant_Call_ID"]: row
            for row in (dict(zip(example_header, fields)) for fields in example_rows)
        }
        assert len(example_rows) == 5
        assert "human_del1" in example_text
        assert example_calls["human_del1"]["BND_Source_VCF_IDs"] == "."
        assert example_calls["human_ins1"]["BND_Source_VCF_IDs"] == "."
        assert example_calls["human_dup1"]["BND_Source_VCF_IDs"] == "."
        assert example_calls["human_inv1"]["BND_Source_VCF_IDs"] == "."
        assert (
            example_calls["human_bnd1"]["BND_Source_VCF_IDs"]
            == "human_bnd1,human_bnd2"
        )
        assert "kssv" not in example_text
        example_report = example_errors.read_text(encoding="utf-8")
        assert "Total issues: 0" in example_report
        assert "No issues." in example_report

        report_failure_result = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--vcf",
                str(vcf),
                "--reference",
                str(reference),
                "--call-tsv",
                str(report_failure_call_tsv),
                "--error-report",
                str(missing_report),
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
