#!/usr/bin/env python3
"""Smoke tests for the public SNP VCF validation tools."""

import gzip
import os
import subprocess
import shutil
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "src" / "kvar_snp_tools" / "Sub_validator.py"
EXAMPLES = ROOT / "examples"


def run_command(args, *, env=None):
    """Run a command and fail with useful output."""
    result = subprocess.run(args, text=True, capture_output=True, check=False, env=env)
    if result.returncode != 0:
        raise AssertionError(
            f"Command failed: {' '.join(str(arg) for arg in args)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result


def run_command_expect_failure(args, *, env=None):
    """Run a command and require a non-zero exit status."""
    result = subprocess.run(args, text=True, capture_output=True, check=False, env=env)
    if result.returncode == 0:
        raise AssertionError(
            f"Command unexpectedly succeeded: {' '.join(str(arg) for arg in args)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result


def without_python_site_packages():
    """Return an environment that cannot inherit a caller-provided PYTHONPATH."""
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    return env


def test_cli_exposes_only_unified_report_option(tmp_dir):
    """Both commands expose one report option and reject the legacy reference report."""
    for command in ("generic-to-dbsnp", "validate-dbsnp"):
        result = run_command([sys.executable, str(CLI), command, "--help"])
        assert "--error-report" in result.stdout
        assert "--reference" in result.stdout
        assert "--reference-report" not in result.stdout

    result = run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(EXAMPLES / "toy.dbsnp.vcf"),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-r",
            str(EXAMPLES / "toy.reference.fa"),
            "-rr",
            str(tmp_dir / "legacy.reference.txt"),
            "-o",
            str(tmp_dir / "legacy.cleaned.vcf"),
            "-e",
            str(tmp_dir / "legacy.errors.txt"),
        ]
    )
    assert "unrecognized arguments" in result.stderr


def test_output_and_report_paths_must_be_distinct(tmp_dir):
    """Generated files must never overwrite each other or an input file."""
    command_inputs = {
        "generic-to-dbsnp": EXAMPLES / "toy.generic.vcf",
        "validate-dbsnp": EXAMPLES / "toy.dbsnp.vcf",
    }
    for command, input_vcf in command_inputs.items():
        shared_path = tmp_dir / f"{command}.same-path.vcf"
        result = run_command_expect_failure(
            [
                sys.executable,
                str(CLI),
                command,
                "-v",
                str(input_vcf),
                "-m",
                str(EXAMPLES / "toy.metadata.txt"),
                "-o",
                str(shared_path),
                "-e",
                str(shared_path),
            ]
        )

        assert "--output and --error-report must use different paths" in result.stderr
        assert not shared_path.exists()

    input_vcf = tmp_dir / "preserved-input.vcf"
    shutil.copyfile(EXAMPLES / "toy.dbsnp.vcf", input_vcf)
    original_input = input_vcf.read_text(encoding="utf-8")
    result = run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(input_vcf),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-o",
            str(tmp_dir / "preserved-input.cleaned.vcf"),
            "-e",
            str(input_vcf),
        ]
    )
    assert "--error-report must use a different path from --vcf" in result.stderr
    assert input_vcf.read_text(encoding="utf-8") == original_input

    default_output = tmp_dir / "default-collision.vcf"
    colliding_metadata = tmp_dir / "default-collision.errors.txt"
    shutil.copyfile(EXAMPLES / "toy.metadata.txt", colliding_metadata)
    original_metadata = colliding_metadata.read_text(encoding="utf-8")
    result = run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(EXAMPLES / "toy.dbsnp.vcf"),
            "-m",
            str(colliding_metadata),
            "-o",
            str(default_output),
        ]
    )
    assert "--error-report must use a different path from --metadata" in result.stderr
    assert colliding_metadata.read_text(encoding="utf-8") == original_metadata

    hard_link_input = tmp_dir / "hard-link-input.vcf"
    hard_link_output = tmp_dir / "hard-link-output.vcf"
    shutil.copyfile(EXAMPLES / "toy.generic.vcf", hard_link_input)
    original_hard_link_input = hard_link_input.read_text(encoding="utf-8")
    os.link(hard_link_input, hard_link_output)
    result = run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "generic-to-dbsnp",
            "-v",
            str(hard_link_input),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-o",
            str(hard_link_output),
            "-e",
            str(tmp_dir / "hard-link.errors.txt"),
        ]
    )
    assert "--output must use a different path from --vcf" in result.stderr
    assert hard_link_input.read_text(encoding="utf-8") == original_hard_link_input


def test_default_report_name_is_consistent(tmp_dir):
    """Both SNP commands derive the same .errors.txt report suffix."""
    command_inputs = {
        "generic-to-dbsnp": EXAMPLES / "toy.generic.vcf",
        "validate-dbsnp": EXAMPLES / "toy.dbsnp.vcf",
    }
    for command, input_vcf in command_inputs.items():
        output_vcf = tmp_dir / f"{command}.default.vcf"
        expected_report = tmp_dir / f"{command}.default.errors.txt"
        legacy_report = tmp_dir / f"{command}.default_errors.txt"
        run_command(
            [
                sys.executable,
                str(CLI),
                command,
                "-v",
                str(input_vcf),
                "-m",
                str(EXAMPLES / "toy.metadata.txt"),
                "-o",
                str(output_vcf),
            ]
        )

        assert output_vcf.exists()
        assert expected_report.exists()
        assert not legacy_report.exists()


def test_default_output_name_is_command_specific(tmp_dir):
    """Each SNP command derives a distinct VCF name and preserves gzip output."""
    generic_input = tmp_dir / "cohort.vcf"
    shutil.copyfile(EXAMPLES / "toy.generic.vcf", generic_input)
    run_command(
        [
            sys.executable,
            str(CLI),
            "generic-to-dbsnp",
            "-v",
            str(generic_input),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
        ]
    )

    generic_output = tmp_dir / "cohort.dbsnp.vcf"
    generic_report = tmp_dir / "cohort.dbsnp.errors.txt"
    assert generic_output.exists()
    assert generic_report.exists()

    dbsnp_input = tmp_dir / "submission.dbsnp.vcf.gz"
    with gzip.open(dbsnp_input, "wt", encoding="utf-8") as handle:
        handle.write((EXAMPLES / "toy.dbsnp.vcf").read_text(encoding="utf-8"))

    run_command(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(dbsnp_input),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
        ]
    )

    cleaned_output = tmp_dir / "submission.dbsnp.cleaned.vcf.gz"
    cleaned_report = tmp_dir / "submission.dbsnp.cleaned.errors.txt"
    assert cleaned_output.exists()
    assert cleaned_report.exists()
    with gzip.open(cleaned_output, "rt", encoding="utf-8") as handle:
        assert handle.readline().strip() == "##fileformat=VCFv4.1"


def test_explicit_output_requires_vcf_suffix(tmp_dir):
    """Explicit SNP output names must identify VCF or compressed VCF files."""
    invalid_output = tmp_dir / "invalid-output.txt"
    result = run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "generic-to-dbsnp",
            "-v",
            str(EXAMPLES / "toy.generic.vcf"),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-o",
            str(invalid_output),
        ]
    )

    assert "--output must end with .vcf or .vcf.gz" in result.stderr
    assert not invalid_output.exists()


def test_early_failures_still_write_the_unified_report(tmp_dir):
    """Failures before row validation still produce the requested single report."""
    missing_metadata = tmp_dir / "missing.metadata.txt"

    generic_output = tmp_dir / "early.generic.vcf"
    generic_report = tmp_dir / "early.generic.report.txt"
    run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "generic-to-dbsnp",
            "-v",
            str(EXAMPLES / "toy.generic.vcf"),
            "-m",
            str(missing_metadata),
            "-o",
            str(generic_output),
            "-e",
            str(generic_report),
        ]
    )
    assert not generic_output.exists()
    assert generic_report.exists()
    assert "FILE_NOT_FOUND" in generic_report.read_text(encoding="utf-8")

    cleaned_output = tmp_dir / "early.cleaned.vcf"
    default_report = tmp_dir / "early.cleaned.errors.txt"
    run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(EXAMPLES / "toy.dbsnp.vcf"),
            "-m",
            str(missing_metadata),
            "-o",
            str(cleaned_output),
        ]
    )
    assert not cleaned_output.exists()
    assert default_report.exists()
    assert "FILE_NOT_FOUND" in default_report.read_text(encoding="utf-8")


def test_output_write_failure_is_recorded_in_the_unified_report(tmp_dir):
    """A failed output publish must be visible as a validation issue."""
    output_vcf = tmp_dir / "missing-output-dir" / "cleaned.vcf"
    error_report = tmp_dir / "write-failure.errors.txt"

    run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(EXAMPLES / "toy.dbsnp.vcf"),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    assert not output_vcf.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "FILE_WRITE_ERROR" in report_text
    assert "Total issues: 1" in report_text


def test_metadata_argument_is_required(tmp_dir):
    """Both public SNP commands require the metadata argument."""
    for command in ("generic-to-dbsnp", "validate-dbsnp"):
        result = run_command_expect_failure(
            [
                sys.executable,
                str(CLI),
                command,
                "-v",
                str(EXAMPLES / "toy.generic.vcf"),
                "-o",
                str(tmp_dir / f"{command}.vcf"),
            ]
        )
        assert "-m/--metadata" in result.stderr


def test_required_metadata_ids_block_output(tmp_dir):
    """Missing Experiment_id and SampleSet_id block dbSNP VCF creation."""
    metadata = tmp_dir / "missing-required-ids.metadata.txt"
    output_vcf = tmp_dir / "missing-required-ids.vcf"
    error_report = tmp_dir / "missing-required-ids.errors.txt"
    metadata.write_text("##reference=toy_ref\n", encoding="utf-8")

    run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "generic-to-dbsnp",
            "-v",
            str(EXAMPLES / "toy.generic.vcf"),
            "-m",
            str(metadata),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    assert not output_vcf.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "MISSING_REQUIRED_METADATA" in report_text
    assert "Experiment_id" in report_text
    assert "SampleSet_id" in report_text


def test_metadata_requires_exactly_one_sampleset_id(tmp_dir):
    """Missing, repeated, or multiple SampleSet IDs block conversion."""
    cases = {
        "missing": [],
        "multiple": ["POP1", "POP2"],
        "duplicate": ["POP1", "POP1"],
    }

    for case_name, sampleset_ids in cases.items():
        metadata = tmp_dir / f"{case_name}.metadata.txt"
        output_vcf = tmp_dir / f"{case_name}.dbsnp.vcf"
        error_report = tmp_dir / f"{case_name}.errors.txt"
        metadata_lines = [
            "##Experiment_id=EXP001",
            "##reference=toy_ref",
            *[f"##SampleSet_id={sampleset_id}" for sampleset_id in sampleset_ids],
        ]
        metadata.write_text("\n".join(metadata_lines) + "\n", encoding="utf-8")

        run_command_expect_failure(
            [
                sys.executable,
                str(CLI),
                "generic-to-dbsnp",
                "-v",
                str(EXAMPLES / "toy.generic.vcf"),
                "-m",
                str(metadata),
                "-o",
                str(output_vcf),
                "-e",
                str(error_report),
            ]
        )

        assert not output_vcf.exists()
        report_text = error_report.read_text(encoding="utf-8")
        if case_name == "missing":
            assert "MISSING_REQUIRED_METADATA" in report_text
        else:
            assert "DUPLICATE_METADATA_TAG" in report_text
        assert "SampleSet_id" in report_text


def test_validate_dbsnp_requires_matching_single_population_column(tmp_dir):
    """The sole population column must match the sole population_id."""
    input_vcf = tmp_dir / "mismatched-population-column.dbsnp.vcf"
    output_vcf = tmp_dir / "mismatched-population-column.cleaned.vcf"
    error_report = tmp_dir / "mismatched-population-column.errors.txt"
    input_text = (EXAMPLES / "toy.dbsnp.vcf").read_text(encoding="utf-8")
    input_text = input_text.replace(
        "\tFORMAT\tPOP1\n",
        "\tFORMAT\tWRONG_POPULATION\n",
    )
    input_vcf.write_text(input_text, encoding="utf-8")

    run_command_expect_failure(
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

    assert not output_vcf.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "POPULATION_ID_MISMATCH" in report_text
    assert "WRONG_POPULATION" in report_text


def test_validate_dbsnp_rejects_multiple_populations(tmp_dir):
    """A dbSNP VCF cannot define more than one population."""
    input_vcf = tmp_dir / "multiple-populations.dbsnp.vcf"
    output_vcf = tmp_dir / "multiple-populations.cleaned.vcf"
    error_report = tmp_dir / "multiple-populations.errors.txt"
    input_text = (EXAMPLES / "toy.dbsnp.vcf").read_text(encoding="utf-8")
    input_text = input_text.replace(
        "##population_id=POP1\n",
        "##population_id=POP1\n##population_id=POP2\n",
    )
    input_text = input_text.replace(
        "\tFORMAT\tPOP1\n",
        "\tFORMAT\tPOP1\tPOP2\n",
    )
    input_text = input_text.replace(
        "\tNA:FRQ\t10:0.2\n",
        "\tNA:FRQ\t10:0.2\t10:0.2\n",
    )
    input_vcf.write_text(input_text, encoding="utf-8")

    run_command_expect_failure(
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

    assert not output_vcf.exists()
    assert "POPULATION_ID_MISMATCH" in error_report.read_text(encoding="utf-8")


def test_short_generic_row_blocks_output(tmp_dir):
    """A malformed generic VCF row must be reported instead of silently dropped."""
    input_vcf = tmp_dir / "short-row.generic.vcf"
    output_vcf = tmp_dir / "short-row.dbsnp.vcf"
    error_report = tmp_dir / "short-row.errors.txt"
    input_lines = (EXAMPLES / "toy.generic.vcf").read_text(encoding="utf-8").splitlines()
    valid_row_index = next(
        index for index, line in enumerate(input_lines) if line and not line.startswith("#")
    )
    valid_row = input_lines[valid_row_index]
    input_lines.insert(valid_row_index, "\t".join(valid_row.split("\t")[:7]))
    input_vcf.write_text("\n".join(input_lines) + "\n", encoding="utf-8")

    run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "generic-to-dbsnp",
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

    assert not output_vcf.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "INSUFFICIENT_FIELDS: 1" in report_text


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
            "-r",
            str(EXAMPLES / "toy.reference.fa"),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    output_text = output_vcf.read_text(encoding="utf-8")
    assert output_text.startswith("##fileformat=VCFv4.1\n")
    assert "##batch=EXP001" in output_text
    assert "##population_id=POP1" in output_text
    assert "##contig=<ID=chr1,length=20>" in output_text
    assert "##SampleSet_id=" not in output_text
    assert "VRT=1" in output_text
    assert "AF=0.2" in output_text
    for format_id in ("NA", "NS", "FRQ", "AC"):
        assert output_text.count(f"##FORMAT=<ID={format_id},") == 1
    assert error_report.exists()
    assert "=== Reference Validation ===" in error_report.read_text(encoding="utf-8")


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
        'Description="Number of alleles for the population"',
        r'Description="Number of alleles for the population \tmp"',
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
    assert output_text.startswith("##fileformat=VCFv4.1\n")
    assert "##batch=EXP001" in output_text
    assert "##population_id=POP1" in output_text
    assert "##contig=<ID=chr1,length=20>" in output_text
    assert "##SampleSet_id=" not in output_text
    assert output_text.count("##INFO=<ID=VRT") == 1
    assert r'Description="Number of alleles for the population \\tmp"' in output_text
    assert "\tNA:FRQ\t10:0.2" in output_text
    assert "CSQ" not in output_text
    assert "61KJPN_AC" not in output_text
    assert "JPN61K_AC" not in output_text
    assert "UNSUPPORTED_DBSNP_INFO_TAG" in error_report.read_text(encoding="utf-8")
    assert error_report.exists()


def test_validate_dbsnp_with_reference_check_uses_unified_report(tmp_dir):
    """Reference statistics and ordinary validation share one report."""
    input_vcf = tmp_dir / "toy.dbsnp.vcf"
    reference_fasta = tmp_dir / "toy.reference.fa"
    output_vcf = tmp_dir / "dbsnp.refchecked.cleaned.vcf"
    error_report = tmp_dir / "dbsnp.refchecked.cleaned.errors.txt"
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
            "-o",
            str(output_vcf),
        ]
    )

    assert output_vcf.exists()
    assert error_report.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "=== Reference Validation ===" in report_text
    assert "Status: COMPLETED" in report_text
    assert "Total variants: 1" in report_text
    assert "Matched: 1" in report_text
    assert "Mismatched: 0" in report_text
    assert not list(tmp_dir.glob("*_ref_check_report.txt"))


def test_reference_mismatch_and_general_validation_share_report(tmp_dir):
    """Reference and ordinary validation issues are reported before output is blocked."""
    input_vcf = tmp_dir / "mismatch.dbsnp.vcf"
    reference_fasta = tmp_dir / "toy.reference.fa"
    output_vcf = tmp_dir / "mismatch.cleaned.vcf"
    error_report = tmp_dir / "mismatch.errors.txt"
    input_text = (EXAMPLES / "toy.dbsnp.vcf").read_text(encoding="utf-8")
    input_text = input_text.replace(
        "##FORMAT=<ID=NA",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Unsupported annotation">\n'
        "##FORMAT=<ID=NA",
    )
    input_text = input_text.replace(
        "chr1\t2\tvar1\tC\tT",
        "chr1\t2\tvar1\tA\tT",
    ).replace(
        "AF=0.2\tNA:FRQ",
        "AF=0.2;CSQ=T|missense_variant\tNA:FRQ",
    )
    input_vcf.write_text(input_text, encoding="utf-8")
    shutil.copyfile(EXAMPLES / "toy.reference.fa", reference_fasta)

    run_command_expect_failure(
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
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    assert not output_vcf.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "REF_MISMATCH" in report_text
    assert "UNSUPPORTED_DBSNP_INFO_TAG" in report_text
    assert "=== Reference Validation ===" in report_text


def test_missing_reference_fasta_is_written_to_unified_report(tmp_dir):
    """A requested reference that cannot be loaded blocks output with one report."""
    output_vcf = tmp_dir / "missing-reference.cleaned.vcf"
    error_report = tmp_dir / "missing-reference.errors.txt"

    run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(EXAMPLES / "toy.dbsnp.vcf"),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-r",
            str(tmp_dir / "missing.reference.fa"),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    assert not output_vcf.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "FILE_NOT_FOUND" in report_text
    assert "=== Reference Validation ===" in report_text
    assert "Status: FAILED" in report_text


def test_missing_vcf_is_not_double_counted_by_reference_validation(tmp_dir):
    """The main validator owns a missing input VCF issue exactly once."""
    output_vcf = tmp_dir / "missing-vcf.cleaned.vcf"
    error_report = tmp_dir / "missing-vcf.errors.txt"

    run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(tmp_dir / "missing.vcf"),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-r",
            str(EXAMPLES / "toy.reference.fa"),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    assert not output_vcf.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "Status: FAILED" in report_text
    assert "Total issues: 1" in report_text
    assert "FILE_NOT_FOUND: 1" in report_text


def test_unreadable_vcf_is_not_double_counted_by_reference_validation(tmp_dir):
    """The main validator owns an unreadable input VCF issue exactly once."""
    input_vcf = tmp_dir / "corrupt.vcf.gz"
    output_vcf = tmp_dir / "corrupt.cleaned.vcf"
    error_report = tmp_dir / "corrupt.errors.txt"
    input_vcf.write_bytes(b"not a gzip stream\n")

    run_command_expect_failure(
        [
            sys.executable,
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(input_vcf),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-r",
            str(EXAMPLES / "toy.reference.fa"),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ]
    )

    assert not output_vcf.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "Status: FAILED" in report_text
    assert "Total issues: 1" in report_text
    assert "FILE_READ_ERROR: 1" in report_text


def test_reference_dependency_failure_is_written_to_unified_report(tmp_dir):
    """Missing pyfaidx is reported when reference validation was requested."""
    output_vcf = tmp_dir / "missing-pyfaidx.cleaned.vcf"
    error_report = tmp_dir / "missing-pyfaidx.errors.txt"

    run_command_expect_failure(
        [
            sys.executable,
            "-S",
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(EXAMPLES / "toy.dbsnp.vcf"),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-r",
            str(EXAMPLES / "toy.reference.fa"),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ],
        env=without_python_site_packages(),
    )

    assert not output_vcf.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "FASTA_READ_ERROR" in report_text
    assert "pyfaidx library is required" in report_text


def test_no_reference_does_not_require_pyfaidx(tmp_dir):
    """SNP validation without --reference has no optional dependency."""
    output_vcf = tmp_dir / "no-reference.cleaned.vcf"
    error_report = tmp_dir / "no-reference.errors.txt"

    run_command(
        [
            sys.executable,
            "-S",
            str(CLI),
            "validate-dbsnp",
            "-v",
            str(EXAMPLES / "toy.dbsnp.vcf"),
            "-m",
            str(EXAMPLES / "toy.metadata.txt"),
            "-o",
            str(output_vcf),
            "-e",
            str(error_report),
        ],
        env=without_python_site_packages(),
    )

    assert output_vcf.exists()
    assert error_report.exists()
    report_text = error_report.read_text(encoding="utf-8")
    assert "KVar SNP Validation Report" in report_text
    assert "=== Validation Summary ===" in report_text
    assert "Total issues: 0" in report_text
    assert "No issues." in report_text
    assert "=== Reference Validation ===" not in report_text


def main():
    """Run smoke tests without requiring pytest."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        test_cli_exposes_only_unified_report_option(tmp_dir)
        test_output_and_report_paths_must_be_distinct(tmp_dir)
        test_default_report_name_is_consistent(tmp_dir)
        test_default_output_name_is_command_specific(tmp_dir)
        test_explicit_output_requires_vcf_suffix(tmp_dir)
        test_early_failures_still_write_the_unified_report(tmp_dir)
        test_output_write_failure_is_recorded_in_the_unified_report(tmp_dir)
        test_metadata_argument_is_required(tmp_dir)
        test_required_metadata_ids_block_output(tmp_dir)
        test_metadata_requires_exactly_one_sampleset_id(tmp_dir)
        test_validate_dbsnp_requires_matching_single_population_column(tmp_dir)
        test_validate_dbsnp_rejects_multiple_populations(tmp_dir)
        test_short_generic_row_blocks_output(tmp_dir)
        test_generic_to_dbsnp_writes_population_id_only(tmp_dir)
        test_validate_dbsnp_writes_cleaned_vcf(tmp_dir)
        test_validate_dbsnp_with_reference_check_uses_unified_report(tmp_dir)
        test_reference_mismatch_and_general_validation_share_report(tmp_dir)
        test_missing_reference_fasta_is_written_to_unified_report(tmp_dir)
        test_missing_vcf_is_not_double_counted_by_reference_validation(tmp_dir)
        test_unreadable_vcf_is_not_double_counted_by_reference_validation(tmp_dir)
        test_reference_dependency_failure_is_written_to_unified_report(tmp_dir)
        test_no_reference_does_not_require_pyfaidx(tmp_dir)
    print("Public CLI smoke tests passed.")


if __name__ == "__main__":
    main()
