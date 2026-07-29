#!/usr/bin/env python3
"""Regression tests for SV reference indexing and coordinate bounds."""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "src" / "kvar_sv_tools" / "vcf_to_kvar_tsv.py"
sys.path.insert(0, str(ROOT / "src"))

from kvar_sv_tools.VCF_parser import FastaReference


def write_reference(path: Path) -> None:
    """Write two 240-base contigs without an index."""
    lines = [
        ">chr1",
        *("A" * 60 for _ in range(4)),
        ">chr2",
        *("A" * 60 for _ in range(4)),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_vcf(path: Path, record: str) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=240>\n"
        "##contig=<ID=chr2,length=240>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate">\n'
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + record
        + "\n",
        encoding="utf-8",
    )


def run_cli(
    vcf: Path,
    reference: Path,
    call_tsv: Path,
    report: Path,
) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            str(CLI),
            "--vcf",
            str(vcf),
            "--reference",
            str(reference),
            "--call-tsv",
            str(call_tsv),
            "--error-report",
            str(report),
        ],
        cwd=str(ROOT),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


class ReferenceBoundsAndIndexTests(unittest.TestCase):
    def test_missing_index_is_created_and_used(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.fasta"
            vcf = tmp / "valid.vcf"
            call_tsv = tmp / "valid.call.tsv"
            report = tmp / "valid.errors.txt"
            write_reference(reference)
            write_vcf(
                vcf,
                "chr1\t20\tdel1\tA\t<DEL>\t.\tPASS\t"
                "SVTYPE=DEL;END=30;SVLEN=-10",
            )

            result = run_cli(vcf, reference, call_tsv, report)

            self.assertEqual(0, result.returncode, result.stderr)
            self.assertTrue(Path(f"{reference}.fai").exists())
            self.assertTrue(call_tsv.exists())
            self.assertIn("Total issues: 0", report.read_text(encoding="utf-8"))
            indexed = FastaReference(str(reference))
            self.assertEqual(240, indexed.length("chr1"))
            self.assertEqual("A" * 62, indexed.fetch("chr1", 59, 120))

    def test_repeated_fetch_reuses_one_open_fasta_handle(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "reference.fasta"
            reference.write_text(">chr1\nACGT\nACGT\n", encoding="ascii")
            indexed = FastaReference(str(reference))

            with mock.patch("builtins.open", wraps=open) as open_mock:
                self.assertEqual("ACGT", indexed.fetch("chr1", 1, 4))
                self.assertEqual("ACGT", indexed.fetch("chr1", 5, 8))
                indexed.close()

            self.assertEqual(1, open_mock.call_count)

    def test_stale_index_is_rebuilt(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "reference.fasta"
            write_reference(reference)
            indexed = FastaReference(str(reference))
            self.assertEqual(240, indexed.length("chr2"))

            with reference.open("a", encoding="ascii") as fasta:
                fasta.write("A\n")
            index_path = Path(f"{reference}.fai")
            os.utime(index_path, (1, 1))

            rebuilt = FastaReference(str(reference))
            self.assertEqual(241, rebuilt.length("chr2"))
            self.assertIn("\t241\t", index_path.read_text(encoding="utf-8"))

    def test_newer_but_semantically_wrong_index_is_rebuilt(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "reference.fasta"
            reference.write_text(">chr1\n" + "A" * 60 + "\n", encoding="ascii")
            index_path = Path(f"{reference}.fai")
            index_path.write_text(
                "chr1\t999\t6\t60\t61\n",
                encoding="utf-8",
            )
            newer_mtime = reference.stat().st_mtime_ns + 1_000_000_000
            os.utime(index_path, ns=(newer_mtime, newer_mtime))

            rebuilt = FastaReference(str(reference))

            self.assertEqual(60, rebuilt.length("chr1"))
            self.assertEqual(
                "chr1\t60\t6\t60\t61\n",
                index_path.read_text(encoding="utf-8"),
            )

    def test_newer_index_with_wrong_line_metrics_is_rebuilt(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "reference.fasta"
            reference.write_text(">chr1\nACGT\nACGT\n", encoding="ascii")
            index_path = Path(f"{reference}.fai")
            index_path.write_text(
                "chr1\t8\t6\t5\t6\n",
                encoding="utf-8",
            )
            newer_mtime = reference.stat().st_mtime_ns + 1_000_000_000
            os.utime(index_path, ns=(newer_mtime, newer_mtime))

            rebuilt = FastaReference(str(reference))

            self.assertEqual("ACGTACGT", rebuilt.fetch("chr1", 1, 8))
            self.assertEqual(
                "chr1\t8\t6\t4\t5\n",
                index_path.read_text(encoding="utf-8"),
            )

    def test_newer_index_missing_trailing_contig_is_rebuilt(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "reference.fasta"
            write_reference(reference)
            index_path = Path(f"{reference}.fai")
            index_path.write_text(
                "chr1\t240\t6\t60\t61\n",
                encoding="utf-8",
            )
            newer_mtime = reference.stat().st_mtime_ns + 1_000_000_000
            os.utime(index_path, ns=(newer_mtime, newer_mtime))

            rebuilt = FastaReference(str(reference))

            self.assertEqual(["chr1", "chr2"], list(rebuilt.index))
            self.assertEqual(240, rebuilt.length("chr2"))
            self.assertEqual(
                2,
                len(index_path.read_text(encoding="utf-8").splitlines()),
            )

    def test_newer_index_missing_middle_contig_is_rebuilt(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "reference.fasta"
            reference.write_text(
                ">chr1\nAAAAAAAAAA\n"
                ">chr2\nCCCCCCCCCC\n"
                ">chr3\nGGGGGGGGGG\n",
                encoding="ascii",
            )
            index_path = Path(f"{reference}.fai")
            index_path.write_text(
                "chr1\t10\t6\t10\t11\n"
                "chr3\t10\t40\t10\t11\n",
                encoding="utf-8",
            )
            newer_mtime = reference.stat().st_mtime_ns + 1_000_000_000
            os.utime(index_path, ns=(newer_mtime, newer_mtime))

            rebuilt = FastaReference(str(reference))

            self.assertEqual(["chr1", "chr2", "chr3"], list(rebuilt.index))
            self.assertEqual("CCCC", rebuilt.fetch("chr2", 4, 7))
            self.assertEqual(
                3,
                len(index_path.read_text(encoding="utf-8").splitlines()),
            )

    def test_index_cache_write_failure_uses_memory_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "reference.fasta"
            write_reference(reference)

            with mock.patch(
                "kvar_sv_tools.VCF_parser.os.replace",
                side_effect=PermissionError("read-only directory"),
            ):
                indexed = FastaReference(str(reference))

            self.assertEqual(240, indexed.length("chr2"))
            self.assertEqual("AAAA", indexed.fetch("chr2", 237, 240))
            self.assertFalse(Path(f"{reference}.fai").exists())

    def test_explicit_end_outside_reference_blocks_output_once(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.fasta"
            vcf = tmp / "end_outside.vcf"
            call_tsv = tmp / "end_outside.call.tsv"
            report = tmp / "end_outside.errors.txt"
            write_reference(reference)
            write_vcf(
                vcf,
                "chr1\t20\tdel_outside\tA\t<DEL>\t.\tPASS\t"
                "SVTYPE=DEL;END=9999;SVLEN=-10",
            )

            result = run_cli(vcf, reference, call_tsv, report)

            self.assertNotEqual(0, result.returncode)
            self.assertFalse(call_tsv.exists())
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("POSITION_OUT_OF_RANGE: 1", report_text)
            self.assertIn("field: END", report_text)

    def test_contig_header_alone_blocks_out_of_range_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vcf = tmp / "header_end_outside.vcf"
            call_tsv = tmp / "header_end_outside.call.tsv"
            report = tmp / "header_end_outside.errors.txt"
            write_vcf(
                vcf,
                "chr1\t20\tdel_outside\tA\t<DEL>\t.\tPASS\t"
                "SVTYPE=DEL;END=9999;SVLEN=-10",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(CLI),
                    "--vcf",
                    str(vcf),
                    "--call-tsv",
                    str(call_tsv),
                    "--error-report",
                    str(report),
                ],
                cwd=str(ROOT),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            self.assertNotEqual(0, result.returncode)
            self.assertFalse(call_tsv.exists())
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("POSITION_OUT_OF_RANGE: 1", report_text)
            self.assertIn("source: VCF contig header", report_text)

    def test_contig_header_checks_full_ref_span_without_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vcf = tmp / "ref_span_outside.vcf"
            call_tsv = tmp / "ref_span_outside.call.tsv"
            report = tmp / "ref_span_outside.errors.txt"
            write_vcf(
                vcf,
                "chr1\t239\tref_span_outside\tAAA\t<DEL>\t.\tPASS\t"
                "SVTYPE=DEL;END=240;SVLEN=-1",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(CLI),
                    "--vcf",
                    str(vcf),
                    "--call-tsv",
                    str(call_tsv),
                    "--error-report",
                    str(report),
                ],
                cwd=str(ROOT),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            self.assertNotEqual(0, result.returncode)
            self.assertFalse(call_tsv.exists())
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("POSITION_OUT_OF_RANGE: 1", report_text)
            self.assertIn("field: POS/REF", report_text)
            self.assertIn("actual: chr1:239-241", report_text)

    def test_converter_derived_end_outside_reference_blocks_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.fasta"
            vcf = tmp / "derived_end_outside.vcf"
            call_tsv = tmp / "derived_end_outside.call.tsv"
            report = tmp / "derived_end_outside.errors.txt"
            write_reference(reference)
            write_vcf(
                vcf,
                "chr1\t230\tderived_outside\tA\t<DEL>\t.\tPASS\t"
                "SVTYPE=DEL;SVLEN=-20",
            )

            result = run_cli(vcf, reference, call_tsv, report)

            self.assertNotEqual(0, result.returncode)
            self.assertFalse(call_tsv.exists())
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("POSITION_OUT_OF_RANGE", report_text)
            self.assertIn("field: derived END", report_text)

    def test_insertion_end_equals_pos_with_svlen_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.fasta"
            vcf = tmp / "insertion.vcf"
            call_tsv = tmp / "insertion.call.tsv"
            report = tmp / "insertion.errors.txt"
            write_reference(reference)
            write_vcf(
                vcf,
                "chr1\t90\tins1\tA\t<INS>\t.\tPASS\t"
                "SVTYPE=INS;END=90;SVLEN=12",
            )

            result = run_cli(vcf, reference, call_tsv, report)

            self.assertEqual(0, result.returncode, result.stderr)
            data_row = next(
                line
                for line in call_tsv.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("#")
            ).split("\t")
            self.assertEqual("90", data_row[4])
            self.assertEqual("90", data_row[7])
            self.assertEqual("12", data_row[9])

    def test_breakend_remote_coordinate_outside_length_plus_one_blocks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.fasta"
            vcf = tmp / "bnd_outside.vcf"
            call_tsv = tmp / "bnd_outside.call.tsv"
            report = tmp / "bnd_outside.errors.txt"
            write_reference(reference)
            write_vcf(
                vcf,
                "chr1\t10\tbnd_outside\tA\tA]chr2:242]\t.\tPASS\tSVTYPE=BND",
            )

            result = run_cli(vcf, reference, call_tsv, report)

            self.assertNotEqual(0, result.returncode)
            self.assertFalse(call_tsv.exists())
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("POSITION_OUT_OF_RANGE", report_text)
            self.assertIn("field: ALT remote coordinate", report_text)
            self.assertIn("expected: 1-241", report_text)

    def test_breakend_remote_coordinate_at_length_plus_one_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.fasta"
            vcf = tmp / "bnd_telomere.vcf"
            call_tsv = tmp / "bnd_telomere.call.tsv"
            report = tmp / "bnd_telomere.errors.txt"
            write_reference(reference)
            write_vcf(
                vcf,
                "chr1\t10\tbnd_telomere\tA\tA]chr2:241]\t.\tPASS\tSVTYPE=BND",
            )

            result = run_cli(vcf, reference, call_tsv, report)

            self.assertEqual(0, result.returncode, result.stderr)
            self.assertTrue(call_tsv.exists())
            self.assertIn("Total issues: 0", report.read_text(encoding="utf-8"))

    def test_single_breakend_checks_only_local_position(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.fasta"
            vcf = tmp / "single_bnd.vcf"
            call_tsv = tmp / "single_bnd.call.tsv"
            report = tmp / "single_bnd.errors.txt"
            write_reference(reference)
            write_vcf(
                vcf,
                "chr1\t10\tsingle_bnd\tA\tA.\t.\tPASS\tSVTYPE=BND",
            )

            result = run_cli(vcf, reference, call_tsv, report)

            self.assertEqual(0, result.returncode, result.stderr)
            self.assertTrue(call_tsv.exists())

    def test_bnd_derived_insertion_end_outside_reference_blocks_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "reference.fasta"
            vcf = tmp / "bnd_insertion_outside.vcf"
            call_tsv = tmp / "bnd_insertion_outside.call.tsv"
            report = tmp / "bnd_insertion_outside.errors.txt"
            write_reference(reference)
            write_vcf(
                vcf,
                "chr1\t235\tbnd_insertion_outside\tA\tACCCCCCCCCC.\t.\tPASS\t"
                "SVTYPE=BND",
            )

            result = run_cli(vcf, reference, call_tsv, report)

            self.assertNotEqual(0, result.returncode)
            self.assertFalse(call_tsv.exists())
            self.assertTrue(report.exists())
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("POSITION_OUT_OF_RANGE", report_text)
            self.assertIn("actual: chr1:245", report_text)
            self.assertIn(
                "coordinate: converter-derived BND insertion END",
                report_text,
            )

    def test_invalid_fasta_writes_report_without_traceback(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reference = tmp / "invalid.fasta"
            vcf = tmp / "input.vcf"
            call_tsv = tmp / "invalid.call.tsv"
            report = tmp / "invalid.errors.txt"
            reference.write_text("not-a-fasta\n", encoding="ascii")
            write_vcf(
                vcf,
                "chr1\t20\tdel1\tA\t<DEL>\t.\tPASS\t"
                "SVTYPE=DEL;END=30;SVLEN=-10",
            )

            result = run_cli(vcf, reference, call_tsv, report)

            self.assertNotEqual(0, result.returncode)
            self.assertFalse(call_tsv.exists())
            self.assertTrue(report.exists())
            self.assertIn("FASTA_INDEX_ERROR", report.read_text(encoding="utf-8"))
            self.assertNotIn("Traceback", result.stdout)
            self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main()
