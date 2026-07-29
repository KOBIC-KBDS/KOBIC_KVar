#!/usr/bin/env python3
"""Regression tests for inserted sequence carried by reciprocal BND mates."""

import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Tuple


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kvar_sv_tools.KVar2TSV import KVarTSVConverter
from kvar_sv_tools.error_handler import ErrorCode


VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100>\n"
    "##contig=<ID=chr2,length=100>\n"
    "##contig=<ID=chr:part,length=100>\n"
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
    '##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate breakend ID">\n'
    '##INFO=<ID=EVENT,Number=1,Type=String,Description="Breakend event ID">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def write_pair(path: Path, primary_alt: str, mate_alt: str) -> None:
    path.write_text(
        VCF_HEADER
        + "chr1\t10\tbnd1\tA\t"
        + primary_alt
        + "\t.\tPASS\tSVTYPE=BND;MATEID=bnd2;EVENT=event1\n"
        + "chr1\t20\tbnd2\tT\t"
        + mate_alt
        + "\t.\tPASS\tSVTYPE=BND;MATEID=bnd1;EVENT=event1\n",
        encoding="utf-8",
    )


def read_call_rows(path: Path) -> List[Dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = next(
        line.lstrip("#").split("\t")
        for line in lines
        if line.startswith("#") and not line.startswith("##")
    )
    rows = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != len(header):
            raise AssertionError(
                f"Call TSV row has {len(fields)} fields; expected {len(header)}"
            )
        rows.append(dict(zip(header, fields)))
    return rows


class BndMateInsertionTests(unittest.TestCase):
    def convert(
        self,
        tmp: Path,
        name: str,
        primary_alt: str,
        mate_alt: str,
    ) -> Tuple[Path, Path]:
        vcf = tmp / f"{name}.vcf"
        call_tsv = tmp / f"{name}.call.tsv"
        report = tmp / f"{name}.errors.txt"
        write_pair(vcf, primary_alt, mate_alt)
        KVarTSVConverter().convert_vcf_to_tsv(
            str(vcf),
            str(call_tsv),
            str(report),
        )
        return call_tsv, report

    def test_pair_without_inserted_sequence_produces_one_bnd(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            call_tsv, report = self.convert(
                Path(tmpdir),
                "empty",
                "A]chr1:20]",
                "T]chr1:10]",
            )

            rows = read_call_rows(call_tsv)
            self.assertEqual(1, len(rows))
            self.assertEqual(
                "intrachromosomal translocation",
                rows[0]["Variant_Call_Type"],
            )
            self.assertEqual("bnd1,bnd2", rows[0]["BND_Source_VCF_IDs"])
            report_text = report.read_text(encoding="utf-8")
            self.assertNotIn("MATEID_INSERTION_SEQUENCE_ONE_SIDED", report_text)
            self.assertNotIn("BND_INSERTION_SPLIT", report_text)

    def test_all_reciprocal_orientations_normalize_to_one_sequence(self) -> None:
        cases = (
            ("suffix_open_prefix_close", "AAGT[chr1:20[", "]chr1:10]AGTT"),
            ("suffix_close_suffix_close", "AAGT]chr1:20]", "TACT]chr1:10]"),
            ("prefix_close_suffix_open", "]chr1:20]AGTA", "TAGT[chr1:10["),
            ("prefix_open_prefix_open", "[chr1:20[AGTA", "[chr1:10[ACTT"),
        )

        for name, primary_alt, mate_alt in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                call_tsv, report = self.convert(
                    Path(tmpdir),
                    name,
                    primary_alt,
                    mate_alt,
                )

                rows = read_call_rows(call_tsv)
                self.assertEqual(2, len(rows))
                insertion = next(
                    row for row in rows if row["Variant_Call_Type"] == "insertion"
                )
                translocation = next(
                    row
                    for row in rows
                    if row["Variant_Call_Type"] != "insertion"
                )
                self.assertEqual("AGT", insertion["Sequence"])
                self.assertEqual(
                    "bnd1,bnd2",
                    translocation["BND_Source_VCF_IDs"],
                )
                self.assertEqual(
                    "bnd1,bnd2",
                    insertion["BND_Source_VCF_IDs"],
                )
                report_text = report.read_text(encoding="utf-8")
                self.assertEqual(1, report_text.count("BND_INSERTION_SPLIT"))
                self.assertNotIn(
                    "MATEID_INSERTION_SEQUENCE_ONE_SIDED",
                    report_text,
                )
                self.assertNotIn(
                    "MATEID_INSERTION_SEQUENCE_MISMATCH",
                    report_text,
                )

    def test_mate_only_sequence_is_rescued_at_primary_position(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            call_tsv, report = self.convert(
                Path(tmpdir),
                "mate_only",
                "A]chr1:20]",
                "TACT]chr1:10]",
            )

            rows = read_call_rows(call_tsv)
            self.assertEqual(2, len(rows))
            insertion = next(row for row in rows if row["Variant_Call_Type"] == "insertion")
            self.assertEqual("AGT", insertion["Sequence"])
            self.assertEqual("10", insertion["Start"])
            self.assertEqual("13", insertion["Stop"])
            self.assertEqual("bnd2", insertion["BND_Source_VCF_IDs"])
            report_text = report.read_text(encoding="utf-8")
            self.assertEqual(
                1,
                report_text.count("MATEID_INSERTION_SEQUENCE_ONE_SIDED"),
            )
            self.assertEqual(1, report_text.count("BND_INSERTION_SPLIT"))

    def test_primary_only_sequence_is_preserved_with_one_warning(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            call_tsv, report = self.convert(
                Path(tmpdir),
                "primary_only",
                "AAGT]chr1:20]",
                "T]chr1:10]",
            )

            rows = read_call_rows(call_tsv)
            self.assertEqual(2, len(rows))
            insertion = next(row for row in rows if row["Variant_Call_Type"] == "insertion")
            self.assertEqual("AGT", insertion["Sequence"])
            self.assertEqual("bnd1", insertion["BND_Source_VCF_IDs"])
            report_text = report.read_text(encoding="utf-8")
            self.assertEqual(
                1,
                report_text.count("MATEID_INSERTION_SEQUENCE_ONE_SIDED"),
            )
            self.assertEqual(1, report_text.count("BND_INSERTION_SPLIT"))

    def test_conflicting_sequences_block_output_without_split_warning(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vcf = tmp / "conflict.vcf"
            call_tsv = tmp / "conflict.call.tsv"
            report = tmp / "conflict.errors.txt"
            write_pair(
                vcf,
                "AAGT]chr1:20]",
                "TAAA]chr1:10]",
            )

            with self.assertRaises(RuntimeError):
                KVarTSVConverter().convert_vcf_to_tsv(
                    str(vcf),
                    str(call_tsv),
                    str(report),
                )

            self.assertFalse(call_tsv.exists())
            report_text = report.read_text(encoding="utf-8")
            self.assertEqual(
                1,
                report_text.count("MATEID_INSERTION_SEQUENCE_MISMATCH"),
            )
            self.assertNotIn("BND_INSERTION_SPLIT", report_text)

    def test_unsupported_mate_shapes_block_without_partial_collapse(self) -> None:
        cases = (
            (
                "multiple_mateids",
                "chr1\t10\tbnd1\tA\tAAGT]chr1:20]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd2,bnd3\n"
                "chr1\t20\tbnd2\tT\tT]chr1:10]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd1\n",
                ErrorCode.MULTIPLE_MATEIDS_UNSUPPORTED,
            ),
            (
                "self_mateid",
                "chr1\t10\tbnd1\tA\tAAGT]chr1:10]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd1\n",
                ErrorCode.MATEID_SELF_REFERENCE,
            ),
            (
                "multiple_alts",
                "chr1\t10\tbnd1\tA\tA]chr1:20],A]chr1:30]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd2,bnd3\n",
                ErrorCode.MULTIALLELIC_ALT,
            ),
        )

        for name, body, expected_code in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                tmp = Path(tmpdir)
                vcf = tmp / f"{name}.vcf"
                call_tsv = tmp / f"{name}.call.tsv"
                report = tmp / f"{name}.errors.txt"
                vcf.write_text(VCF_HEADER + body, encoding="utf-8")
                converter = KVarTSVConverter()

                with self.assertRaises(RuntimeError):
                    converter.convert_vcf_to_tsv(
                        str(vcf),
                        str(call_tsv),
                        str(report),
                    )

                self.assertFalse(call_tsv.exists())
                self.assertEqual(
                    1,
                    len(converter.error_handler.get_errors_by_code(expected_code)),
                )
                if expected_code == ErrorCode.MULTIALLELIC_ALT:
                    self.assertFalse(
                        converter.error_handler.get_errors_by_code(
                            ErrorCode.MULTIPLE_MATEIDS_UNSUPPORTED
                        )
                    )
                self.assertNotIn(
                    "BND_INSERTION_SPLIT",
                    report.read_text(encoding="utf-8"),
                )

    def test_invalid_breakend_flanking_sequence_blocks_without_split(self) -> None:
        cases = (
            ("invalid_base", "AX]chr1:20]"),
            ("wrong_ref_anchor", "TT]chr1:20]"),
            ("whitespace_in_target", "AAGT]chr 1:20]"),
            ("forbidden_target_character", "AAGT]chr<bad>:20]"),
            ("mixed_brackets", "AAGT]chr1:20["),
        )

        for name, alt in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                tmp = Path(tmpdir)
                vcf = tmp / f"{name}.vcf"
                call_tsv = tmp / f"{name}.call.tsv"
                report = tmp / f"{name}.errors.txt"
                vcf.write_text(
                    VCF_HEADER
                    + f"chr1\t10\tbnd1\tA\t{alt}\t.\tPASS\t"
                    "SVTYPE=BND;MATEID=bnd2\n",
                    encoding="utf-8",
                )
                converter = KVarTSVConverter()

                with self.assertRaises(RuntimeError):
                    converter.convert_vcf_to_tsv(
                        str(vcf),
                        str(call_tsv),
                        str(report),
                    )

                self.assertFalse(call_tsv.exists())
                self.assertTrue(
                    converter.error_handler.get_errors_by_code(
                        ErrorCode.INVALID_REF_ALT
                    )
                )
                self.assertTrue(
                    converter.error_handler.get_errors_by_code(
                        ErrorCode.INVALID_BREAKEND_FORMAT
                    )
                )
                self.assertNotIn(
                    "BND_INSERTION_SPLIT",
                    report.read_text(encoding="utf-8"),
                )

    def test_reciprocal_pair_output_is_independent_of_record_order(self) -> None:
        bodies = (
            (
                "chr1\t10\tbnd1\tA\tA]chr2:20]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd2;EVENT=event1\n"
                "chr2\t20\tbnd2\tT\tTACT]chr1:10]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd1;EVENT=event1\n"
            ),
            (
                "chr2\t20\tbnd2\tT\tTACT]chr1:10]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd1;EVENT=event1\n"
                "chr1\t10\tbnd1\tA\tA]chr2:20]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd2;EVENT=event1\n"
            ),
        )

        outputs = []
        for index, body in enumerate(bodies):
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp = Path(tmpdir)
                vcf = tmp / f"order_{index}.vcf"
                call_tsv = tmp / f"order_{index}.call.tsv"
                report = tmp / f"order_{index}.errors.txt"
                vcf.write_text(VCF_HEADER + body, encoding="utf-8")
                KVarTSVConverter().convert_vcf_to_tsv(
                    str(vcf),
                    str(call_tsv),
                    str(report),
                )
                outputs.append(read_call_rows(call_tsv))

        self.assertEqual(outputs[0], outputs[1])
        insertion = next(
            row
            for row in outputs[0]
            if row["Variant_Call_Type"] == "insertion"
        )
        self.assertEqual("bnd1_ins", insertion["Variant_Call_ID"])
        self.assertEqual("chr1", insertion["Chr"])
        self.assertEqual("10", insertion["Start"])
        self.assertEqual("AGT", insertion["Sequence"])
        self.assertEqual("bnd2", insertion["BND_Source_VCF_IDs"])

    def test_pair_without_insertion_has_stable_mutation_order(self) -> None:
        bodies = (
            (
                "chr1\t10\tbnd1\tA\tA]chr2:20]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd2;EVENT=event1\n"
                "chr2\t20\tbnd2\tT\tT]chr1:10]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd1;EVENT=event1\n"
            ),
            (
                "chr2\t20\tbnd2\tT\tT]chr1:10]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd1;EVENT=event1\n"
                "chr1\t10\tbnd1\tA\tA]chr2:20]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd2;EVENT=event1\n"
            ),
        )

        outputs = []
        for index, body in enumerate(bodies):
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp = Path(tmpdir)
                vcf = tmp / f"empty_order_{index}.vcf"
                call_tsv = tmp / f"empty_order_{index}.call.tsv"
                report = tmp / f"empty_order_{index}.errors.txt"
                vcf.write_text(VCF_HEADER + body, encoding="utf-8")
                KVarTSVConverter().convert_vcf_to_tsv(
                    str(vcf),
                    str(call_tsv),
                    str(report),
                )
                outputs.append(read_call_rows(call_tsv))

        self.assertEqual(outputs[0], outputs[1])
        self.assertEqual(1, len(outputs[0]))
        self.assertEqual("bnd1", outputs[0][0]["Variant_Call_ID"])
        self.assertEqual("event1", outputs[0][0]["Mutation_ID"])
        self.assertEqual("1", outputs[0][0]["Mutation_Order"])

    def test_colon_in_breakend_contig_name_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vcf = tmp / "colon_contig.vcf"
            call_tsv = tmp / "colon_contig.call.tsv"
            report = tmp / "colon_contig.errors.txt"
            vcf.write_text(
                VCF_HEADER.replace("VCFv4.2", "VCFv4.3", 1)
                + "chr1\t10\tbnd1\tA\tA]chr:part:20]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd2;EVENT=event1\n"
                + "chr:part\t20\tbnd2\tT\tT]chr1:10]\t.\tPASS\t"
                "SVTYPE=BND;MATEID=bnd1;EVENT=event1\n",
                encoding="utf-8",
            )

            KVarTSVConverter().convert_vcf_to_tsv(
                str(vcf),
                str(call_tsv),
                str(report),
            )

            rows = read_call_rows(call_tsv)
            self.assertEqual(1, len(rows))
            self.assertEqual("chr1", rows[0]["From_Chr"])
            self.assertEqual("chr:part", rows[0]["To_Chr"])


if __name__ == "__main__":
    unittest.main()
