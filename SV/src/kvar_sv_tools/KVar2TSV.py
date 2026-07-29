#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KVar SV VCF to TSV converter.
Converts KVar SV VCF files to submission Call TSV format.
"""

import sys
import os
import argparse
import gzip
import re
import tempfile
from typing import Dict, List, Optional, Any, Tuple
from collections import defaultdict

CALL_TSV_HEADER = [
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


# Relative path import support
try:
    from .VCF_parser import KVarVCFParser, VCFDataRow, SVClassifier, BreakendParser
    from .error_handler import ErrorHandler, ErrorCode
except ImportError:
    from VCF_parser import KVarVCFParser, VCFDataRow, SVClassifier, BreakendParser
    from error_handler import ErrorHandler, ErrorCode


class _OutputTransaction:
    """Stage related outputs and publish them only after validation succeeds."""

    def __init__(self) -> None:
        self._staged_paths: Dict[str, str] = {}

    def stage(self, output_path: Optional[str]) -> Optional[str]:
        if not output_path:
            return None

        final_path = os.path.abspath(output_path)
        if final_path in self._staged_paths:
            return self._staged_paths[final_path]

        output_dir = os.path.dirname(final_path) or "."
        suffix = ".tmp.gz" if final_path.endswith(".gz") else ".tmp"
        fd, staged_path = tempfile.mkstemp(
            dir=output_dir,
            prefix=f".{os.path.basename(final_path)}.",
            suffix=suffix,
        )
        os.close(fd)
        self._staged_paths[final_path] = staged_path
        return staged_path

    def publish(self) -> None:
        for final_path, staged_path in list(self._staged_paths.items()):
            os.replace(staged_path, final_path)
            del self._staged_paths[final_path]

    def cleanup(self) -> None:
        for staged_path in self._staged_paths.values():
            if os.path.exists(staged_path):
                os.unlink(staged_path)
        self._staged_paths.clear()

    def __enter__(self) -> "_OutputTransaction":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.cleanup()


IUPAC_SEQUENCE_RE = re.compile(r"^[- .ABCDGHKMNRSTUVWY]+$", re.IGNORECASE)
ARCHIVE_XREF_DBS = ("AE", "dbGaP", "dbSNP", "dbSNP-batch", "DDBJ", "DGV", "EGA", "ENA", "GENBANK", "GENE", "GEO", "SRA", "TRACE", "GEA", "JGA")
PHENOTYPE_XREF_DBS = ("HP", "MedGen", "MeSH", "OMIM", "SNOMED", "UMLS")
OTHER_XREF_DBS = ("CORIELL", "BioProj", "BioSD", "PubMed", "GeneReviews")
ALL_XREF_DBS = ARCHIVE_XREF_DBS + PHENOTYPE_XREF_DBS + OTHER_XREF_DBS
URL_RE = re.compile(r"https?://\S+", re.IGNORECASE)


def _xref_pattern(databases: Tuple[str, ...]) -> re.Pattern:
    escaped = "|".join(re.escape(database) for database in databases)
    return re.compile(rf"(?:^|[,\s])(?:{escaped})\s*:\s*[-A-Za-z0-9]+(?:$|[,\s])", re.IGNORECASE)


PHENOTYPE_XREF_RE = _xref_pattern(PHENOTYPE_XREF_DBS)
ALL_XREF_RE = _xref_pattern(ALL_XREF_DBS)


def _clean_tsv_value(value: Optional[str]) -> Optional[str]:
    text = str(value or "").strip()
    if not text or text == ".":
        return None
    return text


def _validate_sequence_and_reference_fields(
    row: Dict[str, str],
    error_handler: ErrorHandler,
    *,
    submitted_call_id: str,
    line_number: Optional[int] = None,
    line_content: Optional[str] = None,
) -> None:
    sequence = _clean_tsv_value(row.get("Sequence"))
    if sequence is not None and not IUPAC_SEQUENCE_RE.fullmatch(sequence):
        error_handler.create_error(
            ErrorCode.INVALID_SEQUENCE_FIELD,
            line_number=line_number,
            variant_id=submitted_call_id,
            line_content=line_content,
            field_name="Sequence",
            expected_value="IUPAC sequence symbols ABCDGHKMNRSTUVWY plus space, dot, or dash",
            actual_value=sequence,
        )

    phenotype = _clean_tsv_value(row.get("Phenotype"))
    if phenotype is not None and not PHENOTYPE_XREF_RE.search(phenotype):
        error_handler.create_error(
            ErrorCode.INVALID_PHENOTYPE_LINK,
            line_number=line_number,
            variant_id=submitted_call_id,
            line_content=line_content,
            field_name="Phenotype",
            expected_value="phenotype db:id such as HP:0000001, OMIM:123456, or MedGen:C000000",
            actual_value=phenotype,
        )

    evidence = _clean_tsv_value(row.get("Evidence"))
    if evidence is not None and not ALL_XREF_RE.search(evidence):
        error_handler.create_error(
            ErrorCode.INVALID_EVIDENCE_LINK,
            line_number=line_number,
            variant_id=submitted_call_id,
            line_content=line_content,
            field_name="Evidence",
            expected_value="db:id such as SRA:SRR000000 or PubMed:123456",
            actual_value=evidence,
        )

    external_links = _clean_tsv_value(row.get("External Links"))
    if external_links is not None and not (ALL_XREF_RE.search(external_links) or URL_RE.search(external_links)):
        error_handler.create_error(
            ErrorCode.INVALID_EXTERNAL_LINK,
            line_number=line_number,
            variant_id=submitted_call_id,
            line_content=line_content,
            field_name="External Links",
            expected_value="db:id such as GEO:GPL4010 or an http(s) URL",
            actual_value=external_links,
        )


class KVarTSVConverter:
    """Convert KVar SV VCF to submission Call TSV format"""

    def __init__(
        self,
        error_handler: Optional[ErrorHandler] = None,
        reference_fasta_path: Optional[str] = None
    ):
        self.error_handler = error_handler or ErrorHandler()
        self.reference_fasta_path = reference_fasta_path
        self.parser = KVarVCFParser(
            self.error_handler,
            reference_fasta_path=reference_fasta_path,
            strict_kvar_tags=True
        )
        self.sv_classifier = SVClassifier()
        self.breakend_parser = BreakendParser()

    def convert_vcf_to_tsv(
        self,
        vcf_file_path: str,
        output_file_path: str,
        error_report_path: Optional[str] = None,
    ) -> None:
        """Validate a submitted VCF and create an accession-free Call TSV."""
        try:
            # Parse VCF file
            self.parser.parse_file(vcf_file_path)
        except Exception as e:
            # Stop conversion on critical errors
            if self.error_handler.has_critical_errors():
                print("Conversion stopped due to critical errors.")
                if error_report_path:
                    self.error_handler.generate_report(
                        error_report_path,
                        vcf_file_path,
                        output_file_path,
                    )
                raise
        finally:
            if self.parser.reference is not None:
                self.parser.reference.close()

        # Build ID map and BND grouping
        id_map = self._build_id_map()
        mutation_id_map = self._build_mutation_id_map(id_map)
        call_records = self._build_call_records(id_map)

        # Block final output when CRITICAL/ERROR messages were collected.
        self.error_handler.assert_no_blocking_errors(
            stage="SV VCF validation",
            output_file=error_report_path,
            vcf_file_path=vcf_file_path,
            output_tsv_path=output_file_path,
        )

        report_has_errors = self.error_handler.has_errors()
        try:
            with _OutputTransaction() as outputs:
                # Stage the required report first so a report-generation failure
                # cannot leave a newly published Call TSV behind.
                staged_error_report = outputs.stage(error_report_path)
                staged_call_tsv = outputs.stage(output_file_path) or output_file_path

                self._write_tsv_file(staged_call_tsv, mutation_id_map, call_records)
                if staged_error_report:
                    self.error_handler.generate_report(
                        staged_error_report,
                        report_display_path=error_report_path,
                        vcf_file_path=vcf_file_path,
                        output_tsv_path=output_file_path,
                    )
                outputs.publish()
        except RuntimeError:
            raise
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FILE_WRITE_ERROR,
                additional_info={"file_path": output_file_path, "error": str(e)}
            )
            if error_report_path:
                self.error_handler.generate_report(
                    error_report_path,
                    vcf_file_path,
                    output_tsv_path=output_file_path,
                )
            raise

        if error_report_path:
            report_suffix = "" if report_has_errors else " (no issues)"
            print(f"Validation report: {error_report_path}{report_suffix}")

        print(f"Conversion complete: {vcf_file_path} -> {output_file_path}")
        print(f"Total {len(call_records)} variant calls written")

        # Print error summary
        self.error_handler.print_summary()

    def _build_id_map(self) -> Dict[str, int]:
        """Build mapping from variant ID to index"""
        id_map = {}
        for idx, row in enumerate(self.parser.data_rows):
            if row.id and row.id != ".":
                id_map[row.id] = idx
        return id_map

    def _split_info_values(self, value: Any) -> List[str]:
        """Return comma/list INFO values as clean strings."""
        if value is None or value == ".":
            return []
        if isinstance(value, list):
            values = value
        else:
            values = str(value).split(",")
        return [str(item).strip() for item in values if str(item).strip() and str(item).strip() != "."]

    def _canonical_chrom(self, chrom: str) -> str:
        """Return the reference-resolved chromosome name when possible."""
        chrom = str(chrom or "").strip()
        if not chrom or chrom == ".":
            return chrom
        if self.parser.reference:
            resolved = self.parser.reference.resolve_chrom(chrom)
            if resolved:
                return resolved
        return chrom

    def _bnd_endpoints(self, row: VCFDataRow) -> Optional[Tuple[Tuple[str, str, str], Tuple[str, str, str]]]:
        from_strand, to_chr, to_coord, to_strand = self.breakend_parser.parse_breakend_placement(row.alt)
        if "." in {from_strand, to_chr, to_coord, to_strand}:
            return None
        return (
            (self._canonical_chrom(row.chrom), str(row.pos), from_strand),
            (self._canonical_chrom(to_chr), str(to_coord), to_strand),
        )

    def _position_interval(self, row: VCFDataRow) -> Tuple[int, int]:
        """Return the accepted breakpoint interval for a BND row POS."""
        cipos = row.info.get("CIPOS")
        if isinstance(cipos, list) and len(cipos) >= 2:
            try:
                left = int(cipos[0])
                right = int(cipos[1])
                return row.pos + left, row.pos + right
            except (TypeError, ValueError):
                pass
        return row.pos, row.pos

    def _bnd_target_matches_row_position(self, target_coord: str, row: VCFDataRow) -> bool:
        try:
            target = int(target_coord)
        except (TypeError, ValueError):
            return False
        start, stop = self._position_interval(row)
        return start <= target <= stop

    @staticmethod
    def _opposite_strand(strand: str) -> str:
        return {"+": "-", "-": "+"}.get(strand, ".")

    def _bnd_mate_coordinates_are_compatible(self, row: VCFDataRow, mate_row: VCFDataRow) -> bool:
        row_endpoints = self._bnd_endpoints(row)
        mate_endpoints = self._bnd_endpoints(mate_row)
        if not row_endpoints or not mate_endpoints:
            return False
        row_from, row_to = row_endpoints
        mate_from, mate_to = mate_endpoints
        return (
            row_to[0] == mate_from[0]
            and mate_to[0] == row_from[0]
            and self._bnd_target_matches_row_position(row_to[1], mate_row)
            and self._bnd_target_matches_row_position(mate_to[1], row)
        )

    def _bnd_mate_strands_are_compatible(self, row: VCFDataRow, mate_row: VCFDataRow) -> bool:
        row_endpoints = self._bnd_endpoints(row)
        mate_endpoints = self._bnd_endpoints(mate_row)
        if not row_endpoints or not mate_endpoints:
            return False
        row_from, row_to = row_endpoints
        mate_from, mate_to = mate_endpoints
        return (
            mate_from[2] == self._opposite_strand(row_to[2])
            and mate_to[2] == self._opposite_strand(row_from[2])
        )

    def _bnd_mates_are_compatible(self, row: VCFDataRow, mate_row: VCFDataRow) -> bool:
        return (
            self._bnd_mate_coordinates_are_compatible(row, mate_row)
            and self._bnd_mate_strands_are_compatible(row, mate_row)
        )

    def _record_mate_alt_mismatch(self, row: VCFDataRow, mate_row: VCFDataRow) -> None:
        self.error_handler.create_error(
            ErrorCode.MATEID_ALT_MISMATCH,
            variant_id=row.id if row.id and row.id != "." else None,
            field_name="ALT/MATEID",
            expected_value="reciprocal breakend coordinates within each mate's CIPOS interval",
            actual_value=f"{row.id}:{row.alt} / {mate_row.id}:{mate_row.alt}",
        )

    def _record_mate_strand_mismatch(self, row: VCFDataRow, mate_row: VCFDataRow) -> None:
        row_from, row_to = self._bnd_endpoints(row) or ((".", ".", "."), (".", ".", "."))
        mate_from, mate_to = self._bnd_endpoints(mate_row) or ((".", ".", "."), (".", ".", "."))
        self.error_handler.create_error(
            ErrorCode.MATEID_STRAND_MISMATCH,
            variant_id=row.id if row.id and row.id != "." else None,
            field_name="ALT/MATEID strand orientation",
            expected_value=(
                f"mate strands {self._opposite_strand(row_to[2])},"
                f"{self._opposite_strand(row_from[2])}"
            ),
            actual_value=f"{mate_from[2]},{mate_to[2]} ({row.id}:{row.alt} / {mate_row.id}:{mate_row.alt})",
        )

    def _bnd_inserted_sequence(self, row: VCFDataRow) -> str:
        if row.info.get("SVTYPE", "") != "BND":
            return ""
        return self.breakend_parser.inserted_sequence(row.alt, row.ref).upper()

    @staticmethod
    def _reverse_complement(sequence: str) -> str:
        return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

    def _mate_inserted_sequence_in_row_orientation(
        self,
        row: VCFDataRow,
        mate_row: VCFDataRow,
        mate_sequence: str,
    ) -> str:
        if not mate_sequence:
            return ""

        row_from_strand = self.breakend_parser.parse_breakend_placement(row.alt)[0]
        mate_from_strand = self.breakend_parser.parse_breakend_placement(mate_row.alt)[0]
        if row_from_strand == mate_from_strand:
            return self._reverse_complement(mate_sequence)
        return mate_sequence

    @staticmethod
    def _bnd_pair_sort_key(row: VCFDataRow) -> Tuple[str, int, str]:
        """Return the stable key used to select a reciprocal pair's primary row."""
        row_id = row.id if row.id and row.id != "." else ""
        return str(row.chrom), int(row.pos), row_id

    def _record_one_sided_mate_insertion(
        self,
        row: VCFDataRow,
        mate_row: VCFDataRow,
        row_sequence: str,
        mate_sequence: str,
    ) -> None:
        self.error_handler.create_error(
            ErrorCode.MATEID_INSERTION_SEQUENCE_ONE_SIDED,
            variant_id=row.id if row.id and row.id != "." else None,
            field_name="ALT/MATEID inserted sequence",
            expected_value="equivalent inserted sequence in both reciprocal breakend ALTs",
            actual_value=(
                f"{row.id}={'present' if row_sequence else 'missing'}, "
                f"{mate_row.id}={'present' if mate_sequence else 'missing'}"
            ),
        )

    def _record_mate_insertion_mismatch(
        self,
        row: VCFDataRow,
        mate_row: VCFDataRow,
        row_sequence: str,
        mate_sequence: str,
    ) -> None:
        self.error_handler.create_error(
            ErrorCode.MATEID_INSERTION_SEQUENCE_MISMATCH,
            variant_id=row.id if row.id and row.id != "." else None,
            field_name="ALT/MATEID inserted sequence",
            expected_value=f"{row.id} inserted sequence {row_sequence}",
            actual_value=(
                f"{mate_row.id} inserted sequence in {row.id} orientation "
                f"{mate_sequence}"
            ),
        )

    def _bnd_split_mutation_id(self, row: VCFDataRow, inserted_sequence: str) -> Optional[str]:
        if not inserted_sequence:
            return None
        row_id = row.id if row.id and row.id != "." else f"{row.chrom}_{row.pos}"
        if self.breakend_parser.is_single_breakend_alt(row.alt, row.ref):
            return f"{row_id}_sbnd_ins"
        return f"{row_id}_trans_ins"

    def _unique_generated_id(self, base_id: str, used_ids: set) -> str:
        candidate = base_id
        counter = 1
        while candidate in used_ids:
            candidate = f"{base_id}_{counter}"
            counter += 1
        used_ids.add(candidate)
        return candidate

    def _synthetic_insertion_row(self, source_row: VCFDataRow, synthetic_id: str, sequence: str) -> VCFDataRow:
        info = dict(source_row.info)
        info["SVTYPE"] = "INS"
        info["END"] = source_row.pos + len(sequence)
        info["SVLEN"] = len(sequence)
        info["SVINSSEQ"] = sequence
        info.pop("MATEID", None)
        synthetic_row = VCFDataRow(
            chrom=source_row.chrom,
            pos=source_row.pos,
            id=synthetic_id,
            ref=source_row.ref,
            alt=sequence,
            qual=source_row.qual,
            filter=source_row.filter,
            info=info,
        )
        self.parser.validate_effective_end_bounds(
            synthetic_row.chrom,
            synthetic_row.pos,
            synthetic_row.ref,
            synthetic_row.alt,
            synthetic_row.info,
            None,
            synthetic_id,
            coordinate_origin="converter-derived BND insertion END",
        )
        return synthetic_row

    def _build_call_records(self, id_map: Dict[str, int]) -> List[Dict[str, Any]]:
        """Build output call records, collapsing reciprocal BND MATEID pairs."""
        call_records = []
        absorbed_ids = set()
        generated_ids = set(id_map)
        checked_mate_pairs = set()

        for row in self.parser.data_rows:
            row_id = row.id if row.id and row.id != "." else "."
            if row_id in absorbed_ids:
                continue

            bnd_source_ids = []
            source_by_id = {row_id: "submitted"}
            reason_by_id = {row_id: "validated"}
            inserted_sequence = self._bnd_inserted_sequence(row)
            insertion_source_ids = [row_id] if inserted_sequence else []
            insertion_split_allowed = True
            defer_to_canonical_mate = False

            if row.info.get("SVTYPE", "") == "BND" and row_id != ".":
                mate_ids = self._split_info_values(row.info.get("MATEID"))
                if len(mate_ids) > 1 or (mate_ids and mate_ids[0] == row_id):
                    insertion_split_allowed = False
                    mate_ids = []
                for mate_id in mate_ids:
                    if mate_id not in id_map or mate_id in absorbed_ids:
                        continue

                    mate_row = self.parser.data_rows[id_map[mate_id]]
                    if mate_row.info.get("SVTYPE", "") != "BND":
                        continue

                    reciprocal_mates = self._split_info_values(mate_row.info.get("MATEID"))
                    if (
                        len(reciprocal_mates) > 1
                        or (
                            reciprocal_mates
                            and reciprocal_mates[0] == mate_id
                        )
                    ):
                        insertion_split_allowed = False
                        continue
                    if not reciprocal_mates or reciprocal_mates[0] != row_id:
                        continue
                    pair_key = tuple(sorted((row_id, mate_id)))
                    if not self._bnd_mate_coordinates_are_compatible(row, mate_row):
                        if pair_key not in checked_mate_pairs:
                            self._record_mate_alt_mismatch(row, mate_row)
                            checked_mate_pairs.add(pair_key)
                        continue
                    if not self._bnd_mate_strands_are_compatible(row, mate_row):
                        if pair_key not in checked_mate_pairs:
                            self._record_mate_strand_mismatch(row, mate_row)
                            checked_mate_pairs.add(pair_key)
                        continue

                    if self._bnd_pair_sort_key(mate_row) < self._bnd_pair_sort_key(row):
                        defer_to_canonical_mate = True
                        break

                    mate_sequence = self._bnd_inserted_sequence(mate_row)
                    normalized_mate_sequence = self._mate_inserted_sequence_in_row_orientation(
                        row,
                        mate_row,
                        mate_sequence,
                    )
                    if inserted_sequence and normalized_mate_sequence:
                        if inserted_sequence != normalized_mate_sequence:
                            self._record_mate_insertion_mismatch(
                                row,
                                mate_row,
                                inserted_sequence,
                                normalized_mate_sequence,
                            )
                            inserted_sequence = ""
                            insertion_source_ids = []
                        else:
                            insertion_source_ids = [row_id, mate_id]
                    elif inserted_sequence or normalized_mate_sequence:
                        self._record_one_sided_mate_insertion(
                            row,
                            mate_row,
                            inserted_sequence,
                            mate_sequence,
                        )
                        if normalized_mate_sequence:
                            inserted_sequence = normalized_mate_sequence
                            insertion_source_ids = [mate_id]
                        else:
                            insertion_source_ids = [row_id]

                    bnd_source_ids = [row_id, mate_id]
                    source_by_id = {
                        row_id: "submitted_primary",
                        mate_id: "collapsed_mate"
                    }
                    reason_by_id = {
                        row_id: "collapsed_bnd_pair",
                        mate_id: f"absorbed_by={row_id}"
                    }
                    absorbed_ids.add(mate_id)
                    break

            if defer_to_canonical_mate:
                continue

            split_mutation_id = (
                self._bnd_split_mutation_id(row, inserted_sequence)
                if insertion_split_allowed
                else None
            )
            mutation_info_override = None
            if split_mutation_id:
                mutation_info_override = {
                    "mutation_id": split_mutation_id,
                    "mutation_order": "1",
                    "mutation_molecule": ".",
                }

            call_record = {
                "row": row,
                "primary_id": row_id,
                "output_id": row_id,
                "bnd_source_ids": bnd_source_ids,
                "source_by_id": source_by_id,
                "reason_by_id": reason_by_id,
            }
            if mutation_info_override:
                call_record["mutation_info_override"] = mutation_info_override
            call_records.append(call_record)

            if insertion_split_allowed and inserted_sequence:
                synthetic_id = self._unique_generated_id(f"{row_id}_ins", generated_ids)
                synthetic_row = self._synthetic_insertion_row(row, synthetic_id, inserted_sequence)
                call_records.append({
                    "row": synthetic_row,
                    "primary_id": synthetic_id,
                    "output_id": synthetic_id,
                    "bnd_source_ids": insertion_source_ids,
                    "source_by_id": {
                        source_id: "derived_from_submitted"
                        for source_id in insertion_source_ids
                    },
                    "reason_by_id": {
                        source_id: "split_bnd_insertion"
                        for source_id in insertion_source_ids
                    },
                    "call_type_override": "insertion",
                    "mutation_info_override": {
                        "mutation_id": split_mutation_id,
                        "mutation_order": "2",
                        "mutation_molecule": ".",
                    }
                })
                self.error_handler.create_error(
                    ErrorCode.BND_INSERTION_SPLIT,
                    variant_id=row_id if row_id != "." else None,
                    field_name="ALT",
                    actual_value=inserted_sequence,
                    additional_info={"synthetic_call_id": synthetic_id},
                )

        return call_records

    def _build_mutation_id_map(self, id_map: Dict[str, int]) -> Dict[str, Dict[str, str]]:
        """Group BNDs and assign Mutation ID"""
        mutation_id_map = {}
        processed_variants = set()

        for idx, row in enumerate(self.parser.data_rows):
            if idx in processed_variants:
                continue

            svtype = row.info.get("SVTYPE", "")

            if svtype == "BND":
                # Find BND group (DFS)
                group = []
                self._find_bnd_group(idx, group, processed_variants, id_map)

                if group:
                    ordered_group = sorted(
                        group,
                        key=lambda variant_idx: self._bnd_pair_sort_key(
                            self.parser.data_rows[variant_idx]
                        ),
                    )
                    # DDBJ uses EVENT as the submitter/dbVar-style Mutation ID anchor.
                    # Without EVENT, leave Mutation ID blank; related BND calls remain ungrouped.
                    mutation_id = "."
                    for variant_idx in ordered_group:
                        event_ids = self._split_info_values(self.parser.data_rows[variant_idx].info.get("EVENT"))
                        if event_ids:
                            mutation_id = event_ids[0]
                            break
                    if mutation_id == ".":
                        continue

                    for i, variant_idx in enumerate(ordered_group):
                        variant_id = self.parser.data_rows[variant_idx].id
                        mutation_id_map[variant_id] = {
                            "mutation_id": mutation_id,
                            "mutation_order": str(i + 1),
                            "mutation_molecule": ".",
                        }

        return mutation_id_map

    def _find_bnd_group(self, start_idx: int, group: List[int], processed_variants: set, id_map: Dict[str, int]) -> None:
        """Find BND group via DFS"""
        if start_idx in processed_variants:
            return

        row = self.parser.data_rows[start_idx]
        svtype = row.info.get("SVTYPE", "")

        if svtype != "BND":
            return

        group.append(start_idx)
        processed_variants.add(start_idx)

        # Find mate
        mate_str = row.info.get("MATEID", ".")
        if mate_str != ".":
            mates = [mate.strip() for mate in mate_str.split(",")]

            for mate in mates:
                if mate in id_map:
                    mate_idx = id_map[mate]
                    mate_row = self.parser.data_rows[mate_idx]
                    reciprocal_mates = self._split_info_values(mate_row.info.get("MATEID"))
                    if row.id in reciprocal_mates and self._bnd_mates_are_compatible(row, mate_row):
                        self._find_bnd_group(mate_idx, group, processed_variants, id_map)

    def _write_tsv_file(
        self,
        output_file_path: str,
        mutation_id_map: Dict[str, Dict[str, str]],
        call_records: Optional[List[Dict[str, Any]]] = None
    ) -> None:
        """Write TSV file"""
        with open(output_file_path, 'w', encoding='utf-8') as f:
            # Write header
            f.write("##Variant_Call\n")
            f.write("#" + "\t".join(CALL_TSV_HEADER) + "\n")

            # Write data rows
            if call_records is None:
                call_records = [
                    {
                        "row": row,
                        "bnd_source_ids": [],
                    }
                    for row in self.parser.data_rows
                ]

            for call_record in call_records:
                tsv_row = self._convert_row_to_tsv(
                    call_record["row"],
                    mutation_id_map,
                    call_record.get("bnd_source_ids"),
                    call_record.get("output_id", call_record.get("primary_id")),
                    call_record.get("call_type_override"),
                    call_record.get("mutation_info_override"),
                )
                f.write(tsv_row + "\n")

    @staticmethod
    def _format_info_value(value: Any) -> str:
        if value is None:
            return "."
        if isinstance(value, list):
            values = [str(item).strip() for item in value if str(item).strip() and str(item).strip() != "."]
            return ",".join(values) if values else "."
        text = str(value).strip()
        return text if text and text != "." else "."

    def _resolved_end(self, row: VCFDataRow) -> int:
        """Return submitted END or derive it the same way VCF validation accepts it."""
        end = row.info.get("END")
        if isinstance(end, int):
            return end
        derived_end = KVarVCFParser.derive_missing_end(
            row.info.get("SVTYPE", ""),
            row.pos,
            row.ref,
            row.alt,
            row.info.get("SVLEN"),
        )
        return derived_end if derived_end is not None else row.pos

    def _convert_row_to_tsv(
        self,
        row: VCFDataRow,
        mutation_id_map: Dict[str, Dict[str, str]],
        bnd_source_vcf_ids: Optional[List[str]] = None,
        output_call_id: Optional[str] = None,
        call_type_override: Optional[str] = None,
        mutation_info_override: Optional[Dict[str, str]] = None,
    ) -> str:
        """Convert VCF data row to TSV format"""
        info = row.info
        svtype = info.get("SVTYPE", "")
        end = self._resolved_end(row)

        # Variant Call ID
        variant_call_id = output_call_id or (row.id if row.id and row.id != "." else ".")
        chrom = self._canonical_chrom(row.chrom)

        # Variant Call Type
        if call_type_override:
            call_type = call_type_override
        elif svtype == "BND":
            to_chr, to_coord, _ = self.breakend_parser.parse_breakend(row.alt)
            canonical_to_chr = self._canonical_chrom(to_chr)
            if self.breakend_parser.is_single_breakend_alt(row.alt, row.ref):
                call_type = "interchromosomal translocation"
            elif canonical_to_chr != "." and canonical_to_chr != chrom:
                call_type = "interchromosomal translocation"
            elif canonical_to_chr != "." and canonical_to_chr == chrom:
                call_type = "intrachromosomal translocation"
            else:
                call_type = ""
        else:
            call_type = self.sv_classifier.classify_calltype(info, row.alt, ref=row.ref)

        # Coordinate handling
        outer_start, inner_start, inner_stop, outer_stop = self._process_coordinates(info, row.pos, end)

        # Keep POS/END anchors in Start/Stop. CIPOS/CIEND are represented by
        # the optional confidence-bound columns when present.
        start_value = "." if "POSrange" in info else str(row.pos)
        end_value = "." if "ENDrange" in info else str(end)

        # Insertion Length
        insertion_length = self._calculate_insertion_length(svtype, info, row.pos, end)

        # Allele info
        allele_count = str(info.get("AC", "."))
        allele_frequency = str(info.get("AF", "."))
        allele_number = str(info.get("AN", "."))
        copy_number = self._format_info_value(info.get("CN"))

        # Description
        description = info.get("DESC", ".")
        if description != "." and isinstance(description, str) and description.startswith('"') and description.endswith('"'):
            description = description[1:-1]  # Strip quotes

        # Validation
        validation = self._process_validation(row, info)

        # Origin
        origin = info.get("ORIGIN", ".")
        if origin != "." and isinstance(origin, str) and origin.startswith('"') and origin.endswith('"'):
            origin = origin[1:-1]  # Strip quotes

        # Phenotype
        phenotype = info.get("PHENO", ".")
        if phenotype != "." and isinstance(phenotype, str) and phenotype.startswith('"') and phenotype.endswith('"'):
            phenotype = phenotype[1:-1]  # Strip quotes

        # External Links
        external_links = info.get("LINKS", ".")

        # Sequence
        sequence = self._get_sequence_field(svtype, row.ref, row.alt, info)

        _validate_sequence_and_reference_fields(
            {
                "Sequence": sequence,
                "Phenotype": phenotype,
                "External Links": external_links,
                "Evidence": ".",
            },
            self.error_handler,
            submitted_call_id=variant_call_id,
        )

        # Translocation info (BND)
        from_chr = from_coord = from_strand = to_chr = to_coord = to_strand = "."
        mutation_id = mutation_order = mutation_molecule = "."

        if svtype == "BND":
            from_strand, to_chr, to_coord, to_strand = self.breakend_parser.parse_breakend_placement(row.alt)
            to_chr = self._canonical_chrom(to_chr)
            if self.breakend_parser.is_single_breakend_alt(row.alt, row.ref):
                from_strand, _ = self.breakend_parser.parse_single_breakend(row.alt, row.ref)
                to_chr = to_coord = to_strand = "."
            from_chr = chrom
            from_coord = str(row.pos)

            # Get info from mutation_id_map
            if mutation_info_override:
                mutation_info = mutation_info_override
                mutation_id = mutation_info["mutation_id"]
                mutation_order = mutation_info["mutation_order"]
                mutation_molecule = mutation_info["mutation_molecule"]
            elif row.id in mutation_id_map:
                mutation_info = mutation_id_map[row.id]
                mutation_id = mutation_info["mutation_id"]
                mutation_order = mutation_info["mutation_order"]
                mutation_molecule = mutation_info["mutation_molecule"]
        elif mutation_info_override:
            mutation_id = mutation_info_override["mutation_id"]
            mutation_order = mutation_info_override["mutation_order"]
            mutation_molecule = mutation_info_override["mutation_molecule"]

        # TSV field assembly
        tsv_fields = [
            variant_call_id,      # Variant Call ID
            call_type,            # Variant Call Type
            chrom,                # Chr
            outer_start,          # Outer Start
            start_value,          # Start
            inner_start,          # Inner Start
            inner_stop,           # Inner Stop
            end_value,            # Stop
            outer_stop,           # Outer Stop
            insertion_length,     # Insertion Length
            allele_count,         # Allele Count
            allele_frequency,     # Allele Frequency
            allele_number,        # Allele Number
            copy_number,          # Copy Number
            description,          # Description
            validation,           # Validation
            ".",                  # Zygosity
            origin,               # Origin
            phenotype,            # Phenotype
            external_links,       # External Links
            ".",                  # Evidence
            sequence,             # Sequence
            from_chr,             # From Chr
            from_coord,           # From Coord
            from_strand,          # From Strand
            to_chr,               # To Chr
            to_coord,             # To Coord
            to_strand,            # To Strand
            mutation_id,          # Mutation ID
            mutation_order,       # Mutation Order
            mutation_molecule,    # Mutation Molecule
            ",".join(bnd_source_vcf_ids) if bnd_source_vcf_ids else ".",
            # BND_Source_VCF_IDs
        ]

        return "\t".join(tsv_fields)

    def _process_coordinates(self, info: Dict[str, Any], pos: int, end: int) -> Tuple[str, str, str, str]:
        """Coordinate handling (POSrange/ENDrange or CIPOS/CIEND)

        POSrange format: (outer_start, inner_start)
        - outers only: (2500000, .) -> outer_start=2500000, inner_start="."
        - inners only: (., 2501000) -> outer_start=".", inner_start=2501000
        - outers and inners: (2500000, 2501000) -> outer_start=2500000, inner_start=2501000
        - precise (no POSrange): outer_start=".", inner_start="." (Start uses POS)

        ENDrange format: (inner_stop, outer_stop)
        - outers only: (., 3500000) -> inner_stop=".", outer_stop=3500000
        - inners only: (3499000, .) -> inner_stop=3499000, outer_stop="."
        - outers and inners: (3499000, 3500000) -> inner_stop=3499000, outer_stop=3500000
        - precise (no ENDrange): inner_stop=".", outer_stop="." (Stop uses END)
        """
        # Prefer POSrange/ENDrange
        posrange = info.get("POSrange")
        endrange = info.get("ENDrange")

        # POSrange handling
        if posrange is not None:
            if isinstance(posrange, list) and len(posrange) >= 2:
                outer_val = posrange[0]
                inner_val = posrange[1]

                # outer_start: (., inner)->"."; (outer,.)->outer; (outer,inner)->outer
                outer_str = str(outer_val).strip()
                if outer_val == '.' or outer_str == '.':
                    outer_start = "."  # inners only
                else:
                    outer_start = str(outer_val)

                # inner_start: (outer,.)->"."; (.,inner)->inner; (outer,inner)->inner
                inner_str = str(inner_val).strip()
                if inner_val == '.' or inner_str == '.':
                    inner_start = "."  # outers only
                else:
                    inner_start = str(inner_val)
            else:
                outer_start = inner_start = "."
        else:
            # If no POSrange, use CIPOS (else precise)
            cipos = info.get("CIPOS")
            if cipos and isinstance(cipos, list) and len(cipos) >= 2:
                outer_start = self._ci_offset_to_coordinate(pos, cipos[0])
                inner_start = self._ci_offset_to_coordinate(pos, cipos[1])
            else:
                outer_start = inner_start = "."  # precise: Start uses POS

        # ENDrange handling
        if endrange is not None:
            if isinstance(endrange, list) and len(endrange) >= 2:
                inner_val = endrange[0]
                outer_val = endrange[1]

                # inner_stop: (.,outer)->"."; (inner,.)->inner; (inner,outer)->inner
                inner_str = str(inner_val).strip()
                if inner_val == '.' or inner_str == '.':
                    inner_stop = "."  # outers only
                else:
                    inner_stop = str(inner_val)
                # outer_stop: (inner,.)->"."; (.,outer)->outer; (inner,outer)->outer
                outer_str = str(outer_val).strip()
                if outer_val == '.' or outer_str == '.':
                    outer_stop = "."  # inners only
                else:
                    outer_stop = str(outer_val)
            else:
                inner_stop = outer_stop = "."
        else:
            # If CIEND is missing, VCF v4.4 says it is assumed to match CIPOS.
            ciend = info.get("CIEND") or info.get("CIPOS")
            if ciend and isinstance(ciend, list) and len(ciend) >= 2:
                inner_stop = self._ci_offset_to_coordinate(end, ciend[0])
                outer_stop = self._ci_offset_to_coordinate(end, ciend[1])
            else:
                inner_stop = outer_stop = "."  # precise: Stop uses END

        return outer_start, inner_start, inner_stop, outer_stop

    @staticmethod
    def _ci_offset_to_coordinate(anchor: int, offset: Any) -> str:
        """Return an absolute TSV coordinate for a CI offset, or "." if unknown."""
        if offset == "." or str(offset).strip() == ".":
            return "."
        try:
            parsed_offset = int(offset)
        except (TypeError, ValueError):
            return "."
        return str(anchor + parsed_offset)

    def _process_validation(self, row: VCFDataRow, info: Dict[str, Any]) -> str:
        """Process Validation field.

        1. No valEXPERIMENT -> "."
        2. If present: only "ExperimentID:Pass" or "ExperimentID:Fail"; invalid -> report and ".";
           valid -> join multiple with comma (strip spaces).
        """
        if "valEXPERIMENT" not in info:
            return "."

        val_exp = info["valEXPERIMENT"]
        val_exp_str = str(val_exp)
        val_parts = [part.strip() for part in val_exp_str.split(",")]
        valid_format = True
        invalid_parts = []
        pattern = re.compile(r'^[^:]+:(Pass|Fail)$', re.IGNORECASE)

        for part in val_parts:
            if not pattern.match(part):
                valid_format = False
                invalid_parts.append(part)

        if not valid_format:
            variant_id = row.id if row.id and row.id != '.' else None
            self.error_handler.create_error(
                ErrorCode.INVALID_valEXPERIMENT_FORMAT,
                line_number=None,
                variant_id=variant_id,
                field_name="valEXPERIMENT",
                expected_value="ExperimentID:Pass or ExperimentID:Fail (comma-separated for multiple)",
                actual_value=val_exp_str,
                additional_info={
                    "invalid_parts": invalid_parts,
                    "message": f"Invalid format: {', '.join(invalid_parts)}"
                }
            )
            return "."

        cleaned_parts = [part.strip() for part in val_exp_str.split(",")]
        return ",".join(cleaned_parts)

    def _calculate_insertion_length(self, svtype: str, info: Dict[str, Any], pos: int, end: int) -> str:
        """Compute Insertion Length."""
        if svtype in ["INS", "DEL"]:
            svlen = self._first_integer_info_value(info.get("SVLEN"))
            if svlen is not None:
                return str(abs(svlen))
            if svtype == "INS":
                svinslen = self._first_integer_info_value(info.get("SVINSLEN"))
                if svinslen is not None:
                    return str(abs(svinslen))
            return "."

        return "."

    @staticmethod
    def _first_integer_info_value(value: Any) -> Optional[int]:
        """Return the first integer INFO value, or None for missing/non-numeric values."""
        if value is None:
            return None
        values = value if isinstance(value, list) else [value]
        for item in values:
            if item is None:
                continue
            if isinstance(item, int):
                return item
            text = str(item).strip()
            if not text or text == ".":
                continue
            try:
                return int(text)
            except ValueError:
                continue
        return None

    def _get_sequence_field(self, svtype: str, ref: str, alt: str, info: Dict[str, Any]) -> str:
        """Build Sequence field."""
        for key in ("SEQ", "SVINSSEQ"):
            sequence = self._clean_sequence_value(info.get(key))
            if sequence:
                return sequence

        if svtype == "BND":
            return "."

        # Check if actual sequence
        if alt and alt not in [".", "<INS>", "<DEL>", "<DUP>", "<INV>", "<CNV>", "<DUP:TANDEM>"]:
            if not alt.startswith("<") and not alt.endswith(">"):
                if ref and alt[0].upper() == ref[0].upper() and len(alt) > 1:
                    return alt[1:]
                return alt

        return "."

    @staticmethod
    def _clean_sequence_value(value: Any) -> str:
        """Return a non-empty INFO sequence value, preserving submitter casing."""
        if value is None:
            return ""
        if isinstance(value, list):
            value = ",".join(str(item) for item in value if str(item).strip() and str(item).strip() != ".")
        value = str(value).strip()
        if not value or value == ".":
            return ""
        return value
