#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KVar SV VCF to TSV converter.
Converts KVar SV VCF files to Variant_Call.tsv format.
"""

import sys
import os
import argparse
import gzip
import re
import tempfile
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple
from collections import defaultdict

CALL_SUBMITTED_IDS_FIELD = "Submitted_Variant_Call_IDs"
HGVSG_FIELD = "HGVSG"

CALL_INTERNAL_TSV_HEADER = [
    "Variant Call ID",
    "Variant Call Type",
    "Chr",
    "Outer Start",
    "Start",
    "Inner Start",
    "Inner Stop",
    "Stop",
    "Outer Stop",
    "Insertion Length",
    "Allele Count",
    "Allele Frequency",
    "Allele Number",
    "Copy Number",
    "Description",
    "Validation",
    "Zygosity",
    "Origin",
    "Phenotype",
    HGVSG_FIELD,
    "External Links",
    "Evidence",
    "Sequence",
    "From Chr",
    "From Coord",
    "From Strand",
    "To Chr",
    "To Coord",
    "To Strand",
    "Mutation ID",
    "Mutation Order",
    "Mutation Molecule",
    CALL_SUBMITTED_IDS_FIELD,
]
CALL_TSV_HEADER = [field.replace(" ", "_") for field in CALL_INTERNAL_TSV_HEADER]
CALL_FIELD_ALIASES = {
    output_field: internal_field
    for internal_field, output_field in zip(CALL_INTERNAL_TSV_HEADER, CALL_TSV_HEADER)
}
CALL_FIELD_ALIASES.update({field: field for field in CALL_INTERNAL_TSV_HEADER})
CALL_FIELD_ALIASES.update({
    "HGVSg": HGVSG_FIELD,
    "hgvs_name": HGVSG_FIELD,
})


def normalize_call_field_name(field: str) -> str:
    """Return the internal call TSV field name for public or legacy headers."""
    clean_field = str(field or "").strip()
    return CALL_FIELD_ALIASES.get(clean_field, clean_field)


# Relative path import support
try:
    from .VCF_parser import KVarVCFParser, VCFDataRow, SVClassifier, BreakendParser
    from .error_handler import ErrorHandler, ErrorCode
    from .metadata_validator import MetadataValidator
    from .metadata_parser import MetadataParser
    from .sv_type_ontology import normalize_call_type
except ImportError:
    from VCF_parser import KVarVCFParser, VCFDataRow, SVClassifier, BreakendParser
    from error_handler import ErrorHandler, ErrorCode
    from metadata_validator import MetadataValidator
    from metadata_parser import MetadataParser
    from sv_type_ontology import normalize_call_type


def _error_report_kwargs(sanitize_error_report: bool) -> Dict[str, bool]:
    """Return report options for internal vs user-facing error reports."""
    return {
        "sanitize_paths": sanitize_error_report,
        "include_line_content": not sanitize_error_report,
        "include_additional_info": not sanitize_error_report,
    }


def _display_report_path(path: Optional[str], sanitize_error_report: bool) -> str:
    if not path:
        return "(not specified)"
    if sanitize_error_report:
        return os.path.join("(redacted)", os.path.basename(path))
    return os.path.abspath(path)


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

HGVSG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9_]*\.\d+):g\.(\S+)$")
HGVSG_DELETION_TYPES = {
    "alu deletion",
    "copy number loss",
    "deletion",
    "herv deletion",
    "line1 deletion",
    "mobile element deletion",
    "sva deletion",
}
HGVSG_DUPLICATION_TYPES = {"copy number gain", "duplication", "tandem duplication"}
HGVSG_INSERTION_TYPES = {
    "alu insertion",
    "herv insertion",
    "insertion",
    "line1 insertion",
    "mobile element insertion",
    "novel sequence insertion",
    "sva insertion",
}


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


def _clean_hgvsg_value(value: Any) -> str:
    values = value if isinstance(value, list) else [value]
    for item in values:
        text = str(item or "").strip()
        if text and text != ".":
            return text
    return ""


def _positive_coordinate(row: Dict[str, Any], field: str) -> Optional[int]:
    value = str(row.get(field, "") or "").strip()
    if not value or value == ".":
        return None
    try:
        coordinate = int(value)
    except (TypeError, ValueError):
        return None
    return coordinate if coordinate >= 1 else None


def _hgvsg_start_boundary(row: Dict[str, Any]) -> Optional[str]:
    outer = _positive_coordinate(row, "Outer Start")
    inner = _positive_coordinate(row, "Inner Start")
    if outer is not None or inner is not None:
        return f"({outer if outer is not None else '?'}_{inner if inner is not None else '?'})"
    exact = _positive_coordinate(row, "Start")
    return str(exact) if exact is not None else None


def _hgvsg_stop_boundary(row: Dict[str, Any]) -> Optional[str]:
    inner = _positive_coordinate(row, "Inner Stop")
    outer = _positive_coordinate(row, "Outer Stop")
    if inner is not None or outer is not None:
        return f"({inner if inner is not None else '?'}_{outer if outer is not None else '?'})"
    exact = _positive_coordinate(row, "Stop")
    return str(exact) if exact is not None else None


def _hgvsg_interval(row: Dict[str, Any]) -> Optional[str]:
    start = _hgvsg_start_boundary(row)
    stop = _hgvsg_stop_boundary(row)
    if start is None or stop is None:
        return None
    return start if start == stop and not start.startswith("(") else f"{start}_{stop}"


def _hgvsg_insertion_site(row: Dict[str, Any]) -> Optional[str]:
    start = _hgvsg_start_boundary(row)
    stop = _hgvsg_stop_boundary(row)
    if start and start.startswith("("):
        return start
    if stop and stop.startswith("("):
        return stop
    exact = _positive_coordinate(row, "Start") or _positive_coordinate(row, "Stop")
    if exact is None:
        return None
    return f"{exact}_{exact + 1}"


def _hgvsg_insertion_payload(row: Dict[str, Any]) -> str:
    sequence = str(row.get("Sequence", "") or "").strip().upper()
    if sequence and sequence != "." and re.fullmatch(r"[ACGTBDHKMNRSVWY]+", sequence):
        return sequence
    length = _positive_coordinate(row, "Insertion Length")
    return f"N[{length}]" if length is not None else "N[?]"


def build_hgvsg(row: Dict[str, Any], reference: Any) -> str:
    """Build a conservative genomic HGVS expression from a normalized call row."""
    if reference is None:
        return "."
    call_type = normalize_call_type(row.get("Variant Call Type", ""))
    chrom = str(row.get("Chr", "") or "").strip()
    accession = reference.preferred_accession(chrom) if chrom and chrom != "." else None
    if not accession:
        return "."

    if call_type in HGVSG_INSERTION_TYPES:
        site = _hgvsg_insertion_site(row)
        return f"{accession}:g.{site}ins{_hgvsg_insertion_payload(row)}" if site else "."

    interval = _hgvsg_interval(row)
    if not interval:
        return "."
    if call_type in HGVSG_DELETION_TYPES:
        operation = "del"
    elif call_type in HGVSG_DUPLICATION_TYPES:
        operation = "dup"
    elif call_type == "inversion":
        operation = "inv"
    else:
        return "."
    return f"{accession}:g.{interval}{operation}"


def resolve_hgvsg(
    submitted_value: Any,
    row: Dict[str, Any],
    reference: Any,
    error_handler: ErrorHandler,
    *,
    variant_id: Optional[str] = None,
    line_number: Optional[int] = None,
    line_content: Optional[str] = None,
) -> str:
    """Preserve a valid submitted HGVSG or derive one from normalized placement."""
    submitted = _clean_hgvsg_value(submitted_value)
    if submitted:
        match = HGVSG_RE.fullmatch(submitted)
        valid = match is not None
        if valid and reference is not None:
            expected_chrom = reference.resolve_chrom(row.get("Chr", ""))
            submitted_chrom = reference.resolve_chrom(match.group(1))
            valid = expected_chrom is not None and submitted_chrom == expected_chrom
        if valid:
            return submitted
        error_handler.create_error(
            ErrorCode.INVALID_HGVSG,
            line_number=line_number,
            variant_id=variant_id,
            line_content=line_content,
            field_name=HGVSG_FIELD,
            expected_value="<versioned genomic accession>:g.<HGVS description> on the call chromosome",
            actual_value=submitted,
        )
    return build_hgvsg(row, reference)


class KVarTSVConverter:
    """Convert KVar SV VCF to Variant_Call.tsv format"""

    def __init__(
        self,
        error_handler: Optional[ErrorHandler] = None,
        reference_fasta_path: Optional[str] = None
    ):
        self.error_handler = error_handler or ErrorHandler()
        self.reference_fasta_path = reference_fasta_path
        self.organism_taxid = None
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
        metadata_file_path: Optional[str] = None,
        call_accession_start: Optional[int] = None,
        sanitize_error_report: bool = False,
    ) -> None:
        """Main method: convert VCF to TSV"""
        report_kwargs = _error_report_kwargs(sanitize_error_report)
        # Parse metadata file if provided
        metadata_validator = None
        if metadata_file_path:
            try:
                metadata_validator = MetadataValidator(self.error_handler)
                metadata_validator.parse_metadata_file(metadata_file_path)
                print(f"Metadata file loaded: {metadata_file_path}")
            except Exception as e:
                print(f"Warning: Error parsing metadata file: {e}")
            else:
                metadata_info = metadata_validator.metadata_file_info
                if metadata_info:
                    self.organism_taxid = MetadataParser(metadata_file_path).organism_taxid
                    self.parser.set_expected_metadata(
                        sampleset_id=metadata_info.sampleset_id,
                        experiment_id=metadata_info.experiment_id,
                        reference=metadata_info.reference
                    )

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
                        **report_kwargs,
                    )
                raise

        # Metadata validation
        if metadata_validator:
            # Extract SAMPLESET and EXPERIMENT from VCF (must be same across all rows)
            vcf_sampleset = None
            vcf_experiment = None
            # Try INFO tags first
            for row in self.parser.data_rows:
                if 'SAMPLESET' in row.info:
                    if vcf_sampleset is None:
                        vcf_sampleset = row.info['SAMPLESET']
                    elif vcf_sampleset != row.info['SAMPLESET']:
                        self.error_handler.create_error(
                            ErrorCode.METADATA_MISMATCH,
                            field_name="SAMPLESET",
                            expected_value="Same value in all rows",
                            actual_value="Multiple values found",
                            additional_info={"message": "SAMPLESET must be identical across all rows in VCF"}
                        )

                if 'EXPERIMENT' in row.info:
                    if vcf_experiment is None:
                        vcf_experiment = row.info['EXPERIMENT']
                    elif vcf_experiment != row.info['EXPERIMENT']:
                        self.error_handler.create_error(
                            ErrorCode.METADATA_MISMATCH,
                            field_name="EXPERIMENT",
                            expected_value="Same value in all rows",
                            actual_value="Multiple values found",
                            additional_info={"message": "EXPERIMENT must be identical across all rows in VCF"}
                        )

            # If not in INFO, use VCF header batch/population_id (Experiment_id/SampleSet_id mapping)
            if vcf_experiment is None and self.parser.header.metadata.batch:
                vcf_experiment = self.parser.header.metadata.batch
            if vcf_experiment is None and metadata_validator.metadata_file_info:
                vcf_experiment = metadata_validator.metadata_file_info.experiment_id

            if vcf_sampleset is None and self.parser.header.metadata.population_id:
                # population_id is a list; use first value (or check all if multiple)
                if len(self.parser.header.metadata.population_id) == 1:
                    vcf_sampleset = self.parser.header.metadata.population_id[0]
                elif len(self.parser.header.metadata.population_id) > 1:
                    # Multiple values: use first and report warning
                    vcf_sampleset = self.parser.header.metadata.population_id[0]
                    self.error_handler.create_error(
                        ErrorCode.METADATA_MISMATCH,
                        field_name="SampleSet_id/population_id",
                        expected_value="Single value",
                        actual_value=f"Multiple values: {', '.join(self.parser.header.metadata.population_id)}",
                        additional_info={
                            "message": "VCF header has multiple population_id; using first.",
                            "all_population_ids": self.parser.header.metadata.population_id
                        }
                    )
            if vcf_sampleset is None and metadata_validator.metadata_file_info:
                vcf_sampleset = metadata_validator.metadata_file_info.sampleset_id

            validation_passed = metadata_validator.validate_against_vcf(
                vcf_sampleset,
                vcf_experiment,
                self.parser.header.metadata.reference
            )

            if not validation_passed:
                print("Warning: Metadata validation failed - check error report")
            else:
                print("Metadata validation passed")

        # Build ID map and BND grouping
        id_map = self._build_id_map()
        mutation_id_map = self._build_mutation_id_map(id_map)
        call_records = self._build_call_records(id_map)
        if call_accession_start is not None:
            self._assign_call_accessions(call_records, call_accession_start)

        # Block final output when CRITICAL/ERROR messages were collected.
        self.error_handler.assert_no_blocking_errors(
            stage="SV VCF validation",
            output_file=error_report_path,
            vcf_file_path=vcf_file_path,
            output_tsv_path=output_file_path,
            **report_kwargs,
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
                    if report_has_errors:
                        self.error_handler.generate_report(
                            staged_error_report,
                            report_display_path=error_report_path,
                            vcf_file_path=vcf_file_path,
                            output_tsv_path=output_file_path,
                            **report_kwargs,
                        )
                    else:
                        with open(staged_error_report, 'w', encoding='utf-8') as f:
                            f.write("=" * 80 + "\n")
                            f.write("KVar SV VCF Parsing Error Report\n")
                            f.write("=" * 80 + "\n")
                            path_label_suffix = "" if sanitize_error_report else " absolute path"
                            f.write(
                                f"Input VCF{path_label_suffix}: "
                                f"{_display_report_path(vcf_file_path, sanitize_error_report)}\n"
                            )
                            f.write(
                                f"Output TSV{path_label_suffix}: "
                                f"{_display_report_path(output_file_path, sanitize_error_report)}\n"
                            )
                            f.write(
                                f"Error report{path_label_suffix}: "
                                f"{_display_report_path(error_report_path, sanitize_error_report)}\n"
                            )
                            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                            f.write("\nNo errors.\n")
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
                    **report_kwargs,
                )
            raise

        if error_report_path:
            report_suffix = "" if report_has_errors else " (no errors)"
            print(f"Error report: {error_report_path}{report_suffix}")

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
        return self.breakend_parser.inserted_sequence(row.alt, row.ref)

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
        return VCFDataRow(
            chrom=source_row.chrom,
            pos=source_row.pos,
            id=synthetic_id,
            ref=source_row.ref,
            alt=sequence,
            qual=source_row.qual,
            filter=source_row.filter,
            info=info,
        )

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

            submitted_ids = [row_id]
            source_by_id = {row_id: "submitted"}
            reason_by_id = {row_id: "validated"}
            inserted_sequence = self._bnd_inserted_sequence(row)
            split_mutation_id = self._bnd_split_mutation_id(row, inserted_sequence)
            mutation_info_override = None
            if split_mutation_id:
                mutation_info_override = {
                    "mutation_id": split_mutation_id,
                    "mutation_order": "1",
                    "mutation_molecule": ".",
                }

            if row.info.get("SVTYPE", "") == "BND" and row_id != ".":
                for mate_id in self._split_info_values(row.info.get("MATEID")):
                    if mate_id not in id_map or mate_id in absorbed_ids:
                        continue

                    mate_row = self.parser.data_rows[id_map[mate_id]]
                    if mate_row.info.get("SVTYPE", "") != "BND":
                        continue

                    reciprocal_mates = self._split_info_values(mate_row.info.get("MATEID"))
                    if row_id not in reciprocal_mates:
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
                    submitted_ids = [row_id, mate_id]
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

            call_record = {
                "row": row,
                "primary_id": row_id,
                "output_id": row_id,
                "submitted_ids": submitted_ids,
                "source_by_id": source_by_id,
                "reason_by_id": reason_by_id,
            }
            if mutation_info_override:
                call_record["mutation_info_override"] = mutation_info_override
            call_records.append(call_record)

            if inserted_sequence:
                synthetic_id = self._unique_generated_id(f"{row_id}_ins", generated_ids)
                synthetic_row = self._synthetic_insertion_row(row, synthetic_id, inserted_sequence)
                call_records.append({
                    "row": synthetic_row,
                    "primary_id": synthetic_id,
                    "output_id": synthetic_id,
                    "submitted_ids": [row_id],
                    "source_by_id": {row_id: "derived_from_submitted"},
                    "reason_by_id": {row_id: "split_bnd_insertion"},
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

    def _assign_call_accessions(self, call_records: List[Dict[str, Any]], call_accession_start: int) -> None:
        """Assign KVar call accessions to output call records."""
        counter = call_accession_start
        for call_record in call_records:
            call_record["output_id"] = f"kssv{counter}"
            counter += 1

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
                    # DDBJ uses EVENT as the submitter/dbVar-style Mutation ID anchor.
                    # Without EVENT, leave Mutation ID blank; related BND calls remain ungrouped.
                    mutation_id = "."
                    for variant_idx in group:
                        event_ids = self._split_info_values(self.parser.data_rows[variant_idx].info.get("EVENT"))
                        if event_ids:
                            mutation_id = event_ids[0]
                            break
                    if mutation_id == ".":
                        continue

                    for i, variant_idx in enumerate(group):
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
            if self.organism_taxid:
                f.write(f"##organism_taxid={self.organism_taxid}\n")
            f.write("#" + "\t".join(CALL_TSV_HEADER) + "\n")

            # Write data rows
            if call_records is None:
                call_records = [
                    {
                        "row": row,
                        "submitted_ids": [row.id if row.id and row.id != "." else "."]
                    }
                    for row in self.parser.data_rows
                ]

            for call_record in call_records:
                tsv_row = self._convert_row_to_tsv(
                    call_record["row"],
                    mutation_id_map,
                    call_record.get("submitted_ids"),
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
        submitted_call_ids: Optional[List[str]] = None,
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

        hgvsg_row = {
            "Variant Call Type": call_type,
            "Chr": chrom,
            "Outer Start": outer_start,
            "Start": start_value,
            "Inner Start": inner_start,
            "Inner Stop": inner_stop,
            "Stop": end_value,
            "Outer Stop": outer_stop,
            "Insertion Length": insertion_length,
            "Sequence": sequence,
        }
        submitted_hgvsg = next(
            (
                _clean_hgvsg_value(info.get(key))
                for key in ("HGVSG", "HGVSg", "hgvs_name")
                if _clean_hgvsg_value(info.get(key))
            ),
            "",
        )
        hgvsg = resolve_hgvsg(
            submitted_hgvsg,
            hgvsg_row,
            self.parser.reference,
            self.error_handler,
            variant_id=variant_call_id,
        )

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
            hgvsg,                # HGVSG
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
            ",".join(submitted_call_ids or [variant_call_id])  # Submitted_Variant_Call_IDs
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
