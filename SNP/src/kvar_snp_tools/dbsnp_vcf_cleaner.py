#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate and rewrite a cleaned dbSNP VCF file."""

import argparse
import copy
import gzip
import os
import sys
import tempfile
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from kvar_snp_tools.VCF2dbSNP import (
        DBSNP_OUTPUT_INFO_IDS,
        DBSNP_OUTPUT_INFO_EXPECTED_VALUE,
        DBSNP_OUTPUT_INFO_ORDER,
        DBSNP_OUTPUT_INFO_TAG_DEFINITIONS,
        get_vcf_info_output_id,
        is_dbsnp_output_info_id,
    )
    from kvar_snp_tools.dbSNP_parser import FormatTagDefinition, InfoTagDefinition, VCFDataRow, dbSNPVCFParser
    from kvar_snp_tools.error_handler import ErrorCode, ErrorHandler
    from kvar_snp_tools.metadata_validator import MetadataInfo, MetadataValidator
else:
    from .VCF2dbSNP import (
        DBSNP_OUTPUT_INFO_IDS,
        DBSNP_OUTPUT_INFO_EXPECTED_VALUE,
        DBSNP_OUTPUT_INFO_ORDER,
        DBSNP_OUTPUT_INFO_TAG_DEFINITIONS,
        get_vcf_info_output_id,
        is_dbsnp_output_info_id,
    )
    from .dbSNP_parser import FormatTagDefinition, InfoTagDefinition, VCFDataRow, dbSNPVCFParser
    from .error_handler import ErrorCode, ErrorHandler
    from .metadata_validator import MetadataInfo, MetadataValidator


CANONICAL_VRT_HEADER = (
    '##INFO=<ID=VRT,Number=1,Type=Integer,Description="Variation type, '
    '1 - SNV: single nucleotide variation, '
    '2 - DIV: deletion/insertion variation, '
    '3 - HETEROZYGOUS: variable, but undefined at nucleotide level, '
    '4 - STR: short tandem repeat (microsatellite) variation, '
    '5 - NAMED: insertion/deletion variation of named repetitive element, '
    '6 - NO VARIATION: sequence scanned for variation, but none observed, '
    '7 - MIXED: cluster contains submissions from 2 or more allelic classes, '
    '8 - MNV: multiple nucleotide variation with alleles of common length greater than 1">'
)


class DbSNPVCFCleaner:
    """Validate dbSNP VCF input and write a cleaned dbSNP VCF output."""

    def __init__(self, error_handler: Optional[ErrorHandler] = None) -> None:
        self.error_handler = error_handler or ErrorHandler()
        self.parser = dbSNPVCFParser(self.error_handler)
        self.metadata_info: Optional[MetadataInfo] = None
        self.contig_lines: List[str] = []
        self._warned_unsupported_output_info_tags = set()
        self._used_output_info_id_set = set()

    def clean(
        self,
        vcf_file_path: str,
        output_file_path: str,
        metadata_file_path: Optional[str] = None,
        error_report_path: Optional[str] = None,
    ) -> None:
        """Validate and rewrite a dbSNP VCF file."""
        if metadata_file_path:
            metadata_validator = MetadataValidator(self.error_handler)
            self.metadata_info = metadata_validator.parse_metadata_file(metadata_file_path)

        self.parser.parse_header(vcf_file_path)
        self.contig_lines = self._read_contig_lines(vcf_file_path)
        self._record_unsupported_header_info_tags()

        if self.metadata_info:
            metadata_validator = MetadataValidator(self.error_handler)
            metadata_validator.metadata_file_info = self.metadata_info
            metadata_validator.validate_against_vcf(
                self._vcf_metadata_dict(),
                self.parser.header.population_ids,
            )
            self._apply_metadata_values()

        self._validate_streaming_rows(vcf_file_path)

        if error_report_path is None:
            error_report_path = self._default_report_path(output_file_path)

        self.error_handler.assert_no_blocking_errors(
            stage="dbSNP VCF validation",
            output_file=error_report_path,
            vcf_file_path=vcf_file_path,
            output_tsv_path=output_file_path,
        )

        self._write_cleaned_vcf(vcf_file_path, output_file_path)
        self.error_handler.generate_report(
            error_report_path,
            vcf_file_path=vcf_file_path,
            output_tsv_path=output_file_path,
        )

    @staticmethod
    def _default_report_path(output_file_path: str) -> str:
        """Return a default text report path next to the output VCF."""
        for suffix in (".vcf.gz", ".vcf"):
            if output_file_path.endswith(suffix):
                return output_file_path[: -len(suffix)] + ".errors.txt"
        return output_file_path + ".errors.txt"

    @staticmethod
    def _read_contig_lines(vcf_file_path: str) -> List[str]:
        """Return input contig header lines in their original order."""
        opener = gzip.open if vcf_file_path.endswith(".gz") else open
        contig_lines: List[str] = []
        with opener(vcf_file_path, "rt", encoding="utf-8") as handle:
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith("#CHROM"):
                    break
                if line.startswith("##contig="):
                    contig_lines.append(line)
        return contig_lines

    def _vcf_metadata_dict(self) -> Dict[str, Any]:
        """Return parsed VCF metadata as a dictionary."""
        metadata = self.parser.header.metadata
        return {
            "fileformat": metadata.fileformat,
            "filedate": metadata.filedate,
            "handle": metadata.handle,
            "batch": metadata.batch,
            "bioproject_id": metadata.bioproject_id,
            "biosample_id": metadata.biosample_id,
            "reference": metadata.reference,
        }

    def _apply_metadata_values(self) -> None:
        """Apply trusted metadata values to the output header when allowed."""
        if not self.metadata_info:
            return

        header_metadata = self.parser.header.metadata

        if self.metadata_info.experiment_id:
            header_metadata.batch = self.metadata_info.experiment_id
        if self.metadata_info.bioproject_id:
            header_metadata.bioproject_id = self.metadata_info.bioproject_id
        if self.metadata_info.reference and not header_metadata.reference:
            header_metadata.reference = self.metadata_info.reference
        if self.metadata_info.sampleset_ids:
            self.parser.header.population_ids = list(self.metadata_info.sampleset_ids)

    def _validate_streaming_rows(self, vcf_file_path: str) -> None:
        """Validate data rows without retaining them all in memory."""
        self._used_output_info_id_set.clear()
        for row in self.parser.iter_data_rows(vcf_file_path, store_rows=False):
            self._record_unsupported_row_info_tags(row)
            self._record_used_output_info_ids(row.info)
            self._correct_vrt_value(row, record_warning=True)
        self.parser.validate_parsed_data()

    def _correct_vrt_value(self, row: VCFDataRow, record_warning: bool) -> None:
        """Correct one row's VRT value when it differs from REF/ALT classification."""
        expected_vrt = self._calculate_vrt(row.ref, row.alt)
        if expected_vrt is None:
            return
        current_vrt = row.info.get("VRT")
        if current_vrt == expected_vrt:
            return
        if record_warning:
            self.error_handler.create_error(
                ErrorCode.VRT_REF_ALT_MISMATCH,
                field_name="VRT",
                expected_value=str(expected_vrt),
                actual_value=str(current_vrt),
                additional_info={
                    "variant": f"{row.chrom}:{row.pos}:{row.ref}>{row.alt}",
                    "action": "Output VCF uses the REF/ALT-derived VRT value",
                },
            )
        row.info["VRT"] = expected_vrt

    @staticmethod
    def _calculate_vrt(ref: str, alt: str) -> Optional[int]:
        """Calculate dbSNP VRT from REF and ALT."""
        if not ref or not alt or alt == ".":
            return None

        first_alt = alt.split(",", 1)[0].strip()
        if first_alt.startswith("<") or first_alt.startswith("["):
            return 5

        ref_len = len(ref)
        alt_len = len(first_alt)

        if ref_len == 1 and alt_len == 1 and ref.upper() != first_alt.upper():
            return 1
        if alt_len > ref_len and first_alt.upper().startswith(ref.upper()):
            return 2
        if ref_len > alt_len and ref.upper().startswith(first_alt.upper()):
            return 2
        if ref_len == alt_len and ref_len > 1:
            return 8
        if ref_len != alt_len:
            return 2
        return 1

    def _write_cleaned_vcf(self, vcf_file_path: str, output_file_path: str) -> None:
        """Atomically write the cleaned dbSNP VCF file after validation succeeds."""
        output_dir = os.path.dirname(os.path.abspath(output_file_path)) or "."
        output_base = os.path.basename(output_file_path)
        fd, temp_path = tempfile.mkstemp(
            prefix=f".{output_base}.",
            suffix=".tmp",
            dir=output_dir,
        )
        os.close(fd)

        try:
            self._write_cleaned_vcf_to_path(vcf_file_path, temp_path, output_file_path.endswith(".gz"))
            os.replace(temp_path, output_file_path)
        except Exception:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
            raise

    def _write_cleaned_vcf_to_path(self, vcf_file_path: str, path: str, gzip_output: bool) -> None:
        """Stream cleaned VCF content to an already-created temporary path."""
        opener = gzip.open if gzip_output else open
        writer_parser = dbSNPVCFParser(ErrorHandler())
        writer_parser.header = copy.deepcopy(self.parser.header)
        with opener(path, "wt", encoding="utf-8") as handle:
            self._write_metadata(handle)
            self._write_info_tags(handle)
            self._write_format_tags(handle)
            self._write_population_ids(handle)
            self._write_column_header(handle)
            for row in writer_parser.iter_data_rows(vcf_file_path, store_rows=False):
                self._correct_vrt_value(row, record_warning=False)
                self._write_data_row(handle, row)

    def _write_metadata(self, handle) -> None:
        """Write VCF metadata header lines."""
        metadata = self.parser.header.metadata
        handle.write("##fileformat=VCFv4.1\n")
        handle.write(f"##fileDate={metadata.filedate or datetime.now().strftime('%Y%m%d')}\n")

        if metadata.handle:
            handle.write(f"##handle={metadata.handle}\n")
        if metadata.batch:
            handle.write(f"##batch={metadata.batch}\n")
        if metadata.bioproject_id:
            handle.write(f"##bioproject_id={metadata.bioproject_id}\n")
        if metadata.biosample_id:
            handle.write(f"##biosample_id={metadata.biosample_id}\n")
        if metadata.reference:
            handle.write(f"##reference={metadata.reference}\n")
        for contig_line in self.contig_lines:
            handle.write(contig_line + "\n")

    def _write_info_tags(self, handle) -> None:
        """Write INFO tag definitions."""
        handle.write(CANONICAL_VRT_HEADER + "\n")
        written_ids = {"VRT"}
        used_ids = self._used_output_info_ids()
        for tag_def in self.parser.header.info_tags.values():
            output_id = get_vcf_info_output_id(tag_def.id)
            if output_id not in used_ids or output_id not in DBSNP_OUTPUT_INFO_IDS:
                continue
            if output_id in written_ids:
                continue
            written_ids.add(output_id)
            handle.write(self._format_info_definition(output_id, tag_def) + "\n")

        for output_id in DBSNP_OUTPUT_INFO_ORDER:
            if output_id in written_ids or output_id not in used_ids:
                continue
            if output_id not in DBSNP_OUTPUT_INFO_TAG_DEFINITIONS:
                continue
            number, type_str, description = DBSNP_OUTPUT_INFO_TAG_DEFINITIONS[output_id]
            written_ids.add(output_id)
            handle.write(
                f'##INFO=<ID={output_id},Number={number},'
                f'Type={type_str},Description="{description}">\n'
            )

    def _used_output_info_ids(self) -> set:
        """Return retained INFO IDs that are present in parsed rows."""
        if self._used_output_info_id_set:
            return set(self._used_output_info_id_set)

        used_ids = set()
        for row in self.parser.data_rows:
            self._record_used_output_info_ids(row.info, used_ids)
        return used_ids

    def _record_used_output_info_ids(self, info: Dict[str, Any], used_ids: Optional[set] = None) -> None:
        """Record dbSNP output INFO IDs present in one row."""
        target = used_ids if used_ids is not None else self._used_output_info_id_set
        for tag_id in info:
            output_id = get_vcf_info_output_id(tag_id)
            if output_id not in DBSNP_OUTPUT_INFO_IDS:
                continue
            target.add(output_id)

    def _record_unsupported_header_info_tags(self) -> None:
        """Record unsupported INFO header tags that will be excluded from output."""
        for tag_id in self.parser.header.info_tags:
            if not is_dbsnp_output_info_id(tag_id):
                self._warn_unsupported_output_info_tag(tag_id, context="INFO header")

    def _record_unsupported_row_info_tags(self, row: VCFDataRow) -> None:
        """Record unsupported INFO data tags that will be excluded from output."""
        for tag_id in row.info:
            if not is_dbsnp_output_info_id(tag_id):
                self._warn_unsupported_output_info_tag(tag_id, context="INFO field")

    def _warn_unsupported_output_info_tag(self, tag_id: str, context: str) -> None:
        """Record one warning for an INFO tag excluded from dbSNP VCF output."""
        output_id = get_vcf_info_output_id(tag_id)
        warning_key = f"{tag_id}->{output_id}"
        if warning_key in self._warned_unsupported_output_info_tags:
            return
        self._warned_unsupported_output_info_tags.add(warning_key)
        self.error_handler.create_error(
            ErrorCode.UNSUPPORTED_DBSNP_INFO_TAG,
            field_name=tag_id,
            expected_value=DBSNP_OUTPUT_INFO_EXPECTED_VALUE,
            actual_value=output_id,
            additional_info={
                "action": "Excluded from dbSNP VCF output",
                "context": context,
            },
        )

    @staticmethod
    def _format_info_definition(output_id: str, tag_def: InfoTagDefinition) -> str:
        """Return an INFO definition line."""
        description = DbSNPVCFCleaner._escape_vcf_description(tag_def.description)
        return (
            f'##INFO=<ID={output_id},Number={tag_def.number},'
            f'Type={tag_def.type},Description="{description}">'
        )

    def _write_format_tags(self, handle) -> None:
        """Write FORMAT tag definitions."""
        for tag_def in self.parser.header.format_tags.values():
            description = self._escape_vcf_description(tag_def.description)
            handle.write(
                f'##FORMAT=<ID={tag_def.id},Number={tag_def.number},'
                f'Type={tag_def.type},Description="{description}">\n'
            )

    @staticmethod
    def _escape_vcf_description(description: str) -> str:
        """Escape description text for a VCF header line."""
        return description.replace('\\', '\\\\').replace('"', '\\"')

    def _write_population_ids(self, handle) -> None:
        """Write population_id header lines."""
        for population_id in self.parser.header.population_ids:
            handle.write(f"##population_id={population_id}\n")

    def _write_column_header(self, handle) -> None:
        """Write the VCF column header."""
        columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        if self.parser.header.population_ids:
            columns.append("FORMAT")
            columns.extend(self.parser.header.population_ids)
        handle.write("\t".join(columns) + "\n")

    def _write_data_row(self, handle, row: VCFDataRow) -> None:
        """Write one VCF data row."""
        fields = [
            row.chrom,
            str(row.pos),
            row.id,
            row.ref,
            row.alt,
            row.qual or ".",
            row.filter or ".",
            self._format_info_field(row.info),
        ]

        if self.parser.header.population_ids:
            format_keys = self._format_keys_for_row(row)
            fields.append(":".join(format_keys) if format_keys else ".")
            for population_id in self.parser.header.population_ids:
                population_values = row.population_data.get(population_id, {})
                fields.append(self._format_population_field(population_values, format_keys))

        handle.write("\t".join(fields) + "\n")

    @staticmethod
    def _format_keys_for_row(row: VCFDataRow) -> List[str]:
        """Return FORMAT keys used by the row."""
        keys: List[str] = []
        for population_values in row.population_data.values():
            for key in population_values:
                if key not in keys:
                    keys.append(key)
        return keys

    @staticmethod
    def _format_population_field(population_values: Dict[str, Any], format_keys: Iterable[str]) -> str:
        """Return one population data column."""
        keys = list(format_keys)
        if not keys or not population_values:
            return "."
        return ":".join(DbSNPVCFCleaner._format_value(population_values.get(key, ".")) for key in keys)

    @staticmethod
    def _format_info_field(info: Dict[str, Any]) -> str:
        """Return a VCF INFO field string."""
        if not info:
            return "."
        parts: List[str] = []
        for key, value in info.items():
            output_key = get_vcf_info_output_id(key)
            if output_key not in DBSNP_OUTPUT_INFO_IDS:
                continue
            if isinstance(value, bool) and value:
                parts.append(output_key)
            else:
                parts.append(f"{output_key}={DbSNPVCFCleaner._format_value(value)}")
        return ";".join(parts)

    @staticmethod
    def _format_value(value: Any) -> str:
        """Format a scalar or list for VCF output."""
        if isinstance(value, list):
            return ",".join(str(item) for item in value)
        return str(value)
