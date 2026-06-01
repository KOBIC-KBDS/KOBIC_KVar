#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate and rewrite a cleaned dbSNP VCF file."""

import argparse
import gzip
import os
import sys
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from kvar_snp_tools.VCF2dbSNP import get_vcf_info_output_id
    from kvar_snp_tools.dbSNP_parser import FormatTagDefinition, InfoTagDefinition, VCFDataRow, dbSNPVCFParser
    from kvar_snp_tools.error_handler import ErrorCode, ErrorHandler
    from kvar_snp_tools.metadata_validator import MetadataInfo, MetadataValidator
else:
    from .VCF2dbSNP import get_vcf_info_output_id
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

        self.parser.parse_file(vcf_file_path)

        if self.metadata_info:
            metadata_validator = MetadataValidator(self.error_handler)
            metadata_validator.metadata_file_info = self.metadata_info
            metadata_validator.validate_against_vcf(
                self._vcf_metadata_dict(),
                self.parser.header.population_ids,
            )
            self._apply_metadata_values()

        self._correct_vrt_values()

        if error_report_path is None:
            error_report_path = self._default_report_path(output_file_path)

        self.error_handler.assert_no_blocking_errors(
            stage="dbSNP VCF validation",
            output_file=error_report_path,
            vcf_file_path=vcf_file_path,
            output_tsv_path=output_file_path,
        )

        self._write_cleaned_vcf(output_file_path)
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

    def _correct_vrt_values(self) -> None:
        """Correct VRT values that do not match REF/ALT classification."""
        for row in self.parser.data_rows:
            expected_vrt = self._calculate_vrt(row.ref, row.alt)
            if expected_vrt is None:
                continue
            current_vrt = row.info.get("VRT")
            if current_vrt == expected_vrt:
                continue
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

    def _write_cleaned_vcf(self, output_file_path: str) -> None:
        """Write the cleaned dbSNP VCF file."""
        opener = gzip.open if output_file_path.endswith(".gz") else open
        with opener(output_file_path, "wt", encoding="utf-8") as handle:
            self._write_metadata(handle)
            self._write_info_tags(handle)
            self._write_format_tags(handle)
            self._write_population_ids(handle)
            self._write_column_header(handle)
            for row in self.parser.data_rows:
                self._write_data_row(handle, row)

    def _write_metadata(self, handle) -> None:
        """Write VCF metadata header lines."""
        metadata = self.parser.header.metadata
        handle.write(f"##fileformat={metadata.fileformat or 'VCFv4.1'}\n")
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

    def _write_info_tags(self, handle) -> None:
        """Write INFO tag definitions."""
        handle.write(CANONICAL_VRT_HEADER + "\n")
        written_ids = {"VRT"}
        for tag_def in self.parser.header.info_tags.values():
            output_id = get_vcf_info_output_id(tag_def.id)
            if output_id in written_ids:
                continue
            written_ids.add(output_id)
            handle.write(self._format_info_definition(output_id, tag_def) + "\n")

    @staticmethod
    def _format_info_definition(output_id: str, tag_def: InfoTagDefinition) -> str:
        """Return an INFO definition line."""
        description = tag_def.description.replace('"', "'")
        return (
            f'##INFO=<ID={output_id},Number={tag_def.number},'
            f'Type={tag_def.type},Description="{description}">'
        )

    def _write_format_tags(self, handle) -> None:
        """Write FORMAT tag definitions."""
        for tag_def in self.parser.header.format_tags.values():
            description = tag_def.description.replace('"', "'")
            handle.write(
                f'##FORMAT=<ID={tag_def.id},Number={tag_def.number},'
                f'Type={tag_def.type},Description="{description}">\n'
            )

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


def main() -> None:
    """Run the dbSNP VCF cleaner CLI."""
    parser = argparse.ArgumentParser(description="Validate and rewrite a cleaned dbSNP VCF file")
    parser.add_argument("--vcf", required=True, help="Input dbSNP VCF path")
    parser.add_argument("--output", required=True, help="Output cleaned dbSNP VCF path")
    parser.add_argument("--metadata", help="Optional metadata file path")
    parser.add_argument("--error-report", help="Optional validation report path")
    args = parser.parse_args()

    cleaner = DbSNPVCFCleaner()
    cleaner.clean(
        vcf_file_path=args.vcf,
        output_file_path=args.output,
        metadata_file_path=args.metadata,
        error_report_path=args.error_report,
    )
    print(f"Cleaned dbSNP VCF written: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
