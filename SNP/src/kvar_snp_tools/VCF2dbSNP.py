#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF to dbSNP converter.
Converts generic VCF files to dbSNP VCF submission format.
"""

import sys
import os
import gzip
import re
import argparse
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple
from collections import defaultdict

VCF_INFO_ID_RENAMES = {
    '1kGP_AF': 'GP1K_AF',
    '1kGP_AFR_AF': 'GP1K_AFR_AF',
    '1kGP_AMR_AF': 'GP1K_AMR_AF',
    '1kGP_EAS_AF': 'GP1K_EAS_AF',
    '1kGP_EUR_AF': 'GP1K_EUR_AF',
    '1kGP_SAS_AF': 'GP1K_SAS_AF',
    '61KJPN': 'JPN61K',
    '61KJPN_AF': 'JPN61K_AF',
    '61KJPN_AN': 'JPN61K_AN',
    '61KJPN_nhomalt': 'JPN61K_nhomalt',
}

VCF_INFO_ID_PATTERN = re.compile(r'^[A-Za-z_][0-9A-Za-z_.]*$')


def get_vcf_info_output_id(info_id: str) -> str:
    """Return a VCF-safe INFO ID for output."""
    if info_id in VCF_INFO_ID_RENAMES:
        return VCF_INFO_ID_RENAMES[info_id]
    if VCF_INFO_ID_PATTERN.match(info_id):
        return info_id
    cleaned = re.sub(r'[^0-9A-Za-z_.]', '_', info_id)
    if not cleaned or not re.match(r'^[A-Za-z_]', cleaned):
        cleaned = f'X_{cleaned}'
    return cleaned

# Support relative imports.
if __name__ == "__main__":
    from dbSNP_parser import FormatTagDefinition, InfoTagDefinition, dbSNPVCFParser, VCFDataRow
    from error_handler import ErrorHandler, ErrorCode
    from metadata_validator import MetadataValidator, MetadataInfo
else:
    from .dbSNP_parser import FormatTagDefinition, InfoTagDefinition, dbSNPVCFParser, VCFDataRow
    from .error_handler import ErrorHandler, ErrorCode
    from .metadata_validator import MetadataValidator, MetadataInfo


class VCF2dbSNPConverter:
    """Convert generic VCF files to dbSNP VCF format."""
    
    def __init__(self, error_handler: Optional[ErrorHandler] = None):
        self.error_handler = error_handler or ErrorHandler()
        self.parser = dbSNPVCFParser(self.error_handler)
        self.metadata_info: Optional[MetadataInfo] = None
        self.contig_lines: List[str] = []
    
    def convert_vcf_to_dbsnp(
        self,
        vcf_file_path: str,
        output_file_path: str,
        metadata_file_path: Optional[str] = None,
        error_report_path: Optional[str] = None
    ) -> None:
        """Convert a VCF file to dbSNP VCF format."""
        # Parse metadata file.
        if metadata_file_path:
            try:
                metadata_validator = MetadataValidator(self.error_handler)
                self.metadata_info = metadata_validator.parse_metadata_file(metadata_file_path)
                print(f"Metadata file loaded: {metadata_file_path}")
            except Exception as e:
                print(f"Warning: Error parsing metadata file: {e}")
                self.metadata_info = None
        
        # Parse VCF file in generic VCF format.
        try:
            self._parse_vcf_file(vcf_file_path)
        except Exception as e:
            if self.error_handler.has_critical_errors():
                print("Conversion stopped due to critical errors.")
                if error_report_path:
                    self.error_handler.generate_report(
                        error_report_path,
                        vcf_file_path=vcf_file_path,
                        output_tsv_path=output_file_path
                    )
                raise
        
        # Validate and apply metadata.
        self._validate_and_fix_metadata()

        if error_report_path is None:
            error_report_path = output_file_path.replace('.vcf', '_errors.txt')

        self.error_handler.assert_no_blocking_errors(
            stage="VCF2dbSNP conversion",
            output_file=error_report_path,
            vcf_file_path=vcf_file_path,
            output_tsv_path=output_file_path
        )
        
        # Write dbSNP-formatted VCF output.
        try:
            self._write_dbsnp_vcf(output_file_path)
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FILE_WRITE_ERROR,
                additional_info={"file_path": output_file_path, "error": str(e)}
            )
            raise
        
        self.error_handler.generate_report(
            error_report_path,
            vcf_file_path=vcf_file_path,
            output_tsv_path=output_file_path
        )
        print(f"Error report: {error_report_path}")
        
        print(f"Conversion complete: {vcf_file_path} -> {output_file_path}")
        print(f"Total {len(self.parser.data_rows)} variants converted")
        
        # Print error summary.
        self.error_handler.print_summary()
    
    def _parse_vcf_file(self, vcf_file_path: str) -> None:
        """Parse a VCF file with gzip support."""
        self.parser._seen_local_ids.clear()
        self.parser._seen_variant_sites.clear()

        # Split header and data while reading line by line.
        header_lines = []
        data_lines = []
        in_header = True
        
        # Check whether the input is gzip-compressed.
        if vcf_file_path.endswith('.gz'):
            try:
                with gzip.open(vcf_file_path, 'rt', encoding='utf-8') as f:
                    for line_num, line in enumerate(f, 1):
                        original_line = line
                        line = line.strip()
                        if not line:
                            continue
                        
                        if line.startswith('#'):
                            header_lines.append(line)
                        else:
                            if in_header:
                                in_header = False
                            # Buffer data rows for chunked parsing.
                            data_lines.append((line_num, original_line.strip()))
                            
                            # Parse rows in fixed-size chunks.
                            if len(data_lines) >= 10000:
                                self._parse_data(data_lines)
                                data_lines = []
            except Exception as e:
                self.error_handler.create_error(
                    ErrorCode.FILE_READ_ERROR,
                    additional_info={"file_path": vcf_file_path, "error": str(e)}
                )
                raise
        else:
            try:
                with open(vcf_file_path, 'r', encoding='utf-8') as f:
                    for line_num, line in enumerate(f, 1):
                        original_line = line
                        line = line.strip()
                        if not line:
                            continue
                        
                        if line.startswith('#'):
                            header_lines.append(line)
                        else:
                            if in_header:
                                in_header = False
                            data_lines.append((line_num, original_line.strip()))
                            
                            # Parse rows in fixed-size chunks.
                            if len(data_lines) >= 10000:
                                self._parse_data(data_lines)
                                data_lines = []
            except Exception as e:
                self.error_handler.create_error(
                    ErrorCode.FILE_READ_ERROR,
                    additional_info={"file_path": vcf_file_path, "error": str(e)}
                )
                raise
        
        # Parse header lines.
        self._parse_header(header_lines)
        
        # Parse remaining data rows.
        if data_lines:
            self._parse_data(data_lines)
        
        print(f"Parsing complete: {len(header_lines)} header lines, {len(self.parser.data_rows)} variants")
    
    def _parse_header(self, header_lines: List[str]) -> None:
        """Parse VCF header lines."""
        for line_num, line in enumerate(header_lines, 1):
            if line.startswith('##'):
                self._parse_metadata_line(line, line_num)
            elif line.startswith('#'):
                self._parse_column_header(line, line_num)
    
    def _parse_metadata_line(self, line: str, line_number: int) -> None:
        """Parse a metadata line."""
        if '=' not in line:
            return
        
        try:
            key, value = line[2:].split('=', 1)
            key = key.lower().strip()
            value = value.strip()
            
            if key == 'fileformat':
                self.parser.header.metadata.fileformat = value
            elif key == 'filedate':
                self.parser.header.metadata.filedate = value
            elif key == 'handle':
                self.parser.header.metadata.handle = value
            elif key == 'batch':
                self.parser.header.metadata.batch = value
            elif key == 'bioproject_id':
                self.parser.header.metadata.bioproject_id = value
            elif key == 'biosample_id':
                self.parser.header.metadata.biosample_id = value
            elif key == 'reference':
                self.parser.header.metadata.reference = value
            elif key == 'population_id':
                self.parser.header.population_ids.append(value)
            elif key == 'contig':
                self.contig_lines.append(line)
            elif key == 'info':
                self._parse_info_tag_definition(line, line_number)
            elif key == 'format':
                self._parse_format_tag_definition(line, line_number)
        except Exception as e:
            pass  # Let the main parser handle malformed metadata elsewhere.
    
    def _parse_info_tag_definition(self, line: str, line_number: int) -> None:
        """Parse an INFO tag definition."""
        pattern = r'##INFO=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)"'
        match = re.search(pattern, line)
        
        if match:
            tag_id, number, type_str, description = match.groups()
            self.parser.header.info_tags[tag_id] = InfoTagDefinition(
                id=tag_id,
                number=number,
                type=type_str,
                description=description
            )
    
    def _parse_format_tag_definition(self, line: str, line_number: int) -> None:
        """Parse a FORMAT tag definition."""
        pattern = r'##FORMAT=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)"'
        match = re.search(pattern, line)
        
        if match:
            tag_id, number, type_str, description = match.groups()
            self.parser.header.format_tags[tag_id] = FormatTagDefinition(
                id=tag_id,
                number=number,
                type=type_str,
                description=description
            )
    
    def _parse_column_header(self, line: str, line_number: int) -> None:
        """Parse the column header."""
        columns = line[1:].split('\t')
        self.parser.header.column_header = columns
        
        # Ignore individual sample columns here; population IDs come from metadata.
        # Sample columns are aggregated into population data later.
    
    def _parse_data(self, data_lines: List[Tuple[int, str]]) -> None:
        """Parse data rows."""
        for line_number, line in data_lines:
            row = self._parse_data_row(line, line_number)
            if row:
                self.parser.data_rows.append(row)
    
    def _parse_data_row(self, line: str, line_number: int) -> Optional[VCFDataRow]:
        """Parse a single data row"""
        fields = line.split('\t')
        
        if len(fields) < 8:
            return None
        
        # Basic fields
        chrom = fields[0]
        
        try:
            pos = int(fields[1])
        except ValueError:
            return None
        
        id_field = fields[2]
        ref = fields[3]
        alt = fields[4]
        qual = fields[5] if fields[5] != '.' else None
        filter_field = fields[6] if fields[6] != '.' else None

        if not self.parser._validate_local_id(id_field, line, line_number):
            return None
        
        # REF/ALT empty check
        if not ref or ref == '.':
            self.error_handler.create_error(
                ErrorCode.EMPTY_REF_ALT,
                line_number=line_number,
                line_content=line,
                field_name="REF"
            )
            return None
        
        if not alt or alt == '.':
            self.error_handler.create_error(
                ErrorCode.EMPTY_REF_ALT,
                line_number=line_number,
                line_content=line,
                field_name="ALT"
            )
            return None

        if ',' in alt:
            self.error_handler.create_error(
                ErrorCode.MULTIALLELIC_ALT,
                line_number=line_number,
                line_content=line,
                field_name="ALT",
                actual_value=alt,
                expected_value="single ALT allele after normalization"
            )
            return None

        if ref.upper() == alt.upper():
            self.error_handler.create_error(
                ErrorCode.SAME_REF_ALT,
                line_number=line_number,
                line_content=line,
                field_name="REF/ALT",
                actual_value=f"REF={ref}, ALT={alt}",
                expected_value="different REF and ALT alleles"
            )
            return None
        
        # REF base validation: only A, T, G, C allowed (no N, no IUPAC ambiguity codes, no '*')
        valid_bases = set('ATGCatgc')
        if not all(c in valid_bases for c in ref):
            invalid_chars = sorted(set(c for c in ref if c not in valid_bases))
            self.error_handler.create_error(
                ErrorCode.INVALID_REF_BASE,
                line_number=line_number,
                line_content=line,
                field_name="REF",
                actual_value=ref,
                expected_value="Only A, T, G, C characters allowed",
                additional_info={
                    "invalid_characters": invalid_chars,
                    "message": f"REF '{ref}' contains invalid character(s): {invalid_chars}. "
                               f"Only strict ATGC bases are allowed. "
                               f"N, IUPAC ambiguity codes (e.g. R, Y, S, W, K, M, B, D, H, V), "
                               f"and '*' (missing allele) are not permitted."
                }
            )
            return None
        
        # ALT allele validation: each allele must consist of only A, T, G, C
        alt_alleles = alt.split(',')
        for alt_allele in alt_alleles:
            alt_allele_stripped = alt_allele.strip()
            if not alt_allele_stripped:
                continue
            if not all(c in valid_bases for c in alt_allele_stripped):
                invalid_chars = sorted(set(c for c in alt_allele_stripped if c not in valid_bases))
                self.error_handler.create_error(
                    ErrorCode.INVALID_ALT_ALLELE,
                    line_number=line_number,
                    line_content=line,
                    field_name="ALT",
                    actual_value=alt,
                    expected_value="Only A, T, G, C characters allowed per allele",
                    additional_info={
                        "invalid_allele": alt_allele_stripped,
                        "invalid_characters": invalid_chars,
                        "message": f"ALT allele '{alt_allele_stripped}' contains invalid character(s): {invalid_chars}. "
                                   f"Only strict ATGC bases are allowed. "
                                   f"N, IUPAC ambiguity codes (e.g. R, Y, S, W, K, M, B, D, H, V), "
                                   f"and '*' (missing allele) are not permitted."
                    }
                )
                return None

        if len(ref) != len(alt) and ref[0].upper() != alt[0].upper():
            self.error_handler.create_error(
                ErrorCode.INVALID_INDEL_LEADING_BASE,
                line_number=line_number,
                line_content=line,
                field_name="REF/ALT",
                actual_value=f"REF={ref}, ALT={alt}",
                expected_value="same leading base for indels"
            )
            return None
        
        # Variant length check: exclude rows exceeding this limit.
        ref_len = len(ref)
        for alt_allele in alt_alleles:
            alt_allele_stripped = alt_allele.strip()
            alt_len = len(alt_allele_stripped)
            variant_len = max(ref_len, alt_len)
            if variant_len > 50:
                self.error_handler.create_error(
                    ErrorCode.VARIANT_TOO_LONG,
                    line_number=line_number,
                    line_content=line,
                    field_name="REF/ALT",
                    actual_value=f"REF={ref_len}bp, ALT={alt_len}bp (variant_length={variant_len}bp)",
                    expected_value="50bp or less",
                    additional_info={
                        "ref_length": ref_len,
                        "alt_length": alt_len,
                        "variant_length": variant_len,
                        "message": "Variant exceeds 50bp; submit as SV for variants >50bp."
                    }
                )
                return None
        
        # Parse INFO field.
        info_dict = self._parse_info_field(fields[7], line_number, line)
        
        # Calculate VRT when it is missing.
        if 'VRT' not in info_dict:
            vrt = self._calculate_vrt(ref, alt)
            if vrt:
                info_dict['VRT'] = vrt
        elif not isinstance(info_dict['VRT'], int) or info_dict['VRT'] < 1 or info_dict['VRT'] > 8:
            self.error_handler.create_error(
                ErrorCode.INVALID_VRT_VALUE,
                line_number=line_number,
                line_content=line,
                field_name="VRT",
                actual_value=str(info_dict['VRT']),
                expected_value="integer in range 1-8"
            )
            return None

        if not self.parser._validate_ac_an(info_dict, line, line_number):
            return None

        self.parser._derive_af_from_ac_an(info_dict, line, line_number)

        if not self.parser._validate_duplicate_site(chrom, pos, ref, alt, info_dict, line, line_number):
            return None
        
        # Parse FORMAT and sample data.
        format_data = {}
        population_data = {}
        
        if len(fields) > 8:
            format_data, population_data = self._parse_format_and_sample_data(
                fields[8:], line_number, info_dict, line
            )
        
        return VCFDataRow(
            chrom=chrom,
            pos=pos,
            id=id_field,
            ref=ref,
            alt=alt,
            qual=qual,
            filter=filter_field,
            info=info_dict,
            format_data=format_data,
            population_data=population_data
        )
    
    def _parse_info_field(self, info_str: str, line_number: int, line: Optional[str] = None) -> Dict[str, Any]:
        """Parse an INFO field."""
        info_dict = {}
        
        if info_str == '.':
            return info_dict
        
        pairs = info_str.split(';')
        
        for pair in pairs:
            if not pair.strip():
                continue
            if '=' in pair:
                key, value = pair.split('=', 1)
                self.parser._warn_if_undefined_info_tag(key, line_number, line)
                self.parser._validate_target_info_value(key, value, line_number)
                info_dict[key] = self._convert_info_value(key, value, line_number)
            else:
                self.parser._warn_if_undefined_info_tag(pair, line_number, line)
                info_dict[pair] = True
        
        return info_dict
    
    def _convert_info_value(self, key: str, value: str, line_number: int) -> Any:
        """Convert INFO values for known numeric tags."""
        # Convert known integer-like tags.
        if key in ['VRT', 'NIO', 'SAO', 'SSR', 'PMID', 'AN', 'AC']:
            try:
                if ',' in value:
                    return [int(v) for v in value.split(',')]
                return int(value)
            except ValueError:
                self.error_handler.create_error(
                    ErrorCode.INVALID_INFO_TAG_TYPE,
                    line_number=line_number,
                    field_name=key,
                    actual_value=value,
                    expected_value="integer"
                )
                return value
        
        # Convert known float-like tags.
        if key in ['AF', 'MAF']:
            try:
                if ',' in value:
                    return [float(v) for v in value.split(',')]
                return float(value)
            except ValueError:
                self.error_handler.create_error(
                    ErrorCode.INVALID_INFO_TAG_TYPE,
                    line_number=line_number,
                    field_name=key,
                    actual_value=value,
                    expected_value="float"
                )
                return value
        
        return value
    
    def _parse_format_and_sample_data(
        self, fields: List[str], line_number: int, info_dict: Dict[str, Any], line: Optional[str] = None
    ) -> Tuple[Dict[str, Any], Dict[str, Dict[str, Any]]]:
        """Parse FORMAT and sample data."""
        format_data = {}
        population_data = {}
        
        if not fields:
            return format_data, population_data
        
        # FORMAT field.
        format_str = fields[0]
        format_keys = format_str.split(':') if format_str else []
        for format_key in format_keys:
            self.parser._warn_if_undefined_format_tag(format_key, line_number, line)
        
        # Sample data fields from individual sample columns.
        sample_fields = fields[1:] if len(fields) > 1 else []
        
        # Prefer population ID from metadata, then VCF header, then fallback.
        if self.metadata_info and self.metadata_info.sampleset_ids:
            pop_id = self.metadata_info.sampleset_ids[0]
        elif self.parser.header.population_ids:
            pop_id = self.parser.header.population_ids[0]
        else:
            pop_id = 'DEFAULT_POPULATION'
        
        # Prefer NA and FRQ values derived from the INFO field.
        info_na = None
        info_frq = None
        
        # 1. Prefer INFO AN and AF when both are present.
        if 'AN' in info_dict and 'AF' in info_dict:
            try:
                an_value = info_dict['AN']
                info_na = int(an_value) if not isinstance(an_value, int) else an_value
            except (ValueError, TypeError):
                pass
            
            try:
                af_value = info_dict['AF']
                if isinstance(af_value, list):
                    info_frq = float(af_value[0]) if af_value else None
                else:
                    info_frq = float(af_value)
            except (ValueError, TypeError):
                pass
        
        # 2. Calculate FRQ from INFO AN and AC.
        elif 'AN' in info_dict and 'AC' in info_dict:
            try:
                an_value = info_dict['AN']
                info_na = int(an_value) if not isinstance(an_value, int) else an_value
            except (ValueError, TypeError):
                pass
            
            try:
                ac_value = info_dict['AC']
                if isinstance(ac_value, list):
                    total_ac = sum(int(v) if not isinstance(v, int) else v for v in ac_value)
                else:
                    total_ac = int(ac_value) if not isinstance(ac_value, int) else ac_value
                
                if info_na and info_na > 0:
                    info_frq = total_ac / info_na
            except (ValueError, TypeError, ZeroDivisionError):
                pass
        
        # 3. Use INFO AF and estimate AN from sample count when possible.
        elif 'AF' in info_dict:
            try:
                af_value = info_dict['AF']
                if isinstance(af_value, list):
                    info_frq = float(af_value[0]) if af_value else None
                else:
                    info_frq = float(af_value)
            except (ValueError, TypeError):
                pass
            
            # Estimate AN from sample count when AN is missing.
            if sample_fields:
                info_na = len(sample_fields) * 2
        
        # Use values derived from INFO when available.
        if info_na is not None and info_frq is not None:
            population_data[pop_id] = {
                'NA': info_na,
                'FRQ': info_frq
            }
            return format_data, population_data
        
        # 4. Calculate directly from sample data when INFO values are unavailable.
        # Individual genotype fields are used only when INFO AF/AC/AN is unavailable.
        if not sample_fields:
            # Return empty population data when no sample data is available.
            return format_data, population_data
        
        total_samples = len(sample_fields)
        total_alleles = 0  # Total allele count.
        alt_allele_count = 0  # Alternate allele count.
        
        # Locate FORMAT field indexes.
        gt_idx = None
        af_idx = None
        ac_idx = None
        for i, key in enumerate(format_keys):
            if key == 'GT':
                gt_idx = i
            elif key == 'AF':
                af_idx = i
            elif key == 'AC':
                ac_idx = i
        
        # Parse each sample genotype.
        af_values = []
        ac_values = []
        
        for sample_field in sample_fields:
            if sample_field == '.':
                continue
            
            sample_values = sample_field.split(':')
            for i, key in enumerate(format_keys):
                if i < len(sample_values):
                    self.parser._validate_target_format_value(key, sample_values[i], line_number, pop_id)
            
            # 1. Calculate directly from GT when available.
            if gt_idx is not None and gt_idx < len(sample_values):
                gt_value = sample_values[gt_idx]
                if gt_value != '.' and gt_value != './.':
                    # GT examples: 0/0, 0/1, 1/1, 0|1.
                    # Split phased and unphased genotype separators.
                    alleles = re.split(r'[/|]', gt_value)
                    for allele in alleles:
                        if allele != '.':
                            try:
                                allele_num = int(allele)
                                total_alleles += 1
                                if allele_num > 0:  # 0 is REF, values above 0 are ALT.
                                    alt_allele_count += 1
                            except ValueError:
                                pass
            
            # 2. Collect FORMAT AF values as a fallback.
            if af_idx is not None and af_idx < len(sample_values):
                af_value = sample_values[af_idx]
                if af_value != '.':
                    try:
                        af = float(af_value)
                        af_values.append(af)
                    except ValueError:
                        pass
            
            # 3. Collect FORMAT AC values as a fallback.
            if ac_idx is not None and ac_idx < len(sample_values):
                ac_value = sample_values[ac_idx]
                if ac_value != '.':
                    try:
                        if ',' in ac_value:
                            ac = sum(int(v) for v in ac_value.split(','))
                        else:
                            ac = int(ac_value)
                        ac_values.append(ac)
                    except ValueError:
                        pass
        
        # Calculate FRQ.
        frq = None
        na = total_samples * 2 if total_samples > 0 else 0
        
        if total_alleles > 0:
            # Direct GT-based calculation.
            frq = alt_allele_count / total_alleles
            na = total_alleles
        else:
            # Use collected AF or AC values when GT is unavailable.
            if af_values:
                # AF is already a frequency, so use the mean.
                frq = sum(af_values) / len(af_values)
            elif ac_values and na > 0:
                # Calculate FRQ from AC.
                total_ac = sum(ac_values)
                frq = total_ac / na
        
        # Store population data.
        if frq is not None:
            population_data[pop_id] = {
                'NA': na,
                'FRQ': frq
            }
        elif total_samples > 0:
            # Use a fallback value when FRQ cannot be calculated.
            population_data[pop_id] = {
                'NA': na,
                'FRQ': 0.0
            }
        
        return format_data, population_data
    
    def _calculate_vrt(self, ref: str, alt: str) -> Optional[int]:
        """Calculate a VRT value from REF and ALT."""
        if not ref or not alt or alt == '.':
            return None
        
        ref_len = len(ref)
        
        # Use the first ALT allele if multiple alleles are present.
        alt_alleles = alt.split(',')
        first_alt = alt_alleles[0].strip()
        
        # Handle symbolic or named alleles such as <ID> or [Alu].
        if first_alt.startswith('<') or first_alt.startswith('['):
            return 5  # NAMED
        
        alt_len = len(first_alt)
        
        # SNV: REF and ALT are both one base and differ.
        if ref_len == 1 and alt_len == 1 and ref != first_alt:
            return 1
        
        # Insertion: ALT is longer than REF and starts with REF.
        if alt_len > ref_len and first_alt.startswith(ref):
            return 2
        
        # Deletion: REF is longer than ALT and starts with ALT.
        if ref_len > alt_len and ref.startswith(first_alt):
            return 2
        
        # MNV: REF and ALT have the same length greater than 1.
        if ref_len == alt_len and ref_len > 1:
            return 8
        
        # Indel: REF and ALT lengths differ.
        if ref_len != alt_len:
            return 2
        
        # Fallback.
        return 1  # Default to SNV.
    
    def _validate_and_fix_metadata(self) -> None:
        """Validate metadata and apply trusted metadata values to output headers."""
        if not self.metadata_info:
            return

        metadata_experiment_id = self.metadata_info.experiment_id
        vcf_batch = self.parser.header.metadata.batch
        if metadata_experiment_id:
            if not vcf_batch or not vcf_batch.strip():
                self.parser.header.metadata.batch = metadata_experiment_id
                print(f"[Fixed] VCF had no batch; set to metadata experiment_id ({metadata_experiment_id}).")
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_FILLED,
                    field_name="batch",
                    expected_value=metadata_experiment_id,
                    actual_value="Not in VCF file",
                    additional_info={
                        "action": "Updated from metadata file",
                        "source": "metadata file"
                    }
                )
            elif vcf_batch.strip() != metadata_experiment_id:
                old_value = vcf_batch
                self.parser.header.metadata.batch = metadata_experiment_id
                print(f"[Fixed] VCF batch ({old_value}) differed from metadata experiment_id ({metadata_experiment_id}); updated from metadata.")
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_CORRECTED,
                    field_name="batch/Experiment_id",
                    expected_value=metadata_experiment_id,
                    actual_value=old_value,
                    additional_info={
                        "action": "Updated from metadata file",
                        "source": "metadata file"
                    }
                )
        
        vcf_bioproject_id = self.parser.header.metadata.bioproject_id
        metadata_bioproject_id = self.metadata_info.bioproject_id
        
        if metadata_bioproject_id:
            if not vcf_bioproject_id or not vcf_bioproject_id.strip():
                self.parser.header.metadata.bioproject_id = metadata_bioproject_id
                print(f"[Fixed] VCF had no bioproject_id; set to metadata value ({metadata_bioproject_id}).")
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_FILLED,
                    field_name="bioproject_id",
                    expected_value=metadata_bioproject_id,
                    actual_value="Not in VCF file",
                    additional_info={
                        "action": "Updated from metadata file",
                        "source": "metadata file"
                    }
                )
            elif vcf_bioproject_id.strip() != metadata_bioproject_id:
                old_value = vcf_bioproject_id
                self.parser.header.metadata.bioproject_id = metadata_bioproject_id
                print(f"[Fixed] VCF bioproject_id ({old_value}) differed from metadata ({metadata_bioproject_id}); updated from metadata.")
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_CORRECTED,
                    field_name="bioproject_id",
                    expected_value=metadata_bioproject_id,
                    actual_value=old_value,
                    additional_info={
                        "action": "Updated from metadata file",
                        "source": "metadata file"
                    }
                )
        
        vcf_reference = self.parser.header.metadata.reference
        metadata_reference = self.metadata_info.reference
        vcf_ref_stripped = (vcf_reference or "").strip()
        need_fill = metadata_reference and (vcf_reference is None or not vcf_ref_stripped or vcf_ref_stripped != metadata_reference.strip())
        if need_fill:
            self.parser.header.metadata.reference = metadata_reference
            vcf_display = vcf_reference if (vcf_reference and vcf_ref_stripped) else "(none)"
            if vcf_ref_stripped:
                self.error_handler.create_error(
                    ErrorCode.METADATA_REFERENCE_MISMATCH,
                    field_name="reference",
                    expected_value=metadata_reference,
                    actual_value=vcf_reference or "",
                    additional_info={"source": "metadata file"}
                )
            else:
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_FILLED,
                    field_name="reference",
                    expected_value=metadata_reference,
                    actual_value="Not in VCF file",
                    additional_info={
                        "action": "Updated from metadata file",
                        "source": "metadata file"
                    }
                )
            print(f"[Fixed] Reference: VCF had '{vcf_display}', replaced with metadata value '{metadata_reference}' in output.")
            self.error_handler.create_error(
                ErrorCode.REFERENCE_CORRECTED,
                field_name="reference",
                expected_value=metadata_reference,
                actual_value=vcf_reference or "",
                additional_info={
                    "vcf_reference": vcf_reference or "",
                    "metadata_reference": metadata_reference,
                    "action": "Output VCF uses metadata reference"
                }
            )
    
    def _write_dbsnp_vcf(self, output_file_path: str) -> None:
        """Write a dbSNP-formatted VCF file with optional gzip compression."""
        # Use gzip compression when the output path ends with .gz.
        if output_file_path.endswith('.gz'):
            with gzip.open(output_file_path, 'wt', encoding='utf-8') as f:
                # Write metadata.
                self._write_metadata(f)
                
                # Write INFO tag definitions.
                self._write_info_tag_definitions(f)
                
                # Write FORMAT tag definitions.
                self._write_format_tag_definitions(f)
                
                # Write population IDs.
                self._write_population_ids(f)
                
                # Write column header.
                self._write_column_header(f)
                
                # Write data rows.
                for row in self.parser.data_rows:
                    self._write_data_row(f, row)
        else:
            with open(output_file_path, 'w', encoding='utf-8') as f:
                # Write metadata.
                self._write_metadata(f)
                
                # Write INFO tag definitions.
                self._write_info_tag_definitions(f)
                
                # Write FORMAT tag definitions.
                self._write_format_tag_definitions(f)
                
                # Write population IDs.
                self._write_population_ids(f)
                
                # Write column header.
                self._write_column_header(f)
                
                # Write data rows.
                for row in self.parser.data_rows:
                    self._write_data_row(f, row)
    
    def _write_metadata(self, f) -> None:
        """Write metadata."""
        # fileformat
        if self.parser.header.metadata.fileformat:
            f.write(f"##fileformat={self.parser.header.metadata.fileformat}\n")
        else:
            f.write("##fileformat=VCFv4.1\n")
        
        # fileDate
        if self.parser.header.metadata.filedate:
            f.write(f"##fileDate={self.parser.header.metadata.filedate}\n")
        else:
            f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        
        # handle is always set to KVar.
        f.write(f"##handle=KVar\n")
        
        # batch comes from metadata when available.
        batch = self.parser.header.metadata.batch
        if not batch and self.metadata_info:
            batch = self.metadata_info.experiment_id
        if batch:
            f.write(f"##batch={batch}\n")
        
        # bioproject_id
        bioproject_id = self.parser.header.metadata.bioproject_id
        if not bioproject_id and self.metadata_info:
            bioproject_id = self.metadata_info.bioproject_id
        if bioproject_id:
            f.write(f"##bioproject_id={bioproject_id}\n")
        
        # biosample_id
        if self.parser.header.metadata.biosample_id:
            f.write(f"##biosample_id={self.parser.header.metadata.biosample_id}\n")
        elif self.metadata_info and self.metadata_info.biosample_id:
            biosample_id = ','.join(self.metadata_info.biosample_id)
            f.write(f"##biosample_id={biosample_id}\n")
        
        # reference
        reference = self.parser.header.metadata.reference
        if not reference and self.metadata_info:
            reference = self.metadata_info.reference
        if reference:
            f.write(f"##reference={reference}\n")
        
        for contig_line in self.contig_lines:
            f.write(contig_line + '\n')
    
    def _write_info_tag_definitions(self, f) -> None:
        """Write INFO tag definitions."""
        f.write('##INFO=<ID=VRT,Number=1,Type=Integer,Description="Variation type, 1 - SNV: single nucleotide variation, 2 - DIV: deletion/insertion variation, 3 - HETEROZYGOUS: variable, but undefined at nucleotide level, 4 - STR: short tandem repeat (microsatellite) variation, 5 - NAMED: insertion/deletion variation of named repetitive element, 6 - NO VARIATION: sequence scanned for variation, but none observed, 7 - MIXED: cluster contains submissions from 2 or more allelic classes, 8 - MNV: multiple nucleotide variation with alleles of common length greater than 1">\n')
        
        written_info_ids = {'VRT'}
        for tag_id, tag_def in self.parser.header.info_tags.items():
            output_info_id = get_vcf_info_output_id(tag_def.id)
            if output_info_id in written_info_ids:
                continue
            written_info_ids.add(output_info_id)
            f.write(f'##INFO=<ID={output_info_id},Number={tag_def.number},Type={tag_def.type},Description="{tag_def.description}">\n')
    
    def _write_format_tag_definitions(self, f) -> None:
        """Write FORMAT tag definitions."""
        # Standard dbSNP FORMAT tags.
        standard_formats = {
            'NA': 'Number of alleles for the population',
            'NS': 'Number of samples for the population',
            'FRQ': 'Frequency of each alternate allele',
            'AC': 'Allele count for each alternate allele'
        }
        
        for format_id, description in standard_formats.items():
            if format_id not in self.parser.header.format_tags:
                if format_id == 'FRQ':
                    f.write(f'##FORMAT=<ID={format_id},Number=.,Type=Float,Description="{description}">\n')
                elif format_id == 'AC':
                    f.write(f'##FORMAT=<ID={format_id},Number=.,Type=Integer,Description="{description}">\n')
                else:
                    f.write(f'##FORMAT=<ID={format_id},Number=1,Type=Integer,Description="{description}">\n')
        
        # Existing FORMAT tags.
        for tag_id, tag_def in self.parser.header.format_tags.items():
            if tag_id not in standard_formats:
                f.write(f'##FORMAT=<ID={tag_def.id},Number={tag_def.number},Type={tag_def.type},Description="{tag_def.description}">\n')
    
    def _write_population_ids(self, f) -> None:
        """Write population IDs."""
        for population_id in self._output_population_ids():
            f.write(f"##population_id={population_id}\n")
    
    def _output_population_ids(self) -> List[str]:
        """Return population IDs for dbSNP VCF output."""
        if self.metadata_info and self.metadata_info.sampleset_ids:
            return self.metadata_info.sampleset_ids
        if self.parser.header.population_ids:
            return self.parser.header.population_ids
        return ['DEFAULT_POPULATION']
    
    def _write_column_header(self, f) -> None:
        """Write the column header."""
        # Basic columns.
        columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        
        # Add FORMAT column.
        population_ids = self._output_population_ids()
        if population_ids:
            columns.append('FORMAT')
            columns.extend(population_ids)
        
        f.write('\t'.join(columns) + '\n')
    
    def _write_data_row(self, f, row: VCFDataRow) -> None:
        """Write a data row."""
        # Basic fields.
        chrom = row.chrom
        pos = str(row.pos)
        id_field = row.id if row.id != '.' else '.'
        ref = row.ref
        alt = row.alt
        qual = row.qual if row.qual else '.'
        filter_field = row.filter if row.filter else '.'
        
        # Write INFO field.
        info_parts = []
        for key, value in row.info.items():
            output_key = get_vcf_info_output_id(key)
            if isinstance(value, bool) and value:
                info_parts.append(output_key)
            else:
                if isinstance(value, list):
                    value_str = ','.join(str(v) for v in value)
                else:
                    value_str = str(value)
                info_parts.append(f"{output_key}={value_str}")
        
        info_str = ';'.join(info_parts) if info_parts else '.'
        
        # Write basic fields.
        fields = [chrom, pos, id_field, ref, alt, qual, filter_field, info_str]
        
        # FORMAT and population data.
        # Add FORMAT whenever INFO contains AF or AC.
        has_af_or_ac = 'AF' in row.info or 'AC' in row.info
        population_ids = self._output_population_ids()
        has_population_ids = bool(population_ids)
        
        if row.population_data or has_af_or_ac or has_population_ids:
            # FORMAT.
            format_str = 'NA:FRQ'
            fields.append(format_str)
            
            for pop_id in population_ids:
                if pop_id in row.population_data:
                    # Use parsed population_data when available.
                    pop_data = row.population_data[pop_id]
                    na = pop_data.get('NA', 0)
                    frq = pop_data.get('FRQ', 0.0)
                    fields.append(f"{na}:{frq}")
                elif has_af_or_ac:
                    # Calculate directly from INFO AF or AC when available.
                    na = 0
                    frq = 0.0
                    
                    # Prefer AN and AF when both are present.
                    if 'AN' in row.info and 'AF' in row.info:
                        try:
                            an_value = row.info['AN']
                            na = int(an_value) if not isinstance(an_value, int) else an_value
                        except (ValueError, TypeError):
                            na = 0
                        
                        try:
                            af_value = row.info['AF']
                            if isinstance(af_value, list):
                                frq = float(af_value[0]) if af_value else 0.0
                            else:
                                frq = float(af_value)
                        except (ValueError, TypeError):
                            frq = 0.0
                    
                    # Calculate FRQ from AN and AC.
                    elif 'AN' in row.info and 'AC' in row.info:
                        try:
                            an_value = row.info['AN']
                            na = int(an_value) if not isinstance(an_value, int) else an_value
                        except (ValueError, TypeError):
                            na = 0
                        
                        try:
                            ac_value = row.info['AC']
                            if isinstance(ac_value, list):
                                total_ac = sum(int(v) if not isinstance(v, int) else v for v in ac_value)
                            else:
                                total_ac = int(ac_value) if not isinstance(ac_value, int) else ac_value
                            if na > 0:
                                frq = total_ac / na
                        except (ValueError, TypeError, ZeroDivisionError):
                            frq = 0.0
                    
                    # AF alone cannot estimate AN.
                    elif 'AF' in row.info:
                        try:
                            af_value = row.info['AF']
                            if isinstance(af_value, list):
                                frq = float(af_value[0]) if af_value else 0.0
                            else:
                                frq = float(af_value)
                        except (ValueError, TypeError):
                            frq = 0.0
                        # Set AN to 0 when it is missing.
                        na = 0
                    
                    if na > 0 or frq > 0:
                        fields.append(f"{na}:{frq}")
                    else:
                        fields.append('.')
                else:
                    fields.append('.')
        
        f.write('\t'.join(fields) + '\n')


def main():
    """Main function - command-line argument handling"""
    parser = argparse.ArgumentParser(
        description='Convert generic VCF to dbSNP VCF format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python VCF2dbSNP.py -v input.vcf -o output.vcf -m metadata.txt
  python VCF2dbSNP.py -v input.vcf.gz -o output.vcf.gz -m metadata.txt
  python VCF2dbSNP.py -v input.vcf.gz -o output.dbsnp.vcf.gz
        """
    )
    
    parser.add_argument(
        '-v', '--vcf',
        dest='vcf_file',
        required=True,
        help='Input VCF file path (required, .gz supported)'
    )
    
    parser.add_argument(
        '-o', '--output',
        dest='output_file',
        required=True,
        help='Output dbSNP VCF file path (required; use .gz suffix for gzip)'
    )
    
    parser.add_argument(
        '-m', '--metadata',
        dest='metadata_file',
        required=False,
        help='Metadata file path (optional)'
    )
    
    parser.add_argument(
        '-e', '--error-report',
        dest='error_report',
        required=False,
        help='Error report file path (optional)'
    )
    
    args = parser.parse_args()
    
    # Check file existence
    if not os.path.exists(args.vcf_file):
        print(f"Error: VCF file not found: {args.vcf_file}")
        sys.exit(1)
    
    if args.metadata_file and not os.path.exists(args.metadata_file):
        print(f"Error: Metadata file not found: {args.metadata_file}")
        sys.exit(1)
    
    # Run conversion
    converter = VCF2dbSNPConverter()
    
    try:
        converter.convert_vcf_to_dbsnp(
            args.vcf_file,
            args.output_file,
            args.metadata_file,
            args.error_report
        )
        
        print("\n=== Conversion complete ===")
        
    except Exception as e:
        print(f"\nError during conversion: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
