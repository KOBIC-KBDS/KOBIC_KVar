#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dbSNP VCF parser.
Parses VCF files conforming to dbSNP VCF submission format and extracts all information.
"""

import re
import os
import gzip
from typing import Callable, Dict, Iterator, List, Optional, Any, Tuple, Set
from dataclasses import dataclass, field
from collections import defaultdict

try:
    from .error_handler import ErrorHandler, ErrorCode, ErrorSeverity
except ImportError:
    from error_handler import ErrorHandler, ErrorCode, ErrorSeverity


@dataclass
class VCFHeaderMetadata:
    """Class to store VCF header metadata"""
    fileformat: Optional[str] = None
    filedate: Optional[str] = None
    handle: Optional[str] = None
    batch: Optional[str] = None
    bioproject_id: Optional[str] = None
    biosample_id: Optional[str] = None
    sampleset_id: Optional[str] = None  # Alternative to biosample_id
    reference: Optional[str] = None


@dataclass
class InfoTagDefinition:
    """Class to store INFO tag definition"""
    id: str
    number: str
    type: str
    description: str


@dataclass
class FormatTagDefinition:
    """Class to store FORMAT tag definition"""
    id: str
    number: str
    type: str
    description: str


@dataclass
class VCFHeader:
    """Class to store VCF header information"""
    metadata: VCFHeaderMetadata = field(default_factory=VCFHeaderMetadata)
    info_tags: Dict[str, InfoTagDefinition] = field(default_factory=dict)
    format_tags: Dict[str, FormatTagDefinition] = field(default_factory=dict)
    population_ids: List[str] = field(default_factory=list)
    column_header: List[str] = field(default_factory=list)


@dataclass
class VCFDataRow:
    """Class to store VCF data row"""
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: Optional[str]
    filter: Optional[str]
    info: Dict[str, Any] = field(default_factory=dict)
    format_data: Dict[str, Any] = field(default_factory=dict)
    population_data: Dict[str, Dict[str, Any]] = field(default_factory=dict)


class dbSNPVCFParser:
    """Main class to parse dbSNP VCF files"""

    def __init__(self, error_handler: Optional[ErrorHandler] = None):
        self.header = VCFHeader()
        self.data_rows: List[VCFDataRow] = []
        self.error_handler = error_handler or ErrorHandler()
        self._target_info_patterns = {
            'AN': r'[0-9.]+',
            'AC': r'[0-9.]+',
            'AF': r'[0-9]*\.?[0-9]*',
            'DESC': r'[^;=]+',
            'LINKS': r'[A-Za-z]+:[-A-Za-z0-9]+',
        }
        self._target_format_patterns = {
            'AN': r'[0-9.]+',
            'AC': r'[0-9.]+',
            'AF': r'[0-9]*\.?[0-9]*',
            'GT': r'[^;=]+',
            'GL': r'[^;=]+',
            'PL': r'[^;=]+',
            'GP': r'[^;=]+',
            'PP': r'[^;=]+',
        }
        self._singleton_metadata_keys = {
            'fileformat',
            'filedate',
            'handle',
            'batch',
            'bioproject_id',
            'biosample_id',
            'sampleset_id',
            'reference',
        }
        self._seen_metadata_keys: Set[str] = set()
        self._warned_undefined_info_tags: Set[str] = set()
        self._warned_undefined_format_tags: Set[str] = set()
        self._seen_local_ids: Set[str] = set()
        self._seen_variant_sites: Set[Tuple[str, int, str, str, int]] = set()
        self._current_chromosome: Optional[str] = None
        self._closed_chromosomes: Set[str] = set()
        self._previous_position: Optional[int] = None
        self._density_chromosome: Optional[str] = None
        self._density_window_start: Optional[int] = None
        self._density_window_count = 0
        self._parsed_row_count = 0
        self._has_population_data = False

    def parse_header(self, file_path: str) -> None:
        """Parse only VCF header lines."""
        self._reset_parse_state()
        self._ensure_file_exists(file_path)

        header_lines = []
        try:
            with self._open_text(file_path) as f:
                for _line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    if not line.startswith('#'):
                        break
                    header_lines.append(line)
        except UnicodeDecodeError as e:
            self.error_handler.create_error(
                ErrorCode.FILE_ENCODING_ERROR,
                additional_info={"file_path": file_path, "error": str(e)}
            )
            raise
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FILE_READ_ERROR,
                additional_info={"file_path": file_path, "error": str(e)}
            )
            raise

        self._parse_header(header_lines)

    def iter_data_rows(
        self,
        file_path: str,
        store_rows: bool = False,
        row_callback: Optional[Callable[[VCFDataRow], None]] = None,
    ) -> Iterator[VCFDataRow]:
        """Yield parsed data rows without retaining the whole VCF in memory."""
        self._reset_row_state(clear_rows=store_rows)
        self._ensure_file_exists(file_path)

        try:
            with self._open_text(file_path) as f:
                for line_num, line in enumerate(f, 1):
                    original_line = line
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    row = self._parse_data_row(original_line.strip(), line_num)
                    if not row:
                        continue

                    self._record_parsed_row(row)
                    if store_rows:
                        self.data_rows.append(row)
                    if row_callback:
                        row_callback(row)
                    yield row
        except UnicodeDecodeError as e:
            self.error_handler.create_error(
                ErrorCode.FILE_ENCODING_ERROR,
                additional_info={"file_path": file_path, "error": str(e)}
            )
            raise
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FILE_READ_ERROR,
                additional_info={"file_path": file_path, "error": str(e)}
            )
            raise

    def validate_parsed_data(self) -> None:
        """Run validations that depend on all parsed data rows."""
        self._validate_parsed_data()

    def _reset_parse_state(self) -> None:
        """Reset parsed header, retained rows, and streaming validation state."""
        self.header = VCFHeader()
        self.data_rows.clear()
        self._seen_metadata_keys.clear()
        self._warned_undefined_info_tags.clear()
        self._warned_undefined_format_tags.clear()
        self._reset_row_state(clear_rows=False)

    def _reset_row_state(self, clear_rows: bool) -> None:
        """Reset row-level validation state before a data pass."""
        if clear_rows:
            self.data_rows.clear()
        self._seen_local_ids.clear()
        self._seen_variant_sites.clear()
        self._current_chromosome = None
        self._closed_chromosomes.clear()
        self._previous_position = None
        self._density_chromosome = None
        self._density_window_start = None
        self._density_window_count = 0
        self._parsed_row_count = 0
        self._has_population_data = False

    @staticmethod
    def _open_text(file_path: str):
        """Open a plain or gzipped VCF path for text reading."""
        if file_path.endswith('.gz'):
            return gzip.open(file_path, 'rt', encoding='utf-8')
        return open(file_path, 'r', encoding='utf-8')

    def _ensure_file_exists(self, file_path: str) -> None:
        """Record and raise for missing input files."""
        if os.path.exists(file_path):
            return
        self.error_handler.create_error(
            ErrorCode.FILE_NOT_FOUND,
            additional_info={"file_path": file_path}
        )
        raise FileNotFoundError(f"File not found: {file_path}")

    def _parse_header(self, header_lines: List[str]) -> None:
        """Parse VCF header"""
        for line_num, line in enumerate(header_lines, 1):
            if line.startswith('##'):
                self._parse_metadata_line(line, line_num)
            elif line.startswith('#'):
                self._parse_column_header(line, line_num)

        # Required header validation
        if not self.header.metadata.fileformat:
            self.error_handler.create_error(
                ErrorCode.MISSING_FILEFORMAT,
                line_number=1
            )
        else:
            # Check VCF format (must start with VCFv); all VCF versions allowed
            if not self.header.metadata.fileformat.startswith('VCFv'):
                self.error_handler.create_error(
                    ErrorCode.INVALID_FILEFORMAT,
                    line_number=1,
                    actual_value=self.header.metadata.fileformat,
                    expected_value="VCFv format (e.g. VCFv4.1, VCFv4.2, VCFv4.3, VCFv4.4, VCFv4.5)"
                )

        if not self.header.column_header:
            self.error_handler.create_error(
                ErrorCode.MISSING_COLUMN_HEADER
            )

    def _parse_metadata_line(self, line: str, line_number: int) -> None:
        """Parse metadata line"""
        # Parse ##key=value format
        if '=' not in line:
            self.error_handler.create_error(
                ErrorCode.INVALID_METADATA_FORMAT,
                line_number=line_number,
                line_content=line
            )
            return

        try:
            key, value = line[2:].split('=', 1)
            key = key.lower()

            if key == 'fileformat':
                if not self._validate_singleton_metadata_tag(key, line, line_number):
                    return
                self.header.metadata.fileformat = value
            elif key == 'filedate':
                if not self._validate_singleton_metadata_tag(key, line, line_number):
                    return
                self.header.metadata.filedate = value
            elif key == 'handle':
                if not self._validate_singleton_metadata_tag(key, line, line_number):
                    return
                self.header.metadata.handle = value
            elif key == 'batch':
                if not self._validate_singleton_metadata_tag(key, line, line_number):
                    return
                self.header.metadata.batch = value
            elif key == 'bioproject_id':
                if not self._validate_singleton_metadata_tag(key, line, line_number):
                    return
                self.header.metadata.bioproject_id = value
            elif key == 'biosample_id':
                if not self._validate_singleton_metadata_tag(key, line, line_number):
                    return
                # Also store as sampleset_id (compatibility)
                self.header.metadata.biosample_id = value
                self.header.metadata.sampleset_id = value
            elif key == 'sampleset_id':
                if not self._validate_singleton_metadata_tag(key, line, line_number):
                    return
                # Parse SampleSet_id directly
                self.header.metadata.sampleset_id = value
            elif key == 'reference':
                if not self._validate_singleton_metadata_tag(key, line, line_number):
                    return
                self.header.metadata.reference = value
            elif key == 'population_id':
                self.header.population_ids.append(value)
            elif key == 'info':
                self._parse_info_tag_definition(line, line_number)
            elif key == 'format':
                self._parse_format_tag_definition(line, line_number)
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.INVALID_METADATA_FORMAT,
                line_number=line_number,
                line_content=line,
                additional_info={"error": str(e)}
            )

    def _validate_singleton_metadata_tag(self, key: str, line: str, line_number: int) -> bool:
        """Validate metadata tags that may appear only once."""
        if key not in self._singleton_metadata_keys:
            return True

        if key in self._seen_metadata_keys:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_METADATA_TAG,
                line_number=line_number,
                line_content=line,
                field_name=key,
                expected_value="single metadata tag definition",
                actual_value=f"duplicated ##{key}"
            )
            return False

        self._seen_metadata_keys.add(key)
        return True

    def _parse_info_tag_definition(self, line: str, line_number: int) -> None:
        """Parse an INFO tag definition header line."""
        pattern = r'##INFO=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="((?:\\.|[^"\\])*)"'
        match = re.search(pattern, line)

        if match:
            tag_id, number, type_str, description = match.groups()
            self._store_info_tag_definition(tag_id, number, type_str, description, line, line_number)
        else:
            parts = line.split(',')
            if len(parts) >= 4:
                try:
                    id_part = parts[0].replace('##INFO=<ID=', '').strip()
                    number_part = parts[1].replace('Number=', '').strip()
                    type_part = parts[2].replace('Type=', '').strip()
                    desc_key = 'Description="'
                    desc_start = line.find(desc_key)
                    desc_end = line.rfind('"')
                    if desc_start != -1 and desc_end > desc_start:
                        desc_start += len(desc_key)
                        description = line[desc_start:desc_end]

                        self._store_info_tag_definition(
                            id_part,
                            number_part,
                            type_part,
                            description,
                            line,
                            line_number
                        )
                    else:
                        self.error_handler.create_error(
                            ErrorCode.INFO_TAG_PARSE_ERROR,
                            line_number=line_number,
                            line_content=line,
                            field_name="Description"
                        )
                except Exception as e:
                    self.error_handler.create_error(
                        ErrorCode.INFO_TAG_PARSE_ERROR,
                        line_number=line_number,
                        line_content=line,
                        additional_info={"error": str(e)}
                    )
            else:
                self.error_handler.create_error(
                    ErrorCode.INFO_TAG_PARSE_ERROR,
                    line_number=line_number,
                    line_content=line
                )

    def _store_info_tag_definition(
        self,
        tag_id: str,
        number: str,
        type_str: str,
        description: str,
        line: str,
        line_number: int
    ) -> None:
        """Store an INFO tag definition after duplicate and VRT checks."""
        tag_id = tag_id.strip()
        number = number.strip()
        type_str = type_str.strip()
        description = self._unescape_vcf_description(description.strip())

        if tag_id in self.header.info_tags:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_INFO_TAG_DEFINITION,
                line_number=line_number,
                line_content=line,
                field_name=tag_id,
                expected_value="unique INFO tag ID",
                actual_value=tag_id
            )
            return

        self.header.info_tags[tag_id] = InfoTagDefinition(
            id=tag_id,
            number=number,
            type=type_str,
            description=description
        )

        if tag_id == 'VRT':
            self._validate_vrt_header_definition(number, type_str, description, line, line_number)

    def _validate_vrt_header_definition(
        self,
        number: str,
        type_str: str,
        description: str,
        line: str,
        line_number: int
    ) -> None:
        """Validate the SNP VRT INFO header definition against the canonical map."""
        expected_vrt_map = {
            1: 'SNV',
            2: 'DIV',
            3: 'HETEROZYGOUS',
            4: 'STR',
            5: 'NAMED',
            6: 'NO VARIATION',
            7: 'MIXED',
            8: 'MNV',
        }
        parsed_vrt_map = {
            int(vrt_number): vrt_string.strip()
            for vrt_number, vrt_string in re.findall(r'([1-8])\s*-\s*([A-Z ]+)', description)
        }

        issues = []
        if number != '1':
            issues.append('Number must be 1')
        if type_str != 'Integer':
            issues.append('Type must be Integer')
        if parsed_vrt_map != expected_vrt_map:
            issues.append('VRT map must define 1-8 as SNV, DIV, HETEROZYGOUS, STR, NAMED, NO VARIATION, MIXED, and MNV')

        if issues:
            self.error_handler.create_error(
                ErrorCode.INVALID_VRT_HEADER_DEFINITION,
                line_number=line_number,
                line_content=line,
                field_name='VRT',
                expected_value='Number=1; Type=Integer; VRT map 1-8',
                actual_value=f"Number={number}; Type={type_str}; VRT map={parsed_vrt_map}",
                additional_info={"issues": issues}
            )

    def _parse_format_tag_definition(self, line: str, line_number: int) -> None:
        """Parse a FORMAT tag definition header line."""
        pattern = r'##FORMAT=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="((?:\\.|[^"\\])*)"'
        match = re.search(pattern, line)

        if match:
            tag_id, number, type_str, description = match.groups()
            self._store_format_tag_definition(tag_id, number, type_str, description, line, line_number)
        else:
            parts = line.split(',')
            if len(parts) >= 4:
                try:
                    id_part = parts[0].replace('##FORMAT=<ID=', '').strip()
                    number_part = parts[1].replace('Number=', '').strip()
                    type_part = parts[2].replace('Type=', '').strip()
                    desc_key = 'Description="'
                    desc_start = line.find(desc_key)
                    desc_end = line.rfind('"')
                    if desc_start != -1 and desc_end > desc_start:
                        desc_start += len(desc_key)
                        description = line[desc_start:desc_end]

                        self._store_format_tag_definition(
                            id_part,
                            number_part,
                            type_part,
                            description,
                            line,
                            line_number
                        )
                    else:
                        self.error_handler.create_error(
                            ErrorCode.FORMAT_TAG_PARSE_ERROR,
                            line_number=line_number,
                            line_content=line,
                            field_name="Description"
                        )
                except Exception as e:
                    self.error_handler.create_error(
                        ErrorCode.FORMAT_TAG_PARSE_ERROR,
                        line_number=line_number,
                        line_content=line,
                        additional_info={"error": str(e)}
                    )
            else:
                self.error_handler.create_error(
                    ErrorCode.FORMAT_TAG_PARSE_ERROR,
                    line_number=line_number,
                    line_content=line
                )

    def _store_format_tag_definition(
        self,
        tag_id: str,
        number: str,
        type_str: str,
        description: str,
        line: str,
        line_number: int
    ) -> None:
        """Store a FORMAT tag definition after duplicate checks."""
        tag_id = tag_id.strip()
        number = number.strip()
        type_str = type_str.strip()
        description = self._unescape_vcf_description(description.strip())

        if tag_id in self.header.format_tags:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_FORMAT_TAG_DEFINITION,
                line_number=line_number,
                line_content=line,
                field_name=tag_id,
                expected_value="unique FORMAT tag ID",
                actual_value=tag_id
            )
            return

        self.header.format_tags[tag_id] = FormatTagDefinition(
            id=tag_id,
            number=number,
            type=type_str,
            description=description
        )

    @staticmethod
    def _unescape_vcf_description(description: str) -> str:
        """Unescape quoted characters from a VCF header description."""
        return re.sub(r'\\(["\\])', r'\1', description)

    def _parse_column_header(self, line: str, line_number: int) -> None:
        """Parse the VCF column header."""
        if self.header.column_header:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_COLUMN_HEADER,
                line_number=line_number,
                line_content=line,
                field_name="#CHROM",
                expected_value="single column header line",
                actual_value="duplicated column header line"
            )
            return

        # #CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT] [POPULATION_COLS]
        columns = line[1:].split('\t')

        # Check the minimum required columns.
        required_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        if len(columns) < len(required_columns):
            self.error_handler.create_error(
                ErrorCode.INVALID_COLUMN_HEADER,
                line_number=line_number,
                line_content=line,
                expected_value=f"at least {len(required_columns)} columns",
                actual_value=f"{len(columns)} columns"
            )

        self.header.column_header = columns

    def _parse_data(self, data_lines: List[Tuple[int, str]]) -> None:
        """Parse VCF data lines."""
        for line_number, line in data_lines:
            row = self._parse_data_row(line, line_number)
            if row:
                self._record_parsed_row(row)
                self.data_rows.append(row)

    def _record_parsed_row(self, row: VCFDataRow) -> None:
        """Track aggregate facts needed after a streaming data pass."""
        self._parsed_row_count += 1
        if row.population_data:
            self._has_population_data = True

    def _parse_data_row(self, line: str, line_number: int) -> Optional[VCFDataRow]:
        """Parse one VCF data row."""
        fields = line.split('\t')

        if len(fields) < 8:
            self.error_handler.create_error(
                ErrorCode.INSUFFICIENT_FIELDS,
                line_number=line_number,
                line_content=line,
                expected_value="at least 8",
                actual_value=f"{len(fields)}"
            )
            return None

        # Basic fields
        chrom = fields[0]

        # Position validation
        try:
            pos = int(fields[1])
            if pos < 1:
                self.error_handler.create_error(
                    ErrorCode.INVALID_POSITION,
                    line_number=line_number,
                    line_content=line,
                    field_name="POS",
                    actual_value=str(pos),
                    expected_value="integer greater than or equal to 1"
                )
                return None
        except ValueError:
            self.error_handler.create_error(
                ErrorCode.INVALID_POSITION,
                line_number=line_number,
                line_content=line,
                field_name="POS",
                actual_value=fields[1]
            )
            return None

        if not self._validate_chromosome_order(chrom, pos, line, line_number):
            return None
        self._check_snp_density(chrom, pos, line, line_number)

        id_field = fields[2]
        ref = fields[3]
        alt = fields[4]

        if not self._validate_local_id(id_field, line, line_number):
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
        # Multi-allelic sites are split by comma; each allele is validated independently
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

        qual = fields[5] if fields[5] != '.' else None
        filter_field = fields[6] if fields[6] != '.' else None

        # Parse INFO field.
        info_dict = self._parse_info_field(fields[7], line_number, line)

        vrt_value = info_dict.get('VRT')
        if not isinstance(vrt_value, int) or vrt_value < 1 or vrt_value > 8:
            return None

        if not self._validate_ac_an(info_dict, line, line_number):
            return None

        self._derive_af_from_ac_an(info_dict, line, line_number)

        if not self._validate_duplicate_site(chrom, pos, ref, alt, info_dict, line, line_number):
            return None

        # Parse FORMAT and population data.
        format_data = {}
        population_data = {}

        if len(fields) > 8:
            format_data, population_data = self._parse_format_and_population_data(
                fields[8:], line_number, line
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

    def _validate_chromosome_order(self, chrom: str, pos: int, line: str, line_number: int) -> bool:
        """Validate chromosome grouping and position sorting."""
        if self._current_chromosome is None:
            self._current_chromosome = chrom
            self._previous_position = pos
            return True

        if chrom != self._current_chromosome:
            self._closed_chromosomes.add(self._current_chromosome)
            self._current_chromosome = chrom
            self._previous_position = pos

            if chrom in self._closed_chromosomes:
                self.error_handler.create_error(
                    ErrorCode.CHROMOSOME_NOT_GROUPED,
                    line_number=line_number,
                    line_content=line,
                    field_name="CHROM",
                    actual_value=chrom,
                    expected_value="all records for each chromosome grouped together"
                )
                return False
            return True

        if self._previous_position is not None and pos < self._previous_position:
            previous_position = self._previous_position
            self._previous_position = pos
            self.error_handler.create_error(
                ErrorCode.POSITION_NOT_SORTED,
                line_number=line_number,
                line_content=line,
                field_name="POS",
                actual_value=str(pos),
                expected_value=f">= previous position {previous_position}"
            )
            return False

        self._previous_position = pos
        return True

    def _check_snp_density(self, chrom: str, pos: int, line: str, line_number: int) -> None:
        """Record a warning when more than 10 SNP records appear in a 50bp window."""
        if (
            self._density_chromosome != chrom
            or self._density_window_start is None
            or pos - self._density_window_start > 50
        ):
            self._density_chromosome = chrom
            self._density_window_start = pos
            self._density_window_count = 1
            return

        self._density_window_count += 1
        if self._density_window_count > 10:
            self.error_handler.create_error(
                ErrorCode.SNP_DENSITY_TOO_HIGH,
                line_number=line_number,
                line_content=line,
                field_name="POS",
                actual_value=f"{self._density_window_count} records from {chrom}:{self._density_window_start}-{pos}",
                expected_value="10 or fewer records in a 50bp window"
            )

    def _validate_local_id(self, local_id: str, line: str, line_number: int) -> bool:
        """Validate local ID existence, length, and uniqueness."""
        if not local_id or local_id == '.':
            self.error_handler.create_error(
                ErrorCode.MISSING_LOCAL_ID,
                line_number=line_number,
                line_content=line,
                field_name="ID",
                expected_value="non-empty local ID"
            )
            return False

        if len(local_id) > 64:
            self.error_handler.create_error(
                ErrorCode.LOCAL_ID_TOO_LONG,
                line_number=line_number,
                line_content=line,
                field_name="ID",
                actual_value=local_id,
                expected_value="64 characters or fewer"
            )
            return False

        if local_id in self._seen_local_ids:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_LOCAL_ID,
                line_number=line_number,
                line_content=line,
                field_name="ID",
                actual_value=local_id,
                expected_value="unique local ID"
            )
            return False

        self._seen_local_ids.add(local_id)
        return True

    @staticmethod
    def _coerce_int_list(value: Any) -> Optional[List[int]]:
        """Return integer values for scalar or comma-separated INFO values."""
        if isinstance(value, bool):
            return None
        if isinstance(value, int):
            return [value]
        if isinstance(value, list):
            values = value
        elif isinstance(value, str):
            if value in {'', '.'}:
                return None
            values = value.split(',')
        else:
            return None

        int_values = []
        for item in values:
            if isinstance(item, bool):
                return None
            if isinstance(item, int):
                int_values.append(item)
                continue
            if isinstance(item, str) and re.fullmatch(r"-?\d+", item):
                int_values.append(int(item))
                continue
            return None
        return int_values

    def _validate_ac_an(self, info_dict: Dict[str, Any], line: str, line_number: int) -> bool:
        """Validate that every AC value is positive and less than or equal to AN."""
        if 'AC' not in info_dict or 'AN' not in info_dict:
            return True

        ac_values = self._coerce_int_list(info_dict.get('AC'))
        an_values = self._coerce_int_list(info_dict.get('AN'))
        if not ac_values or not an_values or len(an_values) != 1:
            return True

        an_value = an_values[0]
        for ac_value in ac_values:
            if ac_value <= 0:
                self.error_handler.create_error(
                    ErrorCode.INVALID_ALLELE_COUNT,
                    line_number=line_number,
                    line_content=line,
                    field_name="AC",
                    actual_value=f"AC={ac_value}",
                    expected_value="AC > 0"
                )
                return False

            if ac_value > an_value:
                self.error_handler.create_error(
                    ErrorCode.ALLELE_COUNT_EXCEEDS_NUMBER,
                    line_number=line_number,
                    line_content=line,
                    field_name="AC/AN",
                    actual_value=f"AC={ac_value}, AN={an_value}",
                    expected_value="AC <= AN"
                )
                return False
        return True

    def _derive_af_from_ac_an(self, info_dict: Dict[str, Any], line: str, line_number: int) -> None:
        """Calculate AF from AC and AN when AF is not provided."""
        if 'AF' in info_dict or 'AC' not in info_dict or 'AN' not in info_dict:
            return

        ac_values = self._coerce_int_list(info_dict.get('AC'))
        an_values = self._coerce_int_list(info_dict.get('AN'))
        if not ac_values or not an_values or len(an_values) != 1:
            return

        an_value = an_values[0]
        if an_value <= 0:
            return

        af_values = [ac_value / an_value for ac_value in ac_values]
        calculated_af: Any = af_values[0] if len(af_values) == 1 else af_values
        info_dict['AF'] = calculated_af
        self._ensure_info_tag_definition(
            'AF',
            'A',
            'Float',
            'Allele frequency for each alternate allele'
        )

        self.error_handler.create_error(
            ErrorCode.ALLELE_FREQUENCY_CALCULATED,
            line_number=line_number,
            line_content=line,
            field_name="AF",
            actual_value=self._format_info_value(calculated_af),
            expected_value="AF calculated from AC/AN",
            additional_info={
                "AC": ac_values,
                "AN": an_value,
                "calculation": "AF = AC / AN"
            }
        )

    def _ensure_info_tag_definition(self, tag_id: str, number: str, type_str: str, description: str) -> None:
        """Add an INFO tag definition if it is not already present."""
        if tag_id in self.header.info_tags:
            return

        self.header.info_tags[tag_id] = InfoTagDefinition(
            id=tag_id,
            number=number,
            type=type_str,
            description=description
        )

    @staticmethod
    def _format_info_value(value: Any) -> str:
        """Format a scalar or list INFO value for reports."""
        if isinstance(value, list):
            return ','.join(str(item) for item in value)
        return str(value)

    def _validate_duplicate_site(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        info_dict: Dict[str, Any],
        line: str,
        line_number: int,
    ) -> bool:
        """Validate duplicate sites using CHROM, POS, REF, ALT, and VRT."""
        vrt = info_dict.get('VRT')
        if not isinstance(vrt, int) or vrt < 1 or vrt > 8:
            return True

        site_key = (chrom, pos, ref.upper(), alt.upper(), vrt)
        if site_key in self._seen_variant_sites:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_VARIANT_SITE,
                line_number=line_number,
                line_content=line,
                field_name="CHROM/POS/REF/ALT/VRT",
                actual_value=f"{chrom}:{pos}:{ref}:{alt}:{vrt}",
                expected_value="unique variant site"
            )
            return False

        self._seen_variant_sites.add(site_key)
        return True

    def _parse_info_field(self, info_str: str, line_number: int, line: Optional[str] = None) -> Dict[str, Any]:
        """Parse an INFO field."""
        info_dict = {}

        if info_str == '.':
            # VRT is required for SNP submission rows.
            self.error_handler.create_error(
                ErrorCode.MISSING_VRT_TAG,
                line_number=line_number,
                field_name="INFO"
            )
            return info_dict

        # Parse semicolon-separated key/value pairs.
        pairs = info_str.split(';')

        for pair in pairs:
            if not pair.strip():
                continue
            if '=' in pair:
                key, value = pair.split('=', 1)
                self._warn_if_undefined_info_tag(key, line_number, line)
                self._validate_target_info_value(key, value, line_number)
                if key == 'AF':
                    self._validate_info_allele_frequency(value, line_number)
                info_dict[key] = self._convert_info_value(key, value, line_number)
            else:
                # Flag-style tag without an explicit value.
                self._warn_if_undefined_info_tag(pair, line_number, line)
                info_dict[pair] = True

        # Validate the required VRT tag.
        if 'VRT' not in info_dict:
            self.error_handler.create_error(
                ErrorCode.MISSING_VRT_TAG,
                line_number=line_number,
                field_name="INFO",
                additional_info={"available_tags": list(info_dict.keys())}
            )
        else:
            vrt_value = info_dict['VRT']
            if not isinstance(vrt_value, int) or vrt_value < 1 or vrt_value > 8:
                self.error_handler.create_error(
                    ErrorCode.INVALID_VRT_VALUE,
                    line_number=line_number,
                    field_name="VRT",
                    actual_value=str(vrt_value),
                    expected_value="integer in range 1-8"
            )

        return info_dict

    def _warn_if_undefined_info_tag(self, key: str, line_number: int, line: Optional[str]) -> None:
        """Report the first data use of an INFO tag missing from the header."""
        if not key or key == '.' or key in self.header.info_tags or key in self._warned_undefined_info_tags:
            return

        self._warned_undefined_info_tags.add(key)
        self.error_handler.create_error(
            ErrorCode.UNDEFINED_INFO_TAG,
            line_number=line_number,
            line_content=line,
            field_name=key,
            expected_value="INFO tag definition in header",
            actual_value=f"{key} used in data row"
        )

    def _validate_target_info_value(self, key: str, value: str, line_number: int) -> None:
        """Validate selected INFO values against accepted value formats."""
        if key in {'AN', 'AC', 'AF'}:
            # Numeric INFO tags are type-checked during conversion.
            return

        pattern = self._target_info_patterns.get(key)
        if pattern is None or value == '':
            return

        if not re.fullmatch(pattern, value):
            self.error_handler.create_error(
                ErrorCode.INVALID_INFO_TAG_VALUE,
                line_number=line_number,
                field_name=key,
                actual_value=value,
                expected_value=f"value matching /{pattern}/"
            )

    def _validate_info_allele_frequency(self, value: str, line_number: int) -> None:
        """Validate each numeric INFO/AF value against the inclusive 0-1 range."""
        invalid_values = []
        for item in value.split(','):
            if item in {'', '.'}:
                continue
            try:
                frequency = float(item)
            except ValueError:
                continue
            if not 0.0 <= frequency <= 1.0:
                invalid_values.append(item)

        if invalid_values:
            self.error_handler.create_error(
                ErrorCode.INVALID_ALLELE_FREQUENCY,
                line_number=line_number,
                field_name='AF',
                actual_value=','.join(invalid_values),
                expected_value='numeric value in range 0-1'
            )

    def _convert_info_value(self, key: str, value: str, line_number: int) -> Any:
        """Convert INFO values according to their declared or known type."""
        # Keep missing values unchanged.
        if value == '.' or value == '':
            return value

        # Convert according to the declared INFO tag type.
        if key in self.header.info_tags:
            tag_def = self.header.info_tags[key]

            if tag_def.type == 'Integer':
                try:
                    # Handle comma-separated values for variable-length tags.
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
            elif tag_def.type == 'Float':
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
            elif tag_def.type == 'String':
                return value
            elif tag_def.type == 'Flag':
                return True

        # Try common numeric INFO tags even when the header definition is missing.
        if key in ['VRT', 'NIO', 'SAO', 'SSR', 'AC', 'AN']:
            try:
                if ',' in value:
                    return [int(v) for v in value.split(',')]
                return int(value)
            except ValueError:
                # Keep missing values unchanged.
                if value == '.' or value == '':
                    return value
                self.error_handler.create_error(
                    ErrorCode.INVALID_INFO_TAG_TYPE,
                    line_number=line_number,
                    field_name=key,
                    actual_value=value,
                    expected_value="integer"
                )
                return value
        elif key in ['AF', 'MAF']:
            try:
                if ',' in value:
                    return [float(v) for v in value.split(',')]
                return float(value)
            except ValueError:
                if value == '.' or value == '':
                    return value
                self.error_handler.create_error(
                    ErrorCode.INVALID_INFO_TAG_TYPE,
                    line_number=line_number,
                    field_name=key,
                    actual_value=value,
                    expected_value="float"
                )
                return value
        elif key == 'PMID':
            try:
                if ',' in value:
                    return [int(v) for v in value.split(',')]
                return int(value)
            except ValueError:
                return value
        elif key in ['AA', 'AD', 'CMT', 'LKO', 'OMIM', 'OMIA']:
            return value  # String-like tags

        # Treat undefined tags as strings.
        return value

    def _parse_format_and_population_data(
        self, fields: List[str], line_number: int, line: Optional[str] = None
    ) -> Tuple[Dict[str, Any], Dict[str, Dict[str, Any]]]:
        """Parse FORMAT and population data."""
        format_data = {}
        population_data = {}

        if not fields:
            return format_data, population_data

        # FORMAT field, stored in the first extra column.
        format_str = fields[0]
        format_keys = format_str.split(':') if format_str else []
        for format_key in format_keys:
            self._warn_if_undefined_format_tag(format_key, line_number, line)

        # Population data fields.
        population_fields = fields[1:] if len(fields) > 1 else []

        # Check population ID count against data column count.
        expected_pop_count = len(self.header.population_ids)
        actual_pop_count = len(population_fields)

        if expected_pop_count != actual_pop_count:
            self.error_handler.create_error(
                ErrorCode.POPULATION_ID_MISMATCH,
                line_number=line_number,
                expected_value=f"{expected_pop_count} Population columns",
                actual_value=f"{actual_pop_count} columns",
                additional_info={
                    "population_ids": self.header.population_ids,
                    "format_keys": format_keys
                }
            )

        # Parse data for each population.
        for i, pop_field in enumerate(population_fields):
            if i < len(self.header.population_ids):
                pop_id = self.header.population_ids[i]
                population_data[pop_id] = self._parse_population_field(
                    pop_field, format_keys, line_number, pop_id
                )
            else:
                # No population ID is defined for this column.
                self.error_handler.create_error(
                    ErrorCode.MISSING_POPULATION_ID,
                    line_number=line_number,
                    additional_info={"column_index": i}
                )

        return format_data, population_data

    def _warn_if_undefined_format_tag(self, key: str, line_number: int, line: Optional[str]) -> None:
        """Report the first data use of a FORMAT tag missing from the header."""
        if not key or key == '.' or key in self.header.format_tags or key in self._warned_undefined_format_tags:
            return

        self._warned_undefined_format_tags.add(key)
        self.error_handler.create_error(
            ErrorCode.UNDEFINED_FORMAT_TAG,
            line_number=line_number,
            line_content=line,
            field_name=key,
            expected_value="FORMAT tag definition in header",
            actual_value=f"{key} used in data row"
        )

    def _parse_population_field(
        self, pop_field: str, format_keys: List[str], line_number: int, pop_id: str
    ) -> Dict[str, Any]:
        """Parse one population field."""
        pop_data = {}

        if pop_field == '.':
            return pop_data

        # Colon-separated FORMAT values.
        values = pop_field.split(':')

        # Check that FORMAT keys and values have matching counts.
        if len(values) != len(format_keys):
            self.error_handler.create_error(
                ErrorCode.INVALID_POPULATION_FORMAT,
                line_number=line_number,
                field_name=f"Population {pop_id}",
                expected_value=f"{len(format_keys)} values",
                actual_value=f"{len(values)} values",
                additional_info={"format_keys": format_keys, "values": values}
            )

        for i, value in enumerate(values):
            if i < len(format_keys):
                key = format_keys[i]
                self._validate_target_format_value(key, value, line_number, pop_id)
                pop_data[key] = self._convert_format_value(key, value, line_number, pop_id)

        # Validate allele frequency.
        if 'FRQ' in pop_data:
            frq_value = pop_data['FRQ']
            if isinstance(frq_value, (int, float)):
                if frq_value < 0 or frq_value > 1:
                    self.error_handler.create_error(
                        ErrorCode.INVALID_ALLELE_FREQUENCY,
                        line_number=line_number,
                        field_name=f"FRQ ({pop_id})",
                        actual_value=str(frq_value),
                        expected_value="0-1 range"
                    )
            elif isinstance(frq_value, list):
                for idx, frq in enumerate(frq_value):
                    if isinstance(frq, (int, float)) and (frq < 0 or frq > 1):
                        self.error_handler.create_error(
                            ErrorCode.INVALID_ALLELE_FREQUENCY,
                            line_number=line_number,
                            field_name=f"FRQ[{idx}] ({pop_id})",
                            actual_value=str(frq),
                            expected_value="0-1 range"
                        )

        return pop_data

    def _validate_target_format_value(self, key: str, value: str, line_number: int, pop_id: str) -> None:
        """Validate selected FORMAT values against accepted value formats."""
        pattern = self._target_format_patterns.get(key)
        if pattern is None or value == '':
            return

        if not re.fullmatch(pattern, value):
            self.error_handler.create_error(
                ErrorCode.INVALID_FORMAT_TAG_VALUE,
                line_number=line_number,
                field_name=f"{key} ({pop_id})",
                actual_value=value,
                expected_value=f"value matching /{pattern}/"
            )

    def _convert_format_value(
        self, key: str, value: str, line_number: int, pop_id: str
    ) -> Any:
        """Convert a FORMAT value according to its declared type."""
        if key in self.header.format_tags:
            tag_def = self.header.format_tags[key]

            if tag_def.type == 'Integer':
                try:
                    return int(value)
                except ValueError:
                    self.error_handler.create_error(
                        ErrorCode.INVALID_POPULATION_VALUE_TYPE,
                        line_number=line_number,
                        field_name=f"{key} ({pop_id})",
                        actual_value=value,
                        expected_value="integer"
                    )
                    return value
            elif tag_def.type == 'Float':
                try:
                    # Number=. values may contain comma-separated entries.
                    if ',' in value:
                        return [float(v) for v in value.split(',')]
                    return float(value)
                except ValueError:
                    self.error_handler.create_error(
                        ErrorCode.INVALID_POPULATION_VALUE_TYPE,
                        line_number=line_number,
                        field_name=f"{key} ({pop_id})",
                        actual_value=value,
                        expected_value="float"
                    )
                    return value
            elif tag_def.type == 'String':
                return value

        # Treat undefined tags as strings.
        return value

    def _validate_parsed_data(self) -> None:
        """Validate parsed data as a whole."""
        # Check whether every data row contains VRT.
        for row in self.data_rows:
            if 'VRT' not in row.info:
                # The row-level parser may already have recorded an error.
                pass

        # Check for population IDs without population data.
        has_population_data = self._has_population_data or any(
            row.population_data for row in self.data_rows
        )
        if self.header.population_ids and not has_population_data:
            self.error_handler.create_error(
                ErrorCode.MISSING_POPULATION_DATA,
                additional_info={
                    "defined_population_ids": self.header.population_ids,
                    "message": "Population ID is defined but no data present"
                }
            )
