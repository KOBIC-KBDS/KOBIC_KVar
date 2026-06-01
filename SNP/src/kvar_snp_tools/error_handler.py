#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dbSNP/KVar VCF parser error handling module.
Provides error code scheme and error collection/report functionality.
"""

from enum import Enum
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from datetime import datetime
import os
import json


class ErrorSeverity(Enum):
    """Error severity level"""
    CRITICAL = "CRITICAL"  # Severe error that requires parsing to stop
    ERROR = "ERROR"        # Error that may cause data loss
    WARNING = "WARNING"    # Warning level, processable but requires attention
    INFO = "INFO"          # Informational message


class ErrorCategory(Enum):
    """Error category"""
    FILE_ERROR = "FILE_ERROR"           # File read/write error
    HEADER_ERROR = "HEADER_ERROR"       # Header parsing error
    METADATA_ERROR = "METADATA_ERROR"    # Metadata error
    DATA_ERROR = "DATA_ERROR"            # Data line parsing error
    VALIDATION_ERROR = "VALIDATION_ERROR" # Data validation error
    CONVERSION_ERROR = "CONVERSION_ERROR" # Conversion error
    FORMAT_ERROR = "FORMAT_ERROR"        # Format error


class ErrorCode(Enum):
    """Error code definition"""
    # File-related errors (1000s)
    FILE_NOT_FOUND = (1001, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "File not found")
    FILE_READ_ERROR = (1002, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "File read error")
    FILE_WRITE_ERROR = (1003, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "File write error")
    FILE_ENCODING_ERROR = (1004, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "File encoding error")
    
    # Header-related errors (2000s)
    MISSING_FILEFORMAT = (2001, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "##fileformat metadata is missing")
    INVALID_FILEFORMAT = (2002, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "Not a VCF format (must start with VCFv)")
    MISSING_COLUMN_HEADER = (2003, ErrorSeverity.CRITICAL, ErrorCategory.HEADER_ERROR, "Column header is missing")
    INVALID_COLUMN_HEADER = (2004, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "Column header format is invalid")
    DUPLICATE_COLUMN_HEADER = (2005, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "Column header is duplicated")
    
    # Metadata-related errors (2100s)
    MISSING_REQUIRED_METADATA = (2101, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Required metadata is missing")
    INVALID_METADATA_FORMAT = (2102, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Metadata format is invalid")
    METADATA_MISMATCH = (2103, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Metadata file and VCF file metadata do not match")
    METADATA_VALUE_FILLED = (2104, ErrorSeverity.WARNING, ErrorCategory.METADATA_ERROR, "VCF metadata value was filled from metadata file")
    METADATA_VALUE_CORRECTED = (2105, ErrorSeverity.WARNING, ErrorCategory.METADATA_ERROR, "VCF metadata value was replaced with metadata file value")
    DUPLICATE_METADATA_TAG = (2106, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "Metadata tag is duplicated")
    METADATA_NOT_PARSED = (2107, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Metadata file has not been parsed")
    MISSING_BIOSAMPLE_ID = (2108, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "biosample_id is missing from VCF metadata")
    METADATA_IDENTIFIER_CONFLICT = (2109, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Metadata identifiers conflict with each other")
    METADATA_REFERENCE_MISMATCH = (2110, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Metadata reference and VCF reference do not match")
    METADATA_BIOSAMPLE_MISMATCH = (2111, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Metadata biosample_id and VCF biosample_id do not match")
    METADATA_POPULATION_MISMATCH = (2112, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Metadata SampleSet_id and VCF population_id do not match")
    
    # INFO tag-related errors (2200s)
    INFO_TAG_PARSE_ERROR = (2201, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "INFO tag definition parse error")
    MISSING_VRT_TAG = (2202, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Required VRT tag is missing")
    INVALID_VRT_VALUE = (2203, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "VRT value is invalid (must be in range 1-8)")
    INVALID_INFO_TAG_VALUE = (2204, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "INFO tag value is invalid")
    DUPLICATE_INFO_TAG_DEFINITION = (2205, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "INFO tag definition is duplicated")
    INVALID_VRT_HEADER_DEFINITION = (2206, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "VRT INFO header definition is invalid")
    UNDEFINED_INFO_TAG = (2207, ErrorSeverity.WARNING, ErrorCategory.HEADER_ERROR, "INFO tag is not defined in the VCF header")
    INVALID_INFO_TAG_TYPE = (2208, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "INFO tag value type is invalid")
    
    # FORMAT tag-related errors (2300s)
    FORMAT_TAG_PARSE_ERROR = (2301, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "FORMAT tag definition parse error")
    MISSING_POPULATION_ID = (2302, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Population ID is not defined")
    POPULATION_ID_MISMATCH = (2303, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Population ID and data columns do not match")
    DUPLICATE_FORMAT_TAG_DEFINITION = (2304, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "FORMAT tag definition is duplicated")
    INVALID_FORMAT_TAG_VALUE = (2305, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "FORMAT tag value is invalid")
    UNDEFINED_FORMAT_TAG = (2306, ErrorSeverity.WARNING, ErrorCategory.HEADER_ERROR, "FORMAT tag is not defined in the VCF header")
    MISSING_POPULATION_DATA = (2307, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "Population ID is defined but population data is missing")
    
    # Data line-related errors (3000s)
    INSUFFICIENT_FIELDS = (3001, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "Insufficient required fields (minimum 8: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)")
    INVALID_POSITION = (3002, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "Position value is invalid (must be integer)")
    INVALID_REF_ALT = (3003, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "REF or ALT value is invalid")
    EMPTY_REF_ALT = (3004, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "REF or ALT is empty")
    VARIANT_TOO_LONG = (3005, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Variant exceeds 50bp; submit as SV for variants >50bp")
    INVALID_REF_BASE = (3006, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "REF contains non-ATGC characters (only A, T, G, C allowed; N, IUPAC ambiguity codes, and '*' are not permitted)")
    INVALID_ALT_ALLELE = (3007, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "ALT contains non-ATGC characters (only A, T, G, C allowed; N, IUPAC ambiguity codes, and '*' are not permitted)")
    SAME_REF_ALT = (3008, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "REF and ALT must not be identical")
    MULTIALLELIC_ALT = (3009, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "ALT must contain a single allele after normalization")
    INVALID_INDEL_LEADING_BASE = (3010, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Indel REF and ALT must share the leading base")
    MISSING_LOCAL_ID = (3011, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "Local ID is missing")
    LOCAL_ID_TOO_LONG = (3012, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "Local ID exceeds 64 characters")
    DUPLICATE_LOCAL_ID = (3013, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Local ID is duplicated")
    
    # Data validation-related errors (3100s)
    VRT_REF_ALT_MISMATCH = (3101, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "VRT value was corrected to match REF/ALT length-based classification")
    INVALID_ALLELE_FREQUENCY = (3102, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Allele frequency value is invalid (must be in range 0-1)")
    MISSING_ALLELE_FREQUENCY = (3103, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "Allele frequency information is missing")
    DUPLICATE_VARIANT_SITE = (3104, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Variant site is duplicated")
    ALLELE_COUNT_EXCEEDS_NUMBER = (3105, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Allele count exceeds allele number")
    CHROMOSOME_NOT_GROUPED = (3106, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Data on the same chromosome are not grouped together")
    POSITION_NOT_SORTED = (3107, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Positions are not sorted within chromosome")
    SNP_DENSITY_TOO_HIGH = (3108, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "More than 10 SNPs are present in a 50bp window")
    ALLELE_FREQUENCY_CALCULATED = (3109, ErrorSeverity.INFO, ErrorCategory.VALIDATION_ERROR, "Allele frequency was calculated from AC and AN")
    INVALID_ALLELE_COUNT = (3110, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Allele count must be greater than zero")
    
    # Conversion-related errors (4000s)
    CONVERSION_FAILED = (4001, ErrorSeverity.ERROR, ErrorCategory.CONVERSION_ERROR, "TSV conversion failed")
    MISSING_REQUIRED_FIELD = (4002, ErrorSeverity.ERROR, ErrorCategory.CONVERSION_ERROR, "Required field is missing; cannot convert")
    
    # Format-related errors (5000s)
    INVALID_TAB_SEPARATOR = (5001, ErrorSeverity.ERROR, ErrorCategory.FORMAT_ERROR, "Tab separator is invalid")
    INVALID_INFO_FORMAT = (5002, ErrorSeverity.WARNING, ErrorCategory.FORMAT_ERROR, "INFO field format is invalid")
    INVALID_POPULATION_FORMAT = (5003, ErrorSeverity.ERROR, ErrorCategory.FORMAT_ERROR, "Population data format is invalid")
    INVALID_POPULATION_VALUE_TYPE = (5004, ErrorSeverity.ERROR, ErrorCategory.FORMAT_ERROR, "Population data value type is invalid")
    
    # Reference genome validation errors (6000s)
    REF_MISMATCH = (6001, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "VCF REF sequence does not match reference genome")
    CHROMOSOME_NOT_FOUND = (6002, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Chromosome not found in reference genome")
    POSITION_OUT_OF_RANGE = (6003, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Position is out of reference genome range")
    FASTA_INDEX_ERROR = (6004, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "FASTA index file create/read error")
    FASTA_READ_ERROR = (6005, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "FASTA file read error")
    REFERENCE_CORRECTED = (6006, ErrorSeverity.INFO, ErrorCategory.METADATA_ERROR, "Reference replaced with metadata file value in output")
    
    def __init__(self, code: int, severity: ErrorSeverity, category: ErrorCategory, message: str):
        self.code = code
        self.severity = severity
        self.category = category
        self.message = message


@dataclass
class ParseError:
    """Class to store parsing error information"""
    error_code: ErrorCode
    line_number: Optional[int] = None
    line_content: Optional[str] = None
    field_name: Optional[str] = None
    expected_value: Optional[str] = None
    actual_value: Optional[str] = None
    additional_info: Dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    
    def __str__(self) -> str:
        """Represent error as string"""
        parts = [
            f"[{self.error_code.code}] {self.error_code.severity.value}",
            f"{self.error_code.category.value}",
            f"{self.error_code.message}"
        ]
        
        if self.line_number:
            parts.append(f"(line {self.line_number})")
        
        if self.field_name:
            parts.append(f"field: {self.field_name}")
        
        if self.expected_value:
            parts.append(f"expected: {self.expected_value}")
        
        if self.actual_value:
            parts.append(f"actual: {self.actual_value}")
        
        if self.line_content:
            parts.append(f"content: {self.line_content[:100]}")
        
        return " | ".join(parts)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert error to dictionary"""
        return {
            "error_code": self.error_code.code,
            "error_name": self.error_code.name,
            "severity": self.error_code.severity.value,
            "category": self.error_code.category.value,
            "message": self.error_code.message,
            "line_number": self.line_number,
            "field_name": self.field_name,
            "expected_value": self.expected_value,
            "actual_value": self.actual_value,
            "additional_info": self.additional_info,
            "timestamp": self.timestamp
        }


class ErrorHandler:
    """Error collection and management class"""
    
    def __init__(self):
        self.errors: List[ParseError] = []
        self.error_counts: Dict[ErrorCode, int] = {}
        self.severity_counts: Dict[ErrorSeverity, int] = {}
        self.category_counts: Dict[ErrorCategory, int] = {}
    
    def add_error(self, error: ParseError) -> None:
        """Add error"""
        self.errors.append(error)
        
        # Update statistics
        self.error_counts[error.error_code] = self.error_counts.get(error.error_code, 0) + 1
        self.severity_counts[error.error_code.severity] = self.severity_counts.get(error.error_code.severity, 0) + 1
        self.category_counts[error.error_code.category] = self.category_counts.get(error.error_code.category, 0) + 1
    
    def create_error(
        self,
        error_code: ErrorCode,
        line_number: Optional[int] = None,
        line_content: Optional[str] = None,
        field_name: Optional[str] = None,
        expected_value: Optional[str] = None,
        actual_value: Optional[str] = None,
        additional_info: Optional[Dict[str, Any]] = None
    ) -> ParseError:
        """Create and add error"""
        error = ParseError(
            error_code=error_code,
            line_number=line_number,
            line_content=line_content,
            field_name=field_name,
            expected_value=expected_value,
            actual_value=actual_value,
            additional_info=additional_info or {}
        )
        self.add_error(error)
        return error
    
    def has_errors(self, severity: Optional[ErrorSeverity] = None) -> bool:
        """Check if any errors exist"""
        if severity is None:
            return len(self.errors) > 0
        return any(e.error_code.severity == severity for e in self.errors)
    
    def has_critical_errors(self) -> bool:
        """Check if any critical errors exist"""
        return self.has_errors(ErrorSeverity.CRITICAL)

    def has_fatal_errors(self) -> bool:
        """Check if processing must stop immediately."""
        return self.has_critical_errors()

    def has_blocking_errors(self) -> bool:
        """Check if final outputs must be blocked."""
        return self.has_errors(ErrorSeverity.CRITICAL) or self.has_errors(ErrorSeverity.ERROR)

    def get_validation_status(self) -> str:
        """Return the final validation status from collected messages."""
        if self.has_blocking_errors():
            return "ERROR"
        if self.has_errors(ErrorSeverity.WARNING):
            return "WARNING"
        return "OK"

    def assert_no_blocking_errors(
        self,
        stage: str = "validation",
        output_file: Optional[str] = None,
        vcf_file_path: Optional[str] = None,
        output_tsv_path: Optional[str] = None
    ) -> None:
        """Raise if CRITICAL or ERROR messages are present."""
        if not self.has_blocking_errors():
            return

        if output_file:
            self.generate_report(
                output_file=output_file,
                vcf_file_path=vcf_file_path,
                output_tsv_path=output_tsv_path
            )

        summary = self.get_summary()
        raise RuntimeError(
            f"{stage} blocked by "
            f"{summary['critical_count']} critical and "
            f"{summary['error_count']} error message(s)."
        )
    
    def get_errors_by_severity(self, severity: ErrorSeverity) -> List[ParseError]:
        """Filter errors by severity"""
        return [e for e in self.errors if e.error_code.severity == severity]
    
    def get_errors_by_category(self, category: ErrorCategory) -> List[ParseError]:
        """Filter errors by category"""
        return [e for e in self.errors if e.error_code.category == category]
    
    def get_errors_by_code(self, error_code: ErrorCode) -> List[ParseError]:
        """Filter errors by error code"""
        return [e for e in self.errors if e.error_code == error_code]
    
    def get_summary(self) -> Dict[str, Any]:
        """Return error summary"""
        return {
            "total_errors": len(self.errors),
            "critical_count": self.severity_counts.get(ErrorSeverity.CRITICAL, 0),
            "error_count": self.severity_counts.get(ErrorSeverity.ERROR, 0),
            "warning_count": self.severity_counts.get(ErrorSeverity.WARNING, 0),
            "info_count": self.severity_counts.get(ErrorSeverity.INFO, 0),
            "error_counts_by_code": {code.name: count for code, count in self.error_counts.items()},
            "errors_by_category": {cat.value: len(self.get_errors_by_category(cat)) for cat in ErrorCategory},
            "errors_by_severity": {sev.value: len(self.get_errors_by_severity(sev)) for sev in ErrorSeverity}
        }
    
    def generate_report(
        self, 
        output_file: Optional[str] = None,
        vcf_file_path: Optional[str] = None,
        output_tsv_path: Optional[str] = None
    ) -> str:
        """Generate error report"""
        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("")
        report_lines.append("")
        report_lines.append("KVar SNP VCF Parsing Error Report")
        report_lines.append("")
        report_lines.append("=" * 80)
        
        # Add absolute path information
        if vcf_file_path:
            abs_vcf_path = os.path.abspath(vcf_file_path)
            report_lines.append(f"VCF file absolute path: {abs_vcf_path}")
        
        if output_file:
            abs_error_path = os.path.abspath(output_file)
            report_lines.append(f"Error report absolute path: {abs_error_path}")
        
        if output_tsv_path:
            abs_tsv_path = os.path.abspath(output_tsv_path)
            report_lines.append(f"Output TSV file absolute path: {abs_tsv_path}")
        
        report_lines.append(f"Generated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append("")
        
        # Summary
        summary = self.get_summary()
        report_lines.append("=== Error Summary ===")
        report_lines.append(f"Total errors: {summary['total_errors']}")
        report_lines.append(f"  Critical: {summary['critical_count']}")
        report_lines.append(f"  Error: {summary['error_count']}")
        report_lines.append(f"  Warning: {summary['warning_count']}")
        report_lines.append(f"  Info: {summary['info_count']}")
        report_lines.append("")
        
        # Statistics by category
        report_lines.append("=== Errors by Category ===")
        for category, count in summary['errors_by_category'].items():
            if count > 0:
                report_lines.append(f"  {category}: {count}")
        report_lines.append("")
        
        # Statistics by error code
        if summary['error_counts_by_code']:
            report_lines.append("=== Statistics by Error Code ===")
            for code_name, count in sorted(summary['error_counts_by_code'].items(), key=lambda x: x[1], reverse=True):
                report_lines.append(f"  {code_name}: {count}")
            report_lines.append("")
        
        # Detailed error list
        if self.errors:
            report_lines.append("=== Detailed Error List ===")
            
            # Group by severity
            for severity in [ErrorSeverity.CRITICAL, ErrorSeverity.ERROR, ErrorSeverity.WARNING, ErrorSeverity.INFO]:
                errors = self.get_errors_by_severity(severity)
                if errors:
                    report_lines.append(f"\n[{severity.value}] ({len(errors)} items)")
                    report_lines.append("-" * 80)
                    for i, error in enumerate(errors, 1):
                        report_lines.append(f"{i}. {error}")
                        if error.additional_info:
                            for key, value in error.additional_info.items():
                                report_lines.append(f"   - {key}: {value}")
        
        report_text = "\n".join(report_lines)
        
        # Save to file
        if output_file:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(report_text)
        
        return report_text
    
    def generate_json_report(self, output_file: Optional[str] = None) -> Dict[str, Any]:
        """Generate error report in JSON format"""
        report = {
            "timestamp": datetime.now().isoformat(),
            "summary": self.get_summary(),
            "errors": [error.to_dict() for error in self.errors]
        }
        
        if output_file:
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(report, f, ensure_ascii=False, indent=2)
        
        return report
    
    def clear(self) -> None:
        """Clear all errors"""
        self.errors.clear()
        self.error_counts.clear()
        self.severity_counts.clear()
        self.category_counts.clear()
    
    def print_summary(self) -> None:
        """Print error summary to console"""
        summary = self.get_summary()
        print("\n=== Parsing Error Summary ===")
        print(f"Total errors: {summary['total_errors']}")
        print(f"  Critical: {summary['critical_count']}")
        print(f"  Error: {summary['error_count']}")
        print(f"  Warning: {summary['warning_count']}")
        print(f"  Info: {summary['info_count']}")
        
        if summary['total_errors'] > 0:
            print("\nTop error codes:")
            for code_name, count in sorted(summary['error_counts_by_code'].items(), key=lambda x: x[1], reverse=True)[:10]:
                print(f"  {code_name}: {count}")
