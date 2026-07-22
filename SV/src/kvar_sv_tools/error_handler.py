#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KVar SV VCF parser error handling module.
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


class ErrorAction(Enum):
    """Processing action taken for an error code."""
    BLOCK = "BLOCK"
    ACCEPT = "ACCEPT"
    SKIP_AND_CONTINUE = "SKIP_AND_CONTINUE"
    REPAIR_AND_CONTINUE = "REPAIR_AND_CONTINUE"
    FILL_AND_CONTINUE = "FILL_AND_CONTINUE"
    INFER_AND_CONTINUE = "INFER_AND_CONTINUE"
    SPLIT_AND_CONTINUE = "SPLIT_AND_CONTINUE"
    GENERATE_AND_CONTINUE = "GENERATE_AND_CONTINUE"


class ErrorCategory(Enum):
    """Error category"""
    FILE_ERROR = "FILE_ERROR"
    HEADER_ERROR = "HEADER_ERROR"
    METADATA_ERROR = "METADATA_ERROR"
    DATA_ERROR = "DATA_ERROR"
    VALIDATION_ERROR = "VALIDATION_ERROR"
    CONVERSION_ERROR = "CONVERSION_ERROR"
    FORMAT_ERROR = "FORMAT_ERROR"


class ErrorCode(Enum):
    """Error code definition"""
    # File-related (1000s)
    FILE_NOT_FOUND = (1001, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "File not found")
    FILE_READ_ERROR = (1002, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "File read error")
    FILE_WRITE_ERROR = (1003, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "File write error")
    FILE_ENCODING_ERROR = (1004, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "File encoding error")

    # Header-related (2000s)
    MISSING_FILEFORMAT = (2001, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "##fileformat metadata is missing")
    INVALID_FILEFORMAT = (2002, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "Not a VCF format (must start with VCFv)")
    MISSING_COLUMN_HEADER = (2003, ErrorSeverity.CRITICAL, ErrorCategory.HEADER_ERROR, "Column header is missing")
    INVALID_COLUMN_HEADER = (2004, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "Column header format is invalid")
    DUPLICATE_COLUMN_HEADER = (2005, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "Column header is duplicated")

    # Metadata-related (2100s)
    MISSING_REQUIRED_METADATA = (2101, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Required metadata is missing")
    INVALID_METADATA_FORMAT = (2102, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Metadata format is invalid")
    METADATA_MISMATCH = (2103, ErrorSeverity.ERROR, ErrorCategory.METADATA_ERROR, "Metadata file and VCF metadata do not match")

    # INFO tag (2200s)
    INFO_TAG_PARSE_ERROR = (2201, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "INFO tag definition parse error")
    INVALID_INFO_TAG_VALUE = (2204, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "INFO tag value is invalid")
    DUPLICATE_INFO_TAG_DEFINITION = (2205, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "INFO tag definition is duplicated")
    UNDEFINED_INFO_TAG = (2207, ErrorSeverity.WARNING, ErrorCategory.HEADER_ERROR, "INFO tag is not defined in the VCF header")

    # FORMAT tag (2300s)
    FORMAT_TAG_PARSE_ERROR = (2301, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "FORMAT tag definition parse error")
    DUPLICATE_FORMAT_TAG_DEFINITION = (2304, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "FORMAT tag definition is duplicated")
    INVALID_FORMAT_TAG_VALUE = (2305, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "FORMAT tag value is invalid")
    UNDEFINED_FORMAT_TAG = (2306, ErrorSeverity.WARNING, ErrorCategory.HEADER_ERROR, "FORMAT tag is not defined in the VCF header")

    # SV-specific INFO (2210s)
    MISSING_SVTYPE_TAG = (2210, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Required SVTYPE tag is missing")
    MISSING_END_TAG = (2211, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Required END tag is missing")
    MISSING_SVLEN_TAG = (2212, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "SVLEN tag is missing; derived length fields may be blank")
    MISSING_SAMPLESET_TAG = (2213, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Required SAMPLESET tag is missing")
    MISSING_EXPERIMENT_TAG = (2214, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Required EXPERIMENT tag is missing")
    INVALID_SVTYPE_VALUE = (2215, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "SVTYPE value is invalid (must be one of DEL, INS, DUP, INV, CNV, BND)")
    INVALID_END_VALUE = (2216, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "END value is invalid (must be >= POS)")

    # Data line (3000s)
    INSUFFICIENT_FIELDS = (3001, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "Insufficient required fields (min 8: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)")
    INVALID_POSITION = (3002, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "Position value is invalid (must be integer >= 1)")
    INVALID_REF_ALT = (3003, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "REF or ALT value is invalid")
    EMPTY_REF_ALT = (3004, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "REF or ALT is empty")
    INVALID_CHROMOSOME_FORMAT = (3005, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "CHROM value is invalid or empty")
    MISSING_REQUIRED_FIELD = (3006, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "Required field is missing")
    MULTIALLELIC_ALT = (3007, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "ALT must contain a single allele")
    LOCAL_ID_TOO_LONG = (3008, ErrorSeverity.ERROR, ErrorCategory.DATA_ERROR, "Local ID exceeds 64 characters")
    DUPLICATE_LOCAL_ID = (3009, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Local ID is duplicated")

    # Data validation (3100s)
    INVALID_BREAKEND_FORMAT = (3101, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Breakend (BND) format is invalid")
    MISSING_MATEID = (3102, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "MATEID missing for BND type")
    INVALID_POSRANGE_FORMAT = (3103, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "POSrange/CIPOS format or anchor is invalid")
    INVALID_ENDRANGE_FORMAT = (3104, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "ENDrange/CIEND format or anchor is invalid")
    INVALID_valEXPERIMENT_FORMAT = (3105, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "valEXPERIMENT format is invalid (use ExperimentID:Pass or ExperimentID:Fail; comma-separate multiple)")
    CHROMOSOME_NOT_GROUPED = (3106, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Data on the same chromosome are not grouped together")
    POSITION_NOT_SORTED = (3107, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Positions are not sorted within chromosome")
    ALLELE_COUNT_EXCEEDS_NUMBER = (3108, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Allele count exceeds allele number")
    MATEID_NOT_FOUND = (3109, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "MATEID target variant is not present in VCF")
    MATEID_NOT_RECIPROCAL = (3110, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "MATEID target does not reciprocally reference the source variant")
    MATEID_ALT_MISMATCH = (3111, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "Reciprocal MATEID breakend ALT coordinates are inconsistent")
    BND_INSERTION_SPLIT = (3112, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "Breakend insertion sequence was split into a separate insertion call", ErrorAction.SPLIT_AND_CONTINUE)
    MATEID_STRAND_MISMATCH = (3113, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "Reciprocal MATEID breakend ALT strand orientations are inconsistent")

    # Conversion (4000s)
    CONVERSION_FAILED = (4001, ErrorSeverity.ERROR, ErrorCategory.CONVERSION_ERROR, "TSV conversion failed")

    # Format (5000s)
    INVALID_TAB_SEPARATOR = (5001, ErrorSeverity.ERROR, ErrorCategory.FORMAT_ERROR, "Tab separator is invalid")
    INVALID_INFO_FORMAT = (5002, ErrorSeverity.WARNING, ErrorCategory.FORMAT_ERROR, "INFO field format is invalid")
    FIELD_WHITESPACE_TRIMMED = (5003, ErrorSeverity.WARNING, ErrorCategory.FORMAT_ERROR, "Leading or trailing whitespace was trimmed from a field", ErrorAction.REPAIR_AND_CONTINUE)
    EMPTY_LINE_REMOVED = (5004, ErrorSeverity.WARNING, ErrorCategory.FORMAT_ERROR, "Empty line was removed", ErrorAction.REPAIR_AND_CONTINUE)

    # Reference genome validation (6000s)
    REF_MISMATCH = (6001, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "VCF REF sequence does not match reference genome")
    CHROMOSOME_NOT_FOUND = (6002, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Chromosome not found in reference genome")
    POSITION_OUT_OF_RANGE = (6003, ErrorSeverity.ERROR, ErrorCategory.VALIDATION_ERROR, "Position is out of reference genome range")
    FASTA_INDEX_ERROR = (6004, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "FASTA index file create/read error")
    FASTA_READ_ERROR = (6005, ErrorSeverity.CRITICAL, ErrorCategory.FILE_ERROR, "FASTA file read error")
    REF_CHECK_SKIPPED = (6006, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "Reference FASTA validation was skipped", ErrorAction.SKIP_AND_CONTINUE)

    # DDBJ-style VCF header validation (6100s)
    FILEFORMAT_NOT_FIRST = (6101, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "##fileformat must be the first non-empty VCF line")
    DUPLICATE_METADATA_TAG = (6102, ErrorSeverity.ERROR, ErrorCategory.HEADER_ERROR, "VCF metadata tag is duplicated")

    # Variant Call output validation (7100s)
    INVALID_SEQUENCE_FIELD = (7106, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "Sequence contains characters outside the accepted IUPAC set")
    INVALID_PHENOTYPE_LINK = (7107, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "Phenotype should reference a valid phenotype db:id")
    INVALID_EVIDENCE_LINK = (7108, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "Evidence should reference a valid db:id")
    INVALID_EXTERNAL_LINK = (7109, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "External Links should reference a valid db:id or URL")
    INVALID_HGVSG = (7124, ErrorSeverity.WARNING, ErrorCategory.VALIDATION_ERROR, "HGVSG is invalid or inconsistent with the call placement", ErrorAction.REPAIR_AND_CONTINUE)

    def __init__(
        self,
        code: int,
        severity: ErrorSeverity,
        category: ErrorCategory,
        message: str,
        action: Optional[ErrorAction] = None
    ):
        self.code = code
        self.severity = severity
        self.category = category
        self.message = message
        if action is None:
            action = ErrorAction.BLOCK if severity in {ErrorSeverity.CRITICAL, ErrorSeverity.ERROR} else ErrorAction.ACCEPT
        self.action = action


@dataclass
class ParseError:
    """Class to store parsing error information"""
    error_code: ErrorCode
    line_number: Optional[int] = None
    variant_id: Optional[str] = None
    line_content: Optional[str] = None
    field_name: Optional[str] = None
    expected_value: Optional[str] = None
    actual_value: Optional[str] = None
    additional_info: Dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    def format(self, include_line_content: bool = True) -> str:
        """Represent error as string, optionally omitting raw input content."""
        parts = [
            f"[{self.error_code.code}] {self.error_code.severity.value}",
            f"action: {self.error_code.action.value}",
            f"{self.error_code.category.value}",
            f"{self.error_code.message}"
        ]

        if self.variant_id:
            parts.append(f"(ID: {self.variant_id})")
        elif self.line_number:
            parts.append(f"(line {self.line_number})")

        if self.field_name:
            parts.append(f"field: {self.field_name}")

        if self.expected_value:
            parts.append(f"expected: {self.expected_value}")

        if self.actual_value:
            parts.append(f"actual: {self.actual_value}")

        if include_line_content and self.line_content:
            parts.append(f"content: {self.line_content[:100]}")

        return " | ".join(parts)

    def __str__(self) -> str:
        """Represent error as string"""
        return self.format()

    def to_dict(self) -> Dict[str, Any]:
        """Convert error to dictionary"""
        return {
            "error_code": self.error_code.code,
            "error_name": self.error_code.name,
            "severity": self.error_code.severity.value,
            "action": self.error_code.action.value,
            "category": self.error_code.category.value,
            "message": self.error_code.message,
            "line_number": self.line_number,
            "variant_id": self.variant_id,
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
        self.action_counts: Dict[ErrorAction, int] = {}
        self.category_counts: Dict[ErrorCategory, int] = {}

    def add_error(self, error: ParseError) -> None:
        """에러를 추가"""
        self.errors.append(error)

        # Update statistics
        self.error_counts[error.error_code] = self.error_counts.get(error.error_code, 0) + 1
        self.severity_counts[error.error_code.severity] = self.severity_counts.get(error.error_code.severity, 0) + 1
        self.action_counts[error.error_code.action] = self.action_counts.get(error.error_code.action, 0) + 1
        self.category_counts[error.error_code.category] = self.category_counts.get(error.error_code.category, 0) + 1

    def create_error(
        self,
        error_code: ErrorCode,
        line_number: Optional[int] = None,
        variant_id: Optional[str] = None,
        line_content: Optional[str] = None,
        field_name: Optional[str] = None,
        expected_value: Optional[str] = None,
        actual_value: Optional[str] = None,
        additional_info: Optional[Dict[str, Any]] = None
    ) -> ParseError:
        """에러를 생성하고 추가"""
        error = ParseError(
            error_code=error_code,
            line_number=line_number,
            variant_id=variant_id,
            line_content=line_content,
            field_name=field_name,
            expected_value=expected_value,
            actual_value=actual_value,
            additional_info=additional_info or {}
        )
        self.add_error(error)
        return error

    def has_errors(self, severity: Optional[ErrorSeverity] = None) -> bool:
        """에러가 있는지 확인"""
        if severity is None:
            return len(self.errors) > 0
        return any(e.error_code.severity == severity for e in self.errors)

    def has_critical_errors(self) -> bool:
        """Critical 에러가 있는지 확인"""
        return self.has_errors(ErrorSeverity.CRITICAL)

    def has_fatal_errors(self) -> bool:
        """Check if processing must stop immediately."""
        return self.has_critical_errors()

    def has_blocking_errors(self) -> bool:
        """Check if final outputs must be blocked."""
        return any(e.error_code.action == ErrorAction.BLOCK for e in self.errors)

    def get_validation_status(self) -> str:
        """Return the final validation status from collected messages."""
        if self.has_blocking_errors():
            return "ERROR"
        if self.has_errors(ErrorSeverity.ERROR):
            return "ERROR"
        if self.has_errors(ErrorSeverity.WARNING):
            return "WARNING"
        return "OK"

    def assert_no_blocking_errors(
        self,
        stage: str = "validation",
        output_file: Optional[str] = None,
        vcf_file_path: Optional[str] = None,
        output_tsv_path: Optional[str] = None,
        tsv_file_path: Optional[str] = None,
        sanitize_paths: bool = False,
        include_line_content: bool = True,
        include_additional_info: bool = True
    ) -> None:
        """Raise if any collected message has a blocking action."""
        if not self.has_blocking_errors():
            return

        report_tsv_path = output_tsv_path or tsv_file_path
        if output_file:
            self.generate_report(
                output_file=output_file,
                vcf_file_path=vcf_file_path,
                output_tsv_path=report_tsv_path,
                sanitize_paths=sanitize_paths,
                include_line_content=include_line_content,
                include_additional_info=include_additional_info
            )

        summary = self.get_summary()
        raise RuntimeError(
            f"{stage} blocked by "
            f"{summary['blocking_count']} blocking message(s) "
            f"({summary['critical_count']} critical, {summary['error_count']} error)."
        )

    def get_errors_by_severity(self, severity: ErrorSeverity) -> List[ParseError]:
        """심각도별로 에러 필터링"""
        return [e for e in self.errors if e.error_code.severity == severity]

    def get_errors_by_action(self, action: ErrorAction) -> List[ParseError]:
        """Filter errors by processing action."""
        return [e for e in self.errors if e.error_code.action == action]

    def get_errors_by_category(self, category: ErrorCategory) -> List[ParseError]:
        """카테고리별로 에러 필터링"""
        return [e for e in self.errors if e.error_code.category == category]

    def get_errors_by_code(self, error_code: ErrorCode) -> List[ParseError]:
        """에러 코드별로 에러 필터링"""
        return [e for e in self.errors if e.error_code == error_code]

    def get_summary(self) -> Dict[str, Any]:
        """에러 요약 정보 반환"""
        return {
            "total_errors": len(self.errors),
            "critical_count": self.severity_counts.get(ErrorSeverity.CRITICAL, 0),
            "error_count": self.severity_counts.get(ErrorSeverity.ERROR, 0),
            "warning_count": self.severity_counts.get(ErrorSeverity.WARNING, 0),
            "info_count": self.severity_counts.get(ErrorSeverity.INFO, 0),
            "blocking_count": len(self.get_errors_by_action(ErrorAction.BLOCK)),
            "error_counts_by_code": {code.name: count for code, count in self.error_counts.items()},
            "errors_by_category": {cat.value: len(self.get_errors_by_category(cat)) for cat in ErrorCategory},
            "errors_by_severity": {sev.value: len(self.get_errors_by_severity(sev)) for sev in ErrorSeverity},
            "errors_by_action": {action.value: len(self.get_errors_by_action(action)) for action in ErrorAction}
        }

    @staticmethod
    def _display_path(path: Optional[str], sanitize_paths: bool) -> str:
        if not path:
            return "(not specified)"
        if sanitize_paths:
            return os.path.join("(redacted)", os.path.basename(path))
        return os.path.abspath(path)

    @classmethod
    def _display_additional_info(cls, key: str, value: Any, sanitize_paths: bool) -> Any:
        if not sanitize_paths:
            return value

        key_lower = key.lower()
        if key_lower.endswith("path") or key_lower in {"path", "file", "filename"}:
            return cls._display_path(str(value), sanitize_paths=True)
        if isinstance(value, (list, tuple)):
            return [
                cls._display_additional_info(key, item, sanitize_paths)
                for item in value
            ]
        if isinstance(value, dict):
            return {
                nested_key: cls._display_additional_info(str(nested_key), nested_value, sanitize_paths)
                for nested_key, nested_value in value.items()
            }
        return value

    def generate_report(
        self,
        output_file: Optional[str] = None,
        vcf_file_path: Optional[str] = None,
        tsv_file_path: Optional[str] = None,
        output_tsv_path: Optional[str] = None,
        sanitize_paths: bool = False,
        include_line_content: bool = True,
        include_additional_info: bool = True,
        report_display_path: Optional[str] = None
    ) -> str:
        """에러 리포트 생성"""
        if output_tsv_path is not None:
            tsv_file_path = output_tsv_path

        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("KVar SV VCF Parsing Error Report")
        report_lines.append("=" * 80)

        path_label_suffix = "" if sanitize_paths else " absolute path"
        report_lines.append(
            f"Input VCF file{path_label_suffix}: {self._display_path(vcf_file_path, sanitize_paths)}"
        )
        report_lines.append(
            f"Output TSV file{path_label_suffix}: {self._display_path(tsv_file_path, sanitize_paths)}"
        )
        report_lines.append(
            f"Error report{path_label_suffix}: "
            f"{self._display_path(report_display_path or output_file, sanitize_paths)}"
        )

        report_lines.append(f"Generated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append("")

        summary = self.get_summary()
        report_lines.append("=== Error Summary ===")
        report_lines.append(f"Total errors: {summary['total_errors']}")
        report_lines.append(f"  Critical: {summary['critical_count']}")
        report_lines.append(f"  Error: {summary['error_count']}")
        report_lines.append(f"  Warning: {summary['warning_count']}")
        report_lines.append(f"  Info: {summary['info_count']}")
        report_lines.append(f"  Blocking: {summary['blocking_count']}")
        report_lines.append("")

        report_lines.append("=== Actions ===")
        for action, count in summary['errors_by_action'].items():
            if count > 0:
                report_lines.append(f"  {action}: {count}")
        report_lines.append("")

        report_lines.append("=== Errors by Category ===")
        for category, count in summary['errors_by_category'].items():
            if count > 0:
                report_lines.append(f"  {category}: {count}")
        report_lines.append("")

        if summary['error_counts_by_code']:
            report_lines.append("=== Statistics by Error Code ===")
            for code_name, count in sorted(summary['error_counts_by_code'].items(), key=lambda x: x[1], reverse=True):
                report_lines.append(f"  {code_name}: {count}")
            report_lines.append("")

        if self.errors:
            report_lines.append("=== Detailed Error List ===")

            for severity in [ErrorSeverity.CRITICAL, ErrorSeverity.ERROR, ErrorSeverity.WARNING, ErrorSeverity.INFO]:
                errors = self.get_errors_by_severity(severity)
                if errors:
                    report_lines.append(f"\n[{severity.value}] ({len(errors)} items)")
                    report_lines.append("-" * 80)
                    for i, error in enumerate(errors, 1):
                        report_lines.append(f"{i}. {error.format(include_line_content=include_line_content)}")
                        if include_additional_info and error.additional_info:
                            for key, value in error.additional_info.items():
                                display_value = self._display_additional_info(key, value, sanitize_paths)
                                report_lines.append(f"   - {key}: {display_value}")

        report_text = "\n".join(report_lines)

        # Save to file
        if output_file:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(report_text)

        return report_text

    def generate_json_report(
        self,
        output_file: Optional[str] = None,
        include_additional_info: bool = True
    ) -> Dict[str, Any]:
        """Generate error report in JSON format"""
        report = {
            "timestamp": datetime.now().isoformat(),
            "summary": self.get_summary(),
            "errors": [error.to_dict() for error in self.errors]
        }
        if not include_additional_info:
            for error in report["errors"]:
                error["additional_info"] = {}

        if output_file:
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(report, f, ensure_ascii=False, indent=2)

        return report

    def clear(self) -> None:
        """Clear all errors"""
        self.errors.clear()
        self.error_counts.clear()
        self.severity_counts.clear()
        self.action_counts.clear()
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
        print(f"  Blocking: {summary['blocking_count']}")

        if summary['total_errors'] > 0:
            print("\nTop error codes:")
            for code_name, count in sorted(summary['error_counts_by_code'].items(), key=lambda x: x[1], reverse=True)[:10]:
                print(f"  {code_name}: {count}")
