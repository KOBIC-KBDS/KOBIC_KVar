#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metadata validation module.
Validates by comparing metadata file with VCF file metadata.
"""

import os
from typing import Dict, Optional, List, Set
from dataclasses import dataclass

try:
    from .error_handler import ErrorHandler, ErrorCode, ErrorSeverity, ErrorCategory
except ImportError:
    from error_handler import ErrorHandler, ErrorCode, ErrorSeverity, ErrorCategory


@dataclass
class MetadataInfo:
    """Class to store metadata information"""
    sampleset_id: Optional[str] = None
    experiment_id: Optional[str] = None
    reference: Optional[str] = None


class MetadataValidator:
    """Class to compare and validate metadata file with VCF file metadata"""

    def __init__(self, error_handler: Optional[ErrorHandler] = None):
        self.error_handler = error_handler or ErrorHandler()
        self.metadata_file_info: Optional[MetadataInfo] = None

    def parse_metadata_file(self, metadata_file_path: str) -> MetadataInfo:
        """Parse metadata file"""
        if not os.path.exists(metadata_file_path):
            self.error_handler.create_error(
                ErrorCode.FILE_NOT_FOUND,
                additional_info={"file_path": metadata_file_path, "file_type": "metadata"}
            )
            raise FileNotFoundError(f"Metadata file not found: {metadata_file_path}")

        try:
            with open(metadata_file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FILE_READ_ERROR,
                additional_info={"file_path": metadata_file_path, "error": str(e)}
            )
            raise

        metadata_info = MetadataInfo()

        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or not line.startswith('##'):
                continue

            # Parse ##key=value format
            if '=' not in line:
                self.error_handler.create_error(
                    ErrorCode.INVALID_METADATA_FORMAT,
                    line_number=line_num,
                    line_content=line,
                    additional_info={"file_type": "metadata"}
                )
                continue

            key, value = line[2:].split('=', 1)
            key = key.lower().strip()
            value = value.strip()

            if key == 'sampleset_id':
                metadata_info.sampleset_id = value
            elif key == 'experiment_id':
                metadata_info.experiment_id = value
            elif key == 'reference':
                metadata_info.reference = value

        self.metadata_file_info = metadata_info
        return metadata_info

    def validate_against_vcf(
        self,
        vcf_sampleset: Optional[str],
        vcf_experiment: Optional[str],
        vcf_reference: Optional[str]
    ) -> bool:
        """Compare and validate VCF file metadata with metadata file"""
        if self.metadata_file_info is None:
            self.error_handler.create_error(
                ErrorCode.MISSING_REQUIRED_METADATA,
                additional_info={"message": "Metadata file has not been parsed"}
            )
            return False

        validation_passed = True

        # Validate SAMPLESET
        if self.metadata_file_info.sampleset_id is not None:
            if vcf_sampleset is None:
                self.error_handler.create_error(
                    ErrorCode.MISSING_REQUIRED_METADATA,
                    field_name="SAMPLESET",
                    expected_value=self.metadata_file_info.sampleset_id,
                    actual_value="Not in VCF file",
                    additional_info={"source": "VCF file"}
                )
                validation_passed = False
            elif self._normalize_metadata_value(vcf_sampleset) != self._normalize_metadata_value(self.metadata_file_info.sampleset_id):
                self.error_handler.create_error(
                    ErrorCode.METADATA_MISMATCH,
                    field_name="SAMPLESET",
                    expected_value=self.metadata_file_info.sampleset_id,
                    actual_value=vcf_sampleset,
                    additional_info={"source": "VCF file"}
                )
                validation_passed = False

        # Validate EXPERIMENT
        if self.metadata_file_info.experiment_id is not None:
            if vcf_experiment is None:
                self.error_handler.create_error(
                    ErrorCode.MISSING_REQUIRED_METADATA,
                    field_name="EXPERIMENT",
                    expected_value=self.metadata_file_info.experiment_id,
                    actual_value="Not in VCF file",
                    additional_info={"source": "VCF file"}
                )
                validation_passed = False
            elif self._normalize_metadata_value(vcf_experiment) != self._normalize_metadata_value(self.metadata_file_info.experiment_id):
                self.error_handler.create_error(
                    ErrorCode.METADATA_MISMATCH,
                    field_name="EXPERIMENT",
                    expected_value=self.metadata_file_info.experiment_id,
                    actual_value=vcf_experiment,
                    additional_info={"source": "VCF file"}
                )
                validation_passed = False

        # Validate reference
        if self.metadata_file_info.reference is not None:
            if vcf_reference is None:
                self.error_handler.create_error(
                    ErrorCode.MISSING_REQUIRED_METADATA,
                    field_name="reference",
                    expected_value=self.metadata_file_info.reference,
                    actual_value="Not in VCF file",
                    additional_info={"source": "VCF file"}
                )
                validation_passed = False
            elif self._normalize_metadata_value(vcf_reference) != self._normalize_metadata_value(self.metadata_file_info.reference):
                self.error_handler.create_error(
                    ErrorCode.METADATA_MISMATCH,
                    field_name="reference",
                    expected_value=self.metadata_file_info.reference,
                    actual_value=vcf_reference,
                    additional_info={"source": "VCF file"}
                )
                validation_passed = False

        return validation_passed

    @staticmethod
    def _normalize_metadata_value(value: Optional[str]) -> Optional[str]:
        if value is None:
            return None
        return str(value).strip()
