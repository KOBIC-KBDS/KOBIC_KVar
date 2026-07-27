#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metadata validation module.
Validates by comparing metadata file with VCF file metadata.
"""

import os
from typing import Dict, Optional, List
from dataclasses import dataclass

try:
    from .error_handler import ErrorHandler, ErrorCode, ErrorSeverity, ErrorCategory
except ImportError:
    from error_handler import ErrorHandler, ErrorCode, ErrorSeverity, ErrorCategory


@dataclass
class MetadataInfo:
    """Class to store metadata information"""
    experiment_id: Optional[str] = None  # Used instead of batch
    reference: Optional[str] = None
    sampleset_ids: Optional[List[str]] = None  # Used instead of population_id


class MetadataValidator:
    """Class to compare and validate metadata file with VCF file metadata"""
    
    def __init__(self, error_handler: Optional[ErrorHandler] = None):
        self.error_handler = error_handler or ErrorHandler()
        self.metadata_file_info: Optional[MetadataInfo] = None

    @staticmethod
    def _is_missing(value: Optional[str]) -> bool:
        """Return True when a metadata value is absent or blank."""
        return value is None or str(value).strip() == ""
    
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
            
            # Parse Experiment_id (instead of batch)
            if key == 'experiment_id':
                metadata_info.experiment_id = value
            elif key == 'reference':
                metadata_info.reference = value
            elif key == 'sampleset_id':
                if metadata_info.sampleset_ids is not None:
                    self.error_handler.create_error(
                        ErrorCode.DUPLICATE_METADATA_TAG,
                        line_number=line_num,
                        line_content=line,
                        field_name="SampleSet_id",
                        expected_value="exactly one ##SampleSet_id entry",
                        actual_value="repeated ##SampleSet_id entry",
                    )
                    continue
                metadata_info.sampleset_ids = [value]

        self._validate_required_fields(metadata_info)
        self.metadata_file_info = metadata_info
        return metadata_info

    def _validate_required_fields(self, metadata_info: MetadataInfo) -> None:
        """Report missing identifiers required for dbSNP VCF headers."""
        if self._is_missing(metadata_info.experiment_id):
            self.error_handler.create_error(
                ErrorCode.MISSING_REQUIRED_METADATA,
                field_name="Experiment_id",
                expected_value="##Experiment_id=<non-empty value>",
                actual_value="Not in metadata file",
            )

        if (
            len(metadata_info.sampleset_ids or []) != 1
            or self._is_missing(metadata_info.sampleset_ids[0])
        ):
            self.error_handler.create_error(
                ErrorCode.MISSING_REQUIRED_METADATA,
                field_name="SampleSet_id",
                expected_value="Exactly one ##SampleSet_id=<non-empty value>",
                actual_value=(
                    "Not in metadata file"
                    if not metadata_info.sampleset_ids
                    else f"{len(metadata_info.sampleset_ids)} entries"
                ),
            )
    
    def validate_against_vcf(
        self,
        vcf_metadata: Dict[str, any],
        vcf_population_ids: List[str]
    ) -> bool:
        """Compare and validate VCF file metadata with metadata file"""
        if self.metadata_file_info is None:
            self.error_handler.create_error(
                ErrorCode.METADATA_NOT_PARSED,
                additional_info={"message": "Metadata file has not been parsed"}
            )
            return False
        
        validation_passed = True
        
        # Validate Experiment_id (maps to VCF batch)
        if self.metadata_file_info.experiment_id is not None:
            vcf_batch = vcf_metadata.get('batch')
            if self._is_missing(vcf_batch):
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_FILLED,
                    field_name="batch",
                    expected_value=self.metadata_file_info.experiment_id,
                    actual_value="Not in VCF file",
                    additional_info={
                        "source": "VCF file",
                        "action": "Output VCF uses metadata experiment_id as batch"
                    }
                )
            elif str(vcf_batch).strip() != self.metadata_file_info.experiment_id:
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_CORRECTED,
                    field_name="batch/Experiment_id",
                    expected_value=self.metadata_file_info.experiment_id,
                    actual_value=vcf_batch,
                    additional_info={
                        "source": "VCF file",
                        "action": "Output VCF uses metadata experiment_id as batch"
                    }
                )
        
        # Validate reference
        if self.metadata_file_info.reference is not None:
            vcf_reference = vcf_metadata.get('reference')
            if self._is_missing(vcf_reference):
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_FILLED,
                    field_name="reference",
                    expected_value=self.metadata_file_info.reference,
                    actual_value="Not in VCF file",
                    additional_info={
                        "source": "VCF file",
                        "action": "Output VCF uses metadata reference"
                    }
                )
            elif str(vcf_reference).strip() != self.metadata_file_info.reference:
                self.error_handler.create_error(
                    ErrorCode.METADATA_REFERENCE_MISMATCH,
                    field_name="reference",
                    expected_value=self.metadata_file_info.reference,
                    actual_value=vcf_reference,
                    additional_info={"source": "VCF file"}
                )
                validation_passed = False
        
        # Validate SampleSet_id (maps to VCF population_id)
        metadata_sampleset_ids = self.metadata_file_info.sampleset_ids or []
        if (
            len(metadata_sampleset_ids) != 1
            or len(vcf_population_ids) != 1
            or self._is_missing(metadata_sampleset_ids[0])
            or self._is_missing(vcf_population_ids[0])
            or metadata_sampleset_ids[0] != vcf_population_ids[0]
        ):
            self.error_handler.create_error(
                ErrorCode.METADATA_POPULATION_MISMATCH,
                field_name="SampleSet_id/population_id",
                expected_value=(
                    metadata_sampleset_ids[0]
                    if len(metadata_sampleset_ids) == 1
                    else "exactly one metadata SampleSet_id"
                ),
                actual_value=(
                    vcf_population_ids[0]
                    if len(vcf_population_ids) == 1
                    else f"{len(vcf_population_ids)} VCF population_id entries"
                ),
                additional_info={
                    "source": "VCF file",
                    "note": (
                        "Metadata SampleSet_id and VCF population_id must each "
                        "contain exactly one identical value"
                    ),
                },
            )
            validation_passed = False
        
        return validation_passed
