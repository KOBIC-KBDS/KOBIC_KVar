#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metadata validation module.
Validates by comparing metadata file with VCF file metadata.
"""

import os
import re
from typing import Dict, Optional, List, Set
from dataclasses import dataclass

try:
    from .error_handler import ErrorHandler, ErrorCode, ErrorSeverity, ErrorCategory
except ImportError:
    from error_handler import ErrorHandler, ErrorCode, ErrorSeverity, ErrorCategory


@dataclass
class MetadataInfo:
    """Class to store metadata information"""
    experiment_id: Optional[str] = None  # Used instead of batch
    bioproject_id: Optional[str] = None
    biosample_id: Optional[List[str]] = None
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
            elif key == 'bioproject_id':
                metadata_info.bioproject_id = value
            elif key == 'biosample_id':
                # Handle comma-, space-, or tab-separated values
                biosample_ids = re.split(r'[,\s]+', value)
                metadata_info.biosample_id = [v.strip() for v in biosample_ids if v.strip()]
            elif key == 'reference':
                metadata_info.reference = value
            elif key == 'sampleset_id':
                # Parse SampleSet_id (instead of population_id)
                if metadata_info.sampleset_ids is None:
                    metadata_info.sampleset_ids = []
                metadata_info.sampleset_ids.append(value)
        
        self._validate_distinct_sample_identifiers(metadata_info)
        self.metadata_file_info = metadata_info
        return metadata_info

    def _validate_distinct_sample_identifiers(self, metadata_info: MetadataInfo) -> None:
        """Report an error when biosample_id and SampleSet_id reuse the same ID."""
        if not metadata_info.biosample_id or not metadata_info.sampleset_ids:
            return

        biosample_ids = set(metadata_info.biosample_id)
        sampleset_ids = set(metadata_info.sampleset_ids)
        duplicated_ids = sorted(biosample_ids & sampleset_ids)
        if not duplicated_ids:
            return

        self.error_handler.create_error(
            ErrorCode.METADATA_IDENTIFIER_CONFLICT,
            field_name="biosample_id/SampleSet_id",
            expected_value="Distinct biosample_id and SampleSet_id values",
            actual_value=", ".join(duplicated_ids),
            additional_info={
                "duplicated_ids": duplicated_ids,
                "message": "biosample_id and SampleSet_id must not use the same identifier"
            }
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
        
        # Validate bioproject_id
        if self.metadata_file_info.bioproject_id is not None:
            vcf_bioproject_id = vcf_metadata.get('bioproject_id')
            if self._is_missing(vcf_bioproject_id):
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_FILLED,
                    field_name="bioproject_id",
                    expected_value=self.metadata_file_info.bioproject_id,
                    actual_value="Not in VCF file",
                    additional_info={
                        "source": "VCF file",
                        "action": "Output VCF uses metadata bioproject_id"
                    }
                )
            elif str(vcf_bioproject_id).strip() != self.metadata_file_info.bioproject_id:
                self.error_handler.create_error(
                    ErrorCode.METADATA_VALUE_CORRECTED,
                    field_name="bioproject_id",
                    expected_value=self.metadata_file_info.bioproject_id,
                    actual_value=vcf_bioproject_id,
                    additional_info={
                        "source": "VCF file",
                        "action": "Output VCF uses metadata bioproject_id"
                    }
                )
        
        # Validate biosample_id
        if self.metadata_file_info.biosample_id is not None:
            vcf_biosample_id_str = vcf_metadata.get('biosample_id')
            if self._is_missing(vcf_biosample_id_str):
                self.error_handler.create_error(
                    ErrorCode.MISSING_BIOSAMPLE_ID,
                    field_name="biosample_id",
                    expected_value=", ".join(self.metadata_file_info.biosample_id),
                    actual_value="Not in VCF file",
                    additional_info={"source": "VCF file"}
                )
                validation_passed = False
            else:
                # Parse comma-, space-, or tab-separated values
                vcf_biosample_ids = re.split(r'[,\s]+', vcf_biosample_id_str)
                vcf_biosample_ids = [v.strip() for v in vcf_biosample_ids if v.strip()]
                metadata_biosample_ids = set(self.metadata_file_info.biosample_id)
                vcf_biosample_ids_set = set(vcf_biosample_ids)
                
                # IDs in metadata file but not in VCF
                missing_in_vcf = metadata_biosample_ids - vcf_biosample_ids_set
                if missing_in_vcf:
                    self.error_handler.create_error(
                        ErrorCode.METADATA_BIOSAMPLE_MISMATCH,
                        field_name="biosample_id",
                        expected_value=", ".join(sorted(metadata_biosample_ids)),
                        actual_value=", ".join(sorted(vcf_biosample_ids_set)),
                        additional_info={
                            "source": "VCF file",
                            "missing_in_vcf": sorted(missing_in_vcf),
                            "extra_in_vcf": sorted(vcf_biosample_ids_set - metadata_biosample_ids)
                        }
                    )
                    validation_passed = False
        
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
        if self.metadata_file_info.sampleset_ids is not None:
            metadata_sampleset_ids = set(self.metadata_file_info.sampleset_ids)
            vcf_pop_ids_set = set(vcf_population_ids)
            
            # SampleSet_id in metadata file but not in VCF population_id
            missing_in_vcf = metadata_sampleset_ids - vcf_pop_ids_set
            if missing_in_vcf:
                self.error_handler.create_error(
                    ErrorCode.METADATA_POPULATION_MISMATCH,
                    field_name="SampleSet_id/population_id",
                    expected_value=", ".join(sorted(metadata_sampleset_ids)),
                    actual_value=", ".join(sorted(vcf_pop_ids_set)) if vcf_pop_ids_set else "Not in VCF file",
                    additional_info={
                        "source": "VCF file",
                        "missing_in_vcf": sorted(missing_in_vcf),
                        "extra_in_vcf": sorted(vcf_pop_ids_set - metadata_sampleset_ids),
                        "note": "Comparing metadata file SampleSet_id with VCF file population_id"
                    }
                )
                validation_passed = False
            
            # population_id in VCF but not in metadata file SampleSet_id
            extra_in_vcf = vcf_pop_ids_set - metadata_sampleset_ids
            if extra_in_vcf:
                self.error_handler.create_error(
                    ErrorCode.METADATA_POPULATION_MISMATCH,
                    field_name="SampleSet_id/population_id",
                    expected_value=", ".join(sorted(metadata_sampleset_ids)),
                    actual_value=", ".join(sorted(vcf_pop_ids_set)),
                    additional_info={
                        "source": "VCF file",
                        "extra_in_vcf": sorted(extra_in_vcf),
                        "note": "Comparing metadata file SampleSet_id with VCF file population_id"
                    }
                )
                validation_passed = False
        
        return validation_passed
