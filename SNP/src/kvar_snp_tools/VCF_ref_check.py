#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF REF sequence validation module.
Validates that VCF REF sequence matches the reference genome (FASTA).
"""

import sys
import os
import gzip
import argparse
from typing import Dict, Optional, Set
from collections import defaultdict

# Use pyfaidx library (FASTA file indexing)
try:
    import pyfaidx
except ImportError:
    print("Warning: pyfaidx library is not installed.")
    print("Install with: pip install pyfaidx")
    pyfaidx = None

# Relative path import support
if __name__ == "__main__":
    from error_handler import ErrorHandler, ErrorCode
else:
    from .error_handler import ErrorHandler, ErrorCode


class VCFRefChecker:
    """Class to validate VCF REF sequence against reference genome"""

    def __init__(self, error_handler: Optional[ErrorHandler] = None):
        self.error_handler = error_handler or ErrorHandler()
        self.fasta_handler = None
        self.stats = {
            'total_variants': 0,
            'matched': 0,
            'mismatched': 0,
            'missing_chrom': 0,
            'out_of_range': 0,
            'skipped': 0
        }
        self.chromosome_mapping: Dict[str, str] = {}  # VCF chromosome name -> FASTA chromosome name mapping

    def check_vcf_against_fasta(
        self,
        vcf_file_path: str,
        fasta_file_path: str,
    ) -> None:
        """Collect optional FASTA validation results in the shared error handler."""
        try:
            self._load_fasta(fasta_file_path)
        except Exception as e:
            if not self.error_handler.has_critical_errors():
                self.error_handler.create_error(
                    ErrorCode.FASTA_READ_ERROR,
                    additional_info={"file_path": fasta_file_path, "error": str(e)}
                )
            self._store_result(
                vcf_file_path,
                fasta_file_path,
                status="FAILED",
            )
            return

        try:
            self._check_vcf_file(vcf_file_path)
        except Exception:
            self._store_result(
                vcf_file_path,
                fasta_file_path,
                status="FAILED",
            )
            return

        self._store_result(
            vcf_file_path,
            fasta_file_path,
            status="COMPLETED",
        )

        print(f"\n=== REF validation complete ===")
        print(f"Total variants: {self.stats['total_variants']}")
        print(f"  Matched: {self.stats['matched']}")
        print(f"  Mismatched: {self.stats['mismatched']}")
        print(f"  Missing chromosome: {self.stats['missing_chrom']}")
        print(f"  Out of range: {self.stats['out_of_range']}")
        print(f"  Skipped: {self.stats['skipped']}")

    def _store_result(
        self,
        vcf_file_path: str,
        fasta_file_path: str,
        *,
        status: str,
    ) -> None:
        """Attach reference statistics to the unified validation report."""
        self.error_handler.set_reference_validation(
            vcf_file_path=vcf_file_path,
            fasta_file_path=fasta_file_path,
            status=status,
            stats=self.stats,
            chromosome_mapping=self.chromosome_mapping,
        )

    def _load_fasta(self, fasta_file_path: str) -> None:
        """Load and index FASTA file"""
        if pyfaidx is None:
            raise ImportError("pyfaidx library is required. Install with: pip install pyfaidx")

        if not os.path.exists(fasta_file_path):
            self.error_handler.create_error(
                ErrorCode.FILE_NOT_FOUND,
                additional_info={"file_path": fasta_file_path, "file_type": "FASTA"}
            )
            raise FileNotFoundError(f"FASTA file not found: {fasta_file_path}")

        try:
            # Index FASTA file (.fai is created automatically)
            self.fasta_handler = pyfaidx.Fasta(fasta_file_path, sequence_always_upper=True)
            print(f"FASTA file loaded: {fasta_file_path}")
            print(f"  Number of chromosomes: {len(self.fasta_handler.keys())}")

            # Print chromosome list (first 10 only)
            chroms = list(self.fasta_handler.keys())[:10]
            print(f"  Chromosome examples: {', '.join(chroms)}")
            if len(self.fasta_handler.keys()) > 10:
                print(f"  ... (total {len(self.fasta_handler.keys())})")
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FASTA_INDEX_ERROR,
                additional_info={"file_path": fasta_file_path, "error": str(e)}
            )
            raise

    def _check_vcf_file(self, vcf_file_path: str) -> None:
        """Parse VCF file and validate REF for each variant"""
        if not os.path.exists(vcf_file_path):
            # The main VCF validator owns input-file errors in the unified report.
            raise FileNotFoundError(f"VCF file not found: {vcf_file_path}")

        # Check if gzipped
        is_gzipped = vcf_file_path.endswith('.gz')

        opener = gzip.open if is_gzipped else open
        with opener(vcf_file_path, 'rt', encoding='utf-8') as f:
            line_number = 0
            batch_size = 10000
            processed_count = 0

            for line in f:
                line_number += 1
                line = line.strip()

                # Skip header lines
                if line.startswith('#'):
                    continue

                if not line:
                    continue

                # Parse and validate data line
                self._check_variant_line(line, line_number)
                self.stats['total_variants'] += 1
                processed_count += 1

                # Progress output
                if processed_count % batch_size == 0:
                    print(f"  {processed_count:,} variants validated...")

        print(f"VCF file parsing complete: {self.stats['total_variants']} variants validated")

    def _check_variant_line(self, line: str, line_number: int) -> None:
        """Validate REF for a single VCF data line"""
        fields = line.split('\t')

        # Check minimum required fields
        if len(fields) < 5:
            self.stats['skipped'] += 1
            return

        chrom = fields[0].strip()
        pos_str = fields[1].strip()
        ref = fields[3].strip()

        # Parse position
        try:
            pos = int(pos_str)
            if pos < 1:
                self.stats['skipped'] += 1
                return
        except ValueError:
            self.stats['skipped'] += 1
            return

        # Skip if REF is empty or '.'
        if not ref or ref == '.':
            self.stats['skipped'] += 1
            return

        # Validate REF
        self._check_ref_sequence(chrom, pos, ref, line_number)

    def _check_ref_sequence(self, chrom: str, pos: int, ref: str, line_number: int) -> None:
        """Compare REF sequence with reference genome"""
        # Normalize and map chromosome name
        fasta_chrom = self._get_fasta_chromosome_name(chrom)

        if fasta_chrom is None:
            # Chromosome not found
            self.stats['missing_chrom'] += 1
            self.error_handler.create_error(
                ErrorCode.CHROMOSOME_NOT_FOUND,
                line_number=line_number,
                field_name="CHROM",
                actual_value=chrom,
                additional_info={
                    "available_chromosomes": list(self.fasta_handler.keys())[:20],
                    "message": f"VCF chromosome name '{chrom}' not found in reference genome"
                }
            )
            return

        # Check position range
        chrom_length = len(self.fasta_handler[fasta_chrom])
        ref_length = len(ref)
        end_pos = pos + ref_length - 1

        if pos < 1 or end_pos > chrom_length:
            self.stats['out_of_range'] += 1
            self.error_handler.create_error(
                ErrorCode.POSITION_OUT_OF_RANGE,
                line_number=line_number,
                field_name=f"{chrom}:{pos}",
                actual_value=f"position {pos}-{end_pos}",
                expected_value=f"1-{chrom_length}",
                additional_info={
                    "chromosome": chrom,
                    "position": pos,
                    "ref_length": ref_length,
                    "chromosome_length": chrom_length,
                    "end_position": end_pos
                }
            )
            return

        # Extract sequence from FASTA (1-based coordinates)
        try:
            # pyfaidx uses 1-based coordinates
            fasta_seq = str(self.fasta_handler[fasta_chrom][pos-1:end_pos])
        except Exception as e:
            self.stats['out_of_range'] += 1
            self.error_handler.create_error(
                ErrorCode.POSITION_OUT_OF_RANGE,
                line_number=line_number,
                field_name=f"{chrom}:{pos}",
                additional_info={"error": str(e)}
            )
            return

        # Compare REF with FASTA sequence
        ref_upper = ref.upper()
        fasta_seq_upper = fasta_seq.upper()

        if ref_upper == fasta_seq_upper:
            # Match
            self.stats['matched'] += 1
        else:
            # Mismatch
            self.stats['mismatched'] += 1
            self.error_handler.create_error(
                ErrorCode.REF_MISMATCH,
                line_number=line_number,
                field_name=f"{chrom}:{pos}",
                expected_value=fasta_seq_upper,
                actual_value=ref_upper,
                additional_info={
                    "chromosome": chrom,
                    "position": pos,
                    "vcf_ref": ref,
                    "fasta_sequence": fasta_seq,
                    "ref_length": len(ref),
                    "fasta_length": len(fasta_seq)
                }
            )

    def _get_fasta_chromosome_name(self, vcf_chrom: str) -> Optional[str]:
        """Convert VCF chromosome name to FASTA chromosome name"""
        # Return if already mapped
        if vcf_chrom in self.chromosome_mapping:
            return self.chromosome_mapping[vcf_chrom]

        # Try direct match
        if vcf_chrom in self.fasta_handler:
            self.chromosome_mapping[vcf_chrom] = vcf_chrom
            return vcf_chrom

        # Try adding/removing chr prefix
        if vcf_chrom.startswith('chr'):
            # chr1 -> 1
            chrom_without_chr = vcf_chrom[3:]
            if chrom_without_chr in self.fasta_handler:
                self.chromosome_mapping[vcf_chrom] = chrom_without_chr
                return chrom_without_chr
        else:
            # 1 -> chr1
            chrom_with_chr = 'chr' + vcf_chrom
            if chrom_with_chr in self.fasta_handler:
                self.chromosome_mapping[vcf_chrom] = chrom_with_chr
                return chrom_with_chr

        # Try case-insensitive match
        vcf_chrom_lower = vcf_chrom.lower()
        for fasta_chrom in self.fasta_handler.keys():
            if fasta_chrom.lower() == vcf_chrom_lower:
                self.chromosome_mapping[vcf_chrom] = fasta_chrom
                return fasta_chrom

        # No match
        return None
