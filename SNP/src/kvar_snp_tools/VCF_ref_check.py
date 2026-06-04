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
from typing import Dict, List, Optional, Any, Set
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
        self.mismatch_cases: List[Dict[str, Any]] = []
        self.chromosome_mapping: Dict[str, str] = {}  # VCF chromosome name -> FASTA chromosome name mapping
    
    def check_vcf_against_fasta(
        self,
        vcf_file_path: str,
        fasta_file_path: str,
        output_report_path: Optional[str] = None
    ) -> None:
        """Main method: compare VCF REF sequence with reference genome"""
        # Load FASTA file
        try:
            self._load_fasta(fasta_file_path)
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FASTA_READ_ERROR,
                additional_info={"file_path": fasta_file_path, "error": str(e)}
            )
            raise
        
        # Parse and validate VCF file
        try:
            self._check_vcf_file(vcf_file_path)
        except Exception as e:
            if self.error_handler.has_critical_errors():
                print("Validation stopped due to critical errors.")
                if output_report_path:
                    self._generate_report(output_report_path, vcf_file_path, fasta_file_path)
                raise
        
        # Generate report
        if output_report_path is None:
            output_report_path = vcf_file_path.replace('.vcf', '_ref_check_report.txt')
            if vcf_file_path.endswith('.gz'):
                output_report_path = vcf_file_path.replace('.vcf.gz', '_ref_check_report.txt')
        
        self._generate_report(output_report_path, vcf_file_path, fasta_file_path)
        
        # Print summary
        print(f"\n=== REF validation complete ===")
        print(f"Total variants: {self.stats['total_variants']}")
        print(f"  Matched: {self.stats['matched']}")
        print(f"  Mismatched: {self.stats['mismatched']}")
        print(f"  Missing chromosome: {self.stats['missing_chrom']}")
        print(f"  Out of range: {self.stats['out_of_range']}")
        print(f"  Skipped: {self.stats['skipped']}")
        print(f"\nReport file: {output_report_path}")
        
        # Print error summary
        self.error_handler.print_summary()

        self.error_handler.assert_no_blocking_errors(
            stage="VCF reference check"
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
            self.error_handler.create_error(
                ErrorCode.FILE_NOT_FOUND,
                additional_info={"file_path": vcf_file_path, "file_type": "VCF"}
            )
            raise FileNotFoundError(f"VCF file not found: {vcf_file_path}")
        
        # Check if gzipped
        is_gzipped = vcf_file_path.endswith('.gz')
        
        try:
            if is_gzipped:
                f = gzip.open(vcf_file_path, 'rt', encoding='utf-8')
            else:
                f = open(vcf_file_path, 'r', encoding='utf-8')
            
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
            
            f.close()
            print(f"VCF file parsing complete: {self.stats['total_variants']} variants validated")
            
        except UnicodeDecodeError as e:
            self.error_handler.create_error(
                ErrorCode.FILE_ENCODING_ERROR,
                additional_info={"file_path": vcf_file_path, "error": str(e)}
            )
            raise
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FILE_READ_ERROR,
                additional_info={"file_path": vcf_file_path, "error": str(e)}
            )
            raise
    
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
            mismatch_info = {
                'chrom': chrom,
                'pos': pos,
                'vcf_ref': ref,
                'fasta_seq': fasta_seq,
                'line_number': line_number
            }
            self.mismatch_cases.append(mismatch_info)
            
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
    
    def _generate_report(
        self,
        output_path: str,
        vcf_file_path: str,
        fasta_file_path: str
    ) -> None:
        """Generate validation result report"""
        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("VCF REF Sequence Validation Report")
        report_lines.append("=" * 80)
        report_lines.append("")
        
        # File information
        report_lines.append("=== Input Files ===")
        report_lines.append(f"VCF file: {os.path.abspath(vcf_file_path)}")
        report_lines.append(f"FASTA file: {os.path.abspath(fasta_file_path)}")
        report_lines.append("")
        
        # Statistics
        report_lines.append("=== Validation Statistics ===")
        report_lines.append(f"Total variants: {self.stats['total_variants']:,}")
        report_lines.append(f"  Matched: {self.stats['matched']:,} ({self.stats['matched']/max(self.stats['total_variants'], 1)*100:.2f}%)")
        report_lines.append(f"  Mismatched: {self.stats['mismatched']:,} ({self.stats['mismatched']/max(self.stats['total_variants'], 1)*100:.2f}%)")
        report_lines.append(f"  Missing chromosome: {self.stats['missing_chrom']:,}")
        report_lines.append(f"  Out of range: {self.stats['out_of_range']:,}")
        report_lines.append(f"  Skipped: {self.stats['skipped']:,}")
        report_lines.append("")
        
        # Chromosome mapping
        if self.chromosome_mapping:
            report_lines.append("=== Chromosome Name Mapping ===")
            for vcf_chrom, fasta_chrom in sorted(self.chromosome_mapping.items()):
                if vcf_chrom != fasta_chrom:
                    report_lines.append(f"  {vcf_chrom} -> {fasta_chrom}")
            report_lines.append("")
        
        # Mismatch details (max 100)
        if self.mismatch_cases:
            report_lines.append("=== REF Mismatch Details ===")
            report_lines.append(f"(Showing up to 100 of {len(self.mismatch_cases)} total)")
            report_lines.append("")
            
            for i, mismatch in enumerate(self.mismatch_cases[:100], 1):
                report_lines.append(f"[{i}] {mismatch['chrom']}:{mismatch['pos']}")
                report_lines.append(f"  Line number: {mismatch['line_number']}")
                report_lines.append(f"  VCF REF: {mismatch['vcf_ref']}")
                report_lines.append(f"  FASTA:   {mismatch['fasta_seq']}")
                report_lines.append("")
            
            if len(self.mismatch_cases) > 100:
                report_lines.append(f"... ({len(self.mismatch_cases) - 100} more omitted)")
                report_lines.append("")
        
        # Error report section
        if self.error_handler.has_errors():
            report_lines.append("=" * 80)
            report_lines.append("Error Details")
            report_lines.append("=" * 80)
            report_lines.append("")
            
            error_report = self.error_handler.generate_report(
                output_file=None,
                vcf_file_path=vcf_file_path,
                output_tsv_path=None
            )
            report_lines.append(error_report)
        
        # Save report file
        report_text = "\n".join(report_lines)
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(report_text)
            print(f"Report file created: {output_path}")
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FILE_WRITE_ERROR,
                additional_info={"file_path": output_path, "error": str(e)}
            )
            print(f"Warning: Failed to save report file: {e}")


def main():
    """Main function - command-line argument handling"""
    parser = argparse.ArgumentParser(
        description='Validate that VCF REF sequence matches reference genome (FASTA)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python VCF_ref_check.py -v input.vcf -f reference.fasta -o report.txt
  python VCF_ref_check.py -v input.vcf.gz -f reference.fasta -o report.txt

Notes:
  - pyfaidx library is required: pip install pyfaidx
  - FASTA index file (.fai) is created automatically
        """
    )
    
    parser.add_argument(
        '-v', '--vcf',
        dest='vcf_file',
        required=True,
        help='Input VCF file path (required, .gz supported)'
    )
    
    parser.add_argument(
        '-f', '--fasta',
        dest='fasta_file',
        required=True,
        help='Reference genome FASTA file path (required)'
    )
    
    parser.add_argument(
        '-o', '--output',
        dest='output_report',
        required=False,
        help='Output report file path (optional, default: VCF_filename_ref_check_report.txt)'
    )
    
    args = parser.parse_args()
    
    # Check file existence
    if not os.path.exists(args.vcf_file):
        print(f"Error: VCF file not found: {args.vcf_file}")
        sys.exit(1)
    
    if not os.path.exists(args.fasta_file):
        print(f"Error: FASTA file not found: {args.fasta_file}")
        sys.exit(1)
    
    # Check pyfaidx library
    if pyfaidx is None:
        print("Error: pyfaidx library is not installed.")
        print("Install with: pip install pyfaidx")
        sys.exit(1)
    
    # Run validation
    checker = VCFRefChecker()
    
    try:
        checker.check_vcf_against_fasta(
            args.vcf_file,
            args.fasta_file,
            args.output_report
        )

        sys.exit(0)

    except RuntimeError as e:
        print(f"\nValidation blocked: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\nError during validation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
