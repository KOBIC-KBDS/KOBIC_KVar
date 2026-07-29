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
import tempfile
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict

# Relative path import support
if __name__ == "__main__":
    from error_handler import ErrorHandler, ErrorCode
else:
    from .error_handler import ErrorHandler, ErrorCode


class IndexedFasta:
    """Small faidx-compatible FASTA reader using only the standard library."""

    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.index_path = f"{fasta_path}.fai"
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        self._fasta_handle = None

        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        if self._index_is_current():
            try:
                self.index = self._load_index()
            except (OSError, ValueError):
                self.index = self._build_index()
                self._cache_index()
        else:
            self.index = self._build_index()
            self._cache_index()

    def _index_is_current(self) -> bool:
        """Return whether the cached index is at least as new as the FASTA."""
        try:
            return (
                os.stat(self.index_path).st_mtime_ns
                >= os.stat(self.fasta_path).st_mtime_ns
            )
        except OSError:
            return False

    def _load_index(self) -> Dict[str, Tuple[int, int, int, int]]:
        """Load an index and verify that it describes the complete FASTA."""
        index: Dict[str, Tuple[int, int, int, int]] = {}
        fasta_size = os.path.getsize(self.fasta_path)
        with open(self.index_path, "r", encoding="utf-8") as fai:
            for line_number, line in enumerate(fai, 1):
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 5 or not parts[0]:
                    raise ValueError(
                        f"Invalid FASTA index line {line_number}: {line.rstrip()}"
                    )
                chrom = parts[0]
                if chrom in index:
                    raise ValueError(f"Duplicate FASTA index sequence name: {chrom}")
                try:
                    length, offset, line_bases, line_width = (
                        int(parts[1]),
                        int(parts[2]),
                        int(parts[3]),
                        int(parts[4]),
                    )
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid FASTA index line {line_number}: {line.rstrip()}"
                    ) from exc
                if (
                    length < 1
                    or offset < 0
                    or offset >= fasta_size
                    or line_bases < 1
                    or line_width < line_bases
                ):
                    raise ValueError(
                        f"Invalid FASTA index values on line {line_number}"
                    )
                index[chrom] = (length, offset, line_bases, line_width)

        if not index:
            raise ValueError(f"FASTA index is empty: {self.index_path}")

        entries = list(index.items())
        next_header_offset = 0
        with open(self.fasta_path, "rb") as fasta:
            for entry_number, (chrom, values) in enumerate(entries):
                fasta.seek(next_header_offset)
                header = fasta.readline()
                if not header.startswith(b">"):
                    raise ValueError(
                        f"FASTA index record does not follow a header: {chrom}"
                    )
                try:
                    header_name = (
                        header[1:].strip().split(None, 1)[0].decode("utf-8")
                    )
                except (IndexError, UnicodeDecodeError) as exc:
                    raise ValueError("Invalid FASTA sequence header") from exc
                if header_name != chrom or fasta.tell() != values[1]:
                    raise ValueError(
                        f"FASTA index order or sequence offset is invalid for {chrom}"
                    )

                next_header_offset = self._validate_cached_record(
                    fasta,
                    fasta_size,
                    chrom,
                    values,
                )
                is_last_entry = entry_number == len(entries) - 1
                if is_last_entry and next_header_offset != fasta_size:
                    raise ValueError(
                        f"FASTA contains records absent from the index after {chrom}"
                    )
                if not is_last_entry and next_header_offset >= fasta_size:
                    raise ValueError(
                        f"FASTA index contains records absent from the FASTA after {chrom}"
                    )
        return index

    @staticmethod
    def _header_name_before_offset(fasta, sequence_offset: int) -> Optional[str]:
        """Return the FASTA sequence name immediately before an indexed offset."""
        if sequence_offset < 2:
            return None

        end = sequence_offset - 1
        fasta.seek(end)
        if fasta.read(1) != b"\n":
            return None
        end -= 1
        if end >= 0:
            fasta.seek(end)
            if fasta.read(1) == b"\r":
                end -= 1

        chunks: List[bytes] = []
        cursor = end + 1
        while cursor > 0:
            chunk_start = max(0, cursor - 4096)
            fasta.seek(chunk_start)
            chunk = fasta.read(cursor - chunk_start)
            newline = chunk.rfind(b"\n")
            if newline >= 0:
                chunks.append(chunk[newline + 1:])
                break
            chunks.append(chunk)
            cursor = chunk_start

        header = b"".join(reversed(chunks)).strip()
        if not header.startswith(b">"):
            return None
        try:
            return header[1:].split(None, 1)[0].decode("utf-8")
        except (IndexError, UnicodeDecodeError):
            return None

    def _validate_cached_record(
        self,
        fasta,
        fasta_size: int,
        chrom: str,
        values: Tuple[int, int, int, int],
    ) -> int:
        """Validate one cached record and return the following header offset."""
        length, offset, line_bases, line_width = values
        if self._header_name_before_offset(fasta, offset) != chrom:
            raise ValueError(
                f"FASTA index header offset does not match sequence {chrom}"
            )

        fasta.seek(offset)
        first_line = fasta.readline()
        first_sequence_line = first_line.rstrip(b"\r\n")
        if (
            not first_sequence_line
            or len(first_sequence_line) != line_bases
            or len(first_line) != line_width
            or any(
                byte <= 32 or byte > 126
                for byte in first_sequence_line
            )
            or b">" in first_sequence_line
        ):
            raise ValueError(
                f"FASTA index line metrics are invalid for {chrom}"
            )

        last_base_offset = (
            offset
            + ((length - 1) // line_bases) * line_width
            + ((length - 1) % line_bases)
        )
        if last_base_offset >= fasta_size:
            raise ValueError(
                f"FASTA index length exceeds the sequence record for {chrom}"
            )

        fasta.seek(last_base_offset)
        last_base = fasta.read(1)
        if (
            not last_base
            or last_base in b" \t\r\n>"
            or last_base[0] > 127
        ):
            raise ValueError(
                f"FASTA index has an invalid last-base offset for {chrom}"
            )

        boundary = fasta.read(1)
        if boundary == b"\r":
            if fasta.read(1) != b"\n":
                raise ValueError(
                    f"FASTA index has an invalid record boundary for {chrom}"
                )
            boundary = fasta.read(1)
        elif boundary == b"\n":
            boundary = fasta.read(1)
        elif boundary:
            raise ValueError(
                f"FASTA index length truncates the sequence record for {chrom}"
            )

        if boundary not in {b"", b">"}:
            raise ValueError(
                f"FASTA index length does not end at the record boundary for {chrom}"
            )
        if boundary == b">":
            return fasta.tell() - 1
        return fasta.tell()

    def _build_index(self) -> Dict[str, Tuple[int, int, int, int]]:
        """Build an faidx-compatible index directly from the FASTA bytes."""
        index: Dict[str, Tuple[int, int, int, int]] = {}
        current_name: Optional[str] = None
        sequence_length = 0
        sequence_offset: Optional[int] = None
        line_bases: Optional[int] = None
        line_width: Optional[int] = None
        terminal_line_seen = False

        def finish_record() -> None:
            if current_name is None:
                return
            if (
                sequence_length < 1
                or sequence_offset is None
                or line_bases is None
                or line_width is None
            ):
                raise ValueError(f"FASTA sequence has no bases: {current_name}")
            index[current_name] = (
                sequence_length,
                sequence_offset,
                line_bases,
                line_width,
            )

        try:
            with open(self.fasta_path, "rb") as fasta:
                while True:
                    line_offset = fasta.tell()
                    raw_line = fasta.readline()
                    if not raw_line:
                        break

                    if raw_line.startswith(b">"):
                        finish_record()
                        header = raw_line[1:].strip()
                        if not header:
                            raise ValueError("FASTA header has no sequence name")
                        try:
                            current_name = header.split(None, 1)[0].decode("utf-8")
                        except UnicodeDecodeError as exc:
                            raise ValueError(
                                "FASTA sequence name is not valid UTF-8"
                            ) from exc
                        if current_name in index:
                            raise ValueError(
                                f"Duplicate FASTA sequence name: {current_name}"
                            )
                        sequence_length = 0
                        sequence_offset = None
                        line_bases = None
                        line_width = None
                        terminal_line_seen = False
                        continue

                    if current_name is None:
                        raise ValueError(
                            "FASTA sequence data appears before the first header"
                        )

                    sequence_line = raw_line.rstrip(b"\r\n")
                    if not sequence_line:
                        raise ValueError(
                            f"Blank FASTA sequence line for {current_name}"
                        )
                    if (
                        any(byte <= 32 or byte > 126 for byte in sequence_line)
                        or b">" in sequence_line
                    ):
                        raise ValueError(
                            f"Invalid FASTA sequence line for {current_name}"
                        )
                    if terminal_line_seen:
                        raise ValueError(
                            f"Inconsistent FASTA line wrapping for {current_name}"
                        )

                    bases_on_line = len(sequence_line)
                    if sequence_offset is None:
                        sequence_offset = line_offset
                        line_bases = bases_on_line
                        line_width = len(raw_line)
                    else:
                        assert line_bases is not None
                        assert line_width is not None
                        if bases_on_line > line_bases:
                            raise ValueError(
                                f"Inconsistent FASTA line wrapping for {current_name}"
                            )
                        if bases_on_line < line_bases:
                            terminal_line_seen = True
                        elif len(raw_line) != line_width:
                            if raw_line.endswith((b"\n", b"\r")):
                                raise ValueError(
                                    f"Inconsistent FASTA line endings for {current_name}"
                                )
                            terminal_line_seen = True

                    sequence_length += bases_on_line

                finish_record()
        except OSError as exc:
            raise OSError(f"Unable to read FASTA file: {self.fasta_path}") from exc

        if not index:
            raise ValueError(f"FASTA contains no sequences: {self.fasta_path}")
        return index

    def _cache_index(self) -> None:
        """Atomically cache the index, or continue with it in memory if unwritable."""
        index_dir = os.path.dirname(os.path.abspath(self.index_path)) or "."
        temporary_path: Optional[str] = None
        try:
            fd, temporary_path = tempfile.mkstemp(
                dir=index_dir,
                prefix=f".{os.path.basename(self.index_path)}.",
                suffix=".tmp",
            )
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fai:
                for chrom, values in self.index.items():
                    fai.write(
                        "\t".join([chrom, *(str(value) for value in values)])
                        + "\n"
                    )
                fai.flush()
                os.fsync(fai.fileno())
            os.replace(temporary_path, self.index_path)
            temporary_path = None
        except OSError:
            pass
        finally:
            if temporary_path:
                try:
                    os.unlink(temporary_path)
                except OSError:
                    pass

    def keys(self):
        return self.index.keys()

    def __contains__(self, chrom: str) -> bool:
        return chrom in self.index

    def length(self, chrom: str) -> int:
        return self.index[chrom][0]

    def close(self) -> None:
        """Close the shared FASTA handle used by repeated interval fetches."""
        if self._fasta_handle is not None:
            self._fasta_handle.close()
            self._fasta_handle = None

    def __del__(self):
        self.close()

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch a 1-based inclusive sequence."""
        length, offset, line_bases, line_width = self.index[chrom]
        if start < 1 or end > length or start > end:
            raise ValueError(f"{chrom}:{start}-{end} is outside reference bounds")

        if self._fasta_handle is None:
            self._fasta_handle = open(self.fasta_path, "rb")
        fasta = self._fasta_handle
        sequence_parts: List[str] = []
        position = start
        while position <= end:
            zero_based = position - 1
            line_index = zero_based // line_bases
            line_offset = zero_based % line_bases
            bases_to_read = min(
                end - position + 1,
                line_bases - line_offset,
            )
            byte_offset = offset + line_index * line_width + line_offset
            fasta.seek(byte_offset)
            sequence_parts.append(
                fasta.read(bases_to_read).decode("ascii")
            )
            position += bases_to_read

        return "".join(sequence_parts).upper()


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
            close_fasta = getattr(self.fasta_handler, "close", None)
            if close_fasta is not None:
                close_fasta()
            self._store_result(
                vcf_file_path,
                fasta_file_path,
                status="FAILED",
            )
            return

        close_fasta = getattr(self.fasta_handler, "close", None)
        if close_fasta is not None:
            close_fasta()
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
        if not os.path.exists(fasta_file_path):
            self.error_handler.create_error(
                ErrorCode.FILE_NOT_FOUND,
                additional_info={"file_path": fasta_file_path, "file_type": "FASTA"}
            )
            raise FileNotFoundError(f"FASTA file not found: {fasta_file_path}")

        try:
            self.fasta_handler = IndexedFasta(fasta_file_path)
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
        chrom_length = self.fasta_handler.length(fasta_chrom)
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

        # Extract sequence from FASTA (1-based inclusive coordinates)
        try:
            fasta_seq = self.fasta_handler.fetch(fasta_chrom, pos, end_pos)
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
