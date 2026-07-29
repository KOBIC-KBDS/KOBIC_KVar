#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KVar SV VCF parser.
Parses VCF files conforming to KVar SV VCF submission format and extracts all information.
"""

import os
import gzip
import re
import tempfile
from typing import Dict, List, Optional, Any, Tuple, Set
from dataclasses import dataclass, field
from collections import defaultdict

# Relative path import support
try:
    from .error_handler import ErrorHandler, ErrorCode, ErrorSeverity
    from .sv_type_ontology import ALT_SYMBOLIC_CALL_TYPES, MOBILE_ELEMENT_TYPES, SVTYPE_CALL_TYPES
except ImportError:
    from error_handler import ErrorHandler, ErrorCode, ErrorSeverity
    from sv_type_ontology import ALT_SYMBOLIC_CALL_TYPES, MOBILE_ELEMENT_TYPES, SVTYPE_CALL_TYPES


@dataclass
class VCFHeaderMetadata:
    """Class to store VCF header metadata"""
    fileformat: Optional[str] = None
    filedate: Optional[str] = None
    source: Optional[str] = None


@dataclass
class InfoTagDefinition:
    """Class to store INFO tag definition"""
    id: str
    number: str
    type: str
    description: str


@dataclass
class AltTagDefinition:
    """Class to store ALT tag definition"""
    id: str
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
    alt_tags: Dict[str, AltTagDefinition] = field(default_factory=dict)
    contig_fields: Dict[str, Dict[str, str]] = field(default_factory=dict)
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


class SVClassifier:
    """Structural Variation classification (VCF v4.1 compliant)"""

    @staticmethod
    def classify_calltype(
        info: Dict[str, Any],
        alt: str,
        ref: Optional[str] = None,
    ) -> str:
        """Classify SV type (VCF v4.1 standard)"""
        svt = info.get("SVTYPE", "")
        alt_call_type = SVClassifier._classify_symbolic_alt(alt)
        call_type = alt_call_type or SVTYPE_CALL_TYPES.get(svt, "")

        if call_type in {"insertion", "deletion"} and SVClassifier._is_indel_sequence(ref, alt, call_type):
            return "indel"

        if call_type in {"alu insertion", "herv insertion", "line1 insertion", "sva insertion"}:
            return call_type
        if call_type in {"alu deletion", "herv deletion", "line1 deletion", "sva deletion"}:
            return call_type
        if call_type in {
            "copy number variation",
            "deletion",
            "duplication",
            "inversion",
            "mobile element deletion",
            "mobile element insertion",
            "novel sequence insertion",
            "tandem duplication",
        }:
            return call_type

        mobile_element = SVClassifier._mobile_element_from_alt(alt)
        is_mobile_element = SVClassifier._is_generic_mobile_element_alt(alt)

        # Insertion (INS) special handling
        if svt == "INS":
            if mobile_element:
                return f"{mobile_element.lower()} insertion"
            elif is_mobile_element:
                return "mobile element insertion"
            return "insertion"

        if svt == "DEL":
            if mobile_element:
                return f"{mobile_element.lower()} deletion"
            elif is_mobile_element:
                return "mobile element deletion"
            return "deletion"
        if svt == "DUP":
            if info.get("SUBTYPE") == "tandem":
                return "tandem duplication"
            return "duplication"
        if svt == "BND":
            return "translocation (to be subclassified)"
        if "STR" in info or "SSR" in info:
            return "short tandem repeat variation"
        if "CPX_TYPE" in info:
            return "sequence alteration"
        return ""

    @staticmethod
    def _classify_symbolic_alt(alt: str) -> Optional[str]:
        if alt:
            alt_upper = alt.upper()
            if alt_upper.startswith("<"):
                alt_id = alt_upper[1:].split(">", 1)[0]
                call_type = None
                for symbol, symbol_call_type in ALT_SYMBOLIC_CALL_TYPES:
                    if alt_id.startswith(symbol):
                        call_type = symbol_call_type
                return call_type
        return None

    @staticmethod
    def _mobile_element_from_alt(alt: str) -> Optional[str]:
        alt_upper = str(alt or "").upper()
        for mobile_element in MOBILE_ELEMENT_TYPES:
            if f":ME:{mobile_element}" in alt_upper:
                return mobile_element
        return None

    @staticmethod
    def _is_generic_mobile_element_alt(alt: str) -> bool:
        alt_upper = str(alt or "").upper()
        return ":ME:" in alt_upper or ":ME>" in alt_upper or "INS:ME" in alt_upper or "DEL:ME" in alt_upper

    @staticmethod
    def _is_indel_sequence(ref: Optional[str], alt: str, call_type: str) -> bool:
        if ref is None or not alt:
            return False
        if "[" in alt or "]" in alt or alt.startswith("<"):
            return False
        ref_upper = str(ref).upper()
        alt_upper = str(alt).upper()
        if len(ref_upper) <= 1 or len(alt_upper) <= 1:
            return False
        if ref_upper[0] != alt_upper[0] or ref_upper[1] == alt_upper[1]:
            return False
        if call_type == "insertion" and len(ref_upper) < len(alt_upper):
            return True
        if call_type == "deletion" and len(ref_upper) > len(alt_upper):
            return True
        return False


class BreakendParser:
    """Breakend parsing (VCF v4.1 compliant)"""

    CONTIG_PATTERN = (
        r'[0-9A-Za-z!#$%&+./:;?@^_|~-]'
        r'[0-9A-Za-z!#$%&*+./:;=?@^_|~-]*'
    )

    @staticmethod
    def _parse_paired_breakend_alt(
        alt: str,
    ) -> Optional[Tuple[str, str, str, str]]:
        """Return local placement, local sequence, remote chromosome, and position."""
        suffix_match = re.fullmatch(
            r'(?P<local>[ACGTNacgtn]+)'
            r'(?P<bracket>[\[\]])'
            rf'(?P<chrom>{BreakendParser.CONTIG_PATTERN}):(?P<pos>[1-9]\d*)'
            r'(?P=bracket)',
            alt or "",
        )
        if suffix_match:
            return (
                "suffix",
                suffix_match.group("local"),
                suffix_match.group("chrom"),
                suffix_match.group("pos"),
            )

        prefix_match = re.fullmatch(
            r'(?P<bracket>[\[\]])'
            rf'(?P<chrom>{BreakendParser.CONTIG_PATTERN}):(?P<pos>[1-9]\d*)'
            r'(?P=bracket)'
            r'(?P<local>[ACGTNacgtn]+)',
            alt or "",
        )
        if prefix_match:
            return (
                "prefix",
                prefix_match.group("local"),
                prefix_match.group("chrom"),
                prefix_match.group("pos"),
            )
        return None

    @staticmethod
    def parse_breakend(alt: str) -> Tuple[str, str, str]:
        """Extract translocation info from Breakend ALT field"""
        _, chr2, pos2, strand = BreakendParser.parse_breakend_placement(alt)
        return chr2, pos2, strand

    @staticmethod
    def parse_breakend_placement(alt: str) -> Tuple[str, str, str, str]:
        """Return from_strand, to_chr, to_pos, to_strand for a paired breakend ALT."""
        parsed = BreakendParser._parse_paired_breakend_alt(alt)
        if not parsed:
            return ".", ".", ".", "."

        placement, _, to_chr, to_pos = parsed
        bracket = (alt or "")[0] if placement == "prefix" else (alt or "")[-1]
        if placement == "prefix" and bracket == "]":
            return "-", to_chr, to_pos, "-"
        if placement == "prefix" and bracket == "[":
            return "-", to_chr, to_pos, "+"
        if bracket == "]":
            return "+", to_chr, to_pos, "-"
        if bracket == "[":
            return "+", to_chr, to_pos, "+"
        return ".", to_chr, to_pos, "."

    @staticmethod
    def is_breakend_alt(alt: str, ref: Optional[str] = None) -> bool:
        """Return True when ALT is one complete, REF-anchored paired breakend."""
        parsed = BreakendParser._parse_paired_breakend_alt(alt)
        if not parsed:
            return False
        if not ref:
            return True

        placement, local_sequence, _, _ = parsed
        ref_upper = str(ref).upper()
        local_upper = local_sequence.upper()
        if placement == "suffix":
            return local_upper.startswith(ref_upper)
        return local_upper.endswith(ref_upper)

    @staticmethod
    def parse_single_breakend(alt: str, ref: str) -> Tuple[str, str]:
        """Return from_strand and inserted sequence for a valid single-breakend ALT."""
        ref_pattern = re.escape(ref or "")
        if not ref_pattern:
            return ".", ""
        right_match = re.fullmatch(rf'{ref_pattern}([ACGTNacgtn]*)\.', alt or "")
        if right_match:
            return "+", right_match.group(1)
        left_match = re.fullmatch(rf'\.([ACGTNacgtn]*){ref_pattern}', alt or "")
        if left_match:
            return "-", left_match.group(1)
        return ".", ""

    @staticmethod
    def is_single_breakend_alt(alt: str, ref: str) -> bool:
        """Return True when ALT is a valid VCF single-breakend allele."""
        from_strand, _ = BreakendParser.parse_single_breakend(alt, ref)
        return from_strand != "."

    @staticmethod
    def inserted_sequence(alt: str, ref: str) -> str:
        """Return sequence inserted within a paired or single breakend ALT."""
        if not ref:
            return ""

        single_strand, single_sequence = BreakendParser.parse_single_breakend(alt, ref)
        if single_strand != ".":
            return single_sequence
        if not BreakendParser.is_breakend_alt(alt, ref):
            return ""

        placement, local_sequence, _, _ = BreakendParser._parse_paired_breakend_alt(alt) or (
            "",
            "",
            "",
            "",
        )
        if placement == "suffix":
            return local_sequence[len(ref):]
        return local_sequence[:-len(ref)]


class FastaReference:
    """Indexed FASTA reader with a reusable, automatically generated .fai."""

    GRCH38_REFSEQ_ALIASES = {
        "NC_000001.11": "chr1",
        "NC_000002.12": "chr2",
        "NC_000003.12": "chr3",
        "NC_000004.12": "chr4",
        "NC_000005.10": "chr5",
        "NC_000006.12": "chr6",
        "NC_000007.14": "chr7",
        "NC_000008.11": "chr8",
        "NC_000009.12": "chr9",
        "NC_000010.11": "chr10",
        "NC_000011.10": "chr11",
        "NC_000012.12": "chr12",
        "NC_000013.11": "chr13",
        "NC_000014.9": "chr14",
        "NC_000015.10": "chr15",
        "NC_000016.10": "chr16",
        "NC_000017.11": "chr17",
        "NC_000018.10": "chr18",
        "NC_000019.10": "chr19",
        "NC_000020.11": "chr20",
        "NC_000021.9": "chr21",
        "NC_000022.11": "chr22",
        "NC_000023.11": "chrX",
        "NC_000024.10": "chrY",
        "NC_012920.1": "chrM",
    }

    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.index_path = f"{fasta_path}.fai"
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        self.alias_to_chrom: Dict[str, str] = {}
        self._fasta_handle = None

        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        if self._index_is_current():
            try:
                self.index = self._load_index(self.index_path)
            except (OSError, ValueError):
                self.index = self._build_index()
                self._cache_index()
        else:
            self.index = self._build_index()
            self._cache_index()

        self._build_aliases()

    def _index_is_current(self) -> bool:
        """Return True when the cached index is at least as new as the FASTA."""
        try:
            return (
                os.stat(self.index_path).st_mtime_ns
                >= os.stat(self.fasta_path).st_mtime_ns
            )
        except OSError:
            return False

    def _load_index(self, index_path: str) -> Dict[str, Tuple[int, int, int, int]]:
        """Load and structurally validate a FASTA index."""
        index: Dict[str, Tuple[int, int, int, int]] = {}
        fasta_size = os.path.getsize(self.fasta_path)
        with open(index_path, "r", encoding="utf-8") as fai:
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
            raise ValueError(f"FASTA index is empty: {index_path}")

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

                next_header_offset = self._validate_cached_index_record(
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
        """Read only the FASTA header immediately preceding a sequence offset."""
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

    def _validate_cached_index_record(
        self,
        fasta,
        fasta_size: int,
        chrom: str,
        values: Tuple[int, int, int, int],
    ) -> int:
        """Run bounded offset/end-of-record checks on one cached index entry."""
        length, offset, line_bases, line_width = values
        header_name = self._header_name_before_offset(fasta, offset)
        if header_name != chrom:
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
            nonlocal current_name
            if current_name is None:
                return
            if (
                sequence_length < 1
                or sequence_offset is None
                or line_bases is None
                or line_width is None
            ):
                raise ValueError(
                    f"FASTA sequence has no bases: {current_name}"
                )
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
        """Atomically cache the in-memory index, continuing if the directory is read-only."""
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
            # The in-memory index remains fully usable when the FASTA directory
            # cannot be written.
            pass
        finally:
            if temporary_path:
                try:
                    os.unlink(temporary_path)
                except OSError:
                    pass

    @staticmethod
    def _alias_key(chrom: str) -> str:
        return str(chrom or "").strip().lower()

    def _register_alias(self, alias: str, chrom: str) -> None:
        alias = str(alias or "").strip()
        if not alias or chrom not in self.index:
            return
        self.alias_to_chrom.setdefault(self._alias_key(alias), chrom)

    def _build_aliases(self) -> None:
        for chrom in self.index:
            self._register_alias(chrom, chrom)
            for alias in self._basic_aliases(chrom):
                self._register_alias(alias, chrom)

        for alias, chrom in self.GRCH38_REFSEQ_ALIASES.items():
            if chrom in self.index:
                self._register_alias(alias, chrom)

        self._load_fasta_header_aliases()

    def _basic_aliases(self, chrom: str) -> List[str]:
        aliases: List[str] = []
        if chrom.startswith("chr") and len(chrom) > 3:
            aliases.append(chrom[3:])
        else:
            aliases.append(f"chr{chrom}")

        if chrom in {"chrM", "chrMT", "M", "MT"}:
            aliases.extend(["chrM", "chrMT", "M", "MT", "MTDNA"])
        return aliases

    def _load_fasta_header_aliases(self) -> None:
        """Register accessions embedded in FASTA headers, e.g. AC:CM000663.2."""
        try:
            with open(self.fasta_path, "rb") as fasta:
                for chrom, (_, offset, _, _) in self.index.items():
                    header = self._read_header_before_offset(fasta, offset)
                    if not header:
                        continue
                    header_text = header.decode("ascii", errors="ignore").strip()
                    if header_text.startswith(">"):
                        header_text = header_text[1:]
                    first_token = header_text.split()[0] if header_text.split() else ""
                    self._register_alias(first_token, chrom)
                    for accession in re.findall(r"\b[A-Z]{1,3}_?\d{5,}(?:\.\d+)?\b", header_text):
                        self._register_alias(accession, chrom)
                    for accession in re.findall(r"\b[A-Za-z][A-Za-z0-9_]*:([A-Z]{1,3}_?\d{5,}(?:\.\d+)?)", header_text):
                        self._register_alias(accession, chrom)
        except OSError:
            return

    @staticmethod
    def _read_header_before_offset(fasta, sequence_offset: int) -> bytes:
        if sequence_offset <= 0:
            return b""
        chunk_size = min(sequence_offset, 4096)
        fasta.seek(sequence_offset - chunk_size)
        chunk = fasta.read(chunk_size)
        newline_pos = chunk.rfind(b"\n", 0, max(0, len(chunk) - 1))
        header = chunk[newline_pos + 1:] if newline_pos >= 0 else chunk
        if header.endswith(b"\n"):
            header = header[:-1]
        return header

    def resolve_chrom(self, chrom: str) -> Optional[str]:
        """Resolve exact and common chr/non-chr chromosome aliases."""
        chrom = str(chrom or "").strip()
        if not chrom:
            return None

        if chrom in self.index:
            return chrom

        mapped = self.alias_to_chrom.get(self._alias_key(chrom))
        if mapped:
            return mapped

        for candidate in self._basic_aliases(chrom):
            if candidate in self.index:
                return candidate
            mapped = self.alias_to_chrom.get(self._alias_key(candidate))
            if mapped:
                return mapped
        return None

    def length(self, chrom: str) -> Optional[int]:
        """Return contig length if present."""
        resolved = self.resolve_chrom(chrom)
        if resolved is None:
            return None
        return self.index[resolved][0]

    def close(self) -> None:
        """Close the shared FASTA handle used by repeated interval fetches."""
        if self._fasta_handle is not None:
            self._fasta_handle.close()
            self._fasta_handle = None

    def __del__(self):
        self.close()

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch 1-based inclusive sequence."""
        resolved = self.resolve_chrom(chrom)
        if resolved is None:
            raise KeyError(chrom)

        length, offset, line_bases, line_width = self.index[resolved]
        if start < 1 or end > length or start > end:
            raise ValueError(f"{chrom}:{start}-{end} is outside reference bounds")

        if self._fasta_handle is None:
            self._fasta_handle = open(self.fasta_path, "rb")
        fasta = self._fasta_handle
        seq_parts: List[str] = []
        pos = start
        while pos <= end:
            zero_based = pos - 1
            line_index = zero_based // line_bases
            line_offset = zero_based % line_bases
            bases_to_read = min(end - pos + 1, line_bases - line_offset)
            byte_offset = offset + line_index * line_width + line_offset
            fasta.seek(byte_offset)
            seq_parts.append(fasta.read(bases_to_read).decode("ascii"))
            pos += bases_to_read

        return "".join(seq_parts).upper()


class KVarVCFParser:
    """Main class to parse KVar SV VCF files"""

    def __init__(
        self,
        error_handler: Optional[ErrorHandler] = None,
        skip_metadata_validation: bool = False,
        reference_fasta_path: Optional[str] = None,
        strict_kvar_tags: bool = True
    ):
        self.error_handler = error_handler or ErrorHandler()
        self.sv_classifier = SVClassifier()
        self.breakend_parser = BreakendParser()
        self.skip_metadata_validation = skip_metadata_validation  # True for generic VCF mode
        self.strict_kvar_tags = strict_kvar_tags
        self.reference_fasta_path = reference_fasta_path
        self.reference: Optional[FastaReference] = None
        self._reset_parse_state()

        if reference_fasta_path:
            try:
                self.reference = FastaReference(reference_fasta_path)
            except ValueError as e:
                self.error_handler.create_error(
                    ErrorCode.FASTA_INDEX_ERROR,
                    additional_info={"file_path": reference_fasta_path, "error": str(e)}
                )
                raise
            except OSError as e:
                self.error_handler.create_error(
                    ErrorCode.FASTA_READ_ERROR,
                    additional_info={"file_path": reference_fasta_path, "error": str(e)}
                )
                raise

    def _reset_parse_state(self) -> None:
        """Reset per-file parser state."""
        self.header = VCFHeader()
        self.data_rows: List[VCFDataRow] = []
        self._metadata_tag_counts: Dict[str, int] = defaultdict(int)
        self._seen_ids: Set[str] = set()
        self._current_chrom: Optional[str] = None
        self._closed_chroms: Set[str] = set()
        self._last_pos_by_chrom: Dict[str, int] = {}
        self._bnd_mate_refs: List[Tuple[int, str, str]] = []

    def parse_file(self, file_path: str) -> None:
        """Main method to parse VCF file"""
        self._reset_parse_state()

        if not os.path.exists(file_path):
            self.error_handler.create_error(
                ErrorCode.FILE_NOT_FOUND,
                additional_info={"file_path": file_path}
            )
            raise FileNotFoundError(f"File not found: {file_path}")

        try:
            open_fn = gzip.open if file_path.endswith(".gz") else open
            with open_fn(file_path, "rt", encoding='utf-8') as f:
                lines = f.readlines()
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

        # Separate header and data
        header_lines = []
        data_lines = []

        first_content_line = None
        first_content_line_number = None
        for line_num, line in enumerate(lines, 1):
            original_line = line
            line = line.strip()
            if not line:
                self.error_handler.create_error(
                    ErrorCode.EMPTY_LINE_REMOVED,
                    line_number=line_num,
                    additional_info={"message": "Empty VCF line was skipped"}
                )
                continue
            if first_content_line is None:
                first_content_line = line
                first_content_line_number = line_num

            if line.startswith('#'):
                header_lines.append(line)
            else:
                data_lines.append((line_num, original_line.rstrip('\n\r')))

        if first_content_line and not first_content_line.startswith("##fileformat="):
            self.error_handler.create_error(
                ErrorCode.FILEFORMAT_NOT_FIRST,
                line_number=first_content_line_number,
                line_content=first_content_line,
                expected_value="##fileformat=VCFv... as the first non-empty line",
                actual_value=first_content_line.split("\t", 1)[0],
            )

        # Parse header
        self._parse_header(header_lines)

        # Parse data
        self._parse_data(data_lines)

        # Required validation
        self._validate_parsed_data()

    def _parse_header(self, header_lines: List[str]) -> None:
        """VCF 헤더를 파싱하는 메서드"""
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
            if not self.header.metadata.fileformat.startswith('VCFv'):
                self.error_handler.create_error(
                    ErrorCode.INVALID_FILEFORMAT,
                    line_number=1,
                    actual_value=self.header.metadata.fileformat,
                    expected_value="VCFv format (e.g. VCFv4.1)"
                )

        if not self.header.column_header:
            self.error_handler.create_error(
                ErrorCode.MISSING_COLUMN_HEADER
            )

        fileformat_count = self._metadata_tag_counts.get("fileformat", 0)
        if fileformat_count > 1:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_METADATA_TAG,
                field_name="fileformat",
                expected_value="unique metadata tag",
                actual_value=f"{fileformat_count} occurrences",
            )

    def _parse_metadata_line(self, line: str, line_number: int) -> None:
        """메타데이터 라인을 파싱하는 메서드"""
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
            self._metadata_tag_counts[key] += 1

            if key == 'fileformat':
                self.header.metadata.fileformat = value
            elif key == 'filedate':
                self.header.metadata.filedate = value
            elif key == 'source':
                self.header.metadata.source = value
            elif key == 'info':
                self._parse_info_tag_definition(line, line_number)
            elif key == 'format':
                self._parse_format_tag_definition(line, line_number)
            elif key == 'alt':
                self._parse_alt_tag_definition(line, line_number)
            elif key == 'contig':
                self._parse_contig_definition(line, line_number)
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.INVALID_METADATA_FORMAT,
                line_number=line_number,
                line_content=line,
                additional_info={"error": str(e)}
            )

    def _parse_info_tag_definition(self, line: str, line_number: int) -> None:
        """INFO 태그 정의를 파싱하는 메서드"""
        pattern = r'##INFO=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)"'
        match = re.search(pattern, line)

        if match:
            tag_id, number, type_str, description = match.groups()
            if tag_id in self.header.info_tags:
                self.error_handler.create_error(
                    ErrorCode.DUPLICATE_INFO_TAG_DEFINITION,
                    line_number=line_number,
                    field_name=tag_id,
                    line_content=line
                )
                return
            self.header.info_tags[tag_id] = InfoTagDefinition(
                id=tag_id,
                number=number,
                type=type_str,
                description=description
            )
        else:
            # More flexible parsing
            parts = line.split(',')
            if len(parts) >= 4:
                try:
                    id_part = parts[0].replace('##INFO=<ID=', '')
                    number_part = parts[1].replace('Number=', '')
                    type_part = parts[2].replace('Type=', '')
                    desc_start = line.find('Description="') + len('Description="')
                    desc_end = line.rfind('"')
                    if desc_start > len('Description="') - 1 and desc_end > desc_start:
                        description = line[desc_start:desc_end]
                        if id_part in self.header.info_tags:
                            self.error_handler.create_error(
                                ErrorCode.DUPLICATE_INFO_TAG_DEFINITION,
                                line_number=line_number,
                                field_name=id_part,
                                line_content=line
                            )
                            return
                        self.header.info_tags[id_part] = InfoTagDefinition(
                            id=id_part,
                            number=number_part,
                            type=type_part,
                            description=description
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

    def _parse_format_tag_definition(self, line: str, line_number: int) -> None:
        """FORMAT 태그 정의를 파싱하는 메서드"""
        pattern = r'##FORMAT=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)"'
        match = re.search(pattern, line)

        if not match:
            self.error_handler.create_error(
                ErrorCode.FORMAT_TAG_PARSE_ERROR,
                line_number=line_number,
                line_content=line
            )
            return

        tag_id, number, type_str, description = match.groups()
        if tag_id in self.header.format_tags:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_FORMAT_TAG_DEFINITION,
                line_number=line_number,
                field_name=tag_id,
                line_content=line
            )
            return

        self.header.format_tags[tag_id] = FormatTagDefinition(
            id=tag_id,
            number=number,
            type=type_str,
            description=description
        )

    def _parse_alt_tag_definition(self, line: str, line_number: int) -> None:
        """ALT 태그 정의를 파싱하는 메서드"""
        pattern = r'##ALT=<ID=([^,]+),Description="([^"]+)"'
        match = re.search(pattern, line)

        if match:
            tag_id, description = match.groups()
            self.header.alt_tags[tag_id] = AltTagDefinition(
                id=tag_id,
                description=description
            )

    def _parse_contig_definition(self, line: str, line_number: int) -> None:
        """CONTIG 정의를 파싱하는 메서드"""
        pattern = r'##contig=<ID=([^,]+),length=(\d+)'
        match = re.search(pattern, line)

        if match:
            contig_id, length = match.groups()
            self.header.contig_fields[contig_id] = {
                "ID": contig_id,
                "length": length
            }
        else:
            self.error_handler.create_error(
                ErrorCode.INVALID_METADATA_FORMAT,
                line_number=line_number,
                field_name="contig",
                expected_value="##contig=<ID=...,length=...>",
                actual_value=line
            )

    def _parse_column_header(self, line: str, line_number: int) -> None:
        """컬럼 헤더를 파싱하는 메서드"""
        if self.header.column_header:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_COLUMN_HEADER,
                line_number=line_number,
                line_content=line
            )
            return

        columns = line[1:].split('\t')

        required_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        if len(columns) < len(required_columns):
            self.error_handler.create_error(
                ErrorCode.INVALID_COLUMN_HEADER,
                line_number=line_number,
                line_content=line,
                expected_value=f"at least {len(required_columns)} columns",
                actual_value=f"{len(columns)} columns"
            )
        elif columns[:len(required_columns)] != required_columns:
            self.error_handler.create_error(
                ErrorCode.INVALID_COLUMN_HEADER,
                line_number=line_number,
                line_content=line,
                expected_value="\t".join(required_columns),
                actual_value="\t".join(columns[:len(required_columns)])
            )

        self.header.column_header = columns

    def _parse_data(self, data_lines: List[Tuple[int, str]]) -> None:
        """데이터 라인들을 파싱하는 메서드"""
        for line_number, line in data_lines:
            row = self._parse_data_row(line, line_number)
            if row:
                self.data_rows.append(row)

    def _parse_data_row(self, line: str, line_number: int) -> Optional[VCFDataRow]:
        """단일 데이터 행을 파싱하는 메서드"""
        raw_fields = line.split('\t')
        fields = [field.strip() for field in raw_fields]
        if fields != raw_fields:
            self.error_handler.create_error(
                ErrorCode.FIELD_WHITESPACE_TRIMMED,
                line_number=line_number,
                line_content=line,
                additional_info={"message": "One or more VCF fields had leading/trailing whitespace"}
            )

        # Parse ID first for variant_id (used in error messages)
        id_field = fields[2] if len(fields) > 2 else None
        variant_id = id_field if id_field and id_field != '.' else None

        if len(fields) < 8:
            self.error_handler.create_error(
                ErrorCode.INSUFFICIENT_FIELDS,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                expected_value="at least 8",
                actual_value=f"{len(fields)}"
            )
            return None

        # Basic fields
        chrom = fields[0]

        # CHROM validation is reference/contig based; only reject empty values here.
        if not chrom or chrom == ".":
            self.error_handler.create_error(
                ErrorCode.INVALID_CHROMOSOME_FORMAT,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="CHROM",
                actual_value=chrom,
                expected_value="non-empty chromosome, contig, or accession matching the reference"
            )

        # Position validation
        try:
            pos = int(fields[1])
            if pos < 1:
                self.error_handler.create_error(
                    ErrorCode.INVALID_POSITION,
                    line_number=line_number,
                    variant_id=variant_id,
                    line_content=line,
                    field_name="POS",
                    actual_value=str(pos)
                )
                return None
        except ValueError:
            self.error_handler.create_error(
                ErrorCode.INVALID_POSITION,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="POS",
                actual_value=fields[1]
            )
            return None

        self._validate_chromosome_order(self._canonical_chrom_for_order(chrom), pos, line_number, variant_id)

        ref = fields[3]
        alt = fields[4]

        # ID required validation
        if not id_field or id_field == '.':
            self.error_handler.create_error(
                ErrorCode.MISSING_REQUIRED_FIELD,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="ID"
            )
        elif len(id_field) > 64:
            self.error_handler.create_error(
                ErrorCode.LOCAL_ID_TOO_LONG,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="ID",
                expected_value="64 characters or fewer",
                actual_value=str(len(id_field))
            )
        elif id_field in self._seen_ids:
            self.error_handler.create_error(
                ErrorCode.DUPLICATE_LOCAL_ID,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="ID"
            )
        else:
            self._seen_ids.add(id_field)

        # REF/ALT validation
        if not ref or ref == '.':
            self.error_handler.create_error(
                ErrorCode.EMPTY_REF_ALT,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="REF"
            )
        if not alt or alt == '.':
            self.error_handler.create_error(
                ErrorCode.EMPTY_REF_ALT,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="ALT"
            )

        # REF validity (A, C, G, T, N only)
        valid_bases = set('ACGTNacgtn')
        if ref and not all(c in valid_bases for c in ref):
            self.error_handler.create_error(
                ErrorCode.INVALID_REF_ALT,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="REF",
                actual_value=ref
            )

        if alt and ',' in alt:
            self.error_handler.create_error(
                ErrorCode.MULTIALLELIC_ALT,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="ALT",
                actual_value=alt
            )
        elif alt and not self._is_valid_alt(alt, ref):
            self.error_handler.create_error(
                ErrorCode.INVALID_REF_ALT,
                line_number=line_number,
                variant_id=variant_id,
                line_content=line,
                field_name="ALT",
                actual_value=alt
            )

        self._validate_reference_allele(chrom, pos, ref, line_number, variant_id)

        qual = fields[5] if len(fields) > 5 and fields[5] != '.' else None
        filter_field = fields[6] if len(fields) > 6 and fields[6] != '.' else None
        if len(fields) > 8:
            self._validate_format_columns(fields[8:], line_number, variant_id)

        # INFO field parsing (VCF standard: 8th field is INFO, index 7)
        if len(fields) < 8:
            info_dict = {}
        else:
            info_str = fields[7]
            if not info_str or info_str.strip() == '' or info_str.strip() == '.':
                info_dict = {}
            else:
                info_dict = self._parse_info_field(info_str, line_number, variant_id)

        self._validate_info_tags(info_dict, chrom, pos, ref, alt, line_number, variant_id)

        return VCFDataRow(
            chrom=chrom,
            pos=pos,
            id=id_field,
            ref=ref,
            alt=alt,
            qual=qual,
            filter=filter_field,
            info=info_dict
        )

    def _is_valid_alt(self, alt: str, ref: Optional[str] = None) -> bool:
        """Validate one ALT allele for SV VCF syntax."""
        if alt.startswith("<") and alt.endswith(">"):
            return True
        if self.breakend_parser.is_breakend_alt(alt, ref):
            return True
        if alt.endswith(".") or alt.startswith("."):
            return bool(re.match(r'^[ACGTNacgtn]*\.$|^\.[ACGTNacgtn]*$', alt))
        return bool(re.match(r'^[ACGTNacgtn]+$', alt))

    def _validate_chromosome_order(
        self,
        chrom: str,
        pos: int,
        line_number: int,
        variant_id: Optional[str]
    ) -> None:
        """Validate chromosome grouping and position sorting."""
        if self._current_chrom is None:
            self._current_chrom = chrom
        elif chrom != self._current_chrom:
            self._closed_chroms.add(self._current_chrom)
            if chrom in self._closed_chroms:
                self.error_handler.create_error(
                    ErrorCode.CHROMOSOME_NOT_GROUPED,
                    line_number=line_number,
                    variant_id=variant_id,
                    field_name="CHROM",
                    actual_value=chrom
                )
            self._current_chrom = chrom

        last_pos = self._last_pos_by_chrom.get(chrom)
        if last_pos is not None and pos < last_pos:
            self.error_handler.create_error(
                ErrorCode.POSITION_NOT_SORTED,
                line_number=line_number,
                variant_id=variant_id,
                field_name="POS",
                expected_value=f">= {last_pos}",
                actual_value=str(pos)
            )
        self._last_pos_by_chrom[chrom] = pos

    def _canonical_chrom_for_order(self, chrom: str) -> str:
        if self.reference:
            resolved = self.reference.resolve_chrom(chrom)
            if resolved:
                return resolved
        if self.header.contig_fields:
            resolved = self._resolve_header_contig(chrom)
            if resolved:
                return resolved
        return chrom

    def _validate_reference_allele(
        self,
        chrom: str,
        pos: int,
        ref: str,
        line_number: int,
        variant_id: Optional[str]
    ) -> None:
        """Validate CHROM/POS/REF against the optional reference FASTA."""
        if not self.reference or not ref or ref == ".":
            return

        resolved_chrom = self.reference.resolve_chrom(chrom)
        if resolved_chrom is None:
            self.error_handler.create_error(
                ErrorCode.CHROMOSOME_NOT_FOUND,
                line_number=line_number,
                variant_id=variant_id,
                field_name="CHROM",
                actual_value=chrom
            )
            return

        ref_end = pos + len(ref) - 1
        chrom_length = self.reference.length(resolved_chrom)
        if chrom_length is None or ref_end > chrom_length:
            self.error_handler.create_error(
                ErrorCode.POSITION_OUT_OF_RANGE,
                line_number=line_number,
                variant_id=variant_id,
                field_name="POS/REF",
                expected_value=f"1-{chrom_length}" if chrom_length else "known contig bounds",
                actual_value=f"{chrom}:{pos}-{ref_end}"
            )
            return

        try:
            expected_ref = self.reference.fetch(resolved_chrom, pos, ref_end)
        except Exception as e:
            self.error_handler.create_error(
                ErrorCode.FASTA_READ_ERROR,
                line_number=line_number,
                variant_id=variant_id,
                additional_info={"error": str(e), "region": f"{chrom}:{pos}-{ref_end}"}
            )
            return

        if expected_ref != ref.upper():
            self.error_handler.create_error(
                ErrorCode.REF_MISMATCH,
                line_number=line_number,
                variant_id=variant_id,
                field_name="REF",
                expected_value=expected_ref,
                actual_value=ref
            )

    def _validate_info_tags(
        self,
        info_dict: Dict[str, Any],
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        line_number: int,
        variant_id: Optional[str]
    ) -> None:
        """Validate SV INFO tag presence, types, and cross-field consistency."""
        svtype = info_dict.get("SVTYPE")
        allowed_svtypes = {"DEL", "INS", "DUP", "INV", "CNV", "BND"}
        has_cpx_type = str(info_dict.get("CPX_TYPE", "")).strip() not in {"", "."}
        if "SVTYPE" not in info_dict:
            self.error_handler.create_error(
                ErrorCode.MISSING_SVTYPE_TAG,
                line_number=line_number,
                variant_id=variant_id,
                field_name="INFO"
            )
        elif svtype not in allowed_svtypes and not has_cpx_type:
            self.error_handler.create_error(
                ErrorCode.INVALID_SVTYPE_VALUE,
                line_number=line_number,
                variant_id=variant_id,
                field_name="SVTYPE",
                expected_value=", ".join(sorted(allowed_svtypes)) + " or CPX_TYPE-backed complex event",
                actual_value=str(svtype)
            )

        is_bnd = svtype == "BND"
        enforce_kvar_tags = self.strict_kvar_tags and not self.skip_metadata_validation

        if enforce_kvar_tags and not is_bnd:
            if "END" not in info_dict and self.derive_missing_end(svtype, pos, ref, alt, info_dict.get("SVLEN")) is None:
                self.error_handler.create_error(
                    ErrorCode.MISSING_END_TAG,
                    line_number=line_number,
                    variant_id=variant_id,
                    field_name="INFO"
                )
            if svtype != "INS" and "SVLEN" not in info_dict:
                self.error_handler.create_error(
                    ErrorCode.MISSING_SVLEN_TAG,
                    line_number=line_number,
                    variant_id=variant_id,
                    field_name="INFO"
                )

        end = info_dict.get("END")
        if end is not None:
            if not isinstance(end, int):
                self.error_handler.create_error(
                    ErrorCode.INVALID_END_VALUE,
                    line_number=line_number,
                    variant_id=variant_id,
                    field_name="END",
                    expected_value="integer",
                    actual_value=str(end)
                )
            elif not is_bnd and end < pos:
                self.error_handler.create_error(
                    ErrorCode.INVALID_END_VALUE,
                    line_number=line_number,
                    variant_id=variant_id,
                    field_name="END",
                    expected_value=f">= POS ({pos})",
                    actual_value=str(end)
                )

        svlen = info_dict.get("SVLEN")
        if svlen is not None and not self._is_integer_or_integer_list(svlen):
            self.error_handler.create_error(
                ErrorCode.INVALID_INFO_TAG_VALUE,
                line_number=line_number,
                variant_id=variant_id,
                field_name="SVLEN",
                expected_value="integer or comma-separated integers",
                actual_value=str(svlen)
            )

        self._validate_range_tag("POSrange", info_dict.get("POSrange"), pos, line_number, variant_id)
        if isinstance(end, int):
            self._validate_range_tag("ENDrange", info_dict.get("ENDrange"), end, line_number, variant_id)
        elif "ENDrange" in info_dict:
            self._validate_range_tag("ENDrange", info_dict.get("ENDrange"), None, line_number, variant_id)

        self._validate_ci_tag("CIPOS", info_dict.get("CIPOS"), line_number, variant_id)
        self._validate_ci_tag("CIEND", info_dict.get("CIEND"), line_number, variant_id)
        self._validate_ac_an(info_dict, line_number, variant_id)

        if is_bnd:
            is_breakend = self.breakend_parser.is_breakend_alt(alt, ref)
            is_single_breakend = self.breakend_parser.is_single_breakend_alt(alt, ref)
            if not is_breakend and not is_single_breakend:
                self.error_handler.create_error(
                    ErrorCode.INVALID_BREAKEND_FORMAT,
                    line_number=line_number,
                    variant_id=variant_id,
                    field_name="ALT",
                    actual_value=alt
                )
            mateid = info_dict.get("MATEID")
            if mateid:
                mate_ids = self._split_info_values(mateid)
                if "," in alt:
                    # MULTIALLELIC_ALT already blocks this unsupported shape.
                    pass
                elif len(mate_ids) > 1:
                    self.error_handler.create_error(
                        ErrorCode.MULTIPLE_MATEIDS_UNSUPPORTED,
                        line_number=line_number,
                        variant_id=variant_id,
                        field_name="MATEID",
                        expected_value="one MATEID for one BND ALT",
                        actual_value=",".join(mate_ids),
                    )
                elif variant_id and mate_ids and mate_ids[0] == variant_id:
                    self.error_handler.create_error(
                        ErrorCode.MATEID_SELF_REFERENCE,
                        line_number=line_number,
                        variant_id=variant_id,
                        field_name="MATEID",
                        expected_value="ID of a different BND record",
                        actual_value=mate_ids[0],
                    )
                else:
                    for mate in mate_ids:
                        self._bnd_mate_refs.append((line_number, variant_id or ".", mate))

    @staticmethod
    def _is_integer_or_integer_list(value: Any) -> bool:
        """Return True when value is an int or list of ints."""
        if isinstance(value, int):
            return True
        if isinstance(value, list):
            return all(isinstance(item, int) for item in value)
        return False

    @staticmethod
    def derive_missing_end(
        svtype: Any,
        pos: int,
        ref: str,
        alt: str,
        svlen: Any,
    ) -> Optional[int]:
        """Derive END for records that omit it, following DDBJ/VCF behavior."""
        if svtype == "BND":
            return None
        if svtype == "INS":
            return pos
        alt_text = str(alt or "")
        if re.fullmatch(r"[ATGCNatgcn]+", alt_text):
            return pos + len(ref or "") + 1
        if svtype in {"DEL", "DUP", "INV", "CNV"}:
            svlen_value = KVarVCFParser._first_integer_value(svlen)
            if svlen_value is not None:
                return pos + abs(svlen_value)
        return None

    @staticmethod
    def _first_integer_value(value: Any) -> Optional[int]:
        if isinstance(value, int):
            return value
        if isinstance(value, list):
            for item in value:
                parsed = KVarVCFParser._first_integer_value(item)
                if parsed is not None:
                    return parsed
            return None
        text = str(value or "").strip()
        if re.fullmatch(r"-?\d+", text):
            return int(text)
        return None

    def _validate_range_tag(
        self,
        tag: str,
        value: Any,
        anchor: Optional[int],
        line_number: int,
        variant_id: Optional[str]
    ) -> None:
        """Validate POSrange/ENDrange shape and anchor inclusion."""
        if value is None:
            return

        code = ErrorCode.INVALID_POSRANGE_FORMAT if tag == "POSrange" else ErrorCode.INVALID_ENDRANGE_FORMAT
        if not isinstance(value, list) or len(value) != 2:
            self.error_handler.create_error(
                code,
                line_number=line_number,
                variant_id=variant_id,
                field_name=tag,
                expected_value="two comma-separated values",
                actual_value=str(value)
            )
            return

        numeric_values = [item for item in value if item != "."]
        if not all(isinstance(item, int) for item in numeric_values):
            self.error_handler.create_error(
                code,
                line_number=line_number,
                variant_id=variant_id,
                field_name=tag,
                expected_value="integer or '.' values",
                actual_value=str(value)
            )
            return

        if anchor is None or len(numeric_values) != 2:
            return

        range_min, range_max = numeric_values
        if range_min > range_max or anchor not in numeric_values:
            self.error_handler.create_error(
                code,
                line_number=line_number,
                variant_id=variant_id,
                field_name=tag,
                expected_value=f"one {tag} value equal to {anchor}",
                actual_value=str(value)
            )

    def _validate_ci_tag(
        self,
        tag: str,
        value: Any,
        line_number: int,
        variant_id: Optional[str]
    ) -> None:
        """Validate CIPOS/CIEND as two offsets, allowing "." for one-sided CI."""
        if value is None:
            return
        code = ErrorCode.INVALID_POSRANGE_FORMAT if tag == "CIPOS" else ErrorCode.INVALID_ENDRANGE_FORMAT
        if (
            not isinstance(value, list)
            or len(value) != 2
            or not all(isinstance(item, int) or item == "." for item in value)
        ):
            self.error_handler.create_error(
                code,
                line_number=line_number,
                variant_id=variant_id,
                field_name=tag,
                expected_value='two comma-separated integer offsets, with "." allowed for an unknown side',
                actual_value=str(value)
            )

    def _validate_ac_an(
        self,
        info_dict: Dict[str, Any],
        line_number: int,
        variant_id: Optional[str]
    ) -> None:
        """Validate AC <= AN when both INFO values are present."""
        if "AC" not in info_dict or "AN" not in info_dict:
            return

        ac = info_dict["AC"]
        an = info_dict["AN"]
        ac_values = ac if isinstance(ac, list) else [ac]
        if not all(isinstance(value, int) for value in ac_values) or not isinstance(an, int):
            return

        if any(value > an for value in ac_values):
            self.error_handler.create_error(
                ErrorCode.ALLELE_COUNT_EXCEEDS_NUMBER,
                line_number=line_number,
                variant_id=variant_id,
                field_name="AC/AN",
                expected_value=f"AC <= AN ({an})",
                actual_value=str(ac)
            )

    def _validate_format_columns(
        self,
        format_and_sample_fields: List[str],
        line_number: int,
        variant_id: Optional[str]
    ) -> None:
        """Validate FORMAT keys and basic sample column arity."""
        if not format_and_sample_fields:
            return

        format_field = format_and_sample_fields[0]
        if not format_field or format_field == ".":
            return

        format_keys = [key.strip() for key in format_field.split(":") if key.strip()]
        for key in format_keys:
            if self.header.format_tags and key not in self.header.format_tags:
                self.error_handler.create_error(
                    ErrorCode.UNDEFINED_FORMAT_TAG,
                    line_number=line_number,
                    variant_id=variant_id,
                    field_name=key
                )

        expected_count = len(format_keys)
        for sample_value in format_and_sample_fields[1:]:
            if sample_value in {"", "."}:
                continue
            sample_count = len(sample_value.split(":"))
            if expected_count and sample_count != expected_count:
                self.error_handler.create_error(
                    ErrorCode.INVALID_FORMAT_TAG_VALUE,
                    line_number=line_number,
                    variant_id=variant_id,
                    field_name="FORMAT",
                    expected_value=f"{expected_count} value(s)",
                    actual_value=f"{sample_count} value(s)"
                )

    def _parse_info_field(self, info_str: str, line_number: int, variant_id: Optional[str] = None) -> Dict[str, Any]:
        """INFO 필드를 파싱하는 메서드"""
        info_dict = {}

        if not info_str or info_str == '.':
            return info_dict

        try:
            pairs = info_str.split(';')

            for pair in pairs:
                if not pair.strip():
                    continue
                if '=' in pair:
                    key, value = pair.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    if key:
                        if self.header.info_tags and key not in self.header.info_tags:
                            self.error_handler.create_error(
                                ErrorCode.UNDEFINED_INFO_TAG,
                                line_number=line_number,
                                variant_id=variant_id,
                                field_name=key
                            )
                        info_dict[key] = self._convert_info_value(key, value, line_number, variant_id)
                else:
                    # Flag-type tag
                    pair = pair.strip()
                    if pair:
                        if self.header.info_tags and pair not in self.header.info_tags:
                            self.error_handler.create_error(
                                ErrorCode.UNDEFINED_INFO_TAG,
                                line_number=line_number,
                                variant_id=variant_id,
                                field_name=pair
                            )
                        info_dict[pair] = True
        except Exception as e:
            # Record error on parse exception
            self.error_handler.create_error(
                ErrorCode.INVALID_INFO_FORMAT,
                line_number=line_number,
                variant_id=variant_id,
                field_name="INFO",
                actual_value=info_str,
                additional_info={"error": str(e)}
            )

        return info_dict

    def _convert_info_value(self, key: str, value: str, line_number: int, variant_id: Optional[str] = None) -> Any:
        """INFO 값의 타입을 변환하는 메서드"""
        # POSrange/ENDrange: Integer but may contain "."; handle first
        # Format: "2500000,2501000" or "2500000,." or ".,2501000"
        if key in ['POSrange', 'ENDrange', 'CIPOS', 'CIEND']:
            if ',' in value:
                parts = value.split(',')
                result = []
                for part in parts:
                    part = part.strip()
                    if part == '.':
                        result.append('.')
                    else:
                        try:
                            result.append(int(part))
                        except ValueError:
                            result.append(part)  # Keep as string on conversion failure
                return result
            else:
                # Single value
                if value.strip() == '.':
                    return '.'
                try:
                    return int(value)
                except ValueError:
                    return value

        if key in self.header.info_tags:
            tag_def = self.header.info_tags[key]

            if tag_def.type == 'Integer':
                try:
                    if ',' in value:
                        return [int(v) for v in value.split(',')]
                    return int(value)
                except ValueError:
                    self.error_handler.create_error(
                        ErrorCode.INVALID_INFO_TAG_VALUE,
                        line_number=line_number,
                        variant_id=variant_id,
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
                        ErrorCode.INVALID_INFO_TAG_VALUE,
                        line_number=line_number,
                        variant_id=variant_id,
                        field_name=key,
                        actual_value=value,
                        expected_value="float"
                    )
                    return value
            elif tag_def.type == 'String':
                return value
            elif tag_def.type == 'Flag':
                return True

        # Undefined tag: try generic type conversion
        if key in ['END', 'SVLEN', 'AC', 'AN', 'SVINSLEN']:
            try:
                if ',' in value:
                    return [int(v) for v in value.split(',')]
                return int(value)
            except ValueError:
                return value

        # Default is string
        return value

    @staticmethod
    def _split_info_values(value: Any) -> List[str]:
        """Return comma/list INFO values as clean strings."""
        if value is None or value == ".":
            return []
        if isinstance(value, list):
            values = value
        else:
            values = str(value).split(",")
        return [str(item).strip() for item in values if str(item).strip() and str(item).strip() != "."]

    def _validate_parsed_data(self) -> None:
        """파싱된 데이터의 전체 검증을 수행하는 메서드"""
        id_map = {row.id: row for row in self.data_rows if row.id and row.id != "."}
        valid_ids = set(id_map)
        for line_number, variant_id, mate_id in self._bnd_mate_refs:
            if mate_id not in valid_ids:
                self.error_handler.create_error(
                    ErrorCode.MATEID_NOT_FOUND,
                    line_number=line_number,
                    variant_id=None if variant_id == "." else variant_id,
                    field_name="MATEID",
                    actual_value=mate_id
                )
                continue

            if variant_id == "." or variant_id not in id_map:
                continue

            source_row = id_map[variant_id]
            mate_row = id_map[mate_id]
            mate_mates = self._split_info_values(mate_row.info.get("MATEID"))
            if variant_id not in mate_mates:
                self.error_handler.create_error(
                    ErrorCode.MATEID_NOT_RECIPROCAL,
                    line_number=line_number,
                    variant_id=variant_id,
                    field_name="MATEID",
                    expected_value=f"{mate_id} INFO/MATEID includes {variant_id}",
                    actual_value=str(mate_row.info.get("MATEID", "."))
                )

        for row in self.data_rows:
            self._validate_row_coordinate_bounds(row)

    def _coordinate_source(
        self,
        chrom: str,
    ) -> Tuple[Optional[str], Optional[str], Optional[int]]:
        """Return source, resolved contig, and length for coordinate checks."""
        if self.reference:
            resolved_chrom = self.reference.resolve_chrom(chrom)
            length = (
                self.reference.length(resolved_chrom)
                if resolved_chrom is not None
                else None
            )
            return "reference FASTA", resolved_chrom, length

        if self.header.contig_fields:
            resolved_chrom = self._resolve_header_contig(chrom)
            if resolved_chrom is None:
                return "VCF contig header", None, None
            raw_length = self.header.contig_fields[resolved_chrom].get("length")
            return (
                "VCF contig header",
                resolved_chrom,
                int(raw_length) if raw_length else None,
            )

        return None, None, None

    def _validate_row_coordinate_bounds(self, row: VCFDataRow) -> None:
        """Validate coordinates that will be represented in the Call TSV."""
        variant_id = row.id if row.id != "." else None
        source, resolved_chrom, chrom_length = self._coordinate_source(row.chrom)
        if source is None:
            return

        # Reference-backed POS/REF validation already reports local chromosome
        # and position errors while parsing the row. Without a reference, use
        # the VCF contig declarations as the local coordinate authority.
        if self.reference is None:
            if resolved_chrom is None:
                self.error_handler.create_error(
                    ErrorCode.CHROMOSOME_NOT_FOUND,
                    variant_id=variant_id,
                    field_name="CHROM",
                    actual_value=row.chrom,
                    additional_info={"source": source},
                )
            else:
                ref_length = (
                    len(row.ref)
                    if row.ref and row.ref != "."
                    else 1
                )
                ref_end = row.pos + ref_length - 1
            if (
                resolved_chrom is not None
                and chrom_length is not None
                and ref_end > chrom_length
            ):
                self.error_handler.create_error(
                    ErrorCode.POSITION_OUT_OF_RANGE,
                    variant_id=variant_id,
                    field_name="POS/REF",
                    expected_value=f"1-{chrom_length}",
                    actual_value=f"{row.chrom}:{row.pos}-{ref_end}",
                    additional_info={
                        "source": source,
                        "resolved_chrom": resolved_chrom,
                    },
                )

        self.validate_effective_end_bounds(
            row.chrom,
            row.pos,
            row.ref,
            row.alt,
            row.info,
            None,
            variant_id,
        )

        if (
            row.info.get("SVTYPE") == "BND"
            and self.breakend_parser.is_breakend_alt(row.alt, row.ref)
        ):
            _, remote_chrom, remote_pos_text, _ = (
                self.breakend_parser.parse_breakend_placement(row.alt)
            )
            remote_source, resolved_remote, remote_length = (
                self._coordinate_source(remote_chrom)
            )
            if remote_source is None:
                return
            if resolved_remote is None:
                self.error_handler.create_error(
                    ErrorCode.CHROMOSOME_NOT_FOUND,
                    variant_id=variant_id,
                    field_name="ALT remote CHROM",
                    actual_value=remote_chrom,
                    additional_info={"source": remote_source},
                )
                return

            remote_pos = int(remote_pos_text)
            # VCF permits a breakend target immediately after the last base to
            # represent the right telomere.
            remote_max = remote_length + 1 if remote_length is not None else None
            if remote_max is not None and remote_pos > remote_max:
                self.error_handler.create_error(
                    ErrorCode.POSITION_OUT_OF_RANGE,
                    variant_id=variant_id,
                    field_name="ALT remote coordinate",
                    expected_value=f"1-{remote_max}",
                    actual_value=f"{remote_chrom}:{remote_pos}",
                    additional_info={
                        "source": remote_source,
                        "resolved_chrom": resolved_remote,
                        "telomere_allowance": "contig length + 1",
                    },
                )

    def validate_effective_end_bounds(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        info: Dict[str, Any],
        line_number: Optional[int],
        variant_id: Optional[str],
        coordinate_origin: Optional[str] = None,
    ) -> None:
        """Validate the END written to TSV without comparing it with SVLEN."""
        svtype = info.get("SVTYPE")
        if svtype == "BND":
            return

        submitted_end = info.get("END")
        if submitted_end is not None and not isinstance(submitted_end, int):
            return

        effective_end = submitted_end
        derived = False
        if effective_end is None:
            effective_end = self.derive_missing_end(
                svtype,
                pos,
                ref,
                alt,
                info.get("SVLEN"),
            )
            derived = effective_end is not None
        if effective_end is None:
            return

        source, resolved_chrom, chrom_length = self._coordinate_source(chrom)
        if source is None or resolved_chrom is None or chrom_length is None:
            return
        if effective_end > chrom_length:
            self.error_handler.create_error(
                ErrorCode.POSITION_OUT_OF_RANGE,
                line_number=line_number,
                variant_id=variant_id,
                field_name="derived END" if derived and not coordinate_origin else "END",
                expected_value=f"1-{chrom_length}",
                actual_value=f"{chrom}:{effective_end}",
                additional_info={
                    "source": source,
                    "resolved_chrom": resolved_chrom,
                    "coordinate": coordinate_origin or (
                        "derived END" if derived else "submitted END"
                    ),
                },
            )

    def _resolve_header_contig(self, chrom: str) -> Optional[str]:
        chrom = str(chrom or "").strip()
        if not chrom:
            return None
        if chrom in self.header.contig_fields:
            return chrom

        candidates = []
        if chrom.startswith("chr"):
            candidates.append(chrom[3:])
        else:
            candidates.append(f"chr{chrom}")
        if chrom in {"chrM", "chrMT", "M", "MT"}:
            candidates.extend(["chrM", "chrMT", "M", "MT"])

        for candidate in candidates:
            if candidate in self.header.contig_fields:
                return candidate

        if self.reference:
            resolved = self.reference.resolve_chrom(chrom)
            if resolved in self.header.contig_fields:
                return resolved
        return None

    def get_summary(self) -> Dict[str, Any]:
        """파싱된 데이터의 요약 정보를 반환하는 메서드"""
        summary = {
            'header_metadata': {
                'fileformat': self.header.metadata.fileformat,
                'filedate': self.header.metadata.filedate,
                'source': self.header.metadata.source
            },
            'info_tags': {tag_id: tag_def.description for tag_id, tag_def in self.header.info_tags.items()},
            'format_tags': {tag_id: tag_def.description for tag_id, tag_def in self.header.format_tags.items()},
            'alt_tags': {tag_id: tag_def.description for tag_id, tag_def in self.header.alt_tags.items()},
            'column_header': self.header.column_header,
            'total_variants': len(self.data_rows),
            'svtype_distribution': defaultdict(int),
            'chromosomes': set(),
            'info_tag_usage': defaultdict(int)
        }

        # SV type statistics
        for row in self.data_rows:
            if 'SVTYPE' in row.info:
                svtype = row.info['SVTYPE']
                summary['svtype_distribution'][svtype] += 1

            summary['chromosomes'].add(row.chrom)

            # INFO tag usage statistics
            for key in row.info.keys():
                summary['info_tag_usage'][key] += 1

        summary['chromosomes'] = sorted(list(summary['chromosomes']))

        return summary


def main():
    """Test main function"""
    parser = KVarVCFParser()
    print("KVar SV VCF parser is ready.")


if __name__ == "__main__":
    main()
