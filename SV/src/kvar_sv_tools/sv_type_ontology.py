#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Canonical type definitions for KVar structural variants."""

MOBILE_ELEMENT_TYPES = ("ALU", "LINE1", "SVA", "HERV")

TYPE_ALIASES = {
    "delins": "indel",
    "short tandem repeat": "short tandem repeat variation",
}

ALT_SYMBOLIC_CALL_TYPES = (
    ("DEL", "deletion"),
    ("INS", "insertion"),
    ("DUP", "duplication"),
    ("INV", "inversion"),
    ("CNV", "copy number variation"),
    ("DUP:TANDEM", "tandem duplication"),
    ("INS:NOVEL", "novel sequence insertion"),
    ("INS:ME", "mobile element insertion"),
    ("INS:ME:ALU", "alu insertion"),
    ("INS:ME:HERV", "herv insertion"),
    ("INS:ME:LINE1", "line1 insertion"),
    ("INS:ME:SVA", "sva insertion"),
    ("DEL:ME", "mobile element deletion"),
    ("DEL:ME:ALU", "alu deletion"),
    ("DEL:ME:HERV", "herv deletion"),
    ("DEL:ME:LINE1", "line1 deletion"),
    ("DEL:ME:SVA", "sva deletion"),
)

SVTYPE_CALL_TYPES = {
    "DEL": "deletion",
    "INS": "insertion",
    "DUP": "duplication",
    "INV": "inversion",
    "CNV": "copy number variation",
}

def normalize_call_type(call_type: str) -> str:
    """Return KVar's canonical Variant Call Type term for a submitted/dbVar term."""
    call_type_lower = str(call_type or "").strip().lower()
    return TYPE_ALIASES.get(call_type_lower, call_type_lower)
