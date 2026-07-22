#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Metadata file parser for the KVar SV VCF-to-TSV converter."""

from typing import Optional


class MetadataParser:
    """Metadata file parser"""

    def __init__(self, metadata_path: str):
        self.metadata_path = metadata_path
        self.sampleset_id: Optional[str] = None
        self.experiment_id: Optional[str] = None
        self.reference: Optional[str] = None
        self.organism_taxid: Optional[str] = None
        self._parse_metadata()

    def _parse_metadata(self):
        """Parse metadata file"""
        try:
            with open(self.metadata_path, "rt", encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    # ##KEY=VALUE (strip ##)
                    if line.startswith("##"):
                        line = line[2:]  # Strip ##

                    # Skip comment lines (# but not ##)
                    if line.startswith("#"):
                        continue

                    # key=value
                    if "=" in line:
                        key, value = line.split("=", 1)
                        key = key.strip()
                        value = value.strip()

                        # Normalize key (case/underscore insensitive)
                        key_normalized = key.upper().replace("_", "").replace("-", "")

                        # SAMPLESET ID
                        if key_normalized in ["SAMPLESET", "SAMPLESETID", "SAMPLESET_ID"]:
                            self.sampleset_id = value
                        # EXPERIMENT ID
                        elif key_normalized in ["EXPERIMENT", "EXPERIMENTID", "EXPERIMENT_ID"]:
                            self.experiment_id = value
                        # Reference assembly/accession name
                        elif key_normalized == "REFERENCE":
                            self.reference = value
                        # Organism NCBI taxonomy id, e.g. "9606 (Homo sapiens)"
                        elif key_normalized == "ORGANISMTAXID":
                            self.organism_taxid = value
                    else:
                        # Tab-separated
                        parts = line.split("\t")
                        if len(parts) >= 2:
                            key = parts[0].strip()
                            value = parts[1].strip()

                            key_normalized = key.upper().replace("_", "").replace("-", "")

                            if key_normalized in ["SAMPLESET", "SAMPLESETID", "SAMPLESET_ID"]:
                                self.sampleset_id = value
                            elif key_normalized in ["EXPERIMENT", "EXPERIMENTID", "EXPERIMENT_ID"]:
                                self.experiment_id = value
                            elif key_normalized == "REFERENCE":
                                self.reference = value
                            elif key_normalized == "ORGANISMTAXID":
                                self.organism_taxid = value

        except Exception as e:
            raise ValueError(f"Failed to parse metadata file: {e}")
