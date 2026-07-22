# KVar-Toolkit SV: SV VCF Validation & Variant Call TSV Creation

## Overview

The **SV** module of KVar-Toolkit validates structural variation VCF input and
produces a KVar-formatted `Variant_Call.tsv` for public variant submission. It
assigns KVar accession-style call IDs, preserves submitted IDs for traceability,
and reports validation issues in a separate error report.

Run the commands below from this `SV/` directory.

## Key Features

- **SV VCF to Variant Call TSV conversion**: Convert a structural variation VCF into `Variant_Call.tsv`.
- **Reference-based validation**: Validate chromosome names, coordinates, and REF alleles against an indexed reference FASTA.
- **SV type classification**: Classify deletion, insertion, duplication, inversion, copy number variation, mobile element insertion/deletion, BND, and complex events.
- **BND/MATEID handling**: Validate reciprocal IDs, target coordinates, and VCF breakend orientations before collapsing two translocation rows into one call.
- **Call accession assignment**: Rewrite call IDs as `kssvN` while retaining submitted IDs in the output.

### Prerequisites

- Python 3.10 or higher

## Quick Start

Convert an SV VCF into a KVar `Variant_Call.tsv`:

```bash
python src/kvar_sv_tools/vcf_to_kvar_tsv.py \
  -v input.sv.vcf.gz \
  -f GRCh38.fa \
  -m metadata.txt \
  -t Variant_Call.tsv \
  -e validation_report.txt \
  -c 1
```

The long-option form is equivalent:

```bash
python src/kvar_sv_tools/vcf_to_kvar_tsv.py \
  --vcf input.sv.vcf.gz \
  --reference-fasta GRCh38.fa \
  --metadata metadata.txt \
  --call-tsv Variant_Call.tsv \
  --error-report validation_report.txt \
  --call-accession-start 1
```

## Common Options

| Option | Description |
| --- | --- |
| `-v`, `--vcf` | Input SV VCF path (`.gz` supported, **required**) |
| `-f`, `-r`, `--reference-fasta` | Reference FASTA path for coordinate and REF validation (`.fai` required, **required**) |
| `-m`, `--metadata` | Metadata file path (**required**) |
| `-t`, `--call-tsv` | Output `Variant_Call.tsv` path (**required**) |
| `-e`, `--error-report` | Validation report path (**required**) |
| `-c`, `--call-accession-start` | Starting number for Variant Call accessions, written as `kssvN` (**required**) |
| `--sanitize-error-report` | Redact absolute paths and raw row content from the validation report |

## Metadata Format

Metadata files use VCF-style lines:

```text
##Experiment_id=EXP001
##reference=GRCh38
##SampleSet_id=POP1
##organism_taxid=9606 (Homo sapiens)
```

The `##reference` value should match the input VCF `##reference` header. The
metadata values are also used as fallback context when `SAMPLESET` or
`EXPERIMENT` are not present in every VCF INFO field. The optional
`##organism_taxid` value is written verbatim to the output TSV headers.

## Outputs

The converter writes:

- `Variant_Call.tsv`: normalized variant call records
- validation report: errors, warnings, and repair actions

Original submitted VCF row IDs are retained in `Submitted_Variant_Call_IDs`.
BND mate rows that collapse into one translocation call are mapped to the same
normalized call accession.
Call TSV column names use underscore-separated headers, such as
`Variant_Call_ID`, `Variant_Call_Type`, and `Outer_Start`.
When the metadata file provides `##organism_taxid`, the call TSV carries an
`##organism_taxid=...` header line after the `##Variant_Call` label line.

`Variant_Call.tsv` also contains `HGVSG` after `Phenotype`. The converter derives
it for deletion/mobile-element deletion/copy-number loss (`del`),
duplication/copy-number gain/tandem duplication (`dup`), inversion (`inv`), and
insertion families (`ins`) when the reference FASTA supplies a versioned genomic
accession. Exact insertions use adjacent flanking positions; unknown inserted
sequence is represented as `N[length]` or `N[?]`. Imprecise outer/inner bounds
are retained as uncertain HGVS coordinates. BND/translocation, complex,
sequence-alteration, STR, indel, and direction-unknown CNV calls remain `.`
rather than receiving a speculative expression. Input aliases `HGVSg` and
dbVar `hgvs_name` are normalized to `HGVSG` when present.

A BND pair is collapsed only when `MATEID` is reciprocal, both ALT target
coordinates resolve to the mate records (including `CIPOS` uncertainty), and the
ALT strand orientations are reciprocal. For internal `(From, To)` strands, the
mate must be `(opposite(To), opposite(From))`. A strand mismatch is reported as
`MATEID_STRAND_MISMATCH`, and the two records remain separate calls.

If a BND ALT carries inserted sequence, the sequence is retained in a derived
insertion call linked by the same Mutation ID. This DDBJ-compatible split uses
`Start=POS` and `Stop=POS+inserted_length`; ordinary point insertions still use
`Start=Stop`.

If neither ALT nor `SVTYPE` identifies a supported type but INFO contains the
GATK-SV-style `CPX_TYPE` tag, the converter follows the DDBJ fallback and writes
the call type as `sequence alteration`. `CPX` means complex structural variant;
an already recognized SV type is not replaced by this fallback.

The call TSV is staged and published only after the full requested conversion
succeeds. Blocking validation errors remove the staged data output while
preserving the validation report.

## Project Structure

```text
SV/
├── README.md            # This file
├── src/kvar_sv_tools/
│   ├── vcf_to_kvar_tsv.py            # Public CLI entry point
│   ├── KVar2TSV.py                   # VCF validation and Variant Call TSV writing
│   ├── VCF_parser.py                 # VCF parser and reference checks
│   ├── sv_type_ontology.py           # SV type constants and mappings
│   ├── metadata_validator.py         # Metadata validation against VCF headers
│   ├── metadata_parser.py            # Metadata file parser (organism_taxid, sampleset, experiment)
│   └── error_handler.py              # Error codes and validation report
└── tests/
    └── test_public_cli_smoke.py
```

## Testing

The smoke test creates synthetic temporary input files at runtime. No real VCF,
reference FASTA, dbVar/DDBJ downloads, or Manta result files are checked into
this public subset.

```bash
python tests/test_public_cli_smoke.py
```

## Notes

- This public CLI is intended for VCF-to-`Variant_Call.tsv` conversion only.
- Reference FASTA validation requires an existing `.fai` index.
- This public subset does not include private datasets, generated full-scale outputs, or internal pipeline reports.
