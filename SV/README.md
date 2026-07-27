# KVar-Toolkit SV: SV VCF Validation & Variant Call TSV Creation

## Overview

The **SV** module of KVar-Toolkit validates structural variation VCF input and
produces a KVar-formatted `Variant_Call.tsv` for public variant submission. It
preserves submitted IDs for traceability and reports validation issues in a
single validation report. Accessions are not assigned during submission conversion.

Run the commands below from this `SV/` directory.

## Key Features

- **SV VCF to Variant Call TSV conversion**: Convert a structural variation VCF into `Variant_Call.tsv`.
- **Optional reference-based validation**: Validate chromosome names, coordinates, and REF alleles when an indexed reference FASTA is supplied.
- **SV type classification**: Classify deletion, insertion, duplication, inversion, copy number variation, mobile element insertion/deletion, BND, and complex events.
- **BND/MATEID handling**: Validate reciprocal IDs, target coordinates, and VCF breakend orientations before collapsing two translocation rows into one call.
- **Submission ID preservation**: Keep submitted VCF IDs in the Call TSV; KVar accessions are assigned only after administrator QC.

### Prerequisites

- Python 3.10 or higher

## Quick Start

Convert an SV VCF into a KVar `Variant_Call.tsv`:

```bash
python src/kvar_sv_tools/vcf_to_kvar_tsv.py \
  -v examples/toy.human.sv.vcf
```

This creates `examples/toy.human.sv.Variant_Call.tsv` and
`examples/toy.human.sv.Variant_Call.errors.txt`.

To enable reference validation, provide the synthetic indexed FASTA:

```bash
python src/kvar_sv_tools/vcf_to_kvar_tsv.py \
  --vcf examples/toy.human.sv.vcf \
  --reference examples/toy.human.GRCh38.fa \
  --call-tsv Variant_Call.tsv \
  --error-report validation_report.txt
```

## Common Options

| Option | Description |
| --- | --- |
| `-v`, `--vcf` | Input SV VCF path (`.gz` supported, **required**) |
| `-r`, `--reference` | Optional reference FASTA for coordinate and REF validation (`.fai` required when used) |
| `-t`, `--call-tsv` | Optional output Call TSV path; must end with `.tsv` |
| `-e`, `--error-report` | Optional validation report path; defaults to the Call TSV basename with `.errors.txt` |

If `--call-tsv` is omitted, `sample.vcf` or `sample.vcf.gz` produces
`sample.Variant_Call.tsv`.
Each run writes one validation report. If `--error-report` is omitted, an output
such as `sample.Variant_Call.tsv` produces `sample.Variant_Call.errors.txt`.
Generated Call TSV and report paths must differ from the input VCF and optional
reference FASTA path, and from each other.

## Input Scope

The SV converter accepts an SV VCF and can optionally use an indexed reference
FASTA. It does not accept a separate metadata file, and it does not require or
validate `SampleSet`, `Experiment`, or the VCF `##reference` header. Chromosome,
coordinate, and REF allele validation runs only when a FASTA is supplied through
`--reference`.

The input VCF must already contain the records selected for submission. The
converter does not remove records based on the VCF `FILTER` column. It also does
not calculate AC, AN, or AF from sample `FORMAT` genotype fields.

## Outputs

The converter writes:

- `Variant_Call.tsv`: normalized variant call records
- validation report: errors, warnings, and repair actions

Original submitted VCF row IDs are retained in `Submitted_Variant_Call_IDs`.
BND mate rows that collapse into one translocation call are mapped to the same
submission call ID. The converter does not accept an accession-start option and
does not create `kssvN` identifiers; those are assigned only after QC.
Call TSV column names use underscore-separated headers, such as
`Variant_Call_ID`, `Variant_Call_Type`, and `Outer_Start`.

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
├── examples/            # Privacy-safe synthetic GRCh38-style VCF and FASTA
├── src/kvar_sv_tools/
│   ├── vcf_to_kvar_tsv.py            # Public CLI entry point
│   ├── KVar2TSV.py                   # VCF validation and Variant Call TSV writing
│   ├── VCF_parser.py                 # VCF parser and reference checks
│   ├── sv_type_ontology.py           # SV type constants and mappings
│   └── error_handler.py              # Error codes and validation report
└── tests/
    └── test_public_cli_smoke.py
```

## Testing

The module includes a small synthetic human example under `examples/`. It uses
GRCh38 chromosome/accession names but contains no real person, sample, genotype,
or cohort data. No real VCF, full reference FASTA, dbVar/DDBJ download, or Manta
result file is checked into this public subset.

```bash
python tests/test_public_cli_smoke.py
```

## Notes

- This public CLI is intended for VCF-to-`Variant_Call.tsv` conversion only.
- Reference validation is optional; when used, the FASTA requires an existing `.fai` index.
- This public subset does not include private datasets, generated full-scale outputs, or internal pipeline reports.
