# KVar-Toolkit SNP: SNP VCF Validation & dbSNP VCF Creation

## Overview

The **SNP** module of KVar-Toolkit validates SNP VCF input and produces cleaned,
dbSNP-formatted VCF output for public archive submission. It can convert a
generic VCF into a dbSNP-formatted VCF, validate and rewrite an existing dbSNP
VCF, and optionally verify REF alleles against a reference FASTA.

Run the commands below from this `SNP/` directory.

## Key Features

- **Generic to dbSNP conversion**: Convert a generic VCF into a dbSNP-formatted VCF.
- **dbSNP validation & cleaning**: Validate an input dbSNP VCF and rewrite it as a cleaned dbSNP VCF.
- **Reference allele validation**: Optionally validate REF alleles against a reference FASTA before writing output.
- **Metadata-driven headers**: Translate VCF-style metadata into output VCF headers.
- **Unified validation reporting**: Write format, metadata, and optional reference findings to one report.

### Prerequisites

- Python 3.10 or higher
- `pyfaidx` >= 0.8 only for optional reference validation (`pip install -r requirements.txt`)

## Quick Start

Convert a generic VCF into a dbSNP-formatted VCF:

```bash
python src/kvar_snp_tools/Sub_validator.py generic-to-dbsnp \
  -v examples/toy.generic.vcf \
  -m examples/toy.metadata.txt
```

This creates `examples/toy.generic.dbsnp.vcf` and
`examples/toy.generic.dbsnp.errors.txt`.

Validate and clean an existing dbSNP VCF:

```bash
python src/kvar_snp_tools/Sub_validator.py validate-dbsnp \
  -v examples/toy.dbsnp.vcf \
  -m examples/toy.metadata.txt
```

This creates `examples/toy.dbsnp.cleaned.vcf` and
`examples/toy.dbsnp.cleaned.errors.txt`.

Reference validation can be added to either command:

```bash
python src/kvar_snp_tools/Sub_validator.py validate-dbsnp \
  -v examples/toy.dbsnp.vcf \
  -m examples/toy.metadata.txt \
  -r examples/toy.reference.fa \
  -o examples/toy.dbsnp.cleaned.vcf \
  -e examples/toy.dbsnp.errors.txt
```

## Common Options

The CLI exposes two commands: `generic-to-dbsnp` and `validate-dbsnp`. Both
share the following options:

| Option | Description |
| --- | --- |
| `-v`, `--vcf` | Input VCF path (**required**) |
| `-o`, `--output` | Optional output dbSNP VCF path; must end with `.vcf` or `.vcf.gz` |
| `-m`, `--metadata` | Metadata file path (**required**) |
| `-e`, `--error-report` | Optional integrated validation report path, including reference results when used |
| `-r`, `--reference` | Optional reference FASTA for REF allele validation |

If `--output` is omitted, `generic-to-dbsnp` changes `sample.vcf.gz` to
`sample.dbsnp.vcf.gz`, while `validate-dbsnp` changes it to
`sample.cleaned.vcf.gz`. The same rule applies to uncompressed `.vcf` input.
Each run writes one validation report. When `--reference` is supplied, REF-check
statistics and findings are included in that same report. If `--error-report` is
omitted, the report uses the output VCF basename with `.errors.txt`. Generated
output and report paths must differ from every input path.

## Metadata Format

Metadata files use VCF-style lines:

```text
##Experiment_id=EXP001
##reference=toy_ref
##SampleSet_id=POP1
```

`SampleSet_id` in the metadata file is written to output VCF headers as
`##population_id=...`. The cleaned VCF output does not emit `##SampleSet_id=...`.
`Experiment_id` is written as `##batch=...`. A non-empty `Experiment_id` and
exactly one non-empty `SampleSet_id` are required.

One SNP VCF represents one SampleSet. During `generic-to-dbsnp`, all individual
sample columns are aggregated into that single SampleSet. During
`validate-dbsnp`, the input must contain exactly one `##population_id` and
exactly one population column after `FORMAT`; both must match the metadata
`SampleSet_id`.

## Project Structure

```
SNP/
├── README.md            # This file
├── requirements.txt
├── src/kvar_snp_tools/
│   ├── Sub_validator.py               # Public CLI entry point
│   ├── VCF2dbSNP.py                   # Generic VCF -> dbSNP VCF conversion
│   ├── dbsnp_vcf_cleaner.py           # dbSNP VCF validation & cleaning
│   ├── dbSNP_parser.py                # Streaming dbSNP VCF parser
│   ├── VCF_ref_check.py               # Optional REF check against a reference FASTA
│   ├── metadata_validator.py          # Metadata -> VCF header translation
│   └── error_handler.py               # Error codes and validation report
├── examples/            # Toy inputs for trying the commands
└── tests/               # CLI smoke tests
```

## Testing

```bash
python tests/test_public_cli_smoke.py
python tests/test_public_dbsnp_cleaner_streaming.py
```

## Notes

- `pyfaidx` is required only for reference FASTA validation.
- This public subset does not include private datasets, generated full-scale outputs, or internal pipeline reports.
