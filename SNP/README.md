# KVar SNP : SNP VCF Validation & dbSNP VCF Creation

## Overview

The **SNP** module of KVar validates SNP VCF input and produces cleaned,
dbSNP-formatted VCF output for public archive submission. It can convert a
generic VCF into a dbSNP-formatted VCF, validate and rewrite an existing dbSNP
VCF, and optionally verify REF alleles against a reference FASTA.

Run the commands below from this `SNP/` directory.

## Key Features

- **Generic to dbSNP conversion**: Convert a generic VCF into a dbSNP-formatted VCF.
- **dbSNP validation & cleaning**: Validate an input dbSNP VCF and rewrite it as a cleaned dbSNP VCF.
- **Reference allele validation**: Optionally validate REF alleles against a reference FASTA before writing output.
- **Metadata-driven headers**: Translate VCF-style metadata into output VCF headers.

### Prerequisites

**Runtime:**

- Python 3.8 or higher

**Python packages:**

- `pyfaidx` (>= 0.8) — required only for reference FASTA validation

**Optional external tools:**

- `bcftools` — required only when using the normalization helper directly

### Installation

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

## Quick Start

Convert a generic VCF into a dbSNP-formatted VCF:

```bash
python src/kvar_snp_tools/Sub_validator.py generic-to-dbsnp \
  -v examples/toy.generic.vcf \
  -m examples/toy.metadata.txt \
  -o examples/toy.generic.cleaned.dbsnp.vcf \
  -e examples/toy.generic.errors.txt
```

Validate and clean an existing dbSNP VCF:

```bash
python src/kvar_snp_tools/Sub_validator.py validate-dbsnp \
  -v examples/toy.dbsnp.vcf \
  -m examples/toy.metadata.txt \
  -o examples/toy.dbsnp.cleaned.vcf \
  -e examples/toy.dbsnp.errors.txt
```

Reference validation can be added to either command:

```bash
python src/kvar_snp_tools/Sub_validator.py validate-dbsnp \
  -v examples/toy.dbsnp.vcf \
  -r examples/toy.reference.fa \
  -o examples/toy.dbsnp.cleaned.vcf
```

## Common Options

The CLI exposes two commands: `generic-to-dbsnp` and `validate-dbsnp`. Both
share the following options:

| Option | Description |
| --- | --- |
| `-v`, `--vcf` | Input VCF path (**required**) |
| `-o`, `--output` | Output dbSNP VCF path (**required**) |
| `-m`, `--metadata` | Metadata file path |
| `-e`, `--error-report` | Validation report path |
| `-r`, `--reference` | Reference FASTA for REF allele validation |
| `-rr`, `--reference-report` | Reference validation report path (used with `--reference`) |

`generic-to-dbsnp` additionally accepts:

| Option | Description |
| --- | --- |
| `-c`, `--preserve-contig` | Preserve input `##contig` header lines in the output VCF |

## Metadata Format

Metadata files use VCF-style lines:

```text
##Experiment_id=EXP001
##bioproject_id=PRJNA000001
##biosample_id=SAMD000001
##reference=toy_ref
##SampleSet_id=POP1
```

`SampleSet_id` in the metadata file is written to output VCF headers as
`##population_id=...`. The cleaned VCF output does not emit `##SampleSet_id=...`.

## Project Structure

```
SNP/
├── README.md            # This file
├── requirements.txt
├── src/kvar_snp_tools/
│   ├── Sub_validator.py               # Public CLI entry point
│   ├── VCF2dbSNP.py                    # Generic VCF -> dbSNP conversion
│   ├── dbsnp_vcf_cleaner.py            # dbSNP VCF validation & cleaning
│   ├── VCF_ref_check.py                # REF allele validation vs FASTA
│   ├── dbSNP_parser.py
│   ├── metadata_validator.py
│   ├── norm_VCF.py
│   └── error_handler.py
├── examples/            # Toy inputs for trying the commands
└── tests/               # CLI smoke tests
```

## Testing

```bash
python -m pytest tests/
```

## Notes

- `bcftools` is required only when using the normalization helper directly.
- `pyfaidx` is required only for reference FASTA validation.
- This public subset does not include private datasets, generated full-scale outputs, or internal pipeline reports.
