# KVar-Toolkit: Variant QC & Validation Toolkit

## Overview

**KVar-Toolkit** is a toolkit for quality control (QC) and submission validation
of genetic variant data managed by KOBIC. It standardizes variant files into
the formats required for public archive submission and verifies that records
are internally consistent.

Reference-genome checks are optional; without a supplied reference, the toolkit
performs the remaining format and submission validations. The toolkit is
organized into the **SNP** and **SV** modules. Each module is self-contained and
can be used independently from its own subdirectory.

## Modules

| Module | Description |
| --- | --- |
| [`SNP/`](SNP/) | SNP VCF validation and dbSNP-formatted VCF creation, including optional REF allele checking and metadata validation. |
| [`SV/`](SV/) | Structural variation (SV) VCF validation and accession-free submission `Variant_Call.tsv` creation, with optional reference-based coordinate/REF validation. |

## Key Features

- **Generic VCF to dbSNP conversion**: Rewrites a generic VCF into a dbSNP-formatted VCF for submission.
- **dbSNP VCF validation & cleaning**: Validates an input dbSNP VCF and emits a cleaned, standardized VCF.
- **Reference allele validation**: Optionally checks REF alleles against a reference FASTA.
- **SNP metadata mapping**: Validates SNP submission metadata and writes the corresponding dbSNP VCF headers.
- **Unified validation reporting**: Produces one report per run, including optional reference-validation findings when enabled.
- **SV VCF to Variant Call TSV conversion**: Rewrites a structural-variation VCF into `Variant_Call.tsv` for submission.

The SNP module supports exactly one SampleSet per VCF. Generic VCF sample
columns are aggregated into that SampleSet, and existing dbSNP VCF input must
contain one matching population ID and population column.

### Prerequisites

**Operating System:**

- Linux (CentOS 7+, Ubuntu 18.04+, Debian 9+)

**Runtime:**

- Python 3.10 or higher (tested with Python 3.10, 3.12, and 3.14)

**Python packages:**

- `pyfaidx` (>= 0.8) — only for optional SNP reference FASTA validation
- The SV module has no third-party dependencies

### Installation

1. Clone the repository

```bash
git clone https://github.com/KOBIC-KBDS/KOBIC_KVar.git
cd KOBIC_KVar
```

2. Optional: install `pyfaidx` only when using SNP reference validation

```bash
pip install -r SNP/requirements.txt   # Only when using SNP reference validation
```

Reference validation is optional in both modules. SNP reference validation uses
`pyfaidx`; SV reference validation requires no additional Python package and only
needs an indexed FASTA (`.fai`) when `--reference` is supplied.

## Quick Start

Run the SNP module from its directory. Convert a generic VCF into a
dbSNP-formatted VCF:

```bash
cd SNP
python src/kvar_snp_tools/Sub_validator.py generic-to-dbsnp \
  -v examples/toy.generic.vcf \
  -m examples/toy.metadata.txt
```

Without `--output`, this command creates `examples/toy.generic.dbsnp.vcf`.

See [SNP/README.md](SNP/README.md) for the full command reference, metadata
format, and reference-validation options.

For the SV module, convert a structural-variation VCF into `Variant_Call.tsv`:

```bash
cd SV
python src/kvar_sv_tools/vcf_to_kvar_tsv.py \
  -v examples/toy.human.sv.vcf
```

Without `--call-tsv`, this command creates
`examples/toy.human.sv.Variant_Call.tsv`.

See [SV/README.md](SV/README.md) for the full SV command reference.

SV submission conversion preserves submitted IDs and does not assign `kssvN`
accessions. It does not accept a separate metadata file or require `SampleSet`,
`Experiment`, or the VCF `##reference` header. Reference validation uses the
optional FASTA supplied through `--reference`. The converter processes
every VCF record and expects users to provide an already filtered submission VCF.

## Project Structure

```
KOBIC_KVar/
├── README.md            # This file
├── LICENSE
├── .gitignore
├── SNP/                 # SNP validation module
│   ├── README.md        # SNP module documentation
│   ├── requirements.txt
│   ├── src/kvar_snp_tools/
│   │   └── Sub_validator.py               # Public CLI entry point
│   ├── examples/        # Toy inputs for trying the commands
│   └── tests/           # CLI smoke tests
└── SV/                  # Structural variation module
    ├── README.md        # SV module documentation
    ├── examples/        # Privacy-safe synthetic human example
    ├── src/kvar_sv_tools/
    │   └── vcf_to_kvar_tsv.py             # Public CLI entry point
    └── tests/           # CLI smoke tests
```

## Documentation

- **[SNP module](SNP/README.md)**: Workflows, command reference, metadata format, and notes.
- **[SV module](SV/README.md)**: SV VCF → `Variant_Call.tsv` conversion, command reference, and input/output details.

## Testing

Each module ships with its own tests. For the SNP module:

```bash
cd SNP
python tests/test_public_cli_smoke.py
python tests/test_public_dbsnp_cleaner_streaming.py
```

For the SV module:

```bash
cd SV
python tests/test_public_cli_smoke.py
```

## Support

- **Issues**: [GitHub Issues](https://github.com/KOBIC-KBDS/KOBIC_KVar/issues)

## Acknowledgments

Developed and maintained by the Korea Bioinformation Center (KOBIC).
