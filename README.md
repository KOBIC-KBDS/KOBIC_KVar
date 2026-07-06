# KVar-Toolkit: Variant QC & Validation Toolkit

## Overview

**KVar-Toolkit** is a toolkit for quality control (QC) and submission validation
of genetic variant data managed by KOBIC. It standardizes variant files into
the formats required for public archive submission and verifies that the records
are internally consistent and consistent with a reference genome.

The toolkit is organized by variant type. Each module is self-contained and can
be used independently from its own subdirectory. The current release provides
the **SNP** module; additional modules can be added under the same structure.

## Modules

| Module | Description |
| --- | --- |
| [`SNP/`](SNP/) | SNP VCF validation and dbSNP-formatted VCF creation, including REF allele checking against a reference FASTA and metadata validation. |

## Key Features

- **Generic VCF to dbSNP conversion**: Rewrites a generic VCF into a dbSNP-formatted VCF for submission.
- **dbSNP VCF validation & cleaning**: Validates an input dbSNP VCF and emits a cleaned, standardized VCF.
- **Reference allele validation**: Optionally checks REF alleles against a reference FASTA.
- **Metadata validation**: Reads VCF-style metadata and writes the corresponding output headers.
- **Validation reporting**: Produces optional error/validation reports for each run.

### Prerequisites

**Operating System:**

- Linux (CentOS 7+, Ubuntu 18.04+, Debian 9+)

**Runtime:**

- Python 3.8 or higher

**Python packages:**

- `pyfaidx` (>= 0.8) — required only for reference FASTA validation

### Installation

1. Clone the repository

```bash
git clone https://github.com/KOBIC-KBDS/KVar-Toolkit.git
cd KVar-Toolkit
```

2. Install Python dependencies for the module you want to use

```bash
pip install -r SNP/requirements.txt
```

## Quick Start

Run the SNP module from its directory. Convert a generic VCF into a
dbSNP-formatted VCF:

```bash
cd SNP
python src/kvar_snp_tools/Sub_validator.py generic-to-dbsnp \
  -v examples/toy.generic.vcf \
  -m examples/toy.metadata.txt \
  -o examples/toy.generic.cleaned.dbsnp.vcf \
  -e examples/toy.generic.errors.txt
```

See [SNP/README.md](SNP/README.md) for the full command reference, metadata
format, and reference-validation options.

## Project Structure

```
KVar-Toolkit/
├── README.md            # This file
├── LICENSE
├── .gitignore
└── SNP/                 # SNP validation module
    ├── README.md        # SNP module documentation
    ├── requirements.txt
    ├── src/kvar_snp_tools/
    │   └── Sub_validator.py               # Public CLI entry point
    ├── examples/        # Toy inputs for trying the commands
    └── tests/           # CLI smoke tests
```

## Documentation

- **[SNP module](SNP/README.md)**: Workflows, command reference, metadata format, and notes.

## Testing

Each module ships with its own tests. For the SNP module:

```bash
cd SNP
python tests/test_public_cli_smoke.py
```

## Support

- **Issues**: [GitHub Issues](https://github.com/KOBIC-KBDS/KVar-Toolkit/issues)

## Acknowledgments

Developed and maintained by the Korea Bioinformation Center (KOBIC).
