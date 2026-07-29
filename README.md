# KVar-Toolkit

Public command-line tools for KOBIC variant submission conversion and
validation.

## Modules

| Module | Main function | Required auxiliary input | Optional reference check |
| --- | --- | --- | --- |
| [SNP](SNP/) | Generic VCF → dbSNP VCF; validate and clean an existing dbSNP VCF | Metadata with a non-empty `Experiment_id` and exactly one `SampleSet_id` | REF allele check with `--reference` |
| [SV](SV/) | SV VCF → pre-QC submission Call TSV | None | Chromosome, coordinate, and REF checks with `--reference` |

The modules run independently. Reference validation is optional in both
modules.

## Requirements

- Linux
- Python 3.10 or later
- No third-party Python packages

## Installation

```bash
git clone https://github.com/KOBIC-KBDS/KOBIC_KVar.git
cd KOBIC_KVar
```

## SNP Quick Start

Convert a generic VCF:

```bash
cd SNP
python src/kvar_snp_tools/Sub_validator.py generic-to-dbsnp \
  -v examples/toy.human.snp.generic.vcf \
  -m examples/toy.human.snp.metadata.txt
```

Automatic outputs:

- `examples/toy.human.snp.generic.dbsnp.vcf`
- `examples/toy.human.snp.generic.dbsnp.errors.txt`

Validate and clean a dbSNP VCF:

```bash
python src/kvar_snp_tools/Sub_validator.py validate-dbsnp \
  -v examples/toy.human.snp.dbsnp.vcf \
  -m examples/toy.human.snp.metadata.txt
```

Automatic outputs:

- `examples/toy.human.snp.dbsnp.cleaned.vcf`
- `examples/toy.human.snp.dbsnp.cleaned.errors.txt`

See [SNP/README.md](SNP/README.md) for metadata, frequency handling, and
validation rules.

## SV Quick Start

Convert an SV VCF:

```bash
cd SV
python src/kvar_sv_tools/vcf_to_kvar_tsv.py \
  -v examples/toy.human.sv.vcf
```

Automatic outputs:

- `examples/toy.human.sv.call.tsv`
- `examples/toy.human.sv.call.errors.txt`

The public SV output is an accession-free submission Call TSV created before
administrator QC. See [SV/README.md](SV/README.md) for the exact schema and
validation rules.

## Optional Reference Validation

Add `--reference` with the FASTA used to create the input VCF:

```bash
# SNP
python SNP/src/kvar_snp_tools/Sub_validator.py validate-dbsnp \
  -v SNP/examples/toy.human.snp.dbsnp.vcf \
  -m SNP/examples/toy.human.snp.metadata.txt \
  -r SNP/examples/toy.human.snp.reference.fasta

# SV
python SV/src/kvar_sv_tools/vcf_to_kvar_tsv.py \
  -v SV/examples/toy.human.sv.vcf \
  -r SV/examples/toy.human.sv.reference.fasta
```

Without `--reference`, all non-reference format and submission checks still
run. When supplied, the FASTA index (`.fai`) is built or refreshed
automatically. If the reference directory is read-only, the generated index is
used in memory.

## Tests

```bash
cd SNP
python tests/test_public_cli_smoke.py
python tests/test_public_dbsnp_cleaner_streaming.py

cd ../SV
python tests/test_bnd_mate_insertions.py
python tests/test_reference_bounds_and_index.py
python tests/test_public_cli_smoke.py
```

The bundled human-style examples are synthetic and contain no real person or
cohort data.

## Support

- [GitHub Issues](https://github.com/KOBIC-KBDS/KOBIC_KVar/issues)

Developed and maintained by the Korea Bioinformation Center (KOBIC).
