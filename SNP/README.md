# KVar-Toolkit SNP

Convert and validate SNP VCF files for dbSNP submission.

Run the commands below from the `SNP/` directory.

## Functions

| Command | Input | Output |
| --- | --- | --- |
| `generic-to-dbsnp` | Generic SNP VCF + metadata | dbSNP-formatted VCF + validation report |
| `validate-dbsnp` | dbSNP-formatted VCF + metadata | Cleaned dbSNP VCF + validation report |

Both commands support optional REF allele validation with `--reference`.

## Format Reference

[dbSNP VCF Submission Format Guidelines (PDF, 2013)](docs/dbSNP_VCF_Submission.pdf)

## Requirements

- Python 3.10 or later
- No third-party Python packages

## Quick Start

### Generic VCF to dbSNP VCF

```bash
python src/kvar_snp_tools/Sub_validator.py generic-to-dbsnp \
  -v examples/toy.human.snp.generic.vcf \
  -m examples/toy.human.snp.metadata.txt
```

Automatic outputs:

- `examples/toy.human.snp.generic.dbsnp.vcf`
- `examples/toy.human.snp.generic.dbsnp.errors.txt`

### Validate and Clean a dbSNP VCF

```bash
python src/kvar_snp_tools/Sub_validator.py validate-dbsnp \
  -v examples/toy.human.snp.dbsnp.vcf \
  -m examples/toy.human.snp.metadata.txt
```

Automatic outputs:

- `examples/toy.human.snp.dbsnp.cleaned.vcf`
- `examples/toy.human.snp.dbsnp.cleaned.errors.txt`

### Enable Reference Validation

```bash
python src/kvar_snp_tools/Sub_validator.py validate-dbsnp \
  -v examples/toy.human.snp.dbsnp.vcf \
  -m examples/toy.human.snp.metadata.txt \
  -r examples/toy.human.snp.reference.fasta
```

Use the FASTA used to create the VCF. Without `--reference`, the remaining
format, metadata, and submission checks still run. The FASTA index (`.fai`) is
built or refreshed automatically when reference validation is requested. If
the reference directory is read-only, the generated index is used in memory.

## Options

| Option | Required | Description |
| --- | --- | --- |
| `-v`, `--vcf` | Yes | Input `.vcf` or `.vcf.gz` |
| `-m`, `--metadata` | Yes | Submission metadata |
| `-o`, `--output` | No | Output `.vcf` or `.vcf.gz`; derived automatically when omitted |
| `-e`, `--error-report` | No | Integrated validation report; derived automatically when omitted |
| `-r`, `--reference` | No | Reference FASTA for REF allele validation |

Generated output and report paths must differ from all input paths.
Compression is preserved when the automatic output name is derived from a
`.vcf.gz` input.

## Metadata

Required format:

```text
##Experiment_id=EXP001
##SampleSet_id=POP1
```

Optional VCF reference header value:

```text
##reference=toy_ref
```

Rules:

- a non-empty `Experiment_id` is required;
- exactly one non-empty `SampleSet_id` is required;
- `Experiment_id` becomes the output `##batch`;
- `SampleSet_id` becomes the output `##population_id`;
- individual sample IDs are not stored in metadata; and
- metadata `reference`, when present, is checked against and written to the VCF
  `##reference` header.

## SampleSet and Frequency Handling

One VCF represents one SampleSet.

For `generic-to-dbsnp`, input INFO or sample genotype values are used to write
SampleSet-level `NA:FRQ` values.

When sample columns are used, all individual samples are aggregated into one
population column named with the metadata `SampleSet_id`. Sample column names
are not compared with metadata.

For `validate-dbsnp`, all three values must agree:

- the single metadata `SampleSet_id`;
- the single VCF `##population_id`; and
- the single population column after `FORMAT`.

A mismatch is a blocking error.

## Record Validation

- Local ID: non-empty, unique, and at most 64 characters
- ALT: exactly one allele
- REF/ALT bases: `A`, `T`, `G`, or `C`
- Indel: REF and ALT must share the leading base
- Allele length: the longer allele must not exceed 50 bp

Only these dbSNP submission INFO tags are retained:

`VRT`, `AF`, `AN`, `AC`, `AD`, `AA`, `CMT`, `LKO`, `NIO`, `OMIM`, `OMIA`,
`PMID`, `SAO`, `SSR`

Unsupported INFO tags are reported and removed from the generated VCF.

## Validation Result

- `Error` or `Critical`: blocks output publication
- `Warning` or `Info`: recorded without blocking a valid output
- One integrated report contains format, metadata, and optional reference
  findings
- Output is published only after the full conversion or validation succeeds

## Tests

```bash
python tests/test_public_cli_smoke.py
python tests/test_public_dbsnp_cleaner_streaming.py
```

The bundled human-style examples are synthetic and contain no real person or
cohort data.
