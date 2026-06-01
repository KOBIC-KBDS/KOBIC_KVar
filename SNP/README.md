# KVar SNP Public Validation Tools

This folder contains the public-facing subset for SNP VCF validation and dbSNP VCF creation.
Run the commands below from this `SNP/` directory.

## Included Workflows

1. Convert a generic VCF into a dbSNP-formatted VCF.
2. Validate an input dbSNP VCF and rewrite it as a cleaned dbSNP VCF.
3. Optionally validate REF alleles against a reference FASTA before writing output.

## Layout

```text
src/kvar_snp_tools/
  VCF2dbSNP.py
  VCF_ref_check.py
  dbSNP_parser.py
  dbsnp_vcf_cleaner.py
  error_handler.py
  metadata_validator.py
  norm_VCF.py
  run_submission_validation.py
examples/
tests/
```

## Metadata Format

Metadata files use VCF-style lines:

```text
##Experiment_id=EXP001
##bioproject_id=PRJNA000001
##biosample_id=SAMD000001
##reference=toy_ref
##SampleSet_id=POP1
```

`SampleSet_id` in the metadata file is written to output VCF headers as `##population_id=...`.
The cleaned VCF output does not emit `##SampleSet_id=...`.

## Commands

Convert a generic VCF:

```bash
python src/kvar_snp_tools/run_submission_validation.py generic-to-dbsnp \
  --vcf examples/toy.generic.vcf \
  --metadata examples/toy.metadata.txt \
  --output examples/toy.generic.cleaned.dbsnp.vcf \
  --error-report examples/toy.generic.errors.txt
```

Validate and clean a dbSNP VCF:

```bash
python src/kvar_snp_tools/run_submission_validation.py validate-dbsnp \
  --vcf examples/toy.dbsnp.vcf \
  --metadata examples/toy.metadata.txt \
  --output examples/toy.dbsnp.cleaned.vcf \
  --error-report examples/toy.dbsnp.errors.txt
```

Reference validation can be added to either command:

```bash
python src/kvar_snp_tools/run_submission_validation.py validate-dbsnp \
  --vcf examples/toy.dbsnp.vcf \
  --reference examples/toy.reference.fa \
  --output examples/toy.dbsnp.cleaned.vcf
```

## Notes

- `bcftools` is required only when using the normalization helper directly.
- `pyfaidx` is required only for reference FASTA validation.
- This public subset does not include private datasets, generated full-scale outputs, or internal pipeline reports.
