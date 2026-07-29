# KVar-Toolkit SV

Convert an SV VCF into the Call TSV used for submission before administrator
QC.

Run the commands below from the `SV/` directory.

## Function

| Input | Output |
| --- | --- |
| SV VCF (`.vcf` or `.vcf.gz`) | Accession-free submission Call TSV + validation report |

## Requirements

- Python 3.10 or later
- Python standard library
- Reference FASTA only when `--reference` is used; its `.fai` index is built
  automatically when missing or stale

## Quick Start

```bash
python src/kvar_sv_tools/vcf_to_kvar_tsv.py \
  -v examples/toy.human.sv.vcf
```

Automatic outputs:

- `examples/toy.human.sv.call.tsv`
- `examples/toy.human.sv.call.errors.txt`

Enable optional reference validation:

```bash
python src/kvar_sv_tools/vcf_to_kvar_tsv.py \
  -v examples/toy.human.sv.vcf \
  -r examples/toy.human.sv.reference.fasta
```

Use the FASTA used to create the VCF. Without `--reference`, all remaining VCF
format and submission checks still run. If the reference directory is
read-only, the generated index is used in memory.

## Options

| Option | Required | Description |
| --- | --- | --- |
| `-v`, `--vcf` | Yes | Input SV VCF; gzip input is supported |
| `-r`, `--reference` | No | Reference FASTA for chromosome, coordinate, and REF checks; `.fai` is built or refreshed automatically |
| `-t`, `--call-tsv` | No | Output Call TSV; must end with `.tsv` |
| `-e`, `--error-report` | No | Validation report path |

Automatic naming:

- `sample.vcf` or `sample.vcf.gz` → `sample.call.tsv`
- `sample.call.tsv` → `sample.call.errors.txt`

Generated output and report paths must differ from the VCF, optional FASTA and
FASTA index, and each other.

## Submission Call TSV

The first line is:

```text
##Variant_Call
```

The second line contains exactly these 32 columns:

```text
#Variant_Call_ID	Variant_Call_Type	Chr	Outer_Start	Start	Inner_Start	Inner_Stop	Stop	Outer_Stop	Insertion_Length	Allele_Count	Allele_Frequency	Allele_Number	Copy_Number	Description	Validation	Zygosity	Origin	Phenotype	External_Links	Evidence	Sequence	From_Chr	From_Coord	From_Strand	To_Chr	To_Coord	To_Strand	Mutation_ID	Mutation_Order	Mutation_Molecule	BND_Source_VCF_IDs
```

Output rules:

- `Variant_Call_ID` retains the submitted VCF ID for one-to-one and collapsed
  calls; a BND-derived insertion receives a deterministic ID such as
  `bndA_ins`.
- Ordinary calls and BND rows that remain separate use `.` in
  `BND_Source_VCF_IDs`.
- A collapsed BND call records both source VCF IDs in
  `BND_Source_VCF_IDs`.
- A BND-derived insertion records only the VCF ID or IDs that supplied its
  inserted sequence.
- Blocking validation errors prevent Call TSV publication.
- The validation report is retained even when Call TSV publication is blocked.

## Core Validation

Without a reference FASTA:

- non-empty chromosome names;
- positive coordinates;
- chromosome and position ordering; and
- chromosome membership and bounds when VCF `##contig` declarations exist.

With `--reference`:

- chromosome-name resolution;
- local POS/REF and non-BND END coordinate bounds;
- derived END bounds, including BND-derived insertion calls;
- paired-BND target contig and coordinate bounds, allowing `length+1` for the
  right telomere; and
- REF allele agreement with the indexed FASTA.

Endpoint bounds are independent of `SVLEN`. For example, an insertion with
`POS=END` and a non-zero `SVLEN` is valid.

## SV and BND Handling

- Supported classifications include deletion, insertion, duplication,
  inversion, CNV, mobile-element events, BND, and complex events.
- A BND pair collapses only when `MATEID`, target coordinates, and strand
  orientations are reciprocal.
- BND ordering remains deterministic across different VCF row orders.
- Inserted BND sequence is retained as a linked insertion call.
- Sequence on only one mate is retained and reported as
  `MATEID_INSERTION_SEQUENCE_ONE_SIDED`.
- Conflicting mate sequences block output with
  `MATEID_INSERTION_SEQUENCE_MISMATCH`.
- Multiple ALT values, multiple MATEID values for one ALT, and self-referencing
  MATEID values are blocking errors.
- A BND-derived insertion uses `Start=POS` and
  `Stop=POS+Insertion_Length`.
- An ordinary insertion preserves `END`; when `END` is omitted, it uses
  `Start=Stop=POS`.
- `CPX_TYPE` supplies the `sequence alteration` fallback for an otherwise
  unclassified record.

## Tests

```bash
python tests/test_bnd_mate_insertions.py
python tests/test_reference_bounds_and_index.py
python tests/test_public_cli_smoke.py
```

The bundled GRCh38-style example is synthetic and contains no real person,
sample, genotype, or cohort data.
