[![Tests](https://github.com/pha4ge/primaschema/actions/workflows/test.yml/badge.svg)](https://github.com/pha4ge/primaschema/actions/workflows/test.yml) [![PyPI version](https://badge.fury.io/py/primaschema.svg)](https://pypi.org/project/primaschema)

# Primaschema

## Dev notes

`uv run pytest -v test/test_all.py::test_dev_scheme`

A toolkit for fetching, validating and interrogating tiled amplicon PCR primer scheme definitions. Provides convenient programmatic accesss to the [PHA4GE primer-schemes repository](https://github.com/pha4ge/primer-schemes), a community repository of tiled amplicons primer schemes.

## Install (Python 3.12+)

```shell
uv tool install primaschema
```

### Development

```shell
git clone https://github.com/pha4ge/primaschema.git
cd primaschema
uv sync --all-extras
uv run primaschema --help
uv run pytest
uv run pre-commit install
uv run pre-commit run --all-files
```

`uv sync --all-extras` installs optional dependencies, including the `dev` extra (e.g. `pytest`, `pre-commit`, `ruff`) defined in `pyproject.toml`.

Some Primaschema commands use components from the [primer-schemes](https://github.com/pha4ge/primer-schemes) repository. To show Primaschema where to find these, create the environment variable `PRIMER_SCHEMES_PATH` pointing to the location of the primer-schemes directory on your machine:

```shell
git clone https://github.com/pha4ge/primer-schemes.git
export PRIMER_SCHEMES_PATH="/path/to/primer-schemes"
```



## Usage

### Scheme creation

```bash
mkdir -p built && rm -rf built/artic && uv run primaschema create \
  --name artic \
  --amplicon-size 400 \
  --version v4.1.0 \
  --contributors "ARTIC network" \
  --target-organisms "sars-cov-2" \
  --status DEPRECATED \
  --bed-path test/data/dev-scheme/primer.bed \
  --reference-path test/data/dev-scheme/reference.fasta \
  --primer-schemes-path built
```

```
% primaschema -h
usage: primaschema [-h] [--version]
                   {validate,build,build-manifest,hash-ref,hash-bed,diff,6to7,7to6,plot,show-intervals,show-discordant-primers,subset,sync} ...

positional arguments:
  {validate,build,build-manifest,hash-ref,hash-bed,diff,6to7,7to6,plot,show-intervals,show-discordant-primers,subset,sync}
    validate            Validate one or more primer scheme definitions comprising info.yml, primer.bed and reference.fasta
    build               Build one or more primer scheme definitions comprising info.yml, primer.bed and reference.fasta
    build-manifest      Build a complete manifest of schemes contained in the specified directory
    hash-ref            Generate reference sequence checksum
    hash-bed            Generate a bed file checksum
    diff                Show the symmetric difference of records in two bed files
    6to7                Convert a 6 column scheme.bed file to a 7 column primer.bed file using reference backfill
    7to6                Convert a 7 column primer.bed file to a 6 column scheme.bed file by removing a column
    plot                Plot amplicon and primer coords from 7 column primer.bed
    show-intervals      Show amplicon start and end coordinates given a BED file of primer coordinates
    show-discordant-primers
                        Show primer records with sequences not matching the reference sequence
    subset              Extract a primer.bed and reference.fasta scheme subset for a single chromosome
    sync                Retrieve/update local copy of remote primer scheme repository

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```
