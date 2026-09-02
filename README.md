[![Tests](https://github.com/pha4ge/primaschema/actions/workflows/test.yml/badge.svg)](https://github.com/pha4ge/primaschema/actions/workflows/test.yml) [![PyPI version](https://badge.fury.io/py/primaschema.svg)](https://pypi.org/project/primaschema)

# Primaschema

A toolkit for fetching, validating and interrogating tiled amplicon PCR primer scheme definitions. Provides convenient programmatic access to the [PHA4GE primer-schemes repository](https://github.com/pha4ge/primer-schemes), a community repository of tiled amplicon primer schemes.

## Production install (Python 3.12+)

Recommended:

```shell
uv tool install primaschema
```

Alternatives:

```shell
pipx install primaschema
# or
python -m pip install primaschema
```

## Quick start

Download a scheme from the default index and validate it:

```shell
primaschema get example-scheme/400/v1.0.0 --output ./schemes
primaschema validate ./schemes/example-scheme/400/v1.0.0/info.json --checks core
```

Show CLI help or version:

```shell
primaschema --help
primaschema --show-version
```

## Common commands

- `create`: Create a new scheme from a `primer.bed` and `reference.fasta`.
- `validate`: Validate a single scheme or all schemes under a directory.
- `rebuild`: Recompute checksums, regenerate scheme README, and optionally reformat `primer.bed`.
- `index`: Build or update an `index.json` for a local schemes directory.
- `get`: Download schemes from an index (default is the PHA4GE primer-schemes index).
- `modify`: Update contributors, vendors, status, license, target organisms, dates and the generator.

## Scheme creation

```shell
mkdir -p ./schemes && rm -rf ./schemes/example-scheme
primaschema create \
  --name example-scheme \
  --amplicon-size 400 \
  --version v1.0.0 \
  --contributor "name=Alice Smith,email=alice@example.org" \
  --target-organism "name=Example organism,ncbi_taxon_id=000001" \
  --development-status VALIDATED \
  --bed-path ./example-scheme.primer.bed \
  --reference-path ./example-scheme.reference.fasta \
  --primer-schemes-path ./schemes
```

The scheme is written to `./schemes/example-scheme/400/v1.0.0/`.

## Environment

Many commands accept `--primer-schemes-path`. You can also set it once:

```shell
export PRIMER_SCHEMES_PATH=./schemes
```

## Development

```shell
git clone https://github.com/pha4ge/primaschema.git
cd primaschema
uv sync --all-extras --all-groups
uv run primaschema --help
uv run pytest
uv run pre-commit install
uv run pre-commit run --all-files
```

`uv sync --all-extras --all-groups` installs optional dependencies, including the `dev` group (e.g. `pytest`, `pre-commit`, `ruff`) and the `docs` group (mkdocs and its plugins) defined in `pyproject.toml`. `--all-extras` alone does not install dependency groups, so omitting `--all-groups` leaves the docs unbuildable.

## Documentation

Documentation is built with [MkDocs](https://www.mkdocs.org):

```shell
uv run mkdocs serve   # live preview on http://127.0.0.1:8000
uv run mkdocs build   # render to ./site
```

The CLI and schema pages are generated from the code, so they track the model automatically.

The Pydantic model (`src/primaschema/schema/info.py`) is generated from the LinkML schema (`src/primaschema/schema/info.yml`). After modifying the schema, regenerate with:

```shell
uv run gen-pydantic src/primaschema/schema/info.yml --meta None > src/primaschema/schema/info.py
```
