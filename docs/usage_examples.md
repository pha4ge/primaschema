# Usage Examples

## create

Create a new primer scheme from an existing `primer.bed` and `reference.fasta`.
Contributors and target organisms are passed as `key=value` pairs.

```bash
primaschema create \
    --name example-scheme \
    --amplicon-size 400 \
    --version v1.0.0 \
    --status validated \
    --contributors "name=Alice Smith,email=alice@example.org" \
    --target-organisms "common_name=Example organism,ncbi_tax_id=000001" \
    --bed-path ./example-scheme.primer.bed \
    --reference-path ./example-scheme.reference.fasta \
    --primer-schemes-path ./schemes
```

The scheme is written to `./schemes/example-scheme/400/v1.0.0/`.

---

## validate

Validate a single scheme:

```bash
primaschema validate ./schemes/example-scheme/400/v1.0.0/info.json
```

Validate all schemes under a directory:

```bash
primaschema validate --all ./schemes
```

Auto-normalise `primer.bed` and `reference.fasta` formatting in place when
hashes differ only by line order:

```bash
primaschema validate --fix ./schemes/example-scheme/400/v1.0.0/info.json
```

---

## rebuild

Recompute checksums and regenerate `README.md` for a single scheme:

```bash
primaschema rebuild ./schemes/example-scheme/400/v1.0.0/info.json
```

Rebuild all schemes, also resorting `primer.bed` lines:

```bash
primaschema rebuild --all --reformat-primer-bed ./schemes
```

---

## index

Build an `index.json` from a local schemes directory:

```bash
primaschema index ./schemes
```

Attach a base URL so the index contains direct download links, and write the
index to a specific output directory:

```bash
primaschema index \
    --base-url https://raw.githubusercontent.com/example-org/primer-schemes/main/schemes \
    --output-path ./dist \
    ./schemes
```

Update an existing index rather than rebuilding from scratch:

```bash
primaschema index \
    --index-path ./dist/index.json \
    --base-url https://raw.githubusercontent.com/example-org/primer-schemes/main/schemes \
    --output-path ./dist \
    ./schemes
```

---

## get

Download a specific scheme by `name/amplicon_size/version`:

```bash
primaschema get example-scheme/400/v1.0.0
```

Download to a specific output directory:

```bash
primaschema get example-scheme/400/v1.0.0 --output ./downloaded-schemes
```

Download all schemes matching a partial identifier:

```bash
primaschema get example-scheme/400 --allow-multiple --output ./downloaded-schemes
```

Download every scheme in the index:

```bash
primaschema get --all --output ./downloaded-schemes
```

Use a local or custom index:

```bash
primaschema get example-scheme/400/v1.0.0 --index ./dist/index.json
```

---

## modify add-contributor

Append a contributor (positional `INFO-PATH` then `NAME`, optional flags):

```bash
primaschema modify add-contributor \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    "Alice Smith" \
    --email alice@example.org
```

Insert at a specific index (0-based):

```bash
primaschema modify add-contributor \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    "Alice Smith" \
    --email alice@example.org \
    --idx 0
```

---

## modify remove-contributor

Remove the contributor at index 1:

```bash
primaschema modify remove-contributor \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    1
```

---

## modify update-contributor

Replace the contributor at index 0 (positional `INFO-PATH IDX NAME`, optional flags):

```bash
primaschema modify update-contributor \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    0 \
    "Alice Smith" \
    --email alice.smith@example.org \
    --orcid-id "0000-0001-2345-6789"
```

---

## modify add-vendor

Add a vendor (positional `INFO-PATH ORGANISATION-NAME`, optional flags):

```bash
primaschema modify add-vendor \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    "ExampleCo Ltd" \
    --home-page "https://example.com" \
    --kit-name "ExamplePanel v1"
```

---

## modify remove-vendor

Remove the vendor at index 0:

```bash
primaschema modify remove-vendor \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    0
```

---

## modify update-vendor

Replace the vendor at index 0 (positional `INFO-PATH IDX ORGANISATION-NAME`, optional flags):

```bash
primaschema modify update-vendor \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    0 \
    "ExampleCo Ltd" \
    --home-page "https://example.com" \
    --kit-name "ExamplePanel v2"
```

---

## modify add-tag

```bash
primaschema modify add-tag \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    clinical
```

---

## modify remove-tag

```bash
primaschema modify remove-tag \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    clinical
```

---

## modify update-license

```bash
primaschema modify update-license \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    cc-by-sa-4-full-stop0
```

---

## modify update-status

```bash
primaschema modify update-status \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    validated
```

---

## modify update-algorithm

Positional `INFO-PATH ALGORITHM.NAME`, optional `--algorithm.version`:

```bash
primaschema modify update-algorithm \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    primalscheme \
    --algorithm.version 3.0.3
```

---

## modify add-target-organism

```bash
primaschema modify add-target-organism \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    --common-name "Example organism" \
    --ncbi-tax-id "000001"
```

---

## modify remove-target-organism

Remove the target organism at index 0:

```bash
primaschema modify remove-target-organism \
    ./schemes/example-scheme/400/v1.0.0/info.json \
    0
```
