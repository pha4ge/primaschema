import pytest

from primaschema.schema.flatten import (
    CSV_FIELDNAMES,
    _pack,
    _unpack,
    flatten_scheme,
    unflatten_scheme,
)
from primaschema.schema.info import (
    PrimerSchemeContributor,
    PrimerSchemeTargetOrganism,
    PrimerSchemeVendor,
)
from primaschema.schema.primer_scheme import PrimerScheme

# --- _pack / _unpack ---


def test_pack_empty_list():
    """Packing an empty list produces an empty string."""
    assert _pack([]) == ""
    assert _unpack("") == []


def test_pack_single_item():
    """A single present value packs and unpacks unchanged."""
    assert _pack(["x"]) == "x"
    assert _unpack("x") == ["x"]


def test_pack_single_empty_item():
    """A single missing value packs to '' (not csv.writer's defensive '""')."""
    assert _pack([None]) == ""


def test_pack_missing_value_start_middle_end():
    """Missing values render as adjacent ';;' wherever they occur positionally."""
    assert _pack(["Alice", None, "Bob"]) == "Alice;;Bob"
    assert _pack([None, "Bob"]) == ";Bob"
    assert _pack(["Alice", None]) == "Alice;"
    assert _pack([None, None]) == ";"


def test_pack_unpack_value_with_embedded_delimiter():
    """A value containing a literal ';' round-trips correctly via CSV quoting."""
    packed = _pack(["Ac;me Inc", "Bob"])
    assert packed == '"Ac;me Inc";Bob'
    assert _unpack(packed) == ["Ac;me Inc", "Bob"]


def test_pack_unpack_value_with_embedded_newline():
    """A value containing a literal newline round-trips correctly, not silently truncated.

    Regression test: _pack previously used lineterminator="" to avoid a
    trailing terminator, which also disabled csv.writer's quoting of
    embedded '\\n'/'\\r' (that quoting is keyed off the dialect's
    lineterminator characters). An unquoted embedded newline meant _unpack's
    next(csv.reader(...)) only read the first physical line, silently
    dropping everything after the newline.
    """
    packed = _pack(["Alice", "Bob\nX", "Carol"])
    assert packed == 'Alice;"Bob\nX";Carol'
    assert _unpack(packed) == ["Alice", "Bob\nX", "Carol"]


def test_pack_unpack_single_item_with_embedded_newline():
    """The single-item packing path also correctly quotes an embedded newline."""
    packed = _pack(["Multi-line\nnote"])
    assert _unpack(packed) == ["Multi-line\nnote"]


def test_pack_unpack_value_with_embedded_carriage_return():
    """A value containing a literal carriage return round-trips instead of crashing.

    Regression test: previously raised an unhandled csv.Error ("new-line
    character seen in unquoted field"), which is not a ValueError and so
    bypassed the CLI's clean-error-message handling entirely.
    """
    packed = _pack(["Old Mac\rDonald Ltd", "Bob"])
    assert _unpack(packed) == ["Old Mac\rDonald Ltd", "Bob"]


# --- flatten_scheme / unflatten_scheme round trip ---


def _sample_scheme(**kwargs) -> PrimerScheme:
    from datetime import date

    from primaschema.schema.info import PrimerSchemeDevelopmentStatus

    defaults = dict(
        schema_version="1.0.0",
        primer_scheme_name="test-scheme",
        amplicon_size=400,
        primer_scheme_version="v1.0.0",
        primer_scheme_contributor=[
            PrimerSchemeContributor(
                primer_scheme_contributor_name="Alice",
                primer_scheme_contributor_email="alice@x.org",
            ),
            PrimerSchemeContributor(primer_scheme_contributor_name="Bob"),
        ],
        primer_scheme_target_organism=[
            PrimerSchemeTargetOrganism(primer_scheme_target_organism_name="SARS-CoV-2"),
        ],
        primer_scheme_vendor=[
            PrimerSchemeVendor(primer_scheme_vendor_name="idt"),
            PrimerSchemeVendor(
                primer_scheme_vendor_name="Ac;me Inc",
                primer_scheme_vendor_kit_name="Panel v1",
            ),
        ],
        primer_scheme_development_status=PrimerSchemeDevelopmentStatus.DRAFT,
        primer_scheme_creation_date=date(2024, 1, 15),
    )
    defaults.update(kwargs)
    return PrimerScheme(**defaults)


def test_flatten_unflatten_round_trip():
    """A scheme with multi-item groups of varying completeness round-trips exactly."""
    ps = _sample_scheme()
    row = flatten_scheme(ps)
    restored = unflatten_scheme(row)
    assert restored.model_dump() == ps.model_dump()


def test_flatten_unflatten_round_trip_with_embedded_newline_in_free_text():
    """A repeatable scalar field with an embedded newline round-trips without data loss."""
    ps = _sample_scheme(
        primer_scheme_details=["Paper A", "Multi-line\nnote", "Paper C"],
    )
    row = flatten_scheme(ps)
    restored = unflatten_scheme(row)
    assert restored.primer_scheme_details == ["Paper A", "Multi-line\nnote", "Paper C"]


def test_flatten_unflatten_round_trip_single_item_groups():
    """A scheme with exactly one (partially empty) item per group round-trips.

    This is the realistic common case that the csv.writer single-empty-field
    quirk broke: one target organism with a name but no ncbi_taxon_id.
    """
    ps = _sample_scheme(
        primer_scheme_contributor=[
            PrimerSchemeContributor(primer_scheme_contributor_name="Alice")
        ],
        primer_scheme_target_organism=[
            PrimerSchemeTargetOrganism(primer_scheme_target_organism_name="SARS-CoV-2")
        ],
        primer_scheme_vendor=[PrimerSchemeVendor(primer_scheme_vendor_name="idt")],
    )
    row = flatten_scheme(ps)
    restored = unflatten_scheme(row)
    assert restored.model_dump() == ps.model_dump()
    assert len(restored.primer_scheme_target_organism) == 1
    assert restored.primer_scheme_target_organism[
        0
    ].primer_scheme_target_organism_name == ("SARS-CoV-2")
    assert (
        restored.primer_scheme_target_organism[
            0
        ].primer_scheme_target_organism_ncbi_taxon_id
        is None
    )


def test_flatten_unflatten_round_trip_no_vendors():
    """An optional group with zero items round-trips to an empty list, not a phantom item."""
    ps = _sample_scheme(primer_scheme_vendor=[])
    row = flatten_scheme(ps)
    restored = unflatten_scheme(row)
    assert restored.primer_scheme_vendor == []


def test_csv_fieldnames_cover_expanded_columns():
    """CSV_FIELDNAMES expands group fields into their per-attribute columns."""
    assert "primer_scheme_contributor_name" in CSV_FIELDNAMES
    assert "primer_scheme_contributor_email" in CSV_FIELDNAMES
    assert "primer_scheme_contributor" not in CSV_FIELDNAMES
    assert len(CSV_FIELDNAMES) == len(set(CSV_FIELDNAMES))


def test_unflatten_raises_on_misaligned_columns():
    """A genuinely mismatched pair of packed columns raises a clear error."""
    ps = _sample_scheme()
    row = flatten_scheme(ps)
    # Corrupt one column to have a different item count than its siblings.
    row["primer_scheme_contributor_name"] = "Alice;Bob;Carol"
    with pytest.raises(ValueError, match="Misaligned columns"):
        unflatten_scheme(row)


def test_unflatten_raises_when_identifier_mismatches_edited_fields():
    """Editing amplicon_size (etc.) without updating the identifier is a hard error.

    A stale primer_scheme_identifier is never silently self-healed or merely
    warned about here — unflatten_scheme must raise, same as validate_identifier
    does for a hand-edited info.json on disk.
    """
    ps = _sample_scheme()
    row = flatten_scheme(ps)
    row["amplicon_size"] = "999"  # edited without updating primer_scheme_identifier
    with pytest.raises(ValueError, match="primer_scheme_identifier mismatch"):
        unflatten_scheme(row)


def test_unflatten_accepts_edited_fields_when_identifier_also_updated():
    """Editing amplicon_size *and* the identifier together is accepted."""
    ps = _sample_scheme()
    row = flatten_scheme(ps)
    row["amplicon_size"] = "999"
    row["primer_scheme_identifier"] = "test-scheme/999/v1.0.0"
    restored = unflatten_scheme(row)
    assert restored.amplicon_size == 999
    assert restored.primer_scheme_identifier == "test-scheme/999/v1.0.0"
