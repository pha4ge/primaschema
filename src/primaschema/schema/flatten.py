import csv
import io

from primaschema.schema.info import (
    PrimerSchemeChecksums,
    PrimerSchemeContributor,
    PrimerSchemeGenerator,
    PrimerSchemeTargetOrganism,
    PrimerSchemeVendor,
)
from primaschema.schema.primer_scheme import (
    PrimerScheme,
    check_primer_scheme_identifier,
)

_DELIMITER = ";"

# Repeatable groups: each attribute of the item class becomes its own
# column; multiple items are packed into that column, one value per item,
# positionally aligned across all of the group's columns.
_REPEATABLE_CLASS_GROUPS = {
    "primer_scheme_contributor": PrimerSchemeContributor,
    "primer_scheme_target_organism": PrimerSchemeTargetOrganism,
    "primer_scheme_vendor": PrimerSchemeVendor,
}

# Singular groups: at most one object, so its attributes promote directly
# to flat columns with no packing.
_SINGULAR_CLASS_GROUPS = {
    "primer_scheme_generator": PrimerSchemeGenerator,
    "primer_scheme_checksums": PrimerSchemeChecksums,
}

# Repeatable scalar fields: a plain list of strings/enum values, packed into
# one column with no cross-column alignment needed.
_REPEATABLE_SCALAR_FIELDS = {
    "primer_scheme_identifier_alias",
    "citation",
    "primer_scheme_details",
}


def _pack(values: list) -> str:
    """Join values with ';', using '' for missing/None so position stays aligned.

    Uses the default (non-empty) lineterminator and strips it afterward,
    rather than passing lineterminator="" to csv.writer. csv.writer only
    quotes an embedded '\\n'/'\\r' in a value when that character appears in
    the dialect's lineterminator - with lineterminator="" that check can
    never trigger, so an embedded newline/CR would be written unquoted and
    silently truncate the value on unpack. Confirmed: with the default
    terminator, csv.writer correctly quotes an embedded '\\n'.
    """
    cells = ["" if v is None else str(v) for v in values]
    if len(cells) == 1:
        # csv.writer defensively quotes a lone empty field as '""' to
        # disambiguate a one-empty-field row from a fully empty row. Work
        # around it by writing a dummy second field and dropping it, so a
        # single empty value packs the same way an empty value in any other
        # position would.
        buf = io.StringIO()
        csv.writer(buf, delimiter=_DELIMITER).writerow([cells[0], ""])
        packed = buf.getvalue().rstrip("\r\n")
        return packed[: -len(_DELIMITER)]
    buf = io.StringIO()
    csv.writer(buf, delimiter=_DELIMITER).writerow(cells)
    return buf.getvalue().rstrip("\r\n")


def _unpack(cell: str) -> list[str]:
    """Split a packed cell back into per-position strings ('' for missing)."""
    if cell == "":
        return []
    return next(csv.reader(io.StringIO(cell), delimiter=_DELIMITER))


def _flatten_group_columns(item_cls: type) -> list[str]:
    return list(item_cls.model_fields)


CSV_FIELDNAMES: list[str] = []
for _field_name in PrimerScheme.model_fields:
    if _field_name in _REPEATABLE_CLASS_GROUPS:
        CSV_FIELDNAMES.extend(
            _flatten_group_columns(_REPEATABLE_CLASS_GROUPS[_field_name])
        )
    elif _field_name in _SINGULAR_CLASS_GROUPS:
        CSV_FIELDNAMES.extend(
            _flatten_group_columns(_SINGULAR_CLASS_GROUPS[_field_name])
        )
    else:
        CSV_FIELDNAMES.append(_field_name)
del _field_name


def flatten_scheme(ps: PrimerScheme) -> dict[str, str]:
    """Flatten a PrimerScheme into a single flat row of str -> str.

    Repeatable groups (PrimerSchemeContributor/PrimerSchemeTargetOrganism/PrimerSchemeVendor lists, and plain
    repeatable scalar fields) are packed into semicolon-delimited cells,
    positionally aligned across a group's columns. Singular groups
    (PrimerSchemeGenerator/PrimerSchemeChecksums) and plain scalar fields map straight to one
    column each.
    """
    data = ps.model_dump(mode="json")
    row: dict[str, str] = {}
    for field_name in PrimerScheme.model_fields:
        if field_name in _REPEATABLE_CLASS_GROUPS:
            item_cls = _REPEATABLE_CLASS_GROUPS[field_name]
            items = data.get(field_name) or []
            for attr in item_cls.model_fields:
                row[attr] = _pack([item.get(attr) for item in items])
        elif field_name in _SINGULAR_CLASS_GROUPS:
            item_cls = _SINGULAR_CLASS_GROUPS[field_name]
            obj = data.get(field_name) or {}
            for attr in item_cls.model_fields:
                value = obj.get(attr)
                row[attr] = "" if value is None else str(value)
        elif field_name in _REPEATABLE_SCALAR_FIELDS:
            values = data.get(field_name) or []
            row[field_name] = _pack(values)
        else:
            value = data.get(field_name)
            row[field_name] = "" if value is None else str(value)
    return row


def unflatten_scheme(row: dict[str, str]) -> PrimerScheme:
    """Reconstruct a PrimerScheme from a flat row produced by flatten_scheme."""
    data: dict = {}
    for field_name in PrimerScheme.model_fields:
        if field_name in _REPEATABLE_CLASS_GROUPS:
            item_cls = _REPEATABLE_CLASS_GROUPS[field_name]
            attrs = list(item_cls.model_fields)
            unpacked = {attr: _unpack(row.get(attr, "")) for attr in attrs}
            # A column that is empty for every item in the group unpacks to
            # length 0, indistinguishable from "the group has zero items" —
            # csv.writer can't tell "no items" from "one item, this
            # attribute unset" apart in a single column. Only lengths that
            # are actually populated (>0) need to agree with each other;
            # any all-empty column is padded out to match once the true
            # count is known.
            lengths_nonzero = {len(v) for v in unpacked.values() if len(v) > 0}
            if len(lengths_nonzero) > 1:
                raise ValueError(
                    f"Misaligned columns for {field_name!r}: "
                    f"{ {attr: len(values) for attr, values in unpacked.items()} }"
                )
            n = lengths_nonzero.pop() if lengths_nonzero else 0
            for attr in attrs:
                if len(unpacked[attr]) < n:
                    unpacked[attr] = unpacked[attr] + [""] * (n - len(unpacked[attr]))
            data[field_name] = [
                {attr: (unpacked[attr][i] or None) for attr in attrs} for i in range(n)
            ]
        elif field_name in _SINGULAR_CLASS_GROUPS:
            item_cls = _SINGULAR_CLASS_GROUPS[field_name]
            attrs = list(item_cls.model_fields)
            obj = {attr: (row.get(attr) or None) for attr in attrs}
            data[field_name] = obj if any(obj.values()) else None
        elif field_name in _REPEATABLE_SCALAR_FIELDS:
            data[field_name] = [v for v in _unpack(row.get(field_name, "")) if v != ""]
        else:
            value = row.get(field_name, "")
            data[field_name] = value if value != "" else None

    check_primer_scheme_identifier(
        data.get("primer_scheme_identifier"),
        data.get("primer_scheme_name"),
        data.get("amplicon_size"),
        data.get("primer_scheme_version"),
        source="CSV row",
    )
    return PrimerScheme.model_validate(data)
