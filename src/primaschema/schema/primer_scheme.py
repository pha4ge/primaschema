from typing import Any

from pydantic import model_validator

from primaschema.schema.info import PrimerScheme as _GeneratedPrimerScheme

# info.py is generated straight from info.yml by gen-pydantic and gets
# overwritten by the pre-commit hook, so it can't hand-encode behaviour that
# gen-pydantic doesn't translate from LinkML into the Pydantic model. One
# such gap: primer_scheme_identifier's LinkML `structured_pattern` (its
# composition from name/amplicon_size/version) is dropped entirely during
# generation — the generated field is just a plain optional string with no
# computation or validation behind it. This subclass adds that behaviour
# back on top of the generated model instead.


def compute_primer_scheme_identifier(
    name: str, amplicon_size: int, version: str
) -> str:
    return f"{name}/{amplicon_size}/{version}"


def check_primer_scheme_identifier(
    provided: str | None, name: str, amplicon_size, version: str, source: str
) -> None:
    """Raise if a provided primer_scheme_identifier disagrees with the computed one.

    Used at every boundary where a primer_scheme_identifier might have been
    supplied by something other than PrimerScheme itself (a hand-edited
    info.json, an edited flattened CSV row) — those are hard errors, always,
    everywhere. This is distinct from PrimerScheme's own `mode="after"`
    validator, which silently self-heals on every construction/attribute
    assignment (required so legitimate multi-step internal updates, e.g.
    `rebuild --sync-metadata`, don't spuriously fail mid-update); this check
    is only ever invoked explicitly, at a genuine load/reconstruction
    boundary, never automatically on every assignment.
    """
    if not provided:
        return
    expected = compute_primer_scheme_identifier(name, amplicon_size, version)
    if provided != expected:
        raise ValueError(
            f"primer_scheme_identifier mismatch in {source}: "
            f"found {provided!r}, expected {expected!r}. "
            "Update primer_scheme_identifier to match, or unset it and it "
            "will be computed automatically."
        )


# primer_scheme_application and primer_scheme_scope became multivalued in
# info.yml. Pydantic won't coerce a bare string into a list, so an info.json
# written against the single-valued schema fails to load outright. Wrap a
# scalar into a one-element list on the way in so older files keep working.
_SCALAR_TO_LIST_FIELDS = (
    "primer_scheme_application",
    "primer_scheme_scope",
)


class PrimerScheme(_GeneratedPrimerScheme):
    @model_validator(mode="before")
    @classmethod
    def _widen_scalars_to_lists(cls, data: Any) -> Any:
        if isinstance(data, dict):
            scalars = {
                field: [data[field]]
                for field in _SCALAR_TO_LIST_FIELDS
                if isinstance(data.get(field), str)
            }
            if scalars:
                # Don't mutate the caller's dict.
                return {**data, **scalars}
        return data

    @model_validator(mode="after")
    def _sync_primer_scheme_identifier(self):
        expected = compute_primer_scheme_identifier(
            self.primer_scheme_name, self.amplicon_size, self.primer_scheme_version
        )
        if self.primer_scheme_identifier != expected:
            self.primer_scheme_identifier = expected
        return self
