"""Griffe extension for primaschema docs.

Handles LinkML-generated field validators that griffe-pydantic misses in static
analysis (typically on models with forward-reference fields).  For each function
whose name starts with ``pattern_`` and lacks the ``pydantic-validator`` label,
this extension:

1. Reads the function source and extracts the first ``re.compile(...)`` pattern.
2. Adds a one-line docstring: "Validates: must match ``<pattern>``".
3. Adds the ``pydantic-validator`` label so mkdocstrings renders it in the
   Validators section alongside validators detected by griffe-pydantic.
"""

import re

import griffe

_COMPILE_RE = re.compile(r're\.compile\(r["\']([^"\']+)["\']')


class InjectValidatorDocs(griffe.Extension):
    """Inject regex docstrings on LinkML pattern validators."""

    def on_class_members(self, *, cls: griffe.Class, **kwargs) -> None:
        for member in cls.members.values():
            if not isinstance(member, griffe.Function):
                continue
            if "pydantic-validator" in member.labels:
                continue  # already handled by griffe_pydantic
            if not member.name.startswith("pattern_"):
                continue
            if not member.filepath:
                continue
            try:
                lines = member.filepath.read_text().splitlines()
                source = "\n".join(lines[member.lineno - 1 : member.endlineno])
                match = _COMPILE_RE.search(source)
                if not match:
                    continue
                pattern = match.group(1)
                member.labels.add("pydantic-validator")
                if member.docstring is None:
                    member.docstring = griffe.Docstring(
                        f"Validates: must match `{pattern}`",
                        parent=member,
                    )
            except Exception:
                pass
