import pytest

from primaschema.cli import (
    parse_algorithm,
    parse_contributor_single,
    parse_contributors_pydantic,
    parse_target_organism_single,
    parse_target_organisms_pydantic,
    parse_vendor_single,
)
from primaschema.schema.info import Algorithm, Contributor, TargetOrganism, Vendor

# --- Contributor Tests ---


def test_parse_contributor_single_object():
    """parse_contributor_single passes through an already-constructed Contributor unchanged."""
    c = Contributor(name="John Doe", email="john@example.com")
    assert parse_contributor_single(c) == c


def test_parse_contributor_single_dict():
    """parse_contributor_single constructs a Contributor from a plain dict."""
    data = {"name": "Jane Doe", "email": "jane@example.com"}
    c = parse_contributor_single(data)
    assert c.name == "Jane Doe"
    assert c.email == "jane@example.com"


def test_parse_contributor_single_json_string():
    """parse_contributor_single parses a JSON-encoded contributor string."""
    json_str = '{"name": "Json User", "email": "json@example.com"}'
    c = parse_contributor_single(json_str)
    assert c.name == "Json User"
    assert c.email == "json@example.com"


def test_parse_contributor_single_kv_string():
    """parse_contributor_single parses a key=value comma-separated contributor string."""
    kv_str = "name=KV User,email=kv@example.com"
    c = parse_contributor_single(kv_str)
    assert c.name == "KV User"
    assert c.email == "kv@example.com"


def test_parse_contributor_single_simple_string():
    """parse_contributor_single treats a bare string as the contributor name."""
    name_str = "Simple User"
    c = parse_contributor_single(name_str)
    assert c.name == "Simple User"
    assert c.email is None


def test_parse_contributor_single_invalid():
    """parse_contributor_single raises ValueError for unsupported input types (e.g. int)."""
    with pytest.raises(ValueError):
        parse_contributor_single(123)


def test_parse_contributors_pydantic_list():
    """parse_contributors_pydantic parses a list of mixed-format contributor strings."""
    input_list = ["User One", "name=User Two"]
    result = parse_contributors_pydantic(input_list)
    assert len(result) == 2
    assert result[0].name == "User One"
    assert result[1].name == "User Two"


# --- Vendor Tests ---


def test_parse_vendor_single_object():
    """parse_vendor_single passes through an already-constructed Vendor unchanged."""
    v = Vendor(organisation_name="Acme Corp")
    assert parse_vendor_single(v) == v


def test_parse_vendor_single_dict():
    """parse_vendor_single constructs a Vendor from a plain dict."""
    data = {"organisation_name": "Beta Inc"}
    v = parse_vendor_single(data)
    assert v.organisation_name == "Beta Inc"


def test_parse_vendor_single_json_string():
    """parse_vendor_single parses a JSON-encoded vendor string."""
    json_str = '{"organisation_name": "Gamma Ltd"}'
    v = parse_vendor_single(json_str)
    assert v.organisation_name == "Gamma Ltd"


def test_parse_vendor_single_kv_string():
    """parse_vendor_single parses a key=value comma-separated vendor string."""
    kv_str = "organisation_name=Delta Co,kit_name=SuperKit"
    v = parse_vendor_single(kv_str)
    assert v.organisation_name == "Delta Co"
    assert v.kit_name == "SuperKit"


def test_parse_vendor_single_simple_string():
    """parse_vendor_single treats a bare string as the organisation_name."""
    name_str = "Epsilon LLC"
    v = parse_vendor_single(name_str)
    assert v.organisation_name == "Epsilon LLC"


def test_parse_vendor_single_invalid():
    """parse_vendor_single raises ValueError for unsupported input types (e.g. int)."""
    with pytest.raises(ValueError):
        parse_vendor_single(123)


# --- Algorithm Tests ---


def test_parse_algorithm_none():
    """parse_algorithm returns None when given None (optional field with no value)."""
    assert parse_algorithm(None) is None


def test_parse_algorithm_object():
    """parse_algorithm passes through an already-constructed Algorithm unchanged."""
    a = Algorithm(name="algo", version="1.0")
    assert parse_algorithm(a) == a


def test_parse_algorithm_dict():
    """parse_algorithm constructs an Algorithm from a plain dict."""
    data = {"name": "dict_algo", "version": "2.0"}
    a = parse_algorithm(data)
    assert a.name == "dict_algo"
    assert a.version == "2.0"


def test_parse_algorithm_string_with_version():
    """parse_algorithm splits 'name:version' strings into separate name and version fields."""
    s = "tool:1.2.3"
    a = parse_algorithm(s)
    assert a.name == "tool"
    assert a.version == "1.2.3"


def test_parse_algorithm_string_name_only():
    """parse_algorithm treats a string without ':' as the algorithm name only."""
    s = "tool_only"
    a = parse_algorithm(s)
    assert a.name == "tool_only"
    assert a.version is None


def test_parse_algorithm_invalid():
    """parse_algorithm raises ValueError for unsupported input types (e.g. int)."""
    with pytest.raises(ValueError):
        parse_algorithm(123)


# --- TargetOrganism Tests ---


def test_parse_target_organism_single_object():
    """parse_target_organism_single passes through an already-constructed TargetOrganism unchanged."""
    to = TargetOrganism(common_name="Virus X")
    assert parse_target_organism_single(to) == to


def test_parse_target_organism_single_dict():
    """parse_target_organism_single constructs a TargetOrganism from a plain dict."""
    data = {"common_name": "Virus Y", "ncbi_tax_id": "12345"}
    to = parse_target_organism_single(data)
    assert to.common_name == "Virus Y"
    assert to.ncbi_tax_id == "12345"


def test_parse_target_organism_single_json_string():
    """parse_target_organism_single parses a JSON-encoded target organism string."""
    json_str = '{"common_name": "Virus Z", "ncbi_tax_id": "67890"}'
    to = parse_target_organism_single(json_str)
    assert to.common_name == "Virus Z"
    assert to.ncbi_tax_id == "67890"


def test_parse_target_organism_single_kv_string():
    """parse_target_organism_single parses a key=value comma-separated target organism string."""
    kv_str = "common_name=Virus A,ncbi_tax_id=11111"
    to = parse_target_organism_single(kv_str)
    assert to.common_name == "Virus A"
    assert to.ncbi_tax_id == "11111"


def test_parse_target_organism_single_tax_id_string():
    """parse_target_organism_single treats an all-digit string as an NCBI tax ID."""
    tax_id = "99999"
    to = parse_target_organism_single(tax_id)
    assert to.ncbi_tax_id == "99999"
    assert to.common_name is None


def test_parse_target_organism_single_common_name_string():
    """parse_target_organism_single treats a non-numeric string as a common name."""
    name = "Virus B"
    to = parse_target_organism_single(name)
    assert to.common_name == "Virus B"
    assert to.ncbi_tax_id is None


def test_parse_target_organism_single_invalid():
    """parse_target_organism_single raises ValueError for unsupported input types (e.g. int)."""
    with pytest.raises(ValueError):
        parse_target_organism_single(123)


def test_parse_target_organisms_pydantic_list():
    """parse_target_organisms_pydantic parses a list of mixed-format target organism strings."""
    input_list = ["Virus C", "12345"]
    result = parse_target_organisms_pydantic(input_list)
    assert len(result) == 2
    assert result[0].common_name == "Virus C"
    assert result[1].ncbi_tax_id == "12345"


def test_parse_target_organisms_pydantic_single_string():
    """parse_target_organisms_pydantic wraps a single string in a list."""
    input_str = "Virus D"
    result = parse_target_organisms_pydantic(input_str)
    assert len(result) == 1
    assert result[0].common_name == "Virus D"
