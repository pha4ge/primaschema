import pytest

from primaschema.create import (
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
    c = Contributor(name="John Doe", email="john@example.com")
    assert parse_contributor_single(c) == c


def test_parse_contributor_single_dict():
    data = {"name": "Jane Doe", "email": "jane@example.com"}
    c = parse_contributor_single(data)
    assert c.name == "Jane Doe"
    assert c.email == "jane@example.com"


def test_parse_contributor_single_json_string():
    json_str = '{"name": "Json User", "email": "json@example.com"}'
    c = parse_contributor_single(json_str)
    assert c.name == "Json User"
    assert c.email == "json@example.com"


def test_parse_contributor_single_kv_string():
    kv_str = "name=KV User,email=kv@example.com"
    c = parse_contributor_single(kv_str)
    assert c.name == "KV User"
    assert c.email == "kv@example.com"


def test_parse_contributor_single_simple_string():
    name_str = "Simple User"
    c = parse_contributor_single(name_str)
    assert c.name == "Simple User"
    assert c.email is None


def test_parse_contributor_single_invalid():
    with pytest.raises(ValueError):
        parse_contributor_single(123)


def test_parse_contributors_pydantic_list():
    input_list = ["User One", "name=User Two"]
    result = parse_contributors_pydantic(input_list)
    assert len(result) == 2
    assert result[0].name == "User One"
    assert result[1].name == "User Two"


# --- Vendor Tests ---


def test_parse_vendor_single_object():
    v = Vendor(organisation_name="Acme Corp")
    assert parse_vendor_single(v) == v


def test_parse_vendor_single_dict():
    data = {"organisation_name": "Beta Inc"}
    v = parse_vendor_single(data)
    assert v.organisation_name == "Beta Inc"


def test_parse_vendor_single_json_string():
    json_str = '{"organisation_name": "Gamma Ltd"}'
    v = parse_vendor_single(json_str)
    assert v.organisation_name == "Gamma Ltd"


def test_parse_vendor_single_kv_string():
    kv_str = "organisation_name=Delta Co,kit_name=SuperKit"
    v = parse_vendor_single(kv_str)
    assert v.organisation_name == "Delta Co"
    assert v.kit_name == "SuperKit"


def test_parse_vendor_single_simple_string():
    name_str = "Epsilon LLC"
    v = parse_vendor_single(name_str)
    assert v.organisation_name == "Epsilon LLC"


def test_parse_vendor_single_invalid():
    with pytest.raises(ValueError):
        parse_vendor_single(123)


# --- Algorithm Tests ---


def test_parse_algorithm_none():
    assert parse_algorithm(None) is None


def test_parse_algorithm_object():
    a = Algorithm(name="algo", version="1.0")
    assert parse_algorithm(a) == a


def test_parse_algorithm_dict():
    data = {"name": "dict_algo", "version": "2.0"}
    a = parse_algorithm(data)
    assert a.name == "dict_algo"
    assert a.version == "2.0"


def test_parse_algorithm_string_with_version():
    s = "tool:1.2.3"
    a = parse_algorithm(s)
    assert a.name == "tool"
    assert a.version == "1.2.3"


def test_parse_algorithm_string_name_only():
    s = "tool_only"
    a = parse_algorithm(s)
    assert a.name == "tool_only"
    assert a.version is None


def test_parse_algorithm_invalid():
    with pytest.raises(ValueError):
        parse_algorithm(123)


# --- TargetOrganism Tests ---


def test_parse_target_organism_single_object():
    to = TargetOrganism(common_name="Virus X")
    assert parse_target_organism_single(to) == to


def test_parse_target_organism_single_dict():
    data = {"common_name": "Virus Y", "ncbi_tax_id": "12345"}
    to = parse_target_organism_single(data)
    assert to.common_name == "Virus Y"
    assert to.ncbi_tax_id == "12345"


def test_parse_target_organism_single_json_string():
    json_str = '{"common_name": "Virus Z", "ncbi_tax_id": "67890"}'
    to = parse_target_organism_single(json_str)
    assert to.common_name == "Virus Z"
    assert to.ncbi_tax_id == "67890"


def test_parse_target_organism_single_kv_string():
    kv_str = "common_name=Virus A,ncbi_tax_id=11111"
    to = parse_target_organism_single(kv_str)
    assert to.common_name == "Virus A"
    assert to.ncbi_tax_id == "11111"


def test_parse_target_organism_single_tax_id_string():
    tax_id = "99999"
    to = parse_target_organism_single(tax_id)
    assert to.ncbi_tax_id == "99999"
    assert to.common_name is None


def test_parse_target_organism_single_common_name_string():
    name = "Virus B"
    to = parse_target_organism_single(name)
    assert to.common_name == "Virus B"
    assert to.ncbi_tax_id is None


def test_parse_target_organism_single_invalid():
    with pytest.raises(ValueError):
        parse_target_organism_single(123)


def test_parse_target_organisms_pydantic_list():
    input_list = ["Virus C", "12345"]
    result = parse_target_organisms_pydantic(input_list)
    assert len(result) == 2
    assert result[0].common_name == "Virus C"
    assert result[1].ncbi_tax_id == "12345"


def test_parse_target_organisms_pydantic_single_string():
    input_str = "Virus D"
    result = parse_target_organisms_pydantic(input_str)
    assert len(result) == 1
    assert result[0].common_name == "Virus D"
