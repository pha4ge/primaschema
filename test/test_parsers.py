import pytest

from primaschema.cli import (
    parse_algorithm,
    parse_contributor_single,
    parse_contributors_pydantic,
    parse_target_organism_single,
    parse_target_organisms_pydantic,
    parse_vendor_single,
)
from primaschema.schema.info import (
    PrimerSchemeGenerator,
    PrimerSchemeContributor,
    PrimerSchemeTargetOrganism,
    PrimerSchemeVendor,
)

# --- PrimerSchemeContributor Tests ---


def test_parse_contributor_single_object():
    """parse_contributor_single passes through an already-constructed PrimerSchemeContributor unchanged."""
    c = PrimerSchemeContributor(
        primer_scheme_contributor_name="John Doe",
        primer_scheme_contributor_email="john@example.com",
    )
    assert parse_contributor_single(c) == c


def test_parse_contributor_single_dict():
    """parse_contributor_single constructs a PrimerSchemeContributor from a plain dict."""
    data = {
        "primer_scheme_contributor_name": "Jane Doe",
        "primer_scheme_contributor_email": "jane@example.com",
    }
    c = parse_contributor_single(data)
    assert c.primer_scheme_contributor_name == "Jane Doe"
    assert c.primer_scheme_contributor_email == "jane@example.com"


def test_parse_contributor_single_json_string():
    """parse_contributor_single parses a JSON-encoded contributor string."""
    json_str = '{"primer_scheme_contributor_name": "Json User", "primer_scheme_contributor_email": "json@example.com"}'
    c = parse_contributor_single(json_str)
    assert c.primer_scheme_contributor_name == "Json User"
    assert c.primer_scheme_contributor_email == "json@example.com"


def test_parse_contributor_single_kv_string():
    """parse_contributor_single parses a key=value comma-separated contributor string."""
    kv_str = "primer_scheme_contributor_name=KV User,primer_scheme_contributor_email=kv@example.com"
    c = parse_contributor_single(kv_str)
    assert c.primer_scheme_contributor_name == "KV User"
    assert c.primer_scheme_contributor_email == "kv@example.com"


def test_parse_contributor_single_kv_string_short_keys():
    """parse_contributor_single accepts short keys (name=, email=, orcid=) as well as the full field names."""
    c = parse_contributor_single(
        "name=KV User,email=kv@example.com,orcid=0000-0001-2345-6789"
    )
    assert c.primer_scheme_contributor_name == "KV User"
    assert c.primer_scheme_contributor_email == "kv@example.com"
    assert c.primer_scheme_contributor_orcid == "0000-0001-2345-6789"


def test_parse_contributor_single_short_and_long_keys_equivalent():
    """Short and fully-qualified keys produce the same PrimerSchemeContributor."""
    short = parse_contributor_single(
        "name=KV User,email=kv@example.com,orcid=0000-0001-2345-6789"
    )
    full = parse_contributor_single(
        "primer_scheme_contributor_name=KV User,primer_scheme_contributor_email=kv@example.com,primer_scheme_contributor_orcid=0000-0001-2345-6789"
    )
    assert short == full


def test_parse_contributor_single_kv_string_unknown_key_raises():
    """An unrecognised key=value key still fails loudly instead of being silently dropped."""
    with pytest.raises(ValueError):
        parse_contributor_single("bogus_key=Someone")


def test_parse_contributor_single_simple_string():
    """parse_contributor_single treats a bare string as the contributor name."""
    name_str = "Simple User"
    c = parse_contributor_single(name_str)
    assert c.primer_scheme_contributor_name == "Simple User"
    assert c.primer_scheme_contributor_email is None


def test_parse_contributor_single_invalid():
    """parse_contributor_single raises ValueError for unsupported input types (e.g. int)."""
    with pytest.raises(ValueError):
        parse_contributor_single(123)


def test_parse_contributors_pydantic_list():
    """parse_contributors_pydantic parses a list of mixed-format contributor strings."""
    input_list = ["User One", "primer_scheme_contributor_name=User Two"]
    result = parse_contributors_pydantic(input_list)
    assert len(result) == 2
    assert result[0].primer_scheme_contributor_name == "User One"
    assert result[1].primer_scheme_contributor_name == "User Two"


# --- PrimerSchemeVendor Tests ---


def test_parse_vendor_single_object():
    """parse_vendor_single passes through an already-constructed PrimerSchemeVendor unchanged."""
    v = PrimerSchemeVendor(primer_scheme_vendor_name="Acme Corp")
    assert parse_vendor_single(v) == v


def test_parse_vendor_single_dict():
    """parse_vendor_single constructs a PrimerSchemeVendor from a plain dict."""
    data = {"primer_scheme_vendor_name": "Beta Inc"}
    v = parse_vendor_single(data)
    assert v.primer_scheme_vendor_name == "Beta Inc"


def test_parse_vendor_single_json_string():
    """parse_vendor_single parses a JSON-encoded vendor string."""
    json_str = '{"primer_scheme_vendor_name": "Gamma Ltd"}'
    v = parse_vendor_single(json_str)
    assert v.primer_scheme_vendor_name == "Gamma Ltd"


def test_parse_vendor_single_kv_string():
    """parse_vendor_single parses a key=value comma-separated vendor string."""
    kv_str = "primer_scheme_vendor_name=Delta Co,primer_scheme_vendor_kit_name=SuperKit"
    v = parse_vendor_single(kv_str)
    assert v.primer_scheme_vendor_name == "Delta Co"
    assert v.primer_scheme_vendor_kit_name == "SuperKit"


def test_parse_vendor_single_kv_string_short_keys():
    """parse_vendor_single accepts short keys (name=, kit_name=, url=) as well as the full field names."""
    v = parse_vendor_single("name=Delta Co,kit_name=SuperKit,url=https://example.com")
    assert v.primer_scheme_vendor_name == "Delta Co"
    assert v.primer_scheme_vendor_kit_name == "SuperKit"
    assert v.primer_scheme_vendor_url == "https://example.com"


def test_parse_vendor_single_short_and_long_keys_equivalent():
    """Short and fully-qualified keys produce the same PrimerSchemeVendor."""
    short = parse_vendor_single(
        "name=Delta Co,kit_name=SuperKit,url=https://example.com"
    )
    full = parse_vendor_single(
        "primer_scheme_vendor_name=Delta Co,primer_scheme_vendor_kit_name=SuperKit,primer_scheme_vendor_url=https://example.com"
    )
    assert short == full


def test_parse_vendor_single_kv_string_unknown_key_raises():
    """An unrecognised key=value key still fails loudly instead of being silently dropped."""
    with pytest.raises(ValueError):
        parse_vendor_single("bogus_key=Someone")


def test_parse_vendor_single_simple_string():
    """parse_vendor_single treats a bare string as the organisation_name."""
    name_str = "Epsilon LLC"
    v = parse_vendor_single(name_str)
    assert v.primer_scheme_vendor_name == "Epsilon LLC"


def test_parse_vendor_single_invalid():
    """parse_vendor_single raises ValueError for unsupported input types (e.g. int)."""
    with pytest.raises(ValueError):
        parse_vendor_single(123)


# --- PrimerSchemeGenerator Tests ---


def test_parse_algorithm_none():
    """parse_algorithm returns None when given None (optional field with no value)."""
    assert parse_algorithm(None) is None


def test_parse_algorithm_object():
    """parse_algorithm passes through an already-constructed PrimerSchemeGenerator unchanged."""
    a = PrimerSchemeGenerator(
        primer_scheme_generator_name="algo", primer_scheme_generator_version="1.0"
    )
    assert parse_algorithm(a) == a


def test_parse_algorithm_dict():
    """parse_algorithm constructs an PrimerSchemeGenerator from a plain dict."""
    data = {
        "primer_scheme_generator_name": "dict_algo",
        "primer_scheme_generator_version": "2.0",
    }
    a = parse_algorithm(data)
    assert a.primer_scheme_generator_name == "dict_algo"
    assert a.primer_scheme_generator_version == "2.0"


def test_parse_algorithm_string_with_version():
    """parse_algorithm splits 'name:version' strings into separate name and version fields."""
    s = "tool:1.2.3"
    a = parse_algorithm(s)
    assert a.primer_scheme_generator_name == "tool"
    assert a.primer_scheme_generator_version == "1.2.3"


def test_parse_algorithm_string_name_only():
    """parse_algorithm treats a string without ':' as the algorithm name only."""
    s = "tool_only"
    a = parse_algorithm(s)
    assert a.primer_scheme_generator_name == "tool_only"
    assert a.primer_scheme_generator_version is None


def test_parse_algorithm_invalid():
    """parse_algorithm raises ValueError for unsupported input types (e.g. int)."""
    with pytest.raises(ValueError):
        parse_algorithm(123)


# --- PrimerSchemeTargetOrganism Tests ---


def test_parse_target_organism_single_object():
    """parse_target_organism_single passes through an already-constructed PrimerSchemeTargetOrganism unchanged."""
    to = PrimerSchemeTargetOrganism(primer_scheme_target_organism_name="Virus X")
    assert parse_target_organism_single(to) == to


def test_parse_target_organism_single_dict():
    """parse_target_organism_single constructs a PrimerSchemeTargetOrganism from a plain dict."""
    data = {
        "primer_scheme_target_organism_name": "Virus Y",
        "primer_scheme_target_organism_ncbi_taxon_id": "12345",
    }
    to = parse_target_organism_single(data)
    assert to.primer_scheme_target_organism_name == "Virus Y"
    assert to.primer_scheme_target_organism_ncbi_taxon_id == "12345"


def test_parse_target_organism_single_json_string():
    """parse_target_organism_single parses a JSON-encoded target organism string."""
    json_str = '{"primer_scheme_target_organism_name": "Virus Z", "primer_scheme_target_organism_ncbi_taxon_id": "67890"}'
    to = parse_target_organism_single(json_str)
    assert to.primer_scheme_target_organism_name == "Virus Z"
    assert to.primer_scheme_target_organism_ncbi_taxon_id == "67890"


def test_parse_target_organism_single_kv_string():
    """parse_target_organism_single parses a key=value comma-separated target organism string."""
    kv_str = "primer_scheme_target_organism_name=Virus A,primer_scheme_target_organism_ncbi_taxon_id=11111"
    to = parse_target_organism_single(kv_str)
    assert to.primer_scheme_target_organism_name == "Virus A"
    assert to.primer_scheme_target_organism_ncbi_taxon_id == "11111"


def test_parse_target_organism_single_kv_string_short_keys():
    """parse_target_organism_single accepts short keys (name=, ncbi_taxon_id=) as well as the full field names."""
    to = parse_target_organism_single("name=Virus A,ncbi_taxon_id=11111")
    assert to.primer_scheme_target_organism_name == "Virus A"
    assert to.primer_scheme_target_organism_ncbi_taxon_id == "11111"


def test_parse_target_organism_single_short_and_long_keys_equivalent():
    """Short and fully-qualified keys produce the same PrimerSchemeTargetOrganism."""
    short = parse_target_organism_single("name=Virus A,ncbi_taxon_id=11111")
    full = parse_target_organism_single(
        "primer_scheme_target_organism_name=Virus A,primer_scheme_target_organism_ncbi_taxon_id=11111"
    )
    assert short == full


def test_parse_target_organism_single_kv_string_unknown_key_raises():
    """An unrecognised key=value key still fails loudly instead of being silently dropped."""
    with pytest.raises(ValueError):
        parse_target_organism_single("bogus_key=Something")


def test_parse_target_organism_single_tax_id_string():
    """parse_target_organism_single treats an all-digit string as an NCBI tax ID."""
    tax_id = "99999"
    to = parse_target_organism_single(tax_id)
    assert to.primer_scheme_target_organism_ncbi_taxon_id == "99999"
    assert to.primer_scheme_target_organism_name is None


def test_parse_target_organism_single_common_name_string():
    """parse_target_organism_single treats a non-numeric string as a common name."""
    name = "Virus B"
    to = parse_target_organism_single(name)
    assert to.primer_scheme_target_organism_name == "Virus B"
    assert to.primer_scheme_target_organism_ncbi_taxon_id is None


def test_parse_target_organism_single_invalid():
    """parse_target_organism_single raises ValueError for unsupported input types (e.g. int)."""
    with pytest.raises(ValueError):
        parse_target_organism_single(123)


def test_parse_target_organisms_pydantic_list():
    """parse_target_organisms_pydantic parses a list of mixed-format target organism strings."""
    input_list = ["Virus C", "12345"]
    result = parse_target_organisms_pydantic(input_list)
    assert len(result) == 2
    assert result[0].primer_scheme_target_organism_name == "Virus C"
    assert result[1].primer_scheme_target_organism_ncbi_taxon_id == "12345"


def test_parse_target_organisms_pydantic_single_string():
    """parse_target_organisms_pydantic wraps a single string in a list."""
    input_str = "Virus D"
    result = parse_target_organisms_pydantic(input_str)
    assert len(result) == 1
    assert result[0].primer_scheme_target_organism_name == "Virus D"
