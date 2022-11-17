import subprocess
from pathlib import Path

import pytest

import primaschema.lib as lib

data_dir = "tests/data"


def run(cmd, cwd="./"):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=data_dir, shell=True, check=True, text=True, capture_output=True
    )


def test_cli_version():
    run("primaschema --version")


def test_hash_ref():
    assert (
        lib.hash_ref("tests/data/primer-schemes/eden/v1/reference.fasta")
        == "7d5621cd3b3e498d0c27fcca9d3d3c5168c7f3d3f9776f3005c7011bd90068ca"
    )


def test_cli_hash_ref():
    run_cmd = run("primaschema hash-ref primer-schemes/eden/v1/reference.fasta")
    assert (
        "7d5621cd3b3e498d0c27fcca9d3d3c5168c7f3d3f9776f3005c7011bd90068ca"
        in run_cmd.stdout
    )


def test_cli_hash_primer_bed():
    run_cmd = run("primaschema hash-primer-bed primer-schemes/artic/v4.1/primer.bed")
    assert (
        "6878eaff17dd3e5815223bba9e9f113075daa7b4c80ff95dfe735c24443969bc"
        in run_cmd.stdout
    )


def test_cli_scheme_bed():
    run_cmd = run(
        "primaschema hash-scheme-bed primer-schemes/artic/v4.1/scheme.bed primer-schemes/artic/v4.1/reference.fasta"
    )
    assert (
        "6878eaff17dd3e5815223bba9e9f113075daa7b4c80ff95dfe735c24443969bc"
        in run_cmd.stdout
    )


def test_artic_v41_scheme_hash_matches_primer_hash():
    scheme_bed_hash = lib.hash_scheme_bed(
        "tests/data/primer-schemes/artic/v4.1/scheme.bed",
        "tests/data/primer-schemes/artic/v4.1/reference.fasta",
    )
    primer_bed_hash = lib.hash_primer_bed(
        "tests/data/primer-schemes/artic/v4.1/primer.bed"
    )
    assert scheme_bed_hash == primer_bed_hash


def test_eden_v1_schema():
    lib.validate_yaml(f"{data_dir}/primer-schemes/eden/v1/scheme.yaml")


def test_artic_v41_schema():
    lib.validate_yaml(f"{data_dir}/primer-schemes/artic/v4.1/scheme.yaml")
