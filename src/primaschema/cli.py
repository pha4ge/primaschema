import sys
import logging

from pathlib import Path

import defopt

from . import lib, logger


def configure_logging(debug: bool):
    if debug:
        logger.setLevel(logging.DEBUG)
        for handler in logger.handlers:
            handler.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        for handler in logger.handlers:
            handler.setLevel(logging.INFO)


def hash_bed(bed_path: Path, debug: bool = False):
    """
    Generate a bed file checksum

    :arg bed_path: path of bed file
    :arg debug: show debug messages
    """
    configure_logging(debug)
    hex_digest = lib.hash_bed(bed_path)
    print("BED checksum:", file=sys.stderr)
    print(hex_digest)


def hash_ref(ref_path: Path, debug: bool = False):
    """
    Generate reference sequence checksum

    :arg ref_path: path of reference sequence
    :arg debug: emit debug messages
    """
    configure_logging(debug)
    hex_digest = lib.hash_ref(ref_path)
    print("Reference checksum:", file=sys.stderr)
    print(hex_digest)


def validate(
    scheme_dir: Path,
    full: bool = False,
    ignore_checksums: bool = False,
    recursive: bool = False,
    rebuild: bool = False,
    debug: bool = False,
):
    """
    Validate one or more primer scheme definitions comprising info.yml, primer.bed and reference.fasta

    :arg scheme_dir: path of scheme.bed file
    :arg full: perform meticulous validation using full model
    :arg ignore_checksums: ignore checksum mismatches
    :arg recursive: recursively find and validate primer scheme definitions
    :arg rebuild: forcibly rebuild the pydantic model from the linkml model
    :arg debug: show debug messages
    """
    configure_logging(debug)
    return lib.validate(
        scheme_dir,
        full=full,
        ignore_checksums=ignore_checksums,
        recursive=recursive,
        rebuild=rebuild,
    )


def build(
    scheme_dir: Path,
    out_dir: Path = Path("built"),
    nested: bool = False,
    plot: bool = False,
    recursive: bool = False,
    ignore_checksums: bool = False,
    debug: bool = False,
):
    """
    Build one or more primer scheme definitions comprising info.yml, primer.bed and reference.fasta

    :arg scheme_dir: path of input scheme directory
    :arg out_dir: path of directory in which to save scheme
    :arg nested: use nested output structure ({organism}/{scheme_name}/{amplicon_length}/{version})
    :arg plot: plot primers in SVG format
    :arg recursive: recursively find, validate and build primer scheme definitions
    :arg ignore_checksums: ignore checksum mismatches
    :arg debug: show debug messages
    """
    configure_logging(debug)
    lib.build(
        scheme_dir=scheme_dir,
        out_dir=out_dir,
        nested=nested,
        plot=plot,
        recursive=recursive,
        ignore_checksums=ignore_checksums,
    )


def build_manifest(root_dir: Path, out_dir: Path = Path()):
    """
    Build a complete manifest of schemes contained in the specified directory

    :arg root_dir: path in which to search for schemes
    :arg out_dir: path of directory in which to save manifest
    """
    lib.build_manifest(root_dir=root_dir, out_dir=out_dir)


def seven_to_six(bed_path: Path):
    """
    Convert a 7 column primer.bed file to a 6 column scheme.bed file by removing a column

    :arg bed_path: path of primer.bed file
    """
    bed_str = lib.convert_primer_bed_to_scheme_bed(bed_path=bed_path)
    print(bed_str)


def six_to_seven(bed_path: Path, fasta_path: Path):
    """
    Convert a 6 column scheme.bed file to a 7 column primer.bed file using reference backfill

    :arg bed_path: path of scheme.bed file
    :arg fasta_path: path of reference sequence
    """
    bed_str = lib.convert_scheme_bed_to_primer_bed(
        bed_path=bed_path, fasta_path=fasta_path
    )
    print(bed_str)


def diff(bed1_path: Path, bed2_path: Path, only_positions: bool = False):
    """
    Show the symmetric difference of records in two bed files

    :arg bed_path1: path of first bed file
    :arg bed_path2: path of second bed file
    :arg only_positions: Use only primer positions when computing differences
    """
    df = lib.diff(bed1_path, bed2_path, only_positions)
    if not df.empty:
        print(df.to_string(index=False))


def discordant_primers(scheme_dir: Path):
    """
    Show primer records with sequences not matching the reference sequence

    :arg scheme_dir: path of input scheme directory
    """
    df = lib.discordant_primers(scheme_dir=scheme_dir)
    if not df.empty:
        print(df.to_string(index=False))


def amplicon_intervals(bed_path: Path):
    """
    Show amplicon start and end coordinates given a BED file of primer coordinates

    :arg bed_path: path of bed file
    """
    all_intervals = lib.amplicon_intervals(bed_path)
    sorted_by_chrom = sorted(all_intervals.items())
    for chrom, intervals in sorted_by_chrom:
        sorted_interval_keys = sorted(intervals, key=lambda x: (x[0], x[1]))
        for name in sorted_interval_keys:
            interval = intervals[name]
            print(f"{chrom}\t{interval[0]}\t{interval[1]}\t{name}")


def subset(scheme_dir: Path, chrom: str, out_dir: Path = Path("built")):
    """
    Extract a primer.bed and reference.fasta scheme subset for a single chromosome

    :arg scheme_dir: path of input scheme directory
    :arg chrom: name of chromosome for which to generate a subset scheme definition
    :arg out_dir: path of directory in which to save subset scheme
    """
    lib.subset(scheme_dir=scheme_dir, chrom=chrom, out_dir=out_dir)


def synchronise():
    """
    Retrieve/update local copy of remote primer scheme repository
    """
    lib.synchronise()


def plot(bed_path: Path, out_path: Path = Path("primer.html")):
    """
    Plot amplicon and primer coords from 7 column primer.bed

    :arg bed_path: path of primer.bed file
    :arg out_path: path of generated plot (with .html, .pdf, .png, or .svg extension)
    """
    lib.plot_primers(bed_path=bed_path, out_path=out_path)


def main():
    defopt.run(
        {
            "validate": validate,
            "build": build,
            "build-manifest": build_manifest,
            "hash-ref": hash_ref,
            "hash-bed": hash_bed,
            "diff": diff,
            "6to7": six_to_seven,
            "7to6": seven_to_six,
            "plot": plot,
            "show-intervals": amplicon_intervals,
            "show-discordant-primers": discordant_primers,
            "subset": subset,
            "sync": synchronise,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
