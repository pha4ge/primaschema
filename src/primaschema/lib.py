import logging
from pathlib import Path

import altair as alt
import pandas as pd
from primalbedtools.scheme import Scheme

logger = logging.getLogger(__name__)

SCHEME_BED_FIELDS = ["chrom", "chromStart", "chromEnd", "name", "poolName", "strand"]
PRIMER_BED_FIELDS = SCHEME_BED_FIELDS + ["sequence"]
POSITION_FIELDS = ["chromStart", "chromEnd"]
MANDATORY_FILES = ("primer.bed", "reference.fasta", "info.yml")


def _scheme_to_primer_df(scheme: Scheme) -> pd.DataFrame:
    rows = [
        {
            "chrom": bedline.chrom,
            "chromStart": bedline.start,
            "chromEnd": bedline.end,
            "name": bedline.primername,
            "poolName": bedline.pool,
            "strand": bedline.strand,
        }
        for bedline in scheme.bedlines
    ]
    return pd.DataFrame.from_records(rows, columns=SCHEME_BED_FIELDS)


def plot_primers(bed_path: Path, out_path: Path = Path("primer.html")) -> None:
    """
    Plot amplicon and primer positions from a 7 column primer.bed file
    Requires primers to be named {$scheme_id}_{$amplicon_id}_{LEFT|RIGHT}_{1|2|3…}
    Plots one vertical panel per pool per reference chromosome
    Supported out_path extensions: html (interactive), pdf, png, svg
    """
    scheme = Scheme.from_file(bed_path)  # type: ignore
    scheme.sort_bedlines()
    primer_df = _scheme_to_primer_df(scheme)
    primer_df["amplicon"] = primer_df["name"].str.split("_").str[1]
    primer_df["poolName"] = primer_df["poolName"].astype(str)

    amp_df = (
        primer_df.groupby(["chrom", "amplicon", "poolName"])
        .agg(min_start=("chromStart", "min"), max_end=("chromEnd", "max"))
        .reset_index()
    )
    amp_df["is_amplicon"] = True
    primer_df["is_amplicon"] = False
    amp_df = amp_df.rename(columns={"min_start": "chromStart", "max_end": "chromEnd"})
    combined_df = pd.concat([primer_df, amp_df], ignore_index=True)

    primer_marks = (
        alt.Chart(combined_df)
        .transform_filter(alt.datum.is_amplicon == False)  # noqa
        .mark_line(size=15)
        .encode(
            x=alt.X("chromStart:Q", title=None),
            x2="chromEnd:Q",
            y=alt.Y("poolName:O", title="pool", scale=alt.Scale(padding=1)),
            color=alt.Color("strand:N").scale(scheme="set2"),
            tooltip=[
                alt.Tooltip("name:N", title="Primer name"),
                alt.Tooltip("chromStart:Q", title="start"),
                alt.Tooltip("chromEnd:Q", title="end"),
            ],
        )
        .properties(
            width=800,
        )
    )
    amplicon_marks = (
        alt.Chart(combined_df)
        .transform_filter(alt.datum.is_amplicon == True)  # noqa
        .mark_rule(strokeWidth=2)
        .encode(
            x=alt.X("chromStart:Q", title=None),
            x2="chromEnd:Q",
            y=alt.Y("poolName:O", title="pool", scale=alt.Scale(padding=1)),
            tooltip=[
                alt.Tooltip("amplicon:N", title="Amplicon name"),
                alt.Tooltip("chromStart:Q", title="Min primer start"),
                alt.Tooltip("chromEnd:Q", title="Max primer end"),
            ],
        )
        .properties(
            width=800,
        )
    )
    combined_chart = (
        alt.layer(primer_marks, amplicon_marks)
        .facet(row=alt.Row("chrom:O", header=alt.Header(labelOrient="top"), title=""))
        .configure_axis(domain=False, ticks=False)
    )

    combined_chart.interactive().save(str(out_path))
    logger.debug(f"Plot saved ({out_path})")
