#!/usr/bin/env python3
"""
prepare_annotation.py

Create segments_atrtd3.csv from an atRTD3 GTF file.
Keeps only exon + CDS rows and extracts gene_id / transcript_id / exon_number.

Usage (CLI):
  python scripts/prepare_annotation.py \
    --gtf annotation/atRTD3_TS_21Feb22_transfix.gtf \
    --out annotation/segments_atrtd3.csv

Importable:
  from scripts.prepare_annotation import build_segments_df, write_segments_csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Optional

import pandas as pd


GTF_COLS = [
    "chrom", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute",
]


_ATTR_CACHE: dict[tuple[str, str], Optional[str]] = {}


def parse_attr(attr_str: str, key: str) -> Optional[str]:
    """
    Parse a GTF attributes string and return the value for a given key.
    Example attribute chunk: gene_id "AT1G01010"; transcript_id "AT1G01010.1";
    """
    s = str(attr_str)
    cache_key = (s, key)
    if cache_key in _ATTR_CACHE:
        return _ATTR_CACHE[cache_key]

    m = re.search(rf'{re.escape(key)} "([^"]+)"', s)
    val = m.group(1) if m else None
    _ATTR_CACHE[cache_key] = val
    return val


def build_segments_df(
    gtf_path: str | Path,
    keep_features: tuple[str, ...] = ("exon", "CDS"),
) -> pd.DataFrame:
    """
    Read a GTF and return a tidy segments DataFrame with:
      chrom, gene_id, transcript_id, feature, start, end, strand, exon_number

    - Keeps only exon + CDS by default.
    - exon_number is nullable Int64 (CDS lines may not always have it).
    """
    gtf_path = Path(gtf_path)

    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=GTF_COLS,
        dtype={"chrom": "string", "source": "string", "feature": "string", "attribute": "string"},
        low_memory=False,
    )

    segments = gtf[gtf["feature"].isin(list(keep_features))].copy()

    # Extract attributes
    segments["gene_id"] = segments["attribute"].map(lambda s: parse_attr(s, "gene_id"))
    segments["transcript_id"] = segments["attribute"].map(lambda s: parse_attr(s, "transcript_id"))
    segments["exon_number"] = segments["attribute"].map(lambda s: parse_attr(s, "exon_number"))

    # Clean types
    segments = segments.dropna(subset=["gene_id", "transcript_id"]).copy()
    segments["start"] = segments["start"].astype(int)
    segments["end"] = segments["end"].astype(int)

    # exon_number can be missing -> nullable integer
    segments["exon_number"] = pd.to_numeric(segments["exon_number"], errors="coerce").astype("Int64")

    # Keep only relevant columns (and in consistent order)
    segments = segments[
        ["chrom", "gene_id", "transcript_id", "feature", "start", "end", "strand", "exon_number"]
    ]

    # sort for nicer diffs / reproducibility
    segments = segments.sort_values(
        by=["gene_id", "transcript_id", "start", "end", "feature"],
        kind="mergesort",
    ).reset_index(drop=True)

    return segments


def write_segments_csv(
    gtf_path: str | Path,
    out_csv: str | Path = "annotation/segments_atrtd3.csv",
) -> Path:
    """
    Convenience wrapper: build segments df and write to CSV.
    Returns output path.
    """
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    df = build_segments_df(gtf_path)
    df.to_csv(out_csv, index=False)

    return out_csv


def _cli() -> int:
    parser = argparse.ArgumentParser(
        description="Extract exon+CDS segments from a GTF and write segments_atrtd3.csv"
    )
    parser.add_argument(
        "--gtf",
        default="annotation/atRTD3_TS_21Feb22_transfix.gtf",
        help="Input GTF path (default: annotation/atRTD3_TS_21Feb22_transfix.gtf)",
    )
    parser.add_argument(
        "--out",
        default="annotation/segments_atrtd3.csv",
        help="Output CSV path (default: annotation/segments_atrtd3.csv)",
    )

    args = parser.parse_args()

    out_path = write_segments_csv(args.gtf, args.out)
    # Print summary for the user
    df = pd.read_csv(out_path)
    print(f"Saved: {out_path}")
    print(f"Rows: {len(df)}")
    print("Columns:", list(df.columns))
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
