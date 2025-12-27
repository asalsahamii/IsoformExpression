#!/usr/bin/env python3
"""
Build transcript-level expression tables from Salmon quant.sf files.

Steps:
1) Collect TPM values from all quant.sf files
2) Build transcript_expression.csv (all replicates)
3) Compute mean TPM per transcript / genotype / timepoint
4) Save transcript_expression_mean.csv

Expected folder structure:
BASE_DIR/
 ├── 7ko_LL18_1/quant.sf
 ├── 7ko_LL18_2/quant.sf
 ├── 7ko_LL18_3/quant.sf
 ├── wt_LL18_1/quant.sf
 └── ...

quant.sf naming pattern:
<genotype>_<timepoint>_<replicate>
"""

import os
import re
import argparse
import pandas as pd


# -------------------------------------------------
# Core logic
# -------------------------------------------------
def build_expression_table(base_dir: str) -> pd.DataFrame:
    """
    Parse all quant.sf files and return transcript-level TPM table
    with replicate information.
    """
    rows = []

    # folder pattern: genotype_timepoint_rep
    pattern = re.compile(r"^(7ko|7ox|8ox|WT|wt)_([A-Za-z0-9]+)_(\d+)$", re.IGNORECASE)

    for folder in sorted(os.listdir(base_dir)):
        folder_path = os.path.join(base_dir, folder)
        if not os.path.isdir(folder_path):
            continue

        quant_path = os.path.join(folder_path, "quant.sf")
        if not os.path.isfile(quant_path):
            continue

        m = pattern.match(folder)
        if not m:
            print(f"  Skipping folder (name mismatch): {folder}")
            continue

        genotype, timepoint, rep = m.groups()
        rep = int(rep)

        df = pd.read_csv(quant_path, sep="\t")

        for _, row in df.iterrows():
            rows.append(
                {
                    "transcript_id": row["Name"],
                    "TPM": row["TPM"],
                    "sample": folder,
                    "genotype": genotype,
                    "timepoint": timepoint,
                    "replicate": rep,
                }
            )

    expr_df = pd.DataFrame(rows)

    if expr_df.empty:
        raise RuntimeError(" No quant.sf files found or parsed!")

    return expr_df


def compute_mean_tpm(expr_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute mean TPM across replicates for each transcript / genotype / timepoint.
    """
    mean_df = (
        expr_df
        .groupby(
            ["transcript_id", "genotype", "timepoint"],
            as_index=False
        )
        .agg(
            mean_TPM=("TPM", "mean")
        )
    )

    return mean_df


# -------------------------------------------------
# CLI
# -------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Build transcript-level TPM tables from Salmon quant.sf files"
    )

    parser.add_argument(
        "--base-dir",
        required=True,
        help="Directory containing sample folders with quant.sf files",
    )

    parser.add_argument(
        "--out-dir",
        default="expression",
        help="Output directory (default: expression/)",
    )

    args = parser.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    out_dir = os.path.abspath(args.out_dir)

    os.makedirs(out_dir, exist_ok=True)

    out_expr = os.path.join(out_dir, "transcript_expression.csv")
    out_mean = os.path.join(out_dir, "transcript_expression_mean.csv")

    print(" Base directory:", base_dir)
    print(" Output directory:", out_dir)

    # ---- Build expression table ----
    expr_df = build_expression_table(base_dir)
    expr_df.to_csv(out_expr, index=False)

    print(" Saved:", out_expr)
    print("  Rows:", len(expr_df))
    print("  Unique transcripts:", expr_df["transcript_id"].nunique())
    print("  Unique samples:", expr_df["sample"].nunique())

    # ---- Compute mean TPM ----
    mean_df = compute_mean_tpm(expr_df)
    mean_df.to_csv(out_mean, index=False)

    print(" Saved:", out_mean)
    print("  Rows:", len(mean_df))
    print("  Conditions:", mean_df[["genotype", "timepoint"]].drop_duplicates().shape[0])


if __name__ == "__main__":
    main()
