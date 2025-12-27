import pandas as pd

# ---------- INPUT (original data) ----------
GENE_ANNOTATION_PATH = "annotation/Thalemine_gene_names.csv"
SEGMENTS_PATH = "annotation/segments_atrtd3.csv"
EXPR_MEAN_PATH = "expression/transcript_expression_mean.csv"

# ---------- OUTPUT (demo data) ----------
OUT_GENE_ANNOTATION = "../IsoformExpression/annotation/testdata/Thalemine_gene_names_demo.csv"
OUT_SEGMENTS = "../IsoformExpression/annotation/testdata/segments_demo.csv"
OUT_EXPR_MEAN = "../IsoformExpression/annotation/testdata/transcript_expression_mean_demo.csv"

# ---------- CHOOSE 2 GENES + 2 CONDITIONS ----------
GENES_KEEP = ["AT1G01010", "AT1G01020"]   # <-- ändere hier deine 2 Gene

CONDITIONS_KEEP = [
    ("7ko", "LL18"),
    ("WT",  "LL18"),
]

def main():
    genes_df = pd.read_csv(GENE_ANNOTATION_PATH, sep=";")
    segments_df = pd.read_csv(SEGMENTS_PATH)
    expr_mean_df = pd.read_csv(EXPR_MEAN_PATH)

    # ---- 1) gene annotation subset ----
    genes_demo = genes_df[genes_df["AGI"].isin(GENES_KEEP)].copy()

    # ---- 2) segments subset (only chosen genes; exon+CDS) ----
    segments_demo = segments_df[
        (segments_df["gene_id"].isin(GENES_KEEP))
        & (segments_df["feature"].isin(["exon", "CDS"]))
    ].copy()

    # if a gene has no segments -> stop early (prevents your dashboard error)
    found_genes = sorted(segments_demo["gene_id"].unique())
    missing = sorted(set(GENES_KEEP) - set(found_genes))
    if missing:
        raise SystemExit(
            f"ERROR: These genes have no exon/CDS segments in segments_atrtd3.csv: {missing}\n"
            f"Pick other genes."
        )

    # ---- 3) get all transcripts from these segments ----
    transcripts_keep = sorted(segments_demo["transcript_id"].dropna().unique())

    # ---- 4) expression mean subset for only 2 conditions + those transcripts ----
    cond_mask = False
    for g, t in CONDITIONS_KEEP:
        cond_mask = cond_mask | ((expr_mean_df["genotype"] == g) & (expr_mean_df["timepoint"] == t))

    expr_demo = expr_mean_df[cond_mask & expr_mean_df["transcript_id"].isin(transcripts_keep)].copy()

    # ---- sanity check: expression exists for both conditions ----
    cond_present = expr_demo[["genotype", "timepoint"]].drop_duplicates()
    print("Conditions present in demo:", cond_present.values.tolist())

    # ---- write demo files ----
    genes_demo.to_csv(OUT_GENE_ANNOTATION, sep=";", index=False)
    segments_demo.to_csv(OUT_SEGMENTS, index=False)
    expr_demo.to_csv(OUT_EXPR_MEAN, index=False)

    print("✅ Demo files created:")
    print(" -", OUT_GENE_ANNOTATION, f"({len(genes_demo)} rows)")
    print(" -", OUT_SEGMENTS, f"({len(segments_demo)} rows)")
    print(" -", OUT_EXPR_MEAN, f"({len(expr_demo)} rows)")

if __name__ == "__main__":
    main()

