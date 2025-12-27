# IsoformExpression

IsoformExpression is a Dash-based interactive dashboard for visualizing
gene isoform structures (exons and CDS) colored by mean transcript expression (TPM).

The tool is designed for exploratory analysis of alternative splicing
and isoform-specific expression patterns using transcript-level RNA-seq data.


---

### IsoformExpression provides:

- An interactive visualization of exon–intron structures alongside transcript expression
- Direct comparison of isoforms within a gene under different experimental conditions
- Rapid exploratory analysis without requiring manual plotting or genome browser inspection

The dashboard is **not intended to replace** existing quantification or statistical tools.
Instead, it is designed as a **complementary visualization layer** that supports:

- Quality control of transcript-level expression results
- Hypothesis generation regarding alternative splicing events
- Teaching and demonstration of isoform-specific regulation

By combining isoform structure and expression in a single interactive view, 
IsoformExpression fills the gap between statistical RNA-seq analysis and genome 
browser–based inspection.

### Complementary tools

IsoformExpression is designed to work alongside:

- Salmon / Kallisto – transcript quantification
- DESeq2 / edgeR – statistical differential expression analysis
- Genome browsers (IGV, UCSC) – detailed genomic inspection

--- 

## Features

- Interactive gene search (AGI identifiers)
- Isoform block model visualization
  - exons (rectangles)
  - UTR (thin Rectangles)
  - CDS (thicker rectangles)
  - introns (connecting lines)
- Color encoding based on **log1p(mean TPM of 3 Replicates)**
- Condition selection based on **genotype × timepoint**
- Color scale based on ColorBrewer (BuPu)

---

## Demo data

This repository contains **demo datasets only** located in:

annotation/testdata/

These files are **synthetic / reduced examples**
These subsets are intended solely for demonstration and testing of the dashboard.
Full experimental datasets are intentionally excluded.

---
## Installation

Create a virtual environment (recommended):

```bash
python -m venv venv
source venv/bin/activate

install dependencies: 

pip install -r requirements.txt


 
Citation

If you use this tool, please cite:

Asal Sahami Moghaddam,
IsoformExpression: A Dash-based dashboard for isoform-level expression visualization,
GitHub, 2025.