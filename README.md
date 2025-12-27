# IsoformExpression

A Dash-based interactive dashboard to visualize gene isoform structures
(exons and CDS) colored by mean transcript expression (TPM).

This tool was developed for exploratory analysis of alternative splicing
and isoform-specific expression patterns.

---

## Features

- Interactive gene search (AGI identifiers)
- Isoform block model visualization
  - exons (thin rectangles)
  - CDS (thicker rectangles)
  - introns (connecting lines)
- Color encoding based on **log1p(mean TPM)**
- Condition selection based on **genotype × timepoint**
- Color scale based on ColorBrewer (BuPu)

---

## Demo data

This repository contains **demo datasets only** located in:

annotation/testdata/
These files are **synthetic / reduced examples**

## Installation

Create a virtual environment (recommended):

```bash
python -m venv venv
source venv/bin/activate

install dependencies: 
requirements.txt

