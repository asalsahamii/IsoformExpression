import dash
from dash import Dash, html, dcc, Input, Output
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
import plotly.colors as pc
import numpy as np

# ---- 1) Paths / Data ----
GENE_ANNOTATION_PATH = "annotation/testdata/Thalemine_gene_names_demo.csv"
SEGMENTS_PATH = "annotation/testdata/segments_demo.csv"   # exon + CDS
EXPR_MEAN_PATH = "annotation/testdata/transcript_expression_mean_demo.csv"

# ColorBrewer-like sequential scale (BuPu)
COLOR_SCALE = px.colors.sequential.BuPu

genes_df = pd.read_csv(GENE_ANNOTATION_PATH, sep=";")
segments_df = pd.read_csv(SEGMENTS_PATH)
expr_mean_df = pd.read_csv(EXPR_MEAN_PATH)

# ---- Dropdown 1: genes ----
gene_options = [
    {"label": f"{row['AGI']} – {row['Name']}", "value": row["AGI"]}
    for _, row in genes_df.iterrows()
]

# ---- Dropdown 2: conditions (16 = 4 genotypes x 4 timepoints) ----
cond_df = (
    expr_mean_df[["genotype", "timepoint"]]
    .dropna()
    .drop_duplicates()
    .copy()
)
cond_df["condition"] = cond_df["genotype"] + "_" + cond_df["timepoint"]
condition_values = sorted(cond_df["condition"].unique())
condition_options = [{"label": c.replace("_", " • "), "value": c} for c in condition_values]
DEFAULT_CONDITION = condition_values[0] if condition_values else None

# ---- 2) Dash App ----
app = Dash(__name__)

app.layout = html.Div(
    [
        html.Div(
            [
                html.H1("Isoform Expression Dashboard", style={"marginBottom": "10px"}),
                html.P(
                    "Isoform structure (exons + CDS) colored by mean TPM (replicates averaged).",
                    style={"color": "#555", "marginTop": "0px", "marginBottom": "30px"},
                ),
            ]
        ),

        # Controls
        html.Div(
            [
                html.Label("Select gene:", style={"fontWeight": "600"}),
                dcc.Dropdown(
                    id="gene-dropdown",
                    options=gene_options,
                    placeholder="Search for a gene (AGI)…",
                    style={"width": "60%"},
                ),

                html.Div(style={"height": "12px"}),

                html.Label("Select condition (genotype • timepoint):", style={"fontWeight": "600"}),
                dcc.Dropdown(
                    id="condition-dropdown",
                    options=condition_options,     # 16 conditions
                    value=DEFAULT_CONDITION,
                    placeholder="Search condition (e.g., 7ko_LL18)…",
                    style={"width": "60%"},
                    searchable=True,
                    clearable=True,
                ),

                html.Div(style={"height": "8px"}),

                html.P(
                    "Coloring is based on log1p(mean TPM) of the selected condition.",
                    style={"color": "#666", "fontSize": "13px"},
                ),
            ],
            style={
                "backgroundColor": "white",
                "padding": "20px",
                "borderRadius": "12px",
                "boxShadow": "0 2px 8px rgba(0,0,0,0.05)",
                "marginBottom": "25px",
            },
        ),

        # Output
        html.Div(
            id="gene-info",
            style={
                "backgroundColor": "white",
                "padding": "20px",
                "borderRadius": "12px",
                "boxShadow": "0 2px 8px rgba(0,0,0,0.05)",
            },
        ),
    ],
    style={
        "fontFamily": "system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif",
        "backgroundColor": "#f5f5f7",
        "minHeight": "100vh",
        "padding": "30px 60px",
    },
)

# ---- Helper: isoform plot (introns + exon thin + CDS thick) ----
def build_isoform_block_plot(segments_sub, gene_id, tpm_map, tpm_max, condition_label):
    segments_sub = segments_sub[segments_sub["feature"].isin(["exon", "CDS"])].copy()

    transcripts = sorted(segments_sub["transcript_id"].dropna().unique())
    y_map = {tid: i for i, tid in enumerate(transcripts)}

    x_min = int(segments_sub["start"].min())
    x_max = int(segments_sub["end"].max())

    shapes = []

    for tid in transcripts:
        y = y_map[tid]
        tdf = segments_sub[segments_sub["transcript_id"] == tid].copy()

        # TPM -> color (normalized)
        tpm = float(tpm_map.get(tid, 0.0))
        norm = (tpm / tpm_max) if tpm_max > 0 else 0.0
        norm = max(0.0, min(1.0, norm))
        iso_color = pc.sample_colorscale(COLOR_SCALE, norm)[0]

        # EXONS (thin, lighter) - black outline
        exons = tdf[tdf["feature"] == "exon"].sort_values("start")
        for _, r in exons.iterrows():
            shapes.append(
                dict(
                    type="rect",
                    xref="x",
                    yref="y",
                    x0=int(r["start"]),
                    x1=int(r["end"]),
                    y0=y - 0.10,
                    y1=y + 0.10,
                    fillcolor=iso_color,
                    opacity=0.95,
                    line=dict(width=1, color="black"),
                )
            )

        # INTRONS (lines between exons)
        exon_coords = exons[["start", "end"]].values
        for i in range(len(exon_coords) - 1):
            intron_start = int(exon_coords[i][1])
            intron_end = int(exon_coords[i + 1][0])
            shapes.append(
                dict(
                    type="line",
                    xref="x",
                    yref="y",
                    x0=intron_start,
                    x1=intron_end,
                    y0=y,
                    y1=y,
                    line=dict(width=2, color="rgba(0,0,0,0.55)"),
                )
            )

        # CDS (thicker, stronger) - black outline
        cds = tdf[tdf["feature"] == "CDS"].sort_values("start")
        for _, r in cds.iterrows():
            shapes.append(
                dict(
                    type="rect",
                    xref="x",
                    yref="y",
                    x0=int(r["start"]),
                    x1=int(r["end"]),
                    y0=y - 0.18,
                    y1=y + 0.18,
                    fillcolor=iso_color,
                    opacity=0.95,
                    line=dict(width=1, color="black"),
                )
            )

    fig = go.Figure()

    # Dummy marker trace to show colorbar
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                colorscale=COLOR_SCALE,
                cmin=0,
                cmax=max(float(tpm_max), 1e-9),
                color=[0],
                showscale=True,
                colorbar=dict(title=f"log1p(mean TPM)<br>({condition_label})"),
            ),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    # Invisible trace to set axes
    fig.add_trace(
        go.Scatter(
            x=[x_min, x_max],
            y=[0, max(0, len(transcripts) - 1)],
            mode="markers",
            opacity=0,
            showlegend=False,
            hoverinfo="skip",
        )
    )

    fig.update_layout(
        title=f"Isoform block model: {gene_id}",
        xaxis_title="Genomic position (bp)",
        yaxis=dict(
            title="Transcript isoform",
            tickmode="array",
            tickvals=list(y_map.values()),
            ticktext=list(y_map.keys()),
            range=[-1, len(transcripts)],
            tickfont=dict(color="#2F3E4E"),
        ),
        shapes=shapes,
        height=280 + 60 * len(transcripts),
        margin=dict(l=90, r=30, t=60, b=60),
        plot_bgcolor="#FAFAF7",
        paper_bgcolor="white",
        font=dict(
            family="system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif",
            size=14,
            color="#1F2D3D",
        ),
    )

    fig.update_xaxes(
        range=[x_min - 50, x_max + 50],
        showgrid=True,
        gridcolor="rgba(0,0,0,0.06)",
        zeroline=False,
    )
    fig.update_yaxes(showgrid=False, zeroline=False)

    return fig


# ---- 3) Callback ----
@app.callback(
    Output("gene-info", "children"),
    Input("gene-dropdown", "value"),
    Input("condition-dropdown", "value"),
)
def show_gene_info(selected_agi, selected_condition):
    if selected_condition is None:
        selected_condition = DEFAULT_CONDITION

    if selected_agi is None:
        return html.P(
            "Please select a gene above.",
            style={"color": "#666", "fontStyle": "italic"},
        )

    # Parse condition -> genotype, timepoint
    genotype, timepoint = selected_condition.split("_", 1)
    condition_label = selected_condition.replace("_", " • ")

    # Gene meta
    row = genes_df.loc[genes_df["AGI"] == selected_agi].iloc[0]
    name = row["Name"]

    # Segments for gene
    sub = segments_df[segments_df["gene_id"] == selected_agi].copy()
    if sub.empty:
        return html.Div(
            [
                html.H3(f"Selected gene: {selected_agi}"),
                html.P(f"Name / description: {name}"),
                html.P("No exon/CDS segments found for this gene.", style={"color": "#aa0000"}),
            ]
        )

    transcripts = sorted(sub["transcript_id"].dropna().unique())

    # Expression subset (mean TPM)
    expr_sub = expr_mean_df[
        (expr_mean_df["genotype"] == genotype)
        & (expr_mean_df["timepoint"] == timepoint)
        & (expr_mean_df["transcript_id"].isin(transcripts))
    ].copy()

    # log1p(mean TPM) for better dynamic range
    expr_sub["TPM_log"] = np.log1p(expr_sub["mean_TPM"])

    tpm_map = dict(zip(expr_sub["transcript_id"], expr_sub["TPM_log"]))
    tpm_max = float(expr_sub["TPM_log"].max()) if not expr_sub.empty else 0.0

    fig = build_isoform_block_plot(sub, selected_agi, tpm_map, tpm_max, condition_label)

    return html.Div(
        [
            html.H3(f"Selected gene: {selected_agi}"),
            html.P(f"Name / description: {name}"),
            html.P(f"Condition: {condition_label}", style={"color": "#555"}),
            html.Br(),
            dcc.Graph(figure=fig),
            html.P(
                "Legend: exon = thin rectangle (lighter), CDS = thicker rectangle, intron = line. "
                "Outline is constant black; fill color reflects log1p(mean TPM).",
                style={"color": "#666", "fontSize": "13px", "marginTop": "8px"},
            ),
        ]
    )


if __name__ == "__main__":
    app.run(debug=True)

