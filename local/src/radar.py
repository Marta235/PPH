#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go



input_dir = snakemake.input.enrichment_dir
output_nt = snakemake.output.output_nt
output_cetux = snakemake.output.output_cetux



categories = [
    "ISC", "TA", "TA2", "Late CC", "Intermediate AE", "Mature AE",
    "AE2", "SI_FAE", "Sec. Progenitor", "BEST4", "EEC", "Goblet", "Paneth"
]

sample_colors = {
    "CRC0542_NT72h_1": "#fc8d62",
    "CRC0322_NT_1_3000": "#8da0cb",
    "CRC0327_NT_2": "#66c2a5",
    "CRC0322_cetux_1": "#8da0cb",
    "CRC0542_cTX72h_1": "#fc8d62",
    "CRC0327_cetux_2": "#66c2a5",
}

fill_missing = 0.0
nt_pattern = "NT"


frames = []

for file in sorted(os.listdir(input_dir)):
    if "tsv" not in file:
        continue

    path = os.path.join(input_dir, file)
    if not os.path.isfile(path):
        continue

    data = pd.read_csv(path, sep="\t", header=0)


    data["sample"] = file.split(".")[1]

    frames.append(data)

gsea_res = pd.concat(frames, ignore_index=True)



def build_sample_vector(data, categories, fill_missing=0.0):
    values = []

    for cat in categories:
        s = data.loc[data["Description"] == cat, "NES"]
        if len(s) > 0:
            values.append(float(s.iloc[0]))
        else:
            values.append(float(fill_missing))

    return values


def get_global_range(gsea_res, categories):
    all_values = []

    for sample in gsea_res["sample"].unique():
        data = gsea_res.loc[gsea_res["sample"] == sample]
        all_values.extend(build_sample_vector(data, categories, fill_missing=fill_missing))

    all_values.append(0.0)
    m = max(abs(np.nanmin(all_values)), abs(np.nanmax(all_values)))

    return -m, m


def make_layout(global_min, global_max, title):
    tickvals = np.linspace(global_min, global_max, 5)
    ticktext = [f"{x:.1f}" for x in tickvals]

    return dict(
        template="plotly_white",
        paper_bgcolor="white",
        plot_bgcolor="white",
        title=dict(text=title, x=0.5, xanchor="center"),
        polar=dict(
            bgcolor="white",
            radialaxis=dict(
                showgrid=True,
                gridcolor="lightgray",
                showline=True,
                linecolor="black",
                range=[global_min, global_max],
                tickmode="array",
                tickvals=list(tickvals),
                ticktext=ticktext,
            ),
            angularaxis=dict(
                showgrid=True,
                gridcolor="lightgray",
                showline=True,
                linecolor="black",
            ),
        ),
        font=dict(color="black"),
    )


def make_radar(gsea_res, samples, title, global_min, global_max):
    fig = go.Figure()

    for sample in samples:
        data = gsea_res.loc[gsea_res["sample"] == sample]
        r = build_sample_vector(data, categories, fill_missing=fill_missing)
        color = sample_colors.get(sample, "#999999")

        fig.add_trace(
            go.Scatterpolar(
                r=r,
                theta=categories,
                name=sample,
                fill="toself",
                opacity=0.4,
                fillcolor="#099963",
                line=dict(color=color, width=2),
            )
        )

    fig.update_layout(**make_layout(global_min, global_max, title))
    return fig



samples = sorted(gsea_res["sample"].unique())

nt_samples = [s for s in samples if nt_pattern in s]
cetux_samples = [s for s in samples if nt_pattern not in s]

global_min, global_max = get_global_range(gsea_res, categories)

fig_nt = make_radar(
    gsea_res=gsea_res,
    samples=nt_samples,
    title="NT samples",
    global_min=global_min,
    global_max=global_max,
)

fig_cetux = make_radar(
    gsea_res=gsea_res,
    samples=cetux_samples,
    title="Treatment samples",
    global_min=global_min,
    global_max=global_max,
)

Path(output_nt).parent.mkdir(parents=True, exist_ok=True)
Path(output_cetux).parent.mkdir(parents=True, exist_ok=True)

fig_nt.write_image(output_nt)
fig_cetux.write_image(output_cetux)

print("Radar plots exported")
print("NT samples:", nt_samples)
print("Treatment samples:", cetux_samples)
print("Output NT:", output_nt)
print("Output treatment:", output_cetux)