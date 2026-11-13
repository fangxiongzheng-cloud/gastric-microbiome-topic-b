#!/usr/bin/env python3
"""Reproduce the standalone Figure 8A v3 presentation candidate.

Run from any directory: python3 -B scripts/plot_figure8a_v3_presentation.py
All writes, including library caches, stay inside this revision directory.
Scientific values come only from its prepared inputs. No inference is run.
Layout coordinates are millimetres; font sizes and line widths are points.
"""

from __future__ import annotations

import csv
from decimal import Decimal
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
INPUTS, OUTPUTS, QA = (ROOT / name for name in ("inputs", "outputs", "qa"))
for directory in (OUTPUTS, QA, QA / ".mplconfig", QA / ".cache", QA / ".tmp"):
    directory.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(QA / ".mplconfig")
os.environ["XDG_CACHE_HOME"] = str(QA / ".cache")
os.environ["TMPDIR"] = str(QA / ".tmp")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path as MplPath
import numpy as np
from PIL import Image

WIDTH_MM, HEIGHT_MM = 170.0, 160.0
# Integer pixel/mm gives exact 170 mm in the PNG pHYs metadata.
PNG_DPI, PREVIEW_DPI = 381, 190.5  # 15 and 7.5 pixels/mm, respectively.
BASE = "Figure8A_v3_presentation"
POSITIVE, NEGATIVE = "#D55E00", "#0072B2"  # Okabe–Ito vermillion and blue.
INK, BORDER = "#222222", "#777777"
LINESTYLES = {"solid": "solid", "dashed": (0, (5.0, 3.0)),
              "dashdot": (0, (7.0, 2.8, 1.2, 2.8))}
JOURNAL_URL = "https://link.springer.com/journal/40168/submission-guidelines"
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.1,
    "svg.fonttype": "none", "svg.hashsalt": "Figure8A-v3-presentation",
    "pdf.fonttype": 42, "ps.fonttype": 42, "pdf.use14corefonts": False,
    "text.usetex": False, "axes.unicode_minus": True,
    "lines.scale_dashes": False, "savefig.facecolor": "white",
})

# Display geometry only: no node, effect, sign or statistical tier is inferred.
NODES = {
    "HP": (5, 85, 39, 23),
    "TISSUE": (120, 134, 45, 17),
    "COMM": (120, 106, 45, 12),
    "M1": (120, 86, 45, 17),
    "M2": (120, 70, 45, 10),
    "M3": (120, 48, 45, 18),
    "GRADE": (120, 27, 45, 16),
    "COMP": (5, 47, 39, 13),
}
ROUTES = {
    "A1": [(44, 106), (64, 146), (120, 146)],
    "F2": [(44, 102), (64, 137), (120, 137)],
    "E1c": [(44, 100), (64, 111), (120, 111)],
    "H-M1": [(44, 98), (64, 97), (120, 97)],
    "A4": [(44, 96), (64, 90), (120, 90)],
    "H-M2": [(44, 94), (64, 75), (120, 75)],
    "H-M3": [(44, 92), (64, 64), (120, 64)],
    "E9 niche": [(44, 90), (64, 53), (120, 53)],
    "E10 IM": [(41, 85), (64, 39), (120, 39)],
    "E10 atrophy": [(35, 85), (64, 31), (120, 31)],
    "E1a": [(162, 134), (162, 118)],
}
EXPECTED_PAIRS = {
    "A1": ("HP", "TISSUE"), "F2": ("HP", "TISSUE"),
    "H-M1": ("HP", "M1"), "H-M2": ("HP", "M2"),
    "H-M3": ("HP", "M3"), "E1c": ("HP", "COMM"),
    "E1a": ("COMM", "TISSUE"), "E10 IM": ("HP", "GRADE"),
    "E10 atrophy": ("HP", "GRADE"), "A4": ("HP", "M1"),
    "E9 niche": ("HP", "M3"),
}
EXPECTED_SOLID = {"A1", "F2", "H-M1", "H-M2", "H-M3"}
CSV_NAMES = ("figure8a_source_associations.csv", "e1a_display_aggregation.csv")


def read_csv(name):
    with (INPUTS / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def v2_hashes():
    """Read only the six paths explicitly listed in the supplied baseline."""
    records = []
    for line in (INPUTS / "v2_candidate_hashes_before_v3.txt").read_text().splitlines():
        expected, filename = line.split(maxsplit=1)
        path = Path(filename.strip())
        actual = sha256(path)
        records.append({"path": str(path), "expected_sha256": expected,
                        "actual_sha256": actual, "matches": actual == expected})
    return records


def width_tier(effect):
    # Point widths: preserved original scripts/s10_e11_schematic.py:39.
    # Thresholds: v0.6d:88–89 and freeze_log.md:68; see input provenance.
    a = abs(effect)
    return ("<0.2", 1.0) if a < 0.2 else (("0.2–0.4", 2.5) if a <= 0.4 else (">0.4", 4.5))


def signed(effect):
    return f"{effect:+.3f}".replace("-", "−")


def edge_label(row):
    a, e = row["analysis_id"], signed(float(row["effect"]))
    if a == "A1":
        return f"A1 · RNA   r = {e}"
    if a == "F2":
        return f"F2 · 16S   r = {e}"
    if a.startswith("H-M"):
        return f"{a} · Pooled ρ = {e}"
    if a == "A4":
        return f"A4 · 16S CLR   ρ = {e}"
    if a == "E1c":
        return f"E1c · ΔH. pylori–ΔGemella\nMedian ρ = {e}"
    if a == "E1a":
        return f"E1a · {e}\n√(median R²), F2 + Korea\n2 cohorts; stored P < 0.05"
    if a == "E10 IM":
        return f"E10 · Intestinal metaplasia\nρ = {e}"
    if a == "E10 atrophy":
        return f"E10 · Atrophy   ρ = {e}"
    if a == "E9 niche":
        return f"E9 · A4GNT · ρ = {e}\nMaximum |ρ| among qualifying E9 genes"
    raise ValueError(a)


def draw_panel(fig, ax, associations):
    ax.set(xlim=(0, WIDTH_MM), ylim=(0, HEIGHT_MM))
    ax.set_axis_off()
    texts, edge_objects = [], []

    def text(x, y, content, *, gid, size=7.1, role="key", **kwargs):
        options = dict(color=INK, va="center", ha="left", zorder=6,
                       linespacing=1.24, clip_on=False)
        options.update(kwargs)
        artist = ax.text(x, y, content, fontsize=size, **options)
        artist.set_gid(gid)
        texts.append((artist, {"id": gid, "text": content, "role": role,
                              "font_size_pt": size}))
        return artist

    text(5, 155.5, "A", gid="panel-letter", size=9.5, weight="bold", role="panel_letter")
    for row in associations:
        a = row["analysis_id"]
        effect = float(row["effect"])
        route = ROUTES[a]
        colour = POSITIVE if effect > 0 else NEGATIVE
        path = MplPath(route, [MplPath.MOVETO] + [MplPath.LINETO] * (len(route) - 1))
        edge = PathPatch(path, fill=False, edgecolor=colour,
                         linewidth=width_tier(effect)[1], linestyle=LINESTYLES[row["style"]],
                         capstyle="butt", joinstyle="round", zorder=1)
        edge.set_gid("association-edge-" + a.replace(" ", "-"))
        ax.add_patch(edge)
        if a == "E1a":
            label = text(121.5, 126, edge_label(row), gid="association-label-E1a",
                         size=7.0, role="edge_label")
        else:
            label = text(67, route[-1][1] + 1.8, edge_label(row),
                         gid="association-label-" + a.replace(" ", "-"),
                         role="edge_label", va="bottom", size=7.1)
        edge_objects.append((edge, label, row))

    for node_id, (x, y, w, h) in NODES.items():
        node = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.9",
                              facecolor="#F7F8F9" if node_id == "HP" else "white",
                              edgecolor=BORDER, linewidth=0.75, zorder=3)
        node.set_gid("node-" + node_id)
        ax.add_patch(node)
    text(24.5, 103, "Helicobacter pylori", gid="node-HP-species", size=7.5,
         style="italic", ha="center", role="node_label")
    text(24.5, 99.1, "(H. pylori) load", gid="node-HP-load", size=7.5,
         weight="bold", ha="center", role="node_label")
    text(24.5, 92.5, "transcript reads /\n16S abundance", gid="node-HP-measure",
         ha="center", size=7.0, role="node_label")
    text(123, 145.3, "Tumour versus\nadjacent normal", gid="node-TISSUE-label",
         size=7.5, weight="bold", role="node_label")
    text(123, 136.3, "(tissue state)", gid="node-TISSUE-description", size=7.0, role="node_label")
    text(123, 112, "Non-Helicobacter\ncommunity (genera)", gid="node-COMM-label",
         size=7.4, weight="bold", role="node_label")
    text(123, 100.1, "M1 · Epithelial identity", gid="node-M1-label",
         size=7.3, weight="bold", role="node_label")
    text(123, 92.8, "intestinal metaplasia −\ngastric differentiation",
         gid="node-M1-description", size=7.0, role="node_label")
    text(123, 75, "M2 · Oxyntic programme", gid="node-M2-label",
         size=7.3, weight="bold", role="node_label")
    text(123, 57, "M3 · Mucin / adhesion\nniche (+E9 genes)", gid="node-M3-label",
         size=7.4, weight="bold", role="node_label")
    text(123, 35, "Molecular metaplasia /\natrophy grade (E10)", gid="node-GRADE-label",
         size=7.4, weight="bold", role="node_label")
    text(8, 55.3, "Cell composition", gid="node-COMP-label", size=7.4,
         weight="bold", role="node_label")
    text(8, 50.8, "(deconvolved)", gid="node-COMP-description", size=7.0, role="node_label")
    text(5, 43.5, "No association edge", gid="node-COMP-note", size=7.0, role="node_note")

    # Compact scientific key. All audit/production information is external.
    ax.add_line(Line2D([5, 165], [23, 23], color="#AAAAAA", lw=0.55))

    def sample(x0, x1, y, *, width=1.5, colour=INK, style="solid", gid):
        line = Line2D([x0, x1], [y, y], linewidth=width, color=colour,
                      linestyle=LINESTYLES[style], solid_capstyle="butt", dash_capstyle="butt")
        line.set_gid(gid)
        ax.add_line(line)

    sample(5, 14, 18.9, gid="key-solid")
    text(17, 18.9, "Primary tests: frozen pass rule met", gid="key-primary", size=6.8)
    sample(69, 79, 18.9, style="dashdot", gid="key-dashdot")
    text(82, 18.9, "A4: local secondary", gid="key-secondary", size=6.8)
    sample(125, 135, 18.9, style="dashed", gid="key-dashed")
    text(138, 18.9, "Exploratory", gid="key-exploratory", size=6.8)
    text(5, 13.1, "Sign", gid="key-sign-label", size=7.0, weight="bold")
    sample(17, 25, 13.1, width=2.0, colour=POSITIVE, gid="key-positive-sample")
    text(28, 13.1, "Positive (+)", gid="key-positive-label", size=7.0)
    sample(56, 64, 13.1, width=2.0, colour=NEGATIVE, gid="key-negative-sample")
    text(67, 13.1, "Negative (−)", gid="key-negative-label", size=7.0)
    text(106, 13.1, "Association ≠ causation", gid="key-association", size=7.0)
    sample(5, 14, 6.8, width=1.0, gid="key-width-small")
    text(17, 6.8, "|e| < 0.2 · 1.0 pt", gid="key-width-small-label", size=7.0)
    sample(57, 66, 6.8, width=2.5, gid="key-width-medium")
    text(69, 6.8, "0.2 ≤ |e| ≤ 0.4 · 2.5 pt", gid="key-width-medium-label", size=7.0)
    sample(125, 134, 6.8, width=4.5, gid="key-width-large")
    text(137, 6.8, "|e| > 0.4 · 4.5 pt", gid="key-width-large-label", size=7.0)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inverse = ax.transData.inverted()
    text_records = []
    for artist, record in texts:
        bounds = artist.get_window_extent(renderer)
        xy = inverse.transform([[bounds.x0, bounds.y0], [bounds.x1, bounds.y1]])
        text_records.append({**record, "bounds_mm": xy.ravel().tolist(),
                             "has_background_mask": artist.get_bbox_patch() is not None})
    edge_records = []
    for edge, label, row in edge_objects:
        edge_records.append({**row, "effect_used": float(row["effect"]),
                             "linewidth_used_pt": float(edge.get_linewidth()),
                             "colour": edge.get_edgecolor(), "svg_id": edge.get_gid(),
                             "route_mm": edge.get_path().vertices.tolist(),
                             "arrowheads": False, "label_text": label.get_text()})
    return {"texts": text_records, "edges": edge_records,
            "nodes": [{"node_id": n, "bounds_mm": [x, y, x+w, y+h]}
                      for n, (x, y, w, h) in NODES.items()]}


def rect_overlap(a, b, pad=0):
    return (a[0] - pad < b[2] and a[2] + pad > b[0]
            and a[1] - pad < b[3] and a[3] + pad > b[1])


def route_samples(route, n=401):
    return np.concatenate([np.linspace(a, b, n) for a, b in zip(route[:-1], route[1:])])


def point_rect_distance(points, rect):
    dx = np.maximum(np.maximum(rect[0] - points[:, 0], 0), points[:, 0] - rect[2])
    dy = np.maximum(np.maximum(rect[1] - points[:, 1], 0), points[:, 1] - rect[3])
    return np.hypot(dx, dy)


def segment_distance(a, b, c, d):
    """Exact Euclidean centreline distance between two nondegenerate segments."""
    a, b, c, d = (np.asarray(p, dtype=float) for p in (a, b, c, d))
    v, w = b-a, d-c
    def cross(p, q):
        return p[0]*q[1] - p[1]*q[0]
    denominator = cross(v, w)
    if abs(denominator) > 1e-12:
        t, u = cross(c-a, w) / denominator, cross(c-a, v) / denominator
        if 0 <= t <= 1 and 0 <= u <= 1:
            return 0.0
    def point_segment(p, q, s):
        delta = s-q
        fraction = np.clip(np.dot(p-q, delta) / np.dot(delta, delta), 0, 1)
        return float(np.linalg.norm(p - (q + fraction*delta)))
    return min(point_segment(a, c, d), point_segment(b, c, d),
               point_segment(c, a, b), point_segment(d, a, b))


def endpoint_nodes(point):
    px, py = point
    matches = []
    for name, (x, y, w, h) in NODES.items():
        on_vertical = (abs(px-x) < 1e-8 or abs(px-x-w) < 1e-8) and y <= py <= y+h
        on_horizontal = (abs(py-y) < 1e-8 or abs(py-y-h) < 1e-8) and x <= px <= x+w
        if on_vertical or on_horizontal:
            matches.append(name)
    return matches


def layout_checks(manifest):
    text_overlaps, clipped, stroke_hits, node_hits, escaped_nodes = [], [], [], [], []
    for a, b in itertools.combinations(manifest["texts"], 2):
        if rect_overlap(a["bounds_mm"], b["bounds_mm"], pad=0.05):
            text_overlaps.append([a["id"], b["id"]])
    for t in manifest["texts"]:
        x0, y0, x1, y1 = t["bounds_mm"]
        if x0 < 0 or y0 < 0 or x1 > WIDTH_MM or y1 > HEIGHT_MM:
            clipped.append(t["id"])
        if t["role"] == "node_label":
            n = t["id"].split("-")[1]
            x, y, w, h = NODES[n]
            if x0 < x + 0.5 or y0 < y + 0.5 or x1 > x+w-0.5 or y1 > y+h-0.5:
                escaped_nodes.append(t["id"])
    minimum_label_clearance = math.inf
    for edge in manifest["edges"]:
        samples = route_samples(edge["route_mm"])
        radius = edge["linewidth_used_pt"] * 25.4 / 72 / 2
        for t in manifest["texts"]:
            clearance = float(point_rect_distance(samples, t["bounds_mm"]).min()) - radius
            minimum_label_clearance = min(minimum_label_clearance, clearance)
            if clearance < 0.15:
                stroke_hits.append({"edge": edge["analysis_id"], "text": t["id"],
                                    "clearance_mm": clearance})
        for node in manifest["nodes"]:
            x0, y0, x1, y1 = node["bounds_mm"]
            inside = ((samples[:, 0] > x0+0.05) & (samples[:, 0] < x1-0.05)
                      & (samples[:, 1] > y0+0.05) & (samples[:, 1] < y1-0.05))
            if inside.any():
                node_hits.append([edge["analysis_id"], node["node_id"]])
    edge_collisions, endpoint_errors, edge_clearances = [], [], []
    for a, b in itertools.combinations(manifest["edges"], 2):
        centre_distance = min(segment_distance(p, q, r, s)
                              for p, q in zip(a["route_mm"][:-1], a["route_mm"][1:])
                              for r, s in zip(b["route_mm"][:-1], b["route_mm"][1:]))
        gap = centre_distance - (a["linewidth_used_pt"] + b["linewidth_used_pt"]) * 25.4 / 72 / 2
        edge_clearances.append({"edges": [a["analysis_id"], b["analysis_id"]], "gap_mm": gap})
        if gap < 0:
            edge_collisions.append(edge_clearances[-1])
    for edge in manifest["edges"]:
        first, last = endpoint_nodes(edge["route_mm"][0]), endpoint_nodes(edge["route_mm"][-1])
        if len(first) != 1 or len(last) != 1 or set(first + last) != set(EXPECTED_PAIRS[edge["analysis_id"]]):
            endpoint_errors.append({"analysis_id": edge["analysis_id"], "first": first, "last": last})
    return {"text_overlaps": text_overlaps, "clipped_text": clipped,
            "text_outside_own_node": escaped_nodes, "stroke_text_intersections": stroke_hits,
            "route_enters_node_interior": node_hits,
            "minimum_stroke_to_text_clearance_mm": minimum_label_clearance,
            "edge_stroke_collisions": edge_collisions,
            "minimum_edge_stroke_gap_mm": min(r["gap_mm"] for r in edge_clearances),
            "endpoint_identity_errors": endpoint_errors,
            "label_background_masks": sum(t["has_background_mask"] for t in manifest["texts"]),
            "all_clear": not any((text_overlaps, clipped, escaped_nodes, stroke_hits,
                                   node_hits, edge_collisions, endpoint_errors))}


TITLE = "Associations linking Helicobacter pylori load, gastric tissue state and host programmes"
LEGEND = """Eleven undirected associations link *Helicobacter pylori* (*H. pylori*) load (transcript reads or 16S abundance), tissue state, non-*Helicobacter* community and host programmes. M1 denotes epithelial identity (intestinal metaplasia minus gastric differentiation), M2 the oxyntic programme, and M3 the mucin/adhesion niche. Deconvolved cell composition is retained without an association edge.

Solid lines denote primary pre-specified tests that met their frozen pass rule (A1, F2, H-M1–M3); they do not imply that all robustness or sensitivity checks passed. A1 is non-blind locally pre-specified consistency; F2 is the fully blind gradient validation. The HKSJ-max1 M1 result (P = 0.159) is unsupported and is not represented as supported by H-M1. H-M2/M3 are Korea-dominated. A4 is locally pre-specified secondary, with prior interim inspection disclosed. Dashed lines denote exploratory associations.

Labels retain signed stored effects: rank-biserial r (A1/F2), DerSimonian–Laird-pooled ρ (H-M1–M3), primary Spearman ρ (A4), median same-sign-cohort ρ for ΔH. pylori–ΔGemella (E1c), and descriptive Spearman ρ for the separate E10 intestinal-metaplasia and atrophy edges. E1a shows descriptive √(median stored R²) = 0.07157835860316486 for F2 and Korea only (two cohorts), selected by the freeze-log stored P < 0.05 rule; PRJNA310127 (P = 0.0714) remains in the source table and is excluded only from this displayed summary. E1a's plus sign indicates non-negative magnitude, not directional correlation. E9 displays the maximum |ρ| among qualifying E9 genes (A4GNT), retaining its signed value.

Colour encodes statistical sign; discrete widths encode |effect| bins (<0.2: 1.0 pt; 0.2–0.4: 2.5 pt; >0.4: 4.5 pt). Effect metrics differ across analyses. Associations do not establish causation. No new statistical analyses were performed."""


def write_notes():
    caption = f"# {TITLE}\n\n{LEGEND}\n"
    (OUTPUTS / "proposed_publication_caption.md").write_text(caption, encoding="utf-8")
    note = f"""# Figure 8A v3 presentation candidate: technical note

This individual panel is a presentation revision of v2. It contains the same
eight nodes and 11 undirected associations, using only this folder's prepared
CSVs. No statistical inference, refit, permutation, correlation, new cohort
selection or causal claim is introduced. Figure S9 and other panels are not
composited. v2 and all official figures, frozen plans and manuscript files
remain untouched.

## Declared physical check size and source

The figure is **170 × 160 mm**. SVG and PDF use fixed dimensions without a
tight bounding-box crop. This leaves 65 mm of the cited 225 mm combined
figure-and-legend height for a manuscript legend; actual typeset placement is
not certified here. This is a **declared check size, not a journal-compliance
claim**.

The [official Microbiome submission guidelines]({JOURNAL_URL}) were checked
on 2026-11-13. The page lists BMC final-PDF widths of 85 or 170 mm, a combined
figure-and-legend height up to 225 mm, approximately 300 dpi, embedded fonts,
and strokes wider than 0.25 pt. It also places titles (≤15 words) and legends
(≤300 words) in the manuscript and keys in the graphic. The prepared
`inputs/journal_size_note.md` records the earlier size provenance. Those
targets informed this revision; no submission or publisher acceptance is
claimed.

The production PNG is **381 dpi** (15 pixels/mm), **2550 × 2400 px** with pHYs
dpi metadata; the convenient preview is **190.5 dpi** (7.5 pixels/mm),
**1275 × 1200 px**. These choices represent 170 × 160 mm exactly in integer
pixels. The physical-size proof is a separate 190 × 188 mm canvas at
190.5 dpi; its 170 mm figure region is aligned to a millimetre ruler. Print
the proof at 100%, disabling fit-to-page. Screen zoom does not establish a
physical viewing size.

On-figure text is ≥6.8 pt (association labels 7.0–7.1 pt; node text 7.0–7.5 pt)
at the declared width. All strokes, including node borders and the key
divider, exceed 0.25 pt. DejaVu Sans is used throughout. SVG text stays
editable (`svg.fonttype = none`); PDF and PostScript font types are set to 42.
The exported PDF's embedding and absence of raster images are checked with
Poppler utilities, with raw reports under `qa/`.

## Presentation changes from v2

- Replace the radial curved routing with ten source-to-destination polylines
  in separate lanes and one direct tissue-to-community E1a connection.
- Put labels beside complete strokes, without white masks, cut-outs or
  arrowheads. Preserve separate A1/F2, A4/H-M1, H-M3/E9 and both E10 paths.
- Align destination labels and reduce box area and grey fill; retain all
  eight nodes, including isolated cell composition.
- Expand visible shorthand into accessible node names, define
  Helicobacter pylori / H. pylori, and spell out M1/M2/M3 programmes.
- Use a modest A, a compact scientific key and Okabe–Ito vermillion/blue
  signed colours. Numeric signs and line patterns also work in greyscale.
- Replace misleading 'Confirmatory passes' with primary-test frozen-pass-rule
  wording. The manuscript caption qualifies A1, F2 and the unsupported
  HKSJ-max1 M1 P = 0.159 result; solid does not assert universal robustness.
- Remove the internal audit footer, script paths, provisional dimensions and
  publisher-pending wording from the scientific graphic. Audit material is
  recorded here and in `qa/`, with the full title/legend in the caption file.

After the first-render inspection, the H-M3 lane was raised by 2 mm to give
the E9 annotation more breathing room, the community/M1 box gap was opened
from 1 to 3 mm, and the tissue box was trimmed by 2 mm. E1a endpoints and its
label were repositioned to match. No locked number or width was changed.
The first preview and geometry are retained in `qa/first_render_*`.

## Scientific lock and width provenance

Five primary-test edges are solid: A1, F2, H-M1, H-M2 and H-M3. A4 remains a
distinct dash-dot locally pre-specified secondary test, with prior interim
inspection disclosed. E1a, E1c, E10 IM, E10 atrophy and E9 remain dashed and
exploratory. A1 is non-blind locally pre-specified consistency, not external
preregistration; F2 is the fully blind gradient validation. H-M2/M3 retain
their disclosed Korea dominance. The HKSJ-max1 M1 result (P = 0.159) is
unsupported and is not implied to have passed by the solid H-M1 edge.

Width **thresholds** come from addendum v0.6d:88–89 and freeze_log.md:68.
The numeric **point widths** come from the preserved original generator
`scripts/s10_e11_schematic.py:39`, as recorded in
`inputs/width_bin_provenance.md`:

```python
def tier(x): x=abs(x); return 1.0 if x<0.2 else (2.5 if x<=0.4 else 4.5)
```

Thus |e| < 0.2 → 1.0 pt; 0.2 ≤ |e| ≤ 0.4 → 2.5 pt; |e| > 0.4 → 4.5 pt.
The point values are not attributed to v0.6d itself. There is no continuous
or gradient width. Effect metrics differ, so width is a frozen descriptive
encoding, not a common inferential-strength or biological-importance scale.
E9 is the maximum |ρ| among qualifying E9 genes (A4GNT), retaining its signed
stored estimate; it is not a unique or causal single-gene effect. The stored
gene set description remains in the source CSV. E1c explicitly concerns
ΔH. pylori–ΔGemella. No Spearman correlations are recomputed.

The tissue-state node retains the pre-result freeze-log resolution of the
v0.6d omitted-node ambiguity: it is needed for A1/F2. Cell composition is
retained with no edge because E7a tissue-beta change is not an association
test. No additional node or edge is introduced.

## E1a: two-cohort display and retained unsupported result

The freeze-log rule is **stored P < 0.05 in at least two cohorts**. The
display uses only F2 SILVA and Korea 16S SILVA, each with stored P = 0.0001:

| Cohort | Stored R² | Stored P | In displayed aggregation |
|---|---:|---:|---|
| F2 SILVA | 0.00546836026298805 | 0.0001 | Yes |
| Korea 16S SILVA | 0.0047785625776584802 | 0.0001 | Yes |
| PRJNA310127 SILVA | 0.0084335046317319494 | 0.071400000000000005 | No |

PRJNA310127 **P = 0.0714** fails the display rule. It is retained in
`outputs/e1a_display_aggregation.csv` and this technical note, and is omitted
only from the displayed summary. The descriptive displayed value is
**√(median stored R²) = 0.07157835860316486** (median R² =
0.0051234614203232651), with two cohorts and the 1.0 pt bin. Its plus sign
marks a non-negative magnitude; PERMANOVA R² is not a directional correlation.
No new test is performed. The historical P < 0.10 three-cohort departure
(0.073948362138644086) stays in its audit row and is not plotted. Both source
CSV copies preserve every original byte and full-precision numeric string,
including restored and historical summary rows.

## Reproduction and audit boundaries

```bash
python3 -B scripts/plot_figure8a_v3_presentation.py
```

New script: `scripts/plot_figure8a_v3_presentation.py` in this v3 directory.
It writes only to this revision's `outputs/` and `qa/`. Inputs are read only.
The old 174 × 190 mm provisional size, pending-specifications wording and
script-path footer belong to the v2 audit history, not the v3 graphic.

Verification is in `outputs/verification.json`; actual drawing geometry and
text bounds are in `qa/render_manifest.json`. Exported SVG path styles are
checked independently of the manifest. The six paths listed in
`inputs/v2_candidate_hashes_before_v3.txt` are re-read for SHA-256 comparison
before and after rendering. That limited check is not a global integrity
claim. The four historical binary hash-exception files and the unrelated
worker directories are neither read nor written by the script.
"""
    (OUTPUTS / "technical_note.md").write_text(note, encoding="utf-8")


def poppler_report(command, filename):
    tool = shutil.which(command[0])
    if not tool:
        return {"available": False, "output": "", "returncode": None}
    result = subprocess.run([tool, *command[1:]], capture_output=True, text=True, check=False)
    output = result.stdout + result.stderr
    (QA / filename).write_text(output, encoding="utf-8")
    return {"available": True, "returncode": result.returncode, "output": output}


def verify(manifest, baseline_before):
    rows = read_csv(CSV_NAMES[0])
    aggregate = read_csv(CSV_NAMES[1])
    svg = ET.parse(OUTPUTS / f"{BASE}.svg").getroot()
    ns = {"s": "http://www.w3.org/2000/svg"}
    checks, details = {}, {}
    def check(name, passed, detail=None):
        checks[name] = bool(passed)
        if detail is not None:
            details[name] = detail
    check("eleven_edges", len(rows) == len(manifest["edges"]) == 11)
    check("association_identities", set(EXPECTED_PAIRS) == {r["analysis_id"] for r in rows}
          and all(set(EXPECTED_PAIRS[e["analysis_id"]]) == {e["node_a"], e["node_b"]}
                  for e in manifest["edges"]))
    check("eight_nodes", set(NODES) == {r["node_id"] for r in read_csv("nodes_and_noedge.csv")})
    check("comp_has_no_edge", all("COMP" not in (r["node_a"], r["node_b"]) for r in rows))
    edge_groups = [g for g in svg.findall(".//s:g", ns) if g.attrib.get("id", "").startswith("association-edge-")]
    check("eleven_exported_SVG_edges", len(edge_groups) == 11)
    effects, styles = [], []
    continuous_paths = True
    for edge in manifest["edges"]:
        row = next(r for r in rows if r["analysis_id"] == edge["analysis_id"])
        error = abs(Decimal(str(edge["effect_used"])) - Decimal(row["effect"]))
        effects.append({"analysis_id": edge["analysis_id"], "effect_csv": row["effect"],
                        "effect_used": edge["effect_used"], "absolute_error": float(error),
                        "sign_matches": (edge["effect_used"] > 0) == (Decimal(row["effect"]) > 0)})
        group = svg.find(f".//s:g[@id='{edge['svg_id']}']", ns)
        path = group.find("s:path", ns)
        continuous_paths = continuous_paths and len(group.findall("s:path", ns)) == 1 and path.attrib["d"].count("M") == 1
        attributes = dict(item.strip().split(":", 1) for item in path.attrib["style"].split(";") if ":" in item)
        attributes = {k: v.strip() for k, v in attributes.items()}
        width = float(attributes.get("stroke-width", 1.0))
        dash = attributes.get("stroke-dasharray")
        style = "solid" if not dash else ("dashdot" if len(dash.split(",")) == 4 else "dashed")
        styles.append({"analysis_id": edge["analysis_id"], "svg_style": style,
                       "svg_linewidth_pt": width, "expected_style": row["style"],
                       "width_bin": row["width_bin"], "stroke_dasharray": dash,
                       "stroke_colour": attributes.get("stroke")})
    check("effects_within_1e_12", all(e["absolute_error"] <= 1e-12 for e in effects), effects)
    check("uninterrupted_edge_paths", continuous_paths)
    check("signs_unchanged", all(e["sign_matches"] for e in effects))
    check("five_solid_A4_dashdot_five_dashed", all(s["svg_style"] == s["expected_style"] for s in styles)
          and {s["analysis_id"] for s in styles if s["svg_style"] == "solid"} == EXPECTED_SOLID
          and [s["analysis_id"] for s in styles if s["svg_style"] == "dashdot"] == ["A4"]
          and sum(s["svg_style"] == "dashed" for s in styles) == 5, styles)
    check("width_bins_original_generator", all(
        width_tier(float(row["effect"])) == (s["width_bin"], s["svg_linewidth_pt"])
        and s["svg_linewidth_pt"] == float(row["linewidth_pt"])
        for row, s in zip(rows, styles)))
    check("width_bin_boundaries", [width_tier(x)[1] for x in (0.1999999, 0.2, 0.4, 0.4000001)]
          == [1.0, 2.5, 2.5, 4.5])
    check("signed_rounded_labels", all(signed(e["effect_used"]) in e["label_text"] for e in manifest["edges"]))
    check("signed_colours", all(s["stroke_colour"].lower() ==
          (POSITIVE if float(r["effect"]) > 0 else NEGATIVE).lower() for r, s in zip(rows, styles)))
    cohort_rows = [r for r in aggregate if not r["cohort_key"].startswith("SUMMARY_")]
    eligible = [r for r in cohort_rows if r["eligible_for_display_aggregation"] == "yes"]
    check("E1a_F2_Korea_only", {r["cohort_key"] for r in eligible} == {"f2_silva", "korea16s_silva"}
          and len(eligible) == 2 and all(Decimal(r["P"]) < Decimal("0.05") for r in eligible),
          {"cohorts": [r["cohort_label"] for r in eligible], "stored_P_rule": "P < 0.05"})
    # Descriptive arithmetic check of supplied estimates only, not inference.
    e1a = next(e for e in manifest["edges"] if e["analysis_id"] == "E1a")
    median_r2 = sum(Decimal(r["R2"]) for r in eligible) / 2
    check("E1a_locked_effect", abs(Decimal(str(e1a["effect_used"])) - Decimal("0.07157835860316486")) <= Decimal("1e-12")
          and abs(median_r2.sqrt() - Decimal(str(e1a["effect_used"]))) <= Decimal("1e-12"),
          {"median_stored_R2": str(median_r2), "descriptive_sqrt": str(median_r2.sqrt())})
    excluded = next(r for r in cohort_rows if r["cohort_key"] == "yu310127_silva")
    check("PRJNA310127_retained_not_displayed", excluded["eligible_for_display_aggregation"] == "no"
          and abs(Decimal(excluded["P"]) - Decimal("0.0714")) < Decimal("1e-15")
          and "0.0714" in (OUTPUTS / "technical_note.md").read_text(), excluded)
    e10 = [e for e in manifest["edges"] if e["analysis_id"].startswith("E10")]
    check("two_distinct_E10_edges", {e["analysis_id"] for e in e10} == {"E10 IM", "E10 atrophy"}
          and e10[0]["route_mm"] != e10[1]["route_mm"])
    check("no_arrowheads", not svg.findall(".//s:marker", ns)
          and "marker-end" not in (OUTPUTS / f"{BASE}.svg").read_text()
          and all(e["arrowheads"] is False for e in manifest["edges"]))
    check("editable_svg_text", len(svg.findall(".//s:text", ns)) >= len(manifest["texts"])
          and matplotlib.rcParams["svg.fonttype"] == "none")
    check("no_SVG_raster_images", not svg.findall(".//s:image", ns))
    check("byte_identical_source_tables", all((INPUTS / n).read_bytes() == (OUTPUTS / n).read_bytes() for n in CSV_NAMES))
    png_details = {}
    for name in (f"{BASE}.png", "Figure8A_v3_preview.png", "Figure8A_v3_physical_size_proof.png"):
        with Image.open(OUTPUTS / name) as im:
            dpi = im.info.get("dpi", (0, 0))
            png_details[name] = {"pixels": list(im.size), "dpi": list(dpi),
                                 "physical_size_mm": [im.size[i] / dpi[i] * 25.4 for i in (0, 1)]}
    prod = png_details[f"{BASE}.png"]
    check("PNG_dpi_at_least_300", min(prod["dpi"]) >= 300, png_details)
    check("PNG_declared_170_mm", abs(prod["physical_size_mm"][0] - WIDTH_MM) < 0.001)
    check("SVG_declared_170_mm", abs(float(svg.attrib["width"].replace("pt", "")) * 25.4 / 72 - WIDTH_MM) < 1e-5)
    check("fonts_configured_42", matplotlib.rcParams["pdf.fonttype"] == matplotlib.rcParams["ps.fonttype"] == 42)
    pdf = str(OUTPUTS / f"{BASE}.pdf")
    fonts = poppler_report(["pdffonts", pdf], "pdf_fonts.txt")
    font_rows = fonts["output"].splitlines()[2:]
    check("PDF_fonts_embedded", fonts["returncode"] == 0 and bool(font_rows)
          and all("TrueType" in line and line.split()[-5] == "yes" for line in font_rows), font_rows)
    images = poppler_report(["pdfimages", "-list", pdf], "pdf_images.txt")
    check("PDF_vector_no_images", images["returncode"] == 0 and len(images["output"].strip().splitlines()) == 2)
    info = poppler_report(["pdfinfo", pdf], "pdf_info.txt")
    size = re.search(r"Page size:\s+([\d.]+) x ([\d.]+) pts", info["output"])
    check("PDF_declared_170_by_160_mm", bool(size) and
          abs(float(size[1]) * 25.4 / 72 - WIDTH_MM) < 0.001 and
          abs(float(size[2]) * 25.4 / 72 - HEIGHT_MM) < 0.001)
    poppler_report(["pdftotext", "-layout", pdf, str(QA / "pdf_text.txt")], "pdftotext_log.txt")
    check("on_figure_text_minimum_6_8_pt", min(t["font_size_pt"] for t in manifest["texts"]) >= 6.8)
    check("no_label_masks", not any(t["has_background_mask"] for t in manifest["texts"]))
    check("no_audit_footer_on_figure", not any(re.search(r"publisher|provisional|scripts/|confirmatory passes", t["text"], re.I)
                                               for t in manifest["texts"]))
    title_count, legend_count = len(TITLE.split()), len(LEGEND.split())
    check("caption_length_limits", title_count <= 15 and legend_count <= 300,
          {"title_words": title_count, "legend_words": legend_count})
    check("geometry_clear", manifest["layout_checks"]["all_clear"], manifest["layout_checks"])
    baseline_after = v2_hashes()
    check("v2_hashes_match_before_and_after", all(r["matches"] for r in baseline_before + baseline_after),
          {"before": baseline_before, "after": baseline_after, "scope": "six supplied baseline paths only"})
    verification = {"all_passed": all(checks.values()), "checks": checks, "details": details,
                    "declared_check_size_mm": [WIDTH_MM, HEIGHT_MM],
                    "journal_compliance_claim": False, "no_new_statistical_inference": True,
                    "residual_issues": [] if all(checks.values()) else [k for k, v in checks.items() if not v]}
    (OUTPUTS / "verification.json").write_text(json.dumps(verification, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return verification


def main():
    before = v2_hashes()
    if not all(r["matches"] for r in before):
        raise RuntimeError("A supplied v2 baseline hash does not match before rendering.")
    associations = read_csv(CSV_NAMES[0])
    assert len(associations) == 11
    assert {r["analysis_id"] for r in associations} == set(ROUTES)
    for row in associations:
        assert width_tier(float(row["effect"])) == (row["width_bin"], float(row["linewidth_pt"]))
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=PREVIEW_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    manifest = draw_panel(fig, ax, associations)
    manifest.update({"declared_check_size_mm": [WIDTH_MM, HEIGHT_MM],
                     "journal_compliance_claim": False, "png_dpi": PNG_DPI,
                     "preview_dpi": PREVIEW_DPI, "width_source": "original scripts/s10_e11_schematic.py:39",
                     "input_sha256": {p.name: sha256(p) for p in INPUTS.iterdir() if p.is_file()}})
    manifest["layout_checks"] = layout_checks(manifest)
    description = ("Eleven undirected stored associations across eight nodes. "
                   "Five solid primary tests meeting frozen pass rules, A4 dash-dot local secondary, "
                   "five exploratory dashed associations. Solid does not imply universal robustness. "
                   "E1a: F2 and Korea only, stored P < 0.05; descriptive effect 0.07157835860316486. "
                   "E9: maximum absolute rho among qualifying E9 genes (A4GNT), retaining its sign. "
                   "No arrows; associations do not establish causation.")
    fig.savefig(OUTPUTS / f"{BASE}.svg", metadata={"Title": "Figure 8A", "Description": description,
                 "Creator": "matplotlib; locked prepared CSVs", "Date": "2026-11-13"})
    fig.savefig(OUTPUTS / f"{BASE}.pdf", metadata={"Title": "Figure 8A", "Subject": description,
                 "Creator": "matplotlib; locked prepared CSVs", "CreationDate": None, "ModDate": None})
    fig.savefig(OUTPUTS / f"{BASE}.png", dpi=PNG_DPI,
                metadata={"Description": description, "Declared check size": "170 x 160 mm"})
    fig.savefig(OUTPUTS / "Figure8A_v3_preview.png", dpi=PREVIEW_DPI,
                metadata={"Description": "Same 170 x 160 mm figure at 190.5 dpi; presentation preview."})
    plt.close(fig)

    # Re-render the scientific panel into an exactly sized sub-region of a
    # separate proof canvas. Ruler and annotations never enter main exports.
    proof_w, proof_h = 190.0, 188.0
    proof = plt.figure(figsize=(proof_w / 25.4, proof_h / 25.4), dpi=PREVIEW_DPI)
    proof_ax = proof.add_axes([10 / proof_w, 8 / proof_h, WIDTH_MM / proof_w, HEIGHT_MM / proof_h])
    draw_panel(proof, proof_ax, associations)
    ruler = proof.add_axes([0, 0, 1, 1], zorder=10, frameon=False)
    ruler.set(xlim=(0, proof_w), ylim=(0, proof_h))
    ruler.set_axis_off()
    ruler.text(95, 184.5, "Physical-size proof · 170 mm figure width · print at 100%", fontsize=8,
               ha="center", va="center", color=INK)
    ruler.plot([10, 180], [173, 173], color=INK, linewidth=0.65)
    for mm in range(0, 171, 10):
        ruler.plot([10+mm, 10+mm], [171.5, 174.5], color=INK, linewidth=0.65)
        ruler.text(10+mm, 178, str(mm), fontsize=6.8, ha="center", va="center", color=INK)
    ruler.text(185, 173, "mm", fontsize=7, ha="center", va="center", color=INK)
    proof.savefig(OUTPUTS / "Figure8A_v3_physical_size_proof.png", dpi=PREVIEW_DPI,
                  metadata={"Description": "190 x 188 mm proof canvas, 190.5 dpi; 170 x 160 mm panel at x=10, y=8 mm. Print 100%, no fit-to-page."})
    plt.close(proof)

    for name in CSV_NAMES:
        shutil.copyfile(INPUTS / name, OUTPUTS / name)
    write_notes()
    (QA / "render_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    verification = verify(manifest, before)
    print(f"Figure: {WIDTH_MM:g} × {HEIGHT_MM:g} mm; 11 undirected associations; 8 nodes.")
    print(f"PNG: {PNG_DPI} dpi; preview/proof: {PREVIEW_DPI} dpi. Outputs: {OUTPUTS}")
    print(f"Verification: {sum(verification['checks'].values())}/{len(verification['checks'])} checks passed.")
    if not verification["all_passed"]:
        print("Needs attention: " + ", ".join(verification["residual_issues"]))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
