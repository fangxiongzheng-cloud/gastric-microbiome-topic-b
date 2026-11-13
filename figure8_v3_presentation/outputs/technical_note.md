# Figure 8A v3 presentation candidate: technical note

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

The [official Microbiome submission guidelines](https://link.springer.com/journal/40168/submission-guidelines) were checked
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
