# Gastric host–microbiome analysis — private code-release candidate

This repository is a **private, unpublished code candidate**. It is **not publicly released**, has no DOI, and is not a public data or code deposit. A license has **not** been added and no copyright statement has been assigned in this upload. Third-party packages remain under their own terms; no dependency code or data are bundled.

## What this is, and is not

- Selected analysis scripts, pipeline templates, a header-only runlist template, documentary/observed partial environment specifications, and Figure 8A v3 **plot source and notes**.
- **Not** a turnkey reproduction package. **No full end-to-end rerun**, installation, statistical recomputation, data download, or figure regeneration was performed for this upload. Isolated syntax/parse checks were applied only to this curated copy where noted.
- **Not** a release of raw sequencing data, results, official assembled figures, sample-level metadata, or sample titles.

Use `RUN_ORDER.md` for the documentation-derived run graph and `IO_ACCESSION_MANIFEST.tsv` for expected inputs, outputs and public accessions. Obtain data from the official repositories under their terms. **TCMA CC0 is a data license, not an ethics waiver.**

## Figure 8

Official assembled Figure 8 assets were replaced with the v3 presentation panel on 2026-11-13 21:15 Asia/Shanghai (UTC+8), per the project receipt `manuscript/qa/figure8_v3_install/canonical_revision_manifest.md`. This repository includes the plot script and technical notes under `figure8_v3_presentation/`. It does **not** include official assembled image files or the derived source CSVs (redistribution rights unresolved). Journal-size or typesetting compliance is **not** claimed.

## Historical result-file hashes

Four large result files were **not** demonstrated hash mismatches. In the reviewed baselines their SHA-256 values were empty and marked `skipped: too_large` (files above the hasher size cutoff). A separate read-only inquiry found that **441 recorded hashes all matched**, including manuscript-used downstream tables included in that recheck. Historical byte identity of the four unhashed files remains **unproven**; R-object internals of one file were not read. This is not labelled as known corruption, and it is not a claim that all source/result integrity has passed. Those result files are not in this repository and were not rebuilt.

## Coverage and environment gaps

- Withheld: the upstream-derived BayesPrism driver (`scripts/s10_e7_bayesprism.R`) pending attribution/redistribution review; figure source CSVs; any LICENSE or copyright file; raw data; results; sample titles.
- Incomplete or missing retained producers, as previously inventoried: corrected A4 ssGSEA integration, selected E8 background-list producer, TCGA A4 normal-score producer. Many main-panel scripts still live under the project `results/` tree and were not relocated here. Older included draft plotters are not asserted to be every final manuscript figure source.
- Environment files are documentary or partial local snapshots. Verified pins mean a project specification or selected installed metadata was read; they do not recover a complete historical, solved environment. Methods documentary versions and current metadata can differ (including Python 3.11.18 vs current shotgun Python 3.14.10, and SpiecEasi 1.1.1/1.1.2 records). R 4.5.3 (rbio) and R 4.4.3 (a5sc) must not be merged into one asserted installation. EPIC’s installed license is limited to academic non-commercial research and is not bundled.

Machine-specific paths in script copies were replaced with `/REPLACE_WITH_...` placeholders in the candidate only. The original project was not modified by this upload packaging.
