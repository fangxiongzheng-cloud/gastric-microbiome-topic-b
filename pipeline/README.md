# R1 host-RNA HP transcript-read quantification: candidate source

This workflow is included for local author review. It has not been released or smoke-tested here. Read the package README and RUN_ORDER before any future use. Four historical result-file hash mismatches remain under investigation, and final updated Figure 8 source is pending.

`hp_hostrna_quant.sh` reads an explicitly selected configuration, processes paired RNA reads with SRA Toolkit/fastp, quantifies host fragments with Salmon, aligns HP candidate reads with Bowtie2, exports paired reads with SAMtools and classifies them with Kraken2. `summarize_hp_load.py` reads saved per-run products and writes a caller-selected table. These source files were not run. Their retained future download commands do not imply downloads in this job.

The metric is sequencing-depth-normalised HP transcript-read load (confirmed clade fragments / host mapped fragments × 10^6). It is a relative transcript-read signal and does not establish viable infection or absolute bacterial abundance. Source documentation uses log1p load and detection at at least 10 confirmed fragments. Library selection, contamination, reference redundancy and classifier database coverage limit interpretation.

`environment.yml` is the unchanged documentary project specification; see `../environments/README.md` for observed metadata, version conflicts and missing dependency closure. Documented references are GENCODE v46, the 12 accession-version HP panel listed in the parent inventory, and Kraken2 Standard-8GB dated 2024-09-04. References/indices are not bundled and their bytes/checksums were not verified.

`config.example.sh` is a template with no credentials. It points to a header-only `../config/runlist.template.tsv` by work-root-relative name `config/runlist.template.tsv`; authors must supply rights-cleared run rows separately. No public-run example rows or tissue labels were copied. All paths and any future output/scratch locations must be reviewed in a separate private work copy. No command is intended to run directly from this package.

The virulence alignment is optional in this standalone script. Historical `core_qc.py` expects `vir.bam`, so its input contract is not satisfied automatically when that option is omitted. DONE markers alone do not demonstrate the full historical QC. Calibration result summaries and their report writer were excluded; E4/E8 analytical source and expected path contracts are documented separately.
