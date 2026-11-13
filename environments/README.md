# Environment specifications for a local candidate

These are documentary specifications and partial local snapshots. No environment was created, updated, solved or tested by running an analysis.

`../pipeline/environment.yml` is a byte-identical copy of the project R1 specification. Its command-line pins are sra-tools 3.4.1, fastp 1.4.0, salmon 2.8.0, bowtie2 2.5.5, samtools 1.25 and kraken2 2.17.2. The corresponding selected local conda metadata contains these versions, distributed across shotgun and salmon environments. The template's Python 3.11 constraint differs from the current shotgun metadata's Python 3.14.10. It is not a validated single combined environment.

`analysis.template.yml` records Methods' Python 3.11.18 statement; `python-observed-requirements.txt` and `singlecell-observed-requirements.txt` record distributions found in the project-associated Python environment. They do not establish the environment used for every historical result. Partial `observed-*.yml` files preserve selected current metadata separately, including the QIIME 2 and A5 environments. Missing tools/versions are not invented.

`VERSION_EVIDENCE.tsv` includes the exact metadata path and its new read-time SHA-256, the declared version/license and evidence scope. R package versions are listed there, separated by environment; the rbio environment has R 4.5.3, while A5's a5sc metadata has R 4.4.3. Do not merge those package catalogs into one asserted working installation. GitHub package revision fields, when present in DESCRIPTION, are current installed-package metadata only.

SpiecEasi's historical installer targets v1.1.1; installed a3net DESCRIPTION reports 1.1.2, and rbio conda metadata reports 1.1.1 while its DESCRIPTION reports 1.1.2. Authors must reconcile these records before deciding reproducibility pins. BayesPrism's copied-function driver is withheld pending notice/redistribution review. EPIC's installed license limits use to academic non-commercial research and restricts redistribution; no EPIC code is bundled here.

Databases, reference sequences, classifier artifacts, all transitive dependencies, OS utilities, BLAS/compiler versions and binary builds remain outside these partial specifications. Their bytes and checksums were not verified. No download or installation command was executed.
