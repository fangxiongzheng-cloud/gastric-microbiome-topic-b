# Candidate example configuration. No sample rows are included; supply a rights-cleared runlist.
RUNLIST=config/runlist.template.tsv          # tab-separated, first column = SRA run accession (header 'run' allowed)
OUTDIR=results/hp_hostrna/per_run
WORKDIR=/REPLACE_WITH_PRIVATE_SCRATCH/hp_hostrna           # needs ~3x the largest run's .sra size free
SALMON_INDEX=refs/gencode/salmon_v46      # salmon index --gencode on gencode.v46.transcripts.fa.gz
HP_BT2_INDEX=refs/hp/hp12                 # bowtie2-build of 12 complete H. pylori genomes (README)
KRAKEN2_DB=db/k2_std8                     # Kraken2 Standard-8GB 2024-09-04
VIR_BT2_INDEX=                            # optional: cagA/vacA/babA index (not used for the load metric)
THREADS=8; K2_THREADS=4; MIN_FREE_GB=15
