#!/bin/bash
set -euo pipefail

# ------------ Activate env ------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ~/Dropbox/SEPAL_AI/data/env

# ------------ Paths ------------
BASE=~/Dropbox/SEPAL_AI/data/0_pericytes_prjna1057204
mkdir -p ${BASE}/{raw,fastq,qc,aligned,counts,logs,hisat2_index,tmp}
export TMPDIR=${BASE}/tmp
cd ${BASE}/raw

# ------------ Runs (single-end) ------------
RUNS=(
SRR27352020 SRR27352021 SRR27352022 SRR27352023 SRR27352024 SRR27352025
SRR27352026 SRR27352027 SRR27352028 SRR27352029 SRR27352030 SRR27352031
SRR27352032 SRR27352033 SRR27352034 SRR27352035 SRR27352036 SRR27352037
)

# ------------ Subset settings ------------
READS=1000000   # reads to keep per sample fastq
SEED=100        # for reproducibility

# ------------ Download + stream-subset (NO full FASTQ on disk) ------------
for r in "${RUNS[@]}"; do
  echo ">>> Prefetch $r"
  prefetch "$r" -O .

  echo ">>> Stream subset $r -> fastq/${r}.subset.fastq.gz"
  fasterq-dump -e 8 -p --split-spot --disk-limit 1TB --stdout -t "${TMPDIR}" "$r" \
    | seqtk sample -s${SEED} - ${READS} \
    | pigz > ${BASE}/fastq/${r}.subset.fastq.gz
done

# ------------ QC ------------
cd ${BASE}/fastq
find . -name '*.subset.fastq.gz' -size -1M -delete || true
fastqc *.subset.fastq.gz -o ${BASE}/qc --threads 8
cd ${BASE}/qc
multiqc . -o .

# ------------ HISAT2 index (GRCh38) ------------
cd ${BASE}/hisat2_index
if [ ! -f genome.1.ht2 ]; then
  echo ">>> Downloading HISAT2 GRCh38 index"
  curl -O https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
  tar -xzf grch38_genome.tar.gz
fi

# ------------ Alignment (single-end) ------------
cd ${BASE}
SAMPLES=($(ls fastq/*.subset.fastq.gz | xargs -n1 basename | sed 's/.subset.fastq.gz//'))

for s in "${SAMPLES[@]}"; do
  echo ">>> Align ${s}"
  hisat2 -p 4 -x hisat2_index/genome \
    -U fastq/${s}.subset.fastq.gz \
    2> logs/${s}_hisat2.log \
    | samtools sort -@ 4 -o aligned/${s}.bam
  samtools index aligned/${s}.bam
done

# ------------ Annotation ------------
cd ${BASE}
if [ ! -f gencode.v43.annotation.gtf ]; then
  echo ">>> Downloading GENCODE v43 annotation"
  curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
  gunzip -f gencode.v43.annotation.gtf.gz
fi

# ------------ featureCounts -> GEO-style matrix ------------
featureCounts -T 8 -t exon -g gene_id \
  -a gencode.v43.annotation.gtf \
  -o counts/tmp_counts.txt aligned/*.bam \
  &> logs/featureCounts.log

tail -n +3 counts/tmp_counts.txt | cut -f1,7- > counts/GSE252046_raw_counts_GRCh38.p13_NCBI.tsv
gzip -f counts/GSE252046_raw_counts_GRCh38.p13_NCBI.tsv

echo "✅ Done: ${BASE}/counts/GSE252046_raw_counts_GRCh38.p13_NCBI.tsv.gz"
