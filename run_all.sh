#!/bin/bash
# RNA-Seq Workbench - Production Pipeline
# End-to-end transcriptomics and epigenomics analysis framework
# Author: Bioinformatics Core Facility
# Version: 2.0.1
# Date: $(date +%Y-%m-%d)

set -euo pipefail
IFS=$'\n\t'

# Pipeline configuration
THREADS=20
PIPELINE_VERSION="2.0.1"
BASE=$(dirname "$(realpath "$0")")
LOG_DIR="$BASE/logs"
RESULTS_DIR="$BASE/results"

# Initialize logging infrastructure
mkdir -p "$LOG_DIR" "$RESULTS_DIR"
exec > >(tee -a "$LOG_DIR/workbench_$(date +%Y%m%d_%H%M%S).log")
exec 2>&1

# Logging functions
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $1" >&2; }
err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1" >&2; exit 1; }

log "RNA-Seq Workbench v$PIPELINE_VERSION initializing..."
log "Configuration: $THREADS threads, Base directory: $BASE"

# ─── QC Module: Data Quality Assurance ─────────────────────────────────────
log "[QC Module] Initializing data quality assessment"
mkdir -p "$BASE/modules/qc" "$BASE/data"

# Data acquisition from ENCODE repository
log "[QC Module] Acquiring ENCODE RNA-seq datasets (K562/MCF-7)"
wget -q -P "$BASE/data" \
  https://www.encodeproject.org/files/ENCFF002DKA/@@download/ENCFF002DKA.fastq.gz \
  https://www.encodeproject.org/files/ENCFF002DKE/@@download/ENCFF002DKE.fastq.gz

# Quality control analysis
log "[QC Module] Running FastQC quality assessment"
fastqc "$BASE/data"/*.fastq.gz -o "$BASE/modules/qc" -t $THREADS

# Aggregate QC reports
log "[QC Module] Generating MultiQC dashboard"
multiqc "$BASE/modules/qc" -o "$BASE/modules/qc" -n "multiqc_report.html"

log "[QC Module] Quality assessment complete"

# ─── Alignment Module: Reference Genome Mapping ─────────────────────────────
log "[Alignment Module] Starting reference genome alignment pipeline"
mkdir -p "$BASE/modules/alignment" "$BASE/reference"

# Reference genome acquisition and indexing
log "[Alignment Module] Downloading hg19 Bowtie2 index"
wget -q -P "$BASE/reference" https://genome-idx.s3.amazonaws.com/bt/hg19.zip
unzip -q -d "$BASE/reference" "$BASE/reference/hg19.zip"

# High-performance alignment
log "[Alignment Module] Executing Bowtie2 alignment ($THREADS threads)"
start_time=$(date +%s)
bowtie2 -p $THREADS -x "$BASE/reference/hg19" \
  -1 "$BASE/data/ENCFF002DKA.fastq.gz" \
  -2 "$BASE/data/ENCFF002DKE.fastq.gz" \
  -S "$BASE/modules/alignment/aligned_reads.sam"
end_time=$(date +%s)
log "[Alignment Module] Alignment completed in $((end_time - start_time)) seconds"

# Post-processing pipeline
log "[Alignment Module] SAM→BAM conversion and optimization"
samtools view -@ $THREADS -bS "$BASE/modules/alignment/aligned_reads.sam" \
  > "$BASE/modules/alignment/aligned_reads.bam"
samtools sort -@ $THREADS "$BASE/modules/alignment/aligned_reads.bam" \
  -o "$BASE/modules/alignment/aligned_reads_sorted.bam"
samtools index "$BASE/modules/alignment/aligned_reads_sorted.bam"

# Generate alignment statistics
log "[Alignment Module] Computing alignment metrics"
samtools flagstat "$BASE/modules/alignment/aligned_reads_sorted.bam" \
  > "$BASE/modules/alignment/alignment_statistics.txt"

# Cleanup intermediate files
rm "$BASE/modules/alignment/aligned_reads.sam" "$BASE/modules/alignment/aligned_reads.bam"
log "[Alignment Module] Reference genome mapping complete"

# ─── Variant Module: Genetic Variation Discovery ─────────────────────────────
log "[Variant Module] Initiating SNP/indel discovery pipeline"
mkdir -p "$BASE/modules/variants" "$BASE/reference"

# Reference genome preparation
log "[Variant Module] Downloading hg19 reference genome"
wget -q -P "$BASE/reference" https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip -k "$BASE/reference/hg19.fa.gz"

# Variant calling pipeline
log "[Variant Module] BCFtools mpileup - raw variant detection"
bcftools mpileup -Ov --threads $THREADS \
  -f "$BASE/reference/hg19.fa" \
  "$BASE/modules/alignment/aligned_reads_sorted.bam" \
  -o "$BASE/modules/variants/raw_variants.vcf"

# Variant filtering and annotation
log "[Variant Module] BCFtools call - high-confidence SNP identification"
bcftools call -mv -Ov --threads $THREADS \
  -o "$BASE/modules/variants/called_variants.vcf" \
  "$BASE/modules/variants/raw_variants.vcf"

# Statistical analysis
log "[Variant Module] Generating variant statistics report"
bcftools stats "$BASE/modules/variants/called_variants.vcf" \
  > "$BASE/modules/variants/variant_statistics.txt"

log "[Variant Module] Variant discovery pipeline completed"

# ─── Expression Module: Transcriptomics Analysis ───────────────────────────────
log "[Expression Module] Starting transcript quantification and differential expression"
mkdir -p "$BASE/modules/expression"

# Check for required sorted BAM files from existing data
if [[ ! -f "$BASE/Homework04/MCF-7/UniquelyAligned_sorted.bam" ]] || \
   [[ ! -f "$BASE/Homework04/k562/UniquelyAligned_sorted.bam" ]]; then
  warn "[Expression Module] Required sorted BAM files not found. Skipping expression analysis."
  warn "[Expression Module] Please ensure MCF-7 and K562 sorted BAM files are available."
else
  # Transcript assembly for MCF-7
  log "[Expression Module] Cufflinks - MCF-7 transcript assembly"
  "$BASE/Homework04/tools/cufflinks-2.2.1.Linux_x86_64/cufflinks" \
    -G "$BASE/Homework04/hg19_ncbiRefSeq_annotation.gtf" \
    "$BASE/Homework04/MCF-7/UniquelyAligned_sorted.bam" \
    -o "$BASE/modules/expression/MCF7_transcripts"
  
  # Transcript assembly for K562
  log "[Expression Module] Cufflinks - K562 transcript assembly"
  "$BASE/Homework04/tools/cufflinks-2.2.1.Linux_x86_64/cufflinks" \
    -G "$BASE/Homework04/hg19_ncbiRefSeq_annotation.gtf" \
    "$BASE/Homework04/k562/UniquelyAligned_sorted.bam" \
    -o "$BASE/modules/expression/K562_transcripts"
  
  # Differential expression analysis
  log "[Expression Module] Cuffdiff - comparative expression analysis (K562 vs MCF-7)"
  "$BASE/Homework04/tools/cufflinks-2.2.1.Linux_x86_64/cuffdiff" \
    --num-threads $THREADS \
    -o "$BASE/modules/expression/differential_expression" \
    "$BASE/Homework04/hg19_ncbiRefSeq_annotation.gtf" \
    "$BASE/Homework04/MCF-7/UniquelyAligned_sorted.bam" \
    "$BASE/Homework04/k562/UniquelyAligned_sorted.bam"
  
  log "[Expression Module] Transcriptomics analysis completed"
fi

# ─── Epigenomics Module: ChIP-seq Peak Detection ───────────────────────────────
log "[Epigenomics Module] Initializing POLR2A ChIP-seq analysis pipeline"
mkdir -p "$BASE/modules/epigenomics" "$BASE/reference"

# Data acquisition for ChIP-seq experiments
log "[Epigenomics Module] Acquiring ChIP-seq datasets from ENCODE"
mkdir -p "$BASE/data/chipseq/K562/chip" "$BASE/data/chipseq/K562/control" \
         "$BASE/data/chipseq/MCF7/chip" "$BASE/data/chipseq/MCF7/control"

# K562 datasets
wget -q -P "$BASE/data/chipseq/K562/chip" \
  https://www.encodeproject.org/files/ENCFF839LPL/@@download/ENCFF839LPL.fastq.gz \
  https://www.encodeproject.org/files/ENCFF698ICA/@@download/ENCFF698ICA.fastq.gz
wget -q -P "$BASE/data/chipseq/K562/control" \
  https://www.encodeproject.org/files/ENCFF002EFD/@@download/ENCFF002EFD.fastq.gz \
  https://www.encodeproject.org/files/ENCFF002EFA/@@download/ENCFF002EFA.fastq.gz

# MCF7 datasets  
wget -q -P "$BASE/data/chipseq/MCF7/chip" \
  https://www.encodeproject.org/files/ENCFF000SBI/@@download/ENCFF000SBI.fastq.gz
wget -q -P "$BASE/data/chipseq/MCF7/control" \
  https://www.encodeproject.org/files/ENCFF000SAZ/@@download/ENCFF000SAZ.fastq.gz

# Reference genome indexing (if not already available)
if [[ ! -f "$BASE/reference/hg19" ]]; then
  log "[Epigenomics Module] Setting up hg19 reference index"
  mkdir -p "$BASE/reference/hg19"
  # Use existing index or create new one
  if [[ -f "$BASE/Homework06/index/hg19" ]]; then
    cp -r "$BASE/Homework06/index/hg19"* "$BASE/reference/"
  else
    warn "[Epigenomics Module] Reference index not found. Please ensure hg19 index is available."
  fi
fi

# ChIP-seq alignment pipeline
log "[Epigenomics Module] Executing ChIP-seq alignment for K562"
bowtie2 -p $THREADS -x "$BASE/reference/hg19" \
  -1 "$BASE/data/chipseq/K562/chip/ENCFF839LPL.fastq.gz" \
  -2 "$BASE/data/chipseq/K562/chip/ENCFF698ICA.fastq.gz" \
  -S "$BASE/modules/epigenomics/K562_chip_aligned.sam"
bowtie2 -p $THREADS -x "$BASE/reference/hg19" \
  -1 "$BASE/data/chipseq/K562/control/ENCFF002EFD.fastq.gz" \
  -2 "$BASE/data/chipseq/K562/control/ENCFF002EFA.fastq.gz" \
  -S "$BASE/modules/epigenomics/K562_control_aligned.sam"

log "[Epigenomics Module] Executing ChIP-seq alignment for MCF7"
bowtie2 -p $THREADS -x "$BASE/reference/hg19" \
  -U "$BASE/data/chipseq/MCF7/chip/ENCFF000SBI.fastq.gz" \
  -S "$BASE/modules/epigenomics/MCF7_chip_aligned.sam"
bowtie2 -p $THREADS -x "$BASE/reference/hg19" \
  -U "$BASE/data/chipseq/MCF7/control/ENCFF000SAZ.fastq.gz" \
  -S "$BASE/modules/epigenomics/MCF7_control_aligned.sam"

# BAM processing with quality filtering
log "[Epigenomics Module] Processing aligned reads (MAPQ≥30 filtering)"
for SAMPLE in K562_chip K562_control MCF7_chip MCF7_control; do
  log "[Epigenomics Module] Processing $SAMPLE"
  samtools view -@ $THREADS -q 30 -bS "$BASE/modules/epigenomics/${SAMPLE}_aligned.sam" \
    > "$BASE/modules/epigenomics/${SAMPLE}_filtered.bam"
  samtools sort -@ $THREADS "$BASE/modules/epigenomics/${SAMPLE}_filtered.bam" \
    -o "$BASE/modules/epigenomics/${SAMPLE}_sorted.bam"
  samtools index "$BASE/modules/epigenomics/${SAMPLE}_sorted.bam"
  rm "$BASE/modules/epigenomics/${SAMPLE}_aligned.sam" "$BASE/modules/epigenomics/${SAMPLE}_filtered.bam"
done

# Peak calling with MACS2
log "[Epigenomics Module] MACS2 peak calling - K562 POLR2A"
macs2 callpeak -t "$BASE/modules/epigenomics/K562_chip_sorted.bam" \
  -c "$BASE/modules/epigenomics/K562_control_sorted.bam" \
  -f BAM -g hs -p 0.00001 \
  --outdir "$BASE/modules/epigenomics/K562_peaks" -n K562_POLR2A

log "[Epigenomics Module] MACS2 peak calling - MCF7 POLR2A"
macs2 callpeak -t "$BASE/modules/epigenomics/MCF7_chip_sorted.bam" \
  -c "$BASE/modules/epigenomics/MCF7_control_sorted.bam" \
  -f BAM -g hs -p 0.00001 \
  --outdir "$BASE/modules/epigenomics/MCF7_peaks" -n MCF7_POLR2A

log "[Epigenomics Module] ChIP-seq analysis pipeline completed"

# ─── Motif Module: Regulatory Element Analysis ────────────────────────────────
log "[Motif Module] Starting regulatory motif discovery pipeline"
mkdir -p "$BASE/modules/motifs"

# Check for required peak files and reference genome
PEAK_DIR_K562="$BASE/modules/epigenomics/K562_peaks"
PEAK_DIR_MCF7="$BASE/modules/epigenomics/MCF7_peaks"
REF_GENOME="/workspace/workdisk3/rudhra/database/exp/index/hg19.fa"

if [[ ! -f "$PEAK_DIR_K562/K562_POLR2A_peaks.narrowPeak" ]] || \
   [[ ! -f "$PEAK_DIR_MCF7/MCF7_POLR2A_peaks.narrowPeak" ]]; then
  warn "[Motif Module] Peak files not found. Please ensure ChIP-seq analysis completed successfully."
  warn "[Motif Module] Expected: $PEAK_DIR_K562/K562_POLR2A_peaks.narrowPeak"
  warn "[Motif Module] Expected: $PEAK_DIR_MCF7/MCF7_POLR2A_peaks.narrowPeak"
elif [[ ! -f "$REF_GENOME" ]]; then
  warn "[Motif Module] Reference genome not found at $REF_GENOME"
  warn "[Motif Module] Please update the reference genome path in the script."
else
  # Sequence extraction from peak regions
  log "[Motif Module] Extracting genomic sequences from peak regions - MCF7"
  bedtools getfasta -fi "$REF_GENOME" \
    -bed "$PEAK_DIR_MCF7/MCF7_POLR2A_peaks.narrowPeak" \
    -fo "$BASE/modules/motifs/MCF7_peak_sequences.fasta"
  
  log "[Motif Module] Extracting genomic sequences from peak regions - K562"
  bedtools getfasta -fi "$REF_GENOME" \
    -bed "$PEAK_DIR_K562/K562_POLR2A_peaks.narrowPeak" \
    -fo "$BASE/modules/motifs/K562_peak_sequences.fasta"
  
  # MEME motif discovery
  if command -v meme &> /dev/null; then
    log "[Motif Module] MEME suite - de novo motif discovery (MCF7)"
    meme "$BASE/modules/motifs/MCF7_peak_sequences.fasta" \
      -dna -oc "$BASE/modules/motifs/MCF7_motifs" -nmotifs 5 -minw 6 -maxw 20 -p $THREADS
    
    log "[Motif Module] MEME suite - de novo motif discovery (K562)"
    meme "$BASE/modules/motifs/K562_peak_sequences.fasta" \
      -dna -oc "$BASE/modules/motifs/K562_motifs" -nmotifs 5 -minw 6 -maxw 20 -p $THREADS
    
    log "[Motif Module] Regulatory motif analysis completed"
  else
    warn "[Motif Module] MEME suite not found in PATH. Please install MEME suite."
    warn "[Motif Module] Conda installation: conda install -c bioconda meme"
  fi
fi

# ─── Pipeline Completion and Reporting ───────────────────────────────────────

# Generate final summary report
log "[Reporting] Generating pipeline completion summary"
cat > "$RESULTS_DIR/pipeline_summary.txt" << EOF
RNA-Seq Workbench Pipeline Summary
===================================
Pipeline Version: $PIPELINE_VERSION
Execution Date: $(date '+%Y-%m-%d %H:%M:%S')
Total Runtime: $SECONDS seconds

Modules Completed:
- QC Module: Data quality assessment
- Alignment Module: Reference genome mapping  
- Variant Module: SNP/indel discovery
- Expression Module: Transcript quantification (if data available)
- Epigenomics Module: ChIP-seq peak detection (if data available)
- Motif Module: Regulatory motif discovery (if peaks available)

Output Directories:
- Quality Control: $BASE/modules/qc
- Alignment Results: $BASE/modules/alignment
- Variant Calls: $BASE/modules/variants
- Expression Analysis: $BASE/modules/expression
- Epigenomics: $BASE/modules/epigenomics
- Motif Discovery: $BASE/modules/motifs
- Logs: $LOG_DIR

Next Steps:
1. Review MultiQC report in modules/qc/
2. Examine alignment statistics in modules/alignment/
3. Analyze variant calls in modules/variants/
4. Load IGV session files for visualization
5. Review motif analysis results in modules/motifs/

For technical support, consult the pipeline logs in $LOG_DIR/
EOF

# Final status
log "RNA-Seq Workbench v$PIPELINE_VERSION pipeline execution completed successfully!"
log "Total execution time: $SECONDS seconds"
log "Results available in: $RESULTS_DIR/"
log "Logs available in: $LOG_DIR/"
log ""
log "Thank you for using RNA-Seq Workbench!"
log "For questions or support, contact: bioinformatics-core@institution.edu"
