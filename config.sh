# config.sh
# Paths and sample names shared by every script. Edit PROJECT_DIR (and the
# reference file names if yours differ), then each job script sources this file.
#
# The analysis was run on the University of Birmingham BlueBEAR cluster (Slurm),
# so the job scripts load software with BlueBEAR's "module load" commands.

PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"

# ---- reference files (GRCh38 / hg38) ----------------------------------------
REF_DIR="$PROJECT_DIR/reference"
REF_FASTA="$REF_DIR/GRCh38.p14.genome.fa"
BOWTIE2_INDEX="$REF_DIR/bowtie2/GRCh38"          # prefix written by bowtie2-build
STAR_INDEX="$REF_DIR/star_index"
GTF="$REF_DIR/Homo_sapiens.GRCh38.115.gtf.gz"     # Ensembl gene models
BLACKLIST="$REF_DIR/hg38-blacklist.v2.bed"        # ENCODE blacklist (Amemiya et al. 2019)
CHROM_SIZES="$REF_DIR/hg38.chrom.sizes"

# ---- DNaseI-seq (single-end) -------------------------------------------------
DNASE_DIR="$PROJECT_DIR/dnase"
DNASE_FASTQ="$DNASE_DIR/fastq"          # <sample>.fastq.gz
DNASE_TRIMMED="$DNASE_DIR/trimmed"
DNASE_BAM="$DNASE_DIR/bam"
DNASE_PEAKS="$DNASE_DIR/peaks"
DNASE_CONSENSUS="$DNASE_DIR/consensus"
DNASE_DA="$DNASE_DIR/differential_accessibility"
FOOTPRINT_DIR="$DNASE_DIR/footprinting"

# Two healthy peripheral blood stem cell (PBSC) controls and five t(8;21) patients.
# (The sixth patient, t821-ITD, had RNA-seq only.)
DNASE_CONTROLS=(PBSC-1 PBSC-2)
DNASE_PATIENTS=(t821-1 t821-1rel t821-3 t821-4 t821-KIT)
DNASE_SAMPLES=("${DNASE_CONTROLS[@]}" "${DNASE_PATIENTS[@]}")

# ---- public ChIP-seq (GEO GSE137670, Nafria et al. 2020) --------------------
CHIP_DIR="$PROJECT_DIR/chipseq"         # *_treat_pileup.bedgraph downloaded from GEO
CHIP_BW="$CHIP_DIR/bigwig"
CHIP_PEAKS="$CHIP_DIR/peaks"
CHIP_PLOTS="$CHIP_DIR/plots"
CHIP_INPUT="$CHIP_DIR/GSM4083965_Input_treat_pileup.bedgraph"

# ---- RNA-seq (paired-end) ----------------------------------------------------
RNA_DIR="$PROJECT_DIR/rnaseq"
RNA_FASTQ="$RNA_DIR/fastq"              # <sample>_R1.fastq.gz and <sample>_R2.fastq.gz
RNA_TRIMMED="$RNA_DIR/trimmed"
RNA_BAM="$RNA_DIR/bam"
RNA_COUNTS="$RNA_DIR/counts"
RNA_SAMPLES=(PBSC-1 PBSC-2 t821-1 t821-1R t821-3 t821-4 t821-ITDx2-5 t821-KIT-2)
