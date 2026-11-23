Here is a complete, production-ready `README.md` tailored specifically to your `WES-ACMG-Pipeline.sh` and `WES-ACMG-Pipeline.config.txt` files. It covers installation (via Conda/Mamba and manual setups), reference bundle preparation, ANNOVAR database downloads, CLI usage, configuration parameters, and output structures.

```markdown
# WES-ACMG-Pipeline

A modular, human-readable Bash pipeline for **Whole-Exome Sequencing (WES) germline short-variant discovery and ACMG/AMP 2015 rule-based classification**.

Designed for small-to-medium single-sample or per-sample clinical research workflows, this pipeline implements **GATK4 Best Practices** for germline short variants using **hard-filtering** (ideal when cohort sizes are too small for Variant Quality Score Recalibration / VQSR) and automatically generates a self-contained HTML summary report.

---

## Table of Contents
1. [Pipeline Workflow](#pipeline-workflow)
2. [Prerequisites & Installation](#prerequisites--installation)
3. [Reference & Database Setup](#reference--database-setup)
4. [Configuration (`WES-ACMG-Pipeline.config.txt`)](#configuration-wes-acmg-pipelineconfigtxt)
5. [Usage & CLI Execution](#usage--cli-execution)
6. [Output Directory Structure](#output-directory-structure)
7. [ACMG Classification Logic & Caveats](#acmg-classification-logic--caveats)

---

## Pipeline Workflow

| Step | Tool | Description |
| :--- | :--- | :--- |
| **0** | Bash Config | Validates software dependencies, reference files, and loads run variables |
| **1** | **FastQC** | Evaluates raw paired-end FASTQ read quality |
| **2** | **cutadapt** | Trims Illumina adapters and filters low-quality/short reads |
| **3** | **BWA-MEM** | Aligns trimmed reads to the reference genome with `@RG` read-group headers |
| **4** | **samtools** | Sorts BAM, creates `.bai` index, and generates `flagstat` & `idxstats` |
| **5** | **Picard** | Marks PCR/optical duplicates (`MarkDuplicates`) and logs library metrics |
| **6** | **GATK4** | Runs BQSR (`BaseRecalibrator` + `ApplyBQSR`), `HaplotypeCaller` (`-ERC GVCF`), `GenotypeGVCFs`, splits SNPs/Indels, applies GATK hard-filters, and merges VCFs |
| **7** | **vcffilter** | Applies site-level quality (`QUAL`) and read-depth (`DP`) thresholds, then compresses (`bgzip`) and indexes (`tabix`) |
| **8** | **ANNOVAR** | Annotates variants against `refGene`, `cytoBand`, `avsnp150`, `gnomad`, `clinvar`, and `dbnsfp` |
| **9** | **awk** | Filters annotated variants against a user-defined target gene panel |
| **10** | **Python 3** | Evaluates functional/population/ClinVar annotations to assign ACMG 2015 evidence codes and tiers (`Pathogenic`, `Likely Pathogenic`, `VUS`, `Likely Benign`, `Benign`) |
| **11** | **Bash / HTML** | Compiles QC metrics, variant counts, and top candidate variants into a standalone HTML report |
| **12** | — | Finalizes logs and reports output file paths |

---

## Prerequisites & Installation

### System Requirements
* **OS:** Linux (Ubuntu/Debian/CentOS) or macOS
* **Shell:** Bash `>= 4.4`
* **Java:** OpenJDK `17+` (required by GATK `>= 4.4` and Picard)
* **Python:** Python `3.7+` (uses standard library only—no `pip` packages required)

### 1. Installing Open-Source Tools via Conda / Mamba (Recommended)
You can install every open-source dependency in a single isolated environment using [Miniforge/Mamba](https://github.com/conda-forge/miniforge) or Miniconda:

```bash
# Create and activate the environment
mamba create -n wes-acmg -c conda-forge -c bioconda \
  fastqc \
  cutadapt>=3.0 \
  bwa \
  samtools>=1.10 \
  picard>=2.20 \
  gatk4>=4.4 \
  vcflib \
  htslib \
  perl \
  python>=3.8

conda activate wes-acmg

```

*(Note: `htslib` provides `bgzip` and `tabix`.)*

### 2. Installing ANNOVAR (Manual Registration Required)

ANNOVAR is free for academic/non-profit users but requires registration and cannot be distributed via Conda:

1. Register and download the archive from the [ANNOVAR download portal](https://www.openbioinformatics.org/annovar/annovar_download_form.php).
2. Unpack it into your tools directory (e.g., `/tools/annovar`) and add it to your `$PATH`:

```bash
tar -xzvf annovar.latest.tar.gz -C /tools/
export PATH="/tools/annovar:$PATH"
# Tip: Add the export line to your ~/.bashrc so table_annovar.pl is always found

```

### 3. Verify All Tools

Confirm all required executables are visible in your `$PATH`:

```bash
for tool in fastqc cutadapt bwa samtools picard gatk vcffilter table_annovar.pl bgzip tabix python3; do
  command -v "$tool" >/dev/null && echo "[OK] $tool" || echo "[MISSING] $tool"
done

```

---

## Reference & Database Setup

### 1. Reference Genome & GATK Known Sites (`hg38`)

Download the reference FASTA and known-sites VCFs from the [GATK Resource Bundle (Google Cloud)](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0):

```bash
mkdir -p /ref/GRCh38 && cd /ref/GRCh38

# Download FASTA, index, dictionary, and known sites
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta GRCh38.fa
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai GRCh38.fa.fai
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict GRCh38.dict
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf dbsnp_146.hg38.vcf
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz .
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi .

# Compress & index dbSNP if downloaded uncompressed
bgzip dbsnp_146.hg38.vcf && tabix -p vcf dbsnp_146.hg38.vcf.gz

# Generate BWA index files (.amb, .ann, .bwt, .pac, .sa)
bwa index GRCh38.fa

```

### 2. ANNOVAR Annotation Databases

Download the databases matching the `ANNOVAR_PROTOCOLS` defined in your config file:

```bash
cd /tools/annovar

perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad41_genome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp47a humandb/

```

### 3. Target Gene List (`GENE_LIST`)

Create a plain-text file containing **one HGNC gene symbol per line with no header**:

```bash
mkdir -p /data/panels
cat << 'EOF' > /data/panels/my_gene_list.txt
BRCA1
BRCA2
TP53
PALB2
CHEK2
EOF

```

---

## Configuration (`WES-ACMG-Pipeline.config.txt`)

Before running the pipeline, copy or edit `WES-ACMG-Pipeline.config.txt`. **All file and directory paths must be absolute.**

```bash
# Sample & Compute
SAMPLE="Sample01"
LIBRARY="Lib01"
R1_FASTQ="/data/fastq/Sample01_R1.fastq.gz"
R2_FASTQ="/data/fastq/Sample01_R2.fastq.gz"
THREADS=8
OUT_DIR="/data/results/WES_${SAMPLE}"

# Reference & Known Sites
REF_FASTA="/ref/GRCh38/GRCh38.fa"
DBSNP_VCF="/ref/GRCh38/dbsnp_146.hg38.vcf.gz"
MILLS_VCF="/ref/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
GENOME_BUILD="hg38"          # hg19 or hg38

# cutadapt Trimming
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
CUTADAPT_QUAL=20
CUTADAPT_MINLEN=36
CUTADAPT_MAXN=0

# vcffilter Thresholds
VCFFILTER_QUAL=30
VCFFILTER_DP=10
VCFFILTER_MAXDP=200

# ANNOVAR Settings
ANNOVAR_DB="/tools/annovar/humandb"
ANNOVAR_PROTOCOLS="refGene,cytoBand,avsnp150,gnomad41_genome,clinvar_20240917,dbnsfp47a"
ANNOVAR_OPERATIONS="gx,r,f,f,f,f"

# Target Gene Panel
GENE_LIST="/data/panels/my_gene_list.txt"

```

> **Note on `ANNOVAR_PROTOCOLS`:** Step 9 filters by gene symbol using column 7 (`Gene.refGene`) of the `_multianno.txt` file. If you alter the order of `refGene` in `ANNOVAR_PROTOCOLS`, update `$7` in Step 9 accordingly.

---

## Usage & CLI Execution

### 1. Make the Script Executable

```bash
chmod +x WES-ACMG-Pipeline.sh

```

### 2. Run on a Single Sample

Pass the path to your configuration file as the first positional argument (if omitted, it defaults to `./WES-ACMG-Pipeline.config.txt`):

```bash
./WES-ACMG-Pipeline.sh /absolute/path/to/WES-ACMG-Pipeline.config.txt

```

### 3. Run in the Background (Recommended for Large Exomes)

```bash
nohup ./WES-ACMG-Pipeline.sh /data/configs/Sample01.config.txt > Sample01.nohup.out 2>&1 &

# Monitor live progress
tail -f /data/results/WES_Sample01/00_logs/pipeline.log

```

### 4. Batch Execution Across Multiple Samples

Create per-sample config files (e.g., `Sample01.config.txt`, `Sample02.config.txt`) and iterate through them sequentially:

```bash
for cfg in /data/configs/*.config.txt; do
  echo "Running pipeline for $cfg..."
  ./WES-ACMG-Pipeline.sh "$cfg"
done

```

---

## Output Directory Structure

All results are organized into numbered subdirectories under `$OUT_DIR`:

```text
/data/results/WES_Sample01/
├── 00_logs/
│   └── pipeline.log                              # Timestamped master execution log
├── 01_fastqc/
│   ├── Sample01_R1_fastqc.html                   # Raw R1 quality report
│   └── Sample01_R2_fastqc.html                   # Raw R2 quality report
├── 02_trimmed/
│   ├── Sample01_R1.trimmed.fastq.gz              # Adapter/quality-trimmed R1 reads
│   ├── Sample01_R2.trimmed.fastq.gz              # Adapter/quality-trimmed R2 reads
│   └── cutadapt.log                              # Trimming summary statistics
├── 03_bam/
│   ├── Sample01.sorted.bam                       # Coordinate-sorted alignment
│   ├── Sample01.sorted.bam.bai                   # BAM index
│   ├── Sample01.flagstat.txt                     # Alignment mapping statistics
│   └── Sample01.idxstats.txt                     # Per-chromosome read counts
├── 04_picard/
│   ├── Sample01.markdup.bam                      # Duplicate-marked BAM
│   ├── Sample01.markdup.bai                      # BAM index
│   └── Sample01.markdup_metrics.txt              # Duplication rate & library metrics
├── 05_gatk/
│   ├── Sample01.recal.table                      # BQSR recalibration table
│   ├── Sample01.bqsr.bam                         # Base-recalibrated BAM
│   ├── Sample01.raw.vcf.gz                       # Raw GVCF from HaplotypeCaller
│   ├── Sample01.genotyped.vcf.gz                 # Joint-genotyped single-sample VCF
│   ├── Sample01.snps.hardfiltered.vcf.gz         # Hard-filtered SNPs
│   ├── Sample01.indels.hardfiltered.vcf.gz       # Hard-filtered Indels
│   └── Sample01.merged.hardfiltered.vcf.gz       # Combined hard-filtered VCF
├── 06_filtered_vcf/
│   ├── Sample01.filtered.vcf.gz                  # Variants passing QUAL & DP filters
│   └── Sample01.filtered.vcf.gz.tbi              # Tabix index
├── 07_annovar/
│   └── Sample01.annovar.hg38_multianno.txt       # Full annotated variant table
├── 08_gene_filtered/
│   └── Sample01.gene_filtered_multianno.txt      # Variants restricted to GENE_LIST
├── 09_acmg/
│   ├── acmg_classify.py                          # Auto-generated classification script
│   ├── Sample01.acmg_classified.txt              # All panel variants with ACMG codes/classes
│   └── Sample01.ACMG_candidates.txt              # Subset classified as Pathogenic / Likely Pathogenic
└── 10_report/
    └── Sample01_WES_report.html                  # Self-contained visual summary report

```

---

## ACMG Classification Logic & Caveats

Step 10 embeds and runs a Python 3 rules engine (`acmg_classify.py`) that parses ANNOVAR annotations and applies the **Richards et al. (2015)** evidence-combining criteria:

* **Functional Impact (`ExonicFunc.refGene`):** Assigns `PVS1` to predicted null variants (`frameshift`, `stopgain`, `stoploss`, `splice`) and `PM1` to `nonsynonymous` variants.
* **Population Allele Frequency (`gnomAD` / `ExAC` / `1000G`):** Assigns `PM2` (`AF < 0.0001`), `PP3` (`0.0001 <= AF < 0.01`), or `BP4` (`AF >= 0.01`).
* **Clinical Assertions (`ClinVar`):** Assigns `PP5` for non-conflicting Pathogenic entries and `BP6` for Benign entries.

> **Disclaimer — Research Use Only:** The automated ACMG classification module is a **screening implementation** based on tabular annotations and does not evaluate segregation data, *de novo* status, functional assay literature, or gene-specific dosage sensitivity. All `Pathogenic`, `Likely Pathogenic`, and `VUS` calls must be manually curated by a certified clinical molecular geneticist prior to clinical reporting.

```

```
