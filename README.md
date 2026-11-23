# WES-ACMG-Pipeline

A  Bash pipeline for **whole-exome sequencing (WES) germline
variant discovery and ACMG 2015 rule-based classification**.

The pipeline is designed for small-to-medium single-sample or per-sample
analysis. It follows GATK4 Best Practices for germline short variants and
uses **hard-filtering** (recommended when VQSR cannot be run on small WES
cohorts)[reference:2].

---

## Pipeline steps

| Step | Tool | Purpose |
|------|------|---------|
| 0 | Bash config | Load paths, genome build, gene list, thresholds |
| 1 | FastQC | Raw read quality control |
| 2 | cutadapt | Adapter trimming + quality/length filtering |
| 3 | BWA-MEM | Paired-end alignment with read-group header |
| 4 | samtools | Sort, index, `flagstat`, `idxstats` |
| 5 | Picard MarkDuplicates | PCR duplicate marking + metrics |
| 6 | GATK4 | BQSR (`BaseRecalibrator` → `ApplyBQSR`), `HaplotypeCaller` (GVCF), `GenotypeGVCFs`, SNP/indel split + `VariantFiltration` |
| 7 | vcffilter (vcflib) | Additional QUAL / DP site-level filtration |
| 8 | ANNOVAR | Functional annotation (`refGene`, `avsnp150`, `gnomad`, `clinvar`, `dbnsfp`)[reference:3] |
| 9 | awk | Restrict to genes in the user-supplied gene list |
| 10 | Python (rule-based) | Assign ACMG 2015 evidence codes and classify P / LP / VUS / LB / B[reference:4] |
| 11 | Bash | Generate a self-contained HTML report with QC, mapping, calling, filtration, and candidate statistics |
| 12 | — | All results organised into numbered output directories |

---

## Requirements

- **Bash** ≥ 4.4
- **FastQC**
- **cutadapt** ≥ 3
- **BWA** (or `bwa-mem2`)
- **samtools** ≥ 1.10
- **Picard** ≥ 2.20
- **GATK** ≥ 4.4 (Java 17+)
- **vcflib** (`vcffilter`)
- **ANNOVAR** (with `table_annovar.pl` and prepared `humandb/`)
- **bgzip / tabix** (htslib)
- **Python 3** (standard library only)

---

## Quick start

```bash
# 1. Clone / copy the three files into the same directory
#    WES-ACMG-Pipeline.sh
#    WES-ACMG-Pipeline.config.txt
#    README.md

# 2. Make the script executable
chmod +x WES-ACMG-Pipeline.sh

# 3. Edit the configuration file — every path must be absolute
nano WES-ACMG-Pipeline.config.txt

# 4. Run
./WES-ACMG-Pipeline.sh /absolute/path/to/WES-ACMG-Pipeline.config.txt
