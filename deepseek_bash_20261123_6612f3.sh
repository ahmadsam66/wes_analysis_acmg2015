#!/usr/bin/env bash
# =============================================================================
#  WES-ACMG-Pipeline.sh
#  Whole-Exome Sequencing → ACMG 2015 Variant Classification
# =============================================================================
#  A modular, human-readable germline WES pipeline:
#    FastQC → cutadapt → BWA-MEM → samtools → Picard → GATK4 (BQSR,
#    HaplotypeCaller, hard-filter) → vcffilter → ANNOVAR → gene-list
#    filtration → ACMG rule-based classification → HTML report.
#
#  Requirements: bash ≥4.4, GATK ≥4.4, BWA, samtools, Picard, FastQC,
#                cutadapt, vcflib, ANNOVAR, bgzip/tabix, python3.
# =============================================================================
set -euo pipefail
IFS=$'\n\t'

# ---------------------------------------------------------------------------
# 0. Load configuration
# ---------------------------------------------------------------------------
CONFIG="${1:-./WES-ACMG-Pipeline.config.txt}"
if [[ ! -f "$CONFIG" ]]; then
  echo "[ERROR] Configuration file not found: $CONFIG" >&2
  echo "Usage: $0 /path/to/WES-ACMG-Pipeline.config.txt" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "$CONFIG"

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG_FILE"; }
die()  { log "FATAL: $*"; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Required tool not found: $1"; }

# ---------------------------------------------------------------------------
# Validate software and references
# ---------------------------------------------------------------------------
for tool in fastqc cutadapt bwa samtools picard gatk vcffilter table_annovar.pl python3; do
  need "$tool"
done
[[ -f "$REF_FASTA" ]]       || die "Reference FASTA missing: $REF_FASTA"
[[ -f "$DBSNP_VCF" ]]       || die "dbSNP VCF missing: $DBSNP_VCF"
[[ -f "$MILLS_VCF" ]]       || die "Mills indels VCF missing: $MILLS_VCF"
[[ -f "$GENE_LIST" ]]       || die "Gene list missing: $GENE_LIST"

# ---------------------------------------------------------------------------
# Create output directories
# ---------------------------------------------------------------------------
mkdir -p \
  "${OUT_DIR}"/{00_logs,01_fastqc,02_trimmed,03_bam,04_picard,05_gatk,06_filtered_vcf,07_annovar,08_gene_filtered,09_acmg,10_report}

LOG_FILE="${OUT_DIR}/00_logs/pipeline.log"
: > "$LOG_FILE"
log "Pipeline started"
log "Config: $CONFIG"
log "Sample: $SAMPLE"
log "Output: $OUT_DIR"

# ===========================================================================
# 1. FastQC — raw read quality
# ===========================================================================
log "STEP 1/12 — FastQC on raw reads"
fastqc -t "$THREADS" -o "${OUT_DIR}/01_fastqc" "$R1_FASTQ" "$R2_FASTQ"

# ===========================================================================
# 2. cutadapt — adapter trimming and quality filtering
# ===========================================================================
log "STEP 2/12 — cutadapt paired-end trimming"
cutadapt \
  -j "$THREADS" \
  -a "$ADAPTER_R1" -A "$ADAPTER_R2" \
  -q "$CUTADAPT_QUAL" \
  -m "$CUTADAPT_MINLEN" \
  --max-n "$CUTADAPT_MAXN" \
  -o "${OUT_DIR}/02_trimmed/${SAMPLE}_R1.trimmed.fastq.gz" \
  -p "${OUT_DIR}/02_trimmed/${SAMPLE}_R2.trimmed.fastq.gz" \
  "$R1_FASTQ" "$R2_FASTQ" \
  > "${OUT_DIR}/02_trimmed/cutadapt.log" 2>&1

# ===========================================================================
# 3. BWA-MEM — alignment
# ===========================================================================
log "STEP 3/12 — BWA-MEM alignment"
bwa mem \
  -t "$THREADS" \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${LIBRARY}" \
  "$REF_FASTA" \
  "${OUT_DIR}/02_trimmed/${SAMPLE}_R1.trimmed.fastq.gz" \
  "${OUT_DIR}/02_trimmed/${SAMPLE}_R2.trimmed.fastq.gz" \
  2> "${OUT_DIR}/03_bam/bwa.log" \
  | samtools sort -@ "$THREADS" -o "${OUT_DIR}/03_bam/${SAMPLE}.sorted.bam" -

samtools index "${OUT_DIR}/03_bam/${SAMPLE}.sorted.bam"

# ===========================================================================
# 4. samtools — mapping statistics
# ===========================================================================
log "STEP 4/12 — samtools flagstat and idxstats"
samtools flagstat "${OUT_DIR}/03_bam/${SAMPLE}.sorted.bam" \
  > "${OUT_DIR}/03_bam/${SAMPLE}.flagstat.txt"
samtools idxstats "${OUT_DIR}/03_bam/${SAMPLE}.sorted.bam" \
  > "${OUT_DIR}/03_bam/${SAMPLE}.idxstats.txt"

# ===========================================================================
# 5. Picard — MarkDuplicates
# ===========================================================================
log "STEP 5/12 — Picard MarkDuplicates"
picard MarkDuplicates \
  I="${OUT_DIR}/03_bam/${SAMPLE}.sorted.bam" \
  O="${OUT_DIR}/04_picard/${SAMPLE}.markdup.bam" \
  M="${OUT_DIR}/04_picard/${SAMPLE}.markdup_metrics.txt" \
  REMOVE_DUPLICATES=false \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  > "${OUT_DIR}/04_picard/picard_markdup.log" 2>&1

# ===========================================================================
# 6. GATK4 — BQSR, HaplotypeCaller, hard-filtering
# ===========================================================================
log "STEP 6/12 — GATK4 BaseRecalibrator"
gatk BaseRecalibrator \
  -R "$REF_FASTA" \
  -I "${OUT_DIR}/04_picard/${SAMPLE}.markdup.bam" \
  --known-sites "$DBSNP_VCF" \
  --known-sites "$MILLS_VCF" \
  -O "${OUT_DIR}/05_gatk/${SAMPLE}.recal.table" \
  > "${OUT_DIR}/05_gatk/bqsr_recal.log" 2>&1

log "STEP 6/12 — GATK4 ApplyBQSR"
gatk ApplyBQSR \
  -R "$REF_FASTA" \
  -I "${OUT_DIR}/04_picard/${SAMPLE}.markdup.bam" \
  --bqsr-recal-file "${OUT_DIR}/05_gatk/${SAMPLE}.recal.table" \
  -O "${OUT_DIR}/05_gatk/${SAMPLE}.bqsr.bam" \
  > "${OUT_DIR}/05_gatk/bqsr_apply.log" 2>&1

log "STEP 6/12 — GATK4 HaplotypeCaller (single-sample direct VCF)"
gatk HaplotypeCaller \
  -R "$REF_FASTA" \
  -I "${OUT_DIR}/05_gatk/${SAMPLE}.bqsr.bam" \
  -O "${OUT_DIR}/05_gatk/${SAMPLE}.raw.vcf.gz" \
  -ERC GVCF \
  --native-pair-hmm-threads "$THREADS" \
  > "${OUT_DIR}/05_gatk/haplotypecaller.log" 2>&1

log "STEP 6/12 — GATK4 GenotypeGVCFs"
gatk GenotypeGVCFs \
  -R "$REF_FASTA" \
  -V "${OUT_DIR}/05_gatk/${SAMPLE}.raw.vcf.gz" \
  -O "${OUT_DIR}/05_gatk/${SAMPLE}.genotyped.vcf.gz" \
  > "${OUT_DIR}/05_gatk/genotypegvcfs.log" 2>&1

# --- Split SNPs and indels, apply GATK-recommended hard filters -------------
log "STEP 6/12 — GATK4 hard-filtering (SNPs and indels)"
# SNPs
gatk SelectVariants \
  -R "$REF_FASTA" \
  -V "${OUT_DIR}/05_gatk/${SAMPLE}.genotyped.vcf.gz" \
  --select-type-to-include SNP \
  -O "${OUT_DIR}/05_gatk/${SAMPLE}.snps.vcf.gz" \
  > "${OUT_DIR}/05_gatk/select_snps.log" 2>&1

gatk VariantFiltration \
  -R "$REF_FASTA" \
  -V "${OUT_DIR}/05_gatk/${SAMPLE}.snps.vcf.gz" \
  -O "${OUT_DIR}/05_gatk/${SAMPLE}.snps.hardfiltered.vcf.gz" \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "SNP_HARD_FILTER" \
  > "${OUT_DIR}/05_gatk/filter_snps.log" 2>&1

# Indels
gatk SelectVariants \
  -R "$REF_FASTA" \
  -V "${OUT_DIR}/05_gatk/${SAMPLE}.genotyped.vcf.gz" \
  --select-type-to-include INDEL \
  -O "${OUT_DIR}/05_gatk/${SAMPLE}.indels.vcf.gz" \
  > "${OUT_DIR}/05_gatk/select_indels.log" 2>&1

gatk VariantFiltration \
  -R "$REF_FASTA" \
  -V "${OUT_DIR}/05_gatk/${SAMPLE}.indels.vcf.gz" \
  -O "${OUT_DIR}/05_gatk/${SAMPLE}.indels.hardfiltered.vcf.gz" \
  --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0" \
  --filter-name "INDEL_HARD_FILTER" \
  > "${OUT_DIR}/05_gatk/filter_indels.log" 2>&1

# Merge PASS SNPs + PASS indels
gatk MergeVcfs \
  -I "${OUT_DIR}/05_gatk/${SAMPLE}.snps.hardfiltered.vcf.gz" \
  -I "${OUT_DIR}/05_gatk/${SAMPLE}.indels.hardfiltered.vcf.gz" \
  -O "${OUT_DIR}/05_gatk/${SAMPLE}.merged.hardfiltered.vcf.gz" \
  > "${OUT_DIR}/05_gatk/merge_vcfs.log" 2>&1

# ===========================================================================
# 7. vcffilter — additional quality / depth filtration
# ===========================================================================
log "STEP 7/12 — vcffilter (QUAL and DP)"
vcffilter \
  -f "QUAL > ${VCFFILTER_QUAL} & DP > ${VCFFILTER_DP} & DP < ${VCFFILTER_MAXDP}" \
  "${OUT_DIR}/05_gatk/${SAMPLE}.merged.hardfiltered.vcf.gz" \
  > "${OUT_DIR}/06_filtered_vcf/${SAMPLE}.filtered.vcf"

bgzip -f "${OUT_DIR}/06_filtered_vcf/${SAMPLE}.filtered.vcf"
tabix -f -p vcf "${OUT_DIR}/06_filtered_vcf/${SAMPLE}.filtered.vcf.gz"

# ===========================================================================
# 8. ANNOVAR — functional annotation
# ===========================================================================
log "STEP 8/12 — ANNOVAR table_annovar.pl"
table_annovar.pl \
  "${OUT_DIR}/06_filtered_vcf/${SAMPLE}.filtered.vcf.gz" \
  "$ANNOVAR_DB" \
  -buildver "$GENOME_BUILD" \
  -out "${OUT_DIR}/07_annovar/${SAMPLE}.annovar" \
  -remove \
  -protocol "$ANNOVAR_PROTOCOLS" \
  -operation "$ANNOVAR_OPERATIONS" \
  -nastring . \
  -polish \
  --vcfinput \
  > "${OUT_DIR}/07_annovar/annovar.log" 2>&1

ANNOVAR_MULTIANNO="${OUT_DIR}/07_annovar/${SAMPLE}.annovar.${GENOME_BUILD}_multianno.txt"
[[ -f "$ANNOVAR_MULTIANNO" ]] || die "ANNOVAR multianno file not found: $ANNOVAR_MULTIANNO"

# ===========================================================================
# 9. Selected-gene filtration (based on input gene list)
# ===========================================================================
log "STEP 9/12 — Gene-list filtration"
# Gene list: one gene symbol per line, no header.
# ANNOVAR 'Gene.refGene' is usually column 7; adjust if your protocol changes.
awk -F'\t' '
  NR==FNR { keep[$1]=1; next }
  FNR==1  { print; next }
  ($7 in keep) { print }
' "$GENE_LIST" "$ANNOVAR_MULTIANNO" \
  > "${OUT_DIR}/08_gene_filtered/${SAMPLE}.gene_filtered_multianno.txt"

# Also keep a VCF-compatible subset for downstream tools (optional)
awk -F'\t' '
  NR==FNR { keep[$1]=1; next }
  FNR==1  { print; next }
  ($7 in keep) { print }
' "$GENE_LIST" "${OUT_DIR}/07_annovar/${SAMPLE}.annovar.${GENOME_BUILD}_multianno.txt" \
  > "${OUT_DIR}/08_gene_filtered/${SAMPLE}.gene_filtered_raw.txt"

# ===========================================================================
# 10. ACMG 2015 rule-based classification
# ===========================================================================
log "STEP 10/12 — ACMG 2015 rule-based classification"

ACMG_PY="${OUT_DIR}/09_acmg/acmg_classify.py"
cat > "$ACMG_PY" <<'PY'
#!/usr/bin/env python3
"""Rule-based ACMG 2015 classification for ANNOVAR multianno output.

NOTE: This is a *screening* implementation. It assigns evidence codes from
available ANNOVAR fields and applies the 2015 combining rules. Clinical
classification MUST be confirmed by a certified variant scientist / MD.
"""
import sys, csv

HEADER = None
rows = []

def parse_exonic_func(func):
    func = func.strip()
    if any(k in func for k in ("frameshift", "stopgain", "stoploss", "splice")):
        return "PVS1"   # very strong — null variant (simplified)
    if "nonsynonymous" in func:
        return "PM1"    # moderate — missense (placeholder; requires gene-specific rules)
    return None

def parse_pop_freq(freq_str):
    try:
        f = float(freq_str)
    except Exception:
        return None
    if f < 0.0001:
        return "PM2"    # moderate — absent/rare in population databases
    if f < 0.01:
        return "PP3"    # supporting — computational / population evidence
    return "BP4"        # benign supporting — common polymorphism

def combine(codes):
    """Apply 2015 ACMG combining rules (simplified, pathogenic focus)."""
    c = set(codes)
    vs = {"PVS1"}; s = {"PS1","PS2","PS3","PS4"}
    m = {"PM1","PM2","PM3","PM4","PM5","PM6"}
    p = {"PP1","PP2","PP3","PP4","PP5"}
    # Pathogenic
    if vs and (len(c & s) >= 1 or len(c & m) >= 2 or (len(c & m) >= 1 and len(c & p) >= 1)):
        return "Pathogenic"
    if len(c & s) >= 2:
        return "Pathogenic"
    if len(c & s) >= 1 and len(c & m) >= 3:
        return "Pathogenic"
    if len(c & s) >= 1 and len(c & m) >= 1 and len(c & p) >= 4:
        return "Pathogenic"
    # Likely Pathogenic
    if vs and len(c & p) >= 1:
        return "Likely Pathogenic"
    if len(c & s) >= 1 and len(c & m) >= 1:
        return "Likely Pathogenic"
    if len(c & s) >= 1 and len(c & p) >= 2:
        return "Likely Pathogenic"
    if len(c & m) >= 3:
        return "Likely Pathogenic"
    # Benign / Likely benign
    if "BA1" in c or len(c & {"BS1","BS2","BS3","BS4"}) >= 2:
        return "Benign"
    if len(c & {"BP1","BP2","BP3","BP4","BP5","BP6","BP7"}) >= 2:
        return "Likely Benign"
    return "Uncertain Significance (VUS)"

# Read ANNOVAR multianno
with open(sys.argv[1], newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        codes = []
        ef = row.get("ExonicFunc.refGene", "") or ""
        code = parse_exonic_func(ef)
        if code: codes.append(code)

        # Population frequency (try gnomAD / 1000G / ExAC columns)
        for col in ("gnomAD_exome_ALL", "gnomAD_genome_ALL",
                    "ExAC_ALL", "1000g2015aug_ALL"):
            if col in row and row[col] not in ("", "."):
                pc = parse_pop_freq(row[col])
                if pc: codes.append(pc)
                break

        # ClinVar — supporting pathogenic evidence
        if "ClinVar" in row:
            cv = (row.get("ClinVar", "") or "").lower()
            if "pathogenic" in cv and "conflict" not in cv:
                codes.append("PP5")
            elif "benign" in cv:
                codes.append("BP6")

        acmg = combine(codes)
        row["ACMG_Codes"] = ";".join(codes) if codes else "."
        row["ACMG_Class"] = acmg
        rows.append(row)

if not rows:
    print("No variants to classify.", file=sys.stderr)
    sys.exit(0)

out_cols = list(rows[0].keys())
w = csv.DictWriter(sys.stdout, fieldnames=out_cols, delimiter="\t",
                   lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar="\\")
w.writeheader()
w.writerows(rows)
PY

python3 "$ACMG_PY" "${OUT_DIR}/08_gene_filtered/${SAMPLE}.gene_filtered_multianno.txt" \
  > "${OUT_DIR}/09_acmg/${SAMPLE}.acmg_classified.txt"

# Candidate list: Pathogenic / Likely Pathogenic
awk -F'\t' 'NR==1 || $NF ~ /Pathogenic/ && $NF !~ /Likely Benign/' \
  "${OUT_DIR}/09_acmg/${SAMPLE}.acmg_classified.txt" \
  > "${OUT_DIR}/09_acmg/${SAMPLE}.ACMG_candidates.txt"

# ===========================================================================
# 11. HTML report
# ===========================================================================
log "STEP 11/12 — Generating HTML report"

FLAGSTAT=$(cat "${OUT_DIR}/03_bam/${SAMPLE}.flagstat.txt")
MAPPING_RATE=$(awk '/mapped \(/ {gsub(/[()%]/,"",$5); print $5}' \
  "${OUT_DIR}/03_bam/${SAMPLE}.flagstat.txt" | head -1)
RAW_VCF_COUNT=$(zcat "${OUT_DIR}/05_gatk/${SAMPLE}.genotyped.vcf.gz" | grep -vc '^#')
FILTERED_VCF_COUNT=$(zcat "${OUT_DIR}/06_filtered_vcf/${SAMPLE}.filtered.vcf.gz" | grep -vc '^#')
GENE_COUNT=$(awk 'NR>1' "${OUT_DIR}/08_gene_filtered/${SAMPLE}.gene_filtered_multianno.txt" | wc -l)
ACMG_TOTAL=$(awk 'NR>1' "${OUT_DIR}/09_acmg/${SAMPLE}.acmg_classified.txt" | wc -l)
CAND_COUNT=$(awk 'NR>1' "${OUT_DIR}/09_acmg/${SAMPLE}.ACMG_candidates.txt" | wc -l)
DUP_METRICS=$(grep -A1 '^LIBRARY' "${OUT_DIR}/04_picard/${SAMPLE}.markdup_metrics.txt" | tail -1 || true)

cat > "${OUT_DIR}/10_report/${SAMPLE}_WES_report.html" <<HTML
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>WES-ACMG Pipeline Report — ${SAMPLE}</title>
<style>
  body{font-family:system-ui,sans-serif;margin:2rem;background:#f7f9fc;color:#1a202c}
  h1{border-bottom:3px solid #2b6cb0;padding-bottom:.4rem}
  h2{color:#2b6cb0;margin-top:2rem}
  .card{background:#fff;border-radius:8px;padding:1.2rem 1.6rem;margin:1rem 0;
        box-shadow:0 1px 4px rgba(0,0,0,.08)}
  .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1rem}
  .stat{background:#ebf4ff;border-left:4px solid #2b6cb0;padding:.8rem 1rem;border-radius:4px}
  .stat b{display:block;font-size:1.4rem;color:#2b6cb0}
  table{width:100%;border-collapse:collapse;font-size:.9rem}
  th,td{padding:.45rem .6rem;border-bottom:1px solid #e2e8f0;text-align:left}
  th{background:#2b6cb0;color:#fff;position:sticky;top:0}
  tr:hover{background:#f0f4ff}
  .warn{background:#fffaf0;border-left:4px solid #dd6b20;padding:.8rem 1rem;border-radius:4px}
  .footer{margin-top:3rem;font-size:.8rem;color:#718096}
</style>
</head>
<body>
<h1>WES-ACMG Pipeline Report</h1>
<p><strong>Sample:</strong> ${SAMPLE} &nbsp;|&nbsp;
   <strong>Genome build:</strong> ${GENOME_BUILD} &nbsp;|&nbsp;
   <strong>Generated:</strong> $(date '+%Y-%m-%d %H:%M:%S')</p>

<h2>1. Read &amp; alignment statistics</h2>
<div class="grid">
  <div class="stat"><b>${MAPPING_RATE:-N/A}%</b>Mapping rate</div>
  <div class="stat"><b>$(echo "$FLAGSTAT" | awk '/in total/ {print $1}')</b>Total reads</div>
  <div class="stat"><b>$(echo "$FLAGSTAT" | awk '/properly paired/ {print $1}')</b>Properly paired</div>
  <div class="stat"><b>$(echo "$FLAGSTAT" | awk '/duplicates/ {print $1}')</b>Duplicates</div>
</div>
<pre class="card">${FLAGSTAT}</pre>

<h2>2. Variant calling &amp; filtration</h2>
<div class="grid">
  <div class="stat"><b>${RAW_VCF_COUNT}</b>Raw variants (GenotypeGVCFs)</div>
  <div class="stat"><b>${FILTERED_VCF_COUNT}</b>PASS after vcffilter</div>
  <div class="stat"><b>${GENE_COUNT}</b>Variants in selected genes</div>
  <div class="stat"><b>${ACMG_TOTAL}</b>Variants classified (ACMG)</div>
  <div class="stat"><b>${CAND_COUNT}</b>ACMG candidates (P / LP)</div>
</div>

<h2>3. ACMG 2015 candidate variants</h2>
<div class="warn">
<strong>Research use only.</strong> ACMG codes are assigned by rule-based screening
from ANNOVAR annotations. Final clinical interpretation requires manual review
according to Richards et al. 2015 (Genet Med 17:405–424).
</div>
<div class="card" style="overflow-x:auto">
<table>
$(head -1 "${OUT_DIR}/09_acmg/${SAMPLE}.ACMG_candidates.txt" | awk -F'\t' '{printf "<tr>"; for(i=1;i<=NF;i++) printf "<th>%s</th>",$i; print "</tr>"}')
$(awk -F'\t' 'NR>1 && NR<=26 {printf "<tr>"; for(i=1;i<=NF;i++) printf "<td>%s</td>",$i; print "</tr>"}' "${OUT_DIR}/09_acmg/${SAMPLE}.ACMG_candidates.txt")
</table>
</div>
<p>Showing at most 25 candidate rows. Full list:
<code>${OUT_DIR}/09_acmg/${SAMPLE}.ACMG_candidates.txt</code></p>

<h2>4. Output directory structure</h2>
<pre class="card">${OUT_DIR}
├── 00_logs/          pipeline.log
├── 01_fastqc/        raw read QC
├── 02_trimmed/       cutadapt outputs
├── 03_bam/           sorted BAM, flagstat, idxstats
├── 04_picard/        markdup BAM + metrics
├── 05_gatk/          BQSR, GVCF, genotyped VCF, hard-filtered VCF
├── 06_filtered_vcf/  vcffilter PASS VCF
├── 07_annovar/       ANNOVAR multianno
├── 08_gene_filtered/ gene-list filtered multianno
├── 09_acmg/          ACMG classified + candidate list
└── 10_report/        this HTML report</pre>

<p class="footer">Generated by WES-ACMG-Pipeline.sh · $(date -u '+%Y-%m-%dT%H:%M:%SZ')</p>
</body>
</html>
HTML

log "STEP 12/12 — Pipeline finished successfully"
log "HTML report: ${OUT_DIR}/10_report/${SAMPLE}_WES_report.html"
log "ACMG candidates: ${OUT_DIR}/09_acmg/${SAMPLE}.ACMG_candidates.txt"