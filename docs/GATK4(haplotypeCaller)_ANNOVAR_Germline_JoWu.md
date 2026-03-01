# Clinical Germline Variant Calling & Annotation Pipeline (Single Sample)
## BWA-MEM/BWA + GATK4 (HaplotypeCaller) + ANNOVAR — hg38

This is a **single-sample germline** pipeline aligned to GATK Best Practices, updated to include the **required BQSR known-sites resources** and cleaner VCF hygiene before annotation.

# Clinical Context — Identifies inherited DNA variants to diagnose genetic disease, assess hereditary risk, and guide medical management. Applications includes: 
Hereditary Cancer: eg., BRCA1/2, TP53, RUNX1, GATA2
Rare Genetic Disease: eg.,SCN1A (epileptic encephalopathy), RAG1 (SCID immunodeficiency)
Carrier Screening: eg., HBB, HBA1/HBA2 (thalassemia)
Pharmacogenomics: eg., CYP2D6, CYP2C19, TPMT (dose optimization, toxicity prevention)



---

## 1) Supported sequencing types

| Sequencing type | Supported | Notes |
|---|---:|---|
| Whole genome (WGS) | ✅ | Best for genome-wide variant discovery |
| Whole exome (WES) | ✅ | Use capture intervals for metrics; may use hard filters |
| Targeted panels | ✅ | Prefer panel intervals + higher depth thresholds |
| RNA-seq | ⚠️ | Not recommended for primary germline calling (expression bias) |

---

## 2) Required inputs

### Per sample
- Paired-end FASTQ: `sample_R1.fastq.gz`, `sample_R2.fastq.gz`
- Sample metadata (at minimum): `sample_id`, `library_id`, platform (ILLUMINA), lane/unit (optional)

### Shared resources
- Reference genome: **hg38/GRCh38 FASTA** (+ BWA, samtools, GATK indices)
- **Known-sites for BQSR**: dbSNP + Mills & 1000G gold standard indels (hg38)
- (Optional but recommended) BED/intervals for WES/panel for coverage/QC

### Annotation
- ANNOVAR installation + `humandb/`
- Databases: refGene, gnomAD, dbNSFP, ClinVar (and others as needed)

---

## 3) Download Broad (GATK) hg38 resource bundle (free/public)

> These links point to the Broad public buckets used in GATK Best Practices examples.
> If your institution uses a local mirror, substitute accordingly.

### Reference FASTA (hg38/assembly38)
```bash
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/resources-broad-hg38-v0/Homo_sapiens_assembly38.fasta
```

### Known-sites for BQSR (required)
```bash
# dbSNP (hg38)
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/resources-broad-hg38-v0/dbsnp_146.hg38.vcf.gz

# Mills & 1000G gold standard indels (hg38)
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/resources-broad-hg38-v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

Index known-sites VCFs:
```bash
gatk IndexFeatureFile -I dbsnp_146.hg38.vcf.gz
gatk IndexFeatureFile -I Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

---

## 4) Phase 0 — Environment + raw read QC

```bash
# (optional) activate env
source ~/myenv1/bin/activate

# QC
fastqc sample_R1.fastq.gz sample_R2.fastq.gz
```

---

## 5) Phase 1 — Reference preparation

```bash
# BWA index
bwa index Homo_sapiens_assembly38.fasta

# FASTA index
samtools faidx Homo_sapiens_assembly38.fasta

# Sequence dictionary (required by GATK)
gatk CreateSequenceDictionary -R Homo_sapiens_assembly38.fasta
```

---

## 6) Phase 2 — Alignment and BAM processing

### 6.1 Align (BWA-MEM; include read groups)
```bash
bwa mem -t 8 \
  -R "@RG\tID:lane1\tLB:lib1\tPL:ILLUMINA\tSM:sample1\tPU:unit1" \
  Homo_sapiens_assembly38.fasta \
  sample_R1.fastq.gz sample_R2.fastq.gz > sample1.sam
```

### 6.2 Sort + index
```bash
samtools sort -@ 8 -o sample1.sorted.bam sample1.sam
samtools index sample1.sorted.bam
```

### 6.3 Mark duplicates + index
```bash
gatk MarkDuplicates \
  -I sample1.sorted.bam \
  -O sample1.markdup.bam \
  -M sample1.markdup.metrics.txt

samtools index sample1.markdup.bam
```

(Optional) quick alignment QC:
```bash
samtools flagstat sample1.markdup.bam
```

---

## 7) Phase 3 — BQSR (required known-sites)

### 7.1 BaseRecalibrator 
```bash
gatk BaseRecalibrator \
  -I sample1.markdup.bam \
  -R Homo_sapiens_assembly38.fasta \
  --known-sites dbsnp_146.hg38.vcf.gz \
  --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O sample1.recal.table
```

### 7.2 ApplyBQSR
```bash
gatk ApplyBQSR \
  -R Homo_sapiens_assembly38.fasta \
  -I sample1.markdup.bam \
  --bqsr-recal-file sample1.recal.table \
  -O sample1.recal.bam
```

Index:
```bash
samtools index sample1.recal.bam
```

---

## 8) Phase 4 — Germline variant calling (HaplotypeCaller)

### Option A (simple single-sample VCF)
```bash
gatk HaplotypeCaller \
  -R Homo_sapiens_assembly38.fasta \
  -I sample1.recal.bam \
  -O sample1.raw.vcf.gz
```

### Option B (recommended pattern even for single sample: GVCF → GenotypeGVCFs)
This makes it easy to scale later to joint genotyping.
```bash
gatk HaplotypeCaller \
  -R Homo_sapiens_assembly38.fasta \
  -I sample1.recal.bam \
  -O sample1.g.vcf.gz \
  -ERC GVCF

gatk GenotypeGVCFs \
  -R Homo_sapiens_assembly38.fasta \
  -V sample1.g.vcf.gz \
  -O sample1.raw.vcf.gz
```

Index the VCF:
```bash
gatk IndexFeatureFile -I sample1.raw.vcf.gz
```

---

## 9) Phase 5 — Single-sample filtering (hard filters)

> For single samples, **VQSR** is often not feasible; use hard filters as a pragmatic alternative.
> If you later have a cohort, switch to VQSR.

### 9.1 Split SNPs and INDELs
```bash
gatk SelectVariants \
  -R Homo_sapiens_assembly38.fasta \
  -V sample1.raw.vcf.gz \
  --select-type-to-include SNP \
  -O sample1.snps.vcf.gz

gatk SelectVariants \
  -R Homo_sapiens_assembly38.fasta \
  -V sample1.raw.vcf.gz \
  --select-type-to-include INDEL \
  -O sample1.indels.vcf.gz
```

### 9.2 Filter SNPs (typical GATK hard filters)
```bash
gatk VariantFiltration \
  -R Homo_sapiens_assembly38.fasta \
  -V sample1.snps.vcf.gz \
  --filter-name "SNP_QD2" --filter-expression "QD < 2.0" \
  --filter-name "SNP_FS60" --filter-expression "FS > 60.0" \
  --filter-name "SNP_SOR3" --filter-expression "SOR > 3.0" \
  --filter-name "SNP_MQ40" --filter-expression "MQ < 40.0" \
  --filter-name "SNP_MQRankSum-12.5" --filter-expression "MQRankSum < -12.5" \
  --filter-name "SNP_ReadPosRankSum-8" --filter-expression "ReadPosRankSum < -8.0" \
  -O sample1.snps.filtered.vcf.gz
```

### 9.3 Filter INDELs (typical GATK hard filters)
```bash
gatk VariantFiltration \
  -R Homo_sapiens_assembly38.fasta \
  -V sample1.indels.vcf.gz \
  --filter-name "INDEL_QD2" --filter-expression "QD < 2.0" \
  --filter-name "INDEL_FS200" --filter-expression "FS > 200.0" \
  --filter-name "INDEL_SOR10" --filter-expression "SOR > 10.0" \
  --filter-name "INDEL_ReadPosRankSum-20" --filter-expression "ReadPosRankSum < -20.0" \
  -O sample1.indels.filtered.vcf.gz
```

### 9.4 Keep PASS variants only
```bash
gatk SelectVariants \
  -V sample1.snps.filtered.vcf.gz \
  --exclude-filtered true \
  -O sample1.snps.PASS.vcf.gz

gatk SelectVariants \
  -V sample1.indels.filtered.vcf.gz \
  --exclude-filtered true \
  -O sample1.indels.PASS.vcf.gz
```

### 9.5 Merge back (PASS-only)
```bash
gatk MergeVcfs \
  -I sample1.snps.PASS.vcf.gz \
  -I sample1.indels.PASS.vcf.gz \
  -O sample1.PASS.merged.vcf.gz
```

---

## 10) Phase 6 — Normalize the VCF (recommended before annotation)

This step improves downstream annotation consistency by splitting multi-allelics and normalizing indels.

```bash
bcftools norm \
  -f Homo_sapiens_assembly38.fasta \
  -m -both \
  sample1.PASS.merged.vcf.gz \
  -Oz -o sample1.PASS.normalized.vcf.gz

tabix -p vcf sample1.PASS.normalized.vcf.gz
```

---

## 11) Phase 7 — ANNOVAR annotation

Set your ANNOVAR DB path:
```bash
export ANNOVAR_DB="/path/to/annovar/humandb"
```

Annotate:
```bash
perl /path/to/annovar/table_annovar.pl \
  sample1.PASS.normalized.vcf.gz $ANNOVAR_DB \
  -buildver hg38 \
  -out sample1.annotated \
  -remove \
  -protocol refGene,gnomad211_exome,dbnsfp42a,clinvar_20221231 \
  -operation g,f,f,f \
  -nastring . \
  -vcfinput
```

---

## 12) Phase 8 — Example clinical filtering (Python/Pandas)

> Adjust thresholds per indication and lab SOP.

```python
import pandas as pd

df = pd.read_csv("sample1.annotated.hg38_multianno.txt", sep="\t")

# Coding + splice
df = df[df["Func.refGene"].isin(["exonic", "splicing"])]

# Rare variants (example: gnomAD popmax < 0.1%)
df["AF_popmax"] = pd.to_numeric(df.get("AF_popmax", 0), errors="coerce").fillna(0)
df = df[df["AF_popmax"] < 0.001]

# Remove synonymous
df = df[df["ExonicFunc.refGene"] != "synonymous SNV"]

# Example damaging threshold (CADD > 20)
df["CADD_phred"] = pd.to_numeric(df.get("CADD_phred", 0), errors="coerce").fillna(0)
df = df[df["CADD_phred"] >= 20]

# ClinVar pathogenic/likely pathogenic (example)
clnsig = df.get("CLNSIG", pd.Series([""] * len(df)))
df = df[clnsig.str.contains("Pathogenic|Likely_pathogenic", case=False, na=False)]

df.to_csv("sample1.final_clinical_candidates.csv", index=False)
print("Final candidates:", len(df))
```

---

## 13) Notes for clinical interpretation (germline)

Germline reporting commonly follows ACMG/AMP classification:
- Pathogenic / Likely pathogenic / VUS / Likely benign / Benign

Recommended downstream additions (outside scope of core calling):
- phenotype-driven prioritization (HPO terms)
- ACMG rules engine (e.g., InterVar) or clinical interpretation platform
- confirmatory orthogonal testing where required

---

## 14) Conceptual workflow

FASTQ  
↓ FastQC  
↓ BWA-MEM (+ RG)  
↓ Sort + MarkDuplicates  
↓ BQSR (dbSNP + Mills)  
↓ HaplotypeCaller (VCF or GVCF→GenotypeGVCFs)  
↓ SNP/INDEL hard filters (single sample)  
↓ PASS-only + normalize  
↓ ANNOVAR  
↓ clinical filtering & interpretation
