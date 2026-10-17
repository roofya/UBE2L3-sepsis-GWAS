#!/bin/bash
# ============================================================
# REGENIE Sepsis Case/Control Association Analysis
# 210 ICU patients | 68 cases | 142 controls
# Covariates: Age, Sex, BMI, SOFA, APACHE II, ICU days
# ============================================================
set -euo pipefail

# ============================================================
# PATHS — adjust if needed
# ============================================================
VCFDIR=/local/storage/Precision_Pharmacogenomics/WGS/Robay_Crystal/allBatch/vcf
OUTDIR=/local/storage/Precision_Pharmacogenomics/WGS/Robay_Crystal/allBatch/regenie
PHENO=${OUTDIR}/sepsis_pheno.txt
COV=${OUTDIR}/sepsis_covariates.txt
VCF=${VCFDIR}/allsample_phased_dbnsfp.vcf.gz
REGENIE=/programs/regenie_v3.2.9/regenie          # adjust path to your REGENIE binary
PLINK2=/programs/plink2_linux_20260504/plink2     # adjust path to your plink2 binary



# ============================================================
# STEP 0 — Copy phenotype and covariate files to OUTDIR
# ============================================================
# (already prepared — upload sepsis_pheno.txt and sepsis_covariates.txt to OUTDIR)

echo "======================================================"
echo "STEP 1: Prepare PLINK2 genotype files for Step 1"
echo "======================================================"
# REGENIE Step 1 needs PLINK2 pgen/psam/pvar format
# Filter: common variants only (MAF>5%), biallelic SNPs, max missingness 10%
# These are used ONLY for the null model (Step 1), not for association (Step 2)

${PLINK2} \
    -vcf ${VCFDIR}/allsample_phased_dbnsfp.vcf.gz \
    --update-sex ${OUTDIR}/sex_file.txt \
    --split-par hg38 \
    --set-all-var-ids '@:#:$r:$a' \
    --rm-dup exclude-mismatch \
    --maf 0.05 \
    --max-alleles 2 \
    --min-alleles 2 \
    --snps-only just-acgt \
    --geno 0.10 \
    --hwe 1e-6 \
    --make-pgen \
    --out ${OUTDIR}/geno_step1 \
    --threads 8

# Check how many variants remain for Step 1
echo "Variants for Step 1 null model:"
wc -l ${OUTDIR}/geno_step1.pvar

echo "======================================================"
echo "STEP 2: LD pruning for Step 1 null model"
echo "======================================================"
# LD pruning reduces collinearity in the null model
${PLINK2} \
    --pfile ${OUTDIR}/geno_step1 \
    --indep-pairwise 1000 100 0.9 \
    --out ${OUTDIR}/geno_step1_pruned \
    --threads 8

${PLINK2} \
    --pfile ${OUTDIR}/geno_step1 \
    --extract ${OUTDIR}/geno_step1_pruned.prune.in \
    --make-pgen \
    --out ${OUTDIR}/geno_step1_pruned_final \
    --threads 8

echo "Variants after LD pruning:"
wc -l ${OUTDIR}/geno_step1_pruned_final.pvar

echo "======================================================"
echo "STEP 3: REGENIE Step 1 — Fit null model"
echo "======================================================"
# Step 1 fits a whole-genome regression model to account for
# population structure and relatedness
# --bt = binary trait (case/control)
# --lowmem = memory-efficient mode for large datasets

${REGENIE} \
    --step 1 \
    --pgen ${OUTDIR}/geno_step1_pruned_final \
    --phenoFile ${PHENO} \
    --covarFile ${COV} \
    --covarColList AGE,SEX,BMI,SOFA,APACHE2,ICU_DAYS \
    --bt \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix ${OUTDIR}/tmp_regenie_step1 \
    --out ${OUTDIR}/regenie_step1 \
    --threads 8

echo "Step 1 complete. Files:"
ls -lh ${OUTDIR}/regenie_step1*

echo "======================================================"
echo "STEP 4: Prepare full VCF for Step 2 (all variants)"
echo "======================================================"
# Step 2 tests ALL variants for association
# We keep rare variants too (for gene-based tests)
# Only filter extreme missingness
${PLINK2} \
    --vcf ${VCF} \
    --max-alleles 2 \
    --geno 0.10 \
    --make-pgen \
    --out ${OUTDIR}/geno_step2_all \
    --threads 8

echo "======================================================"
echo "STEP 5: REGENIE Step 2 — Single variant association"
echo "======================================================"
# Tests each variant individually for association with sepsis
# --firth = Firth logistic regression (handles case/control imbalance, 68/142)
# --approx = approximate Firth for speed (use for genome-wide)
# --pThresh = only apply Firth for variants with p < 0.05

${REGENIE} \
    --step 2 \
    --pgen ${OUTDIR}/geno_step2_all \
    --phenoFile ${PHENO} \
    --covarFile ${COV} \
    --covarColList AGE,SEX,BMI,SOFA,APACHE2,ICU_DAYS \
    --bt \
    --firth --approx \
    --pThresh 0.05 \
    --pred ${OUTDIR}/regenie_step1_pred.list \
    --bsize 400 \
    --minMAC 3 \
    --out ${OUTDIR}/regenie_step2_single \
    --threads 8

echo "Step 2 single-variant complete."
ls -lh ${OUTDIR}/regenie_step2_single*

echo "======================================================"
echo "STEP 6: REGENIE Step 2 — Gene-based burden tests"
echo "======================================================"
# Gene-based tests are critical for rare variants in 210 samples
# Requires annotation file mapping variants to genes
# We use the ANN field from SnpEff to build this

# First build annotation file from your annotated VCF
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' \
    ${VCF} | \
    awk -F'\t' '{
        split($5, ann, "|");
        gene = ann[4];
        if (gene != "" && gene != ".") {
            print $1"_"$2"_"$3"_"$4"\t"gene
        }
    }' > ${OUTDIR}/variant_gene_annotations.txt

echo "Annotation file built:"
wc -l ${OUTDIR}/variant_gene_annotations.txt
head -5 ${OUTDIR}/variant_gene_annotations.txt

# Gene-based burden test (SKAT, ACAT-V, BURDEN, ACAT-O)
${REGENIE} \
    --step 2 \
    --pgen ${OUTDIR}/geno_step2_all \
    --phenoFile ${PHENO} \
    --covarFile ${COV} \
    --covarColList AGE,SEX,BMI,SOFA,APACHE2,ICU_DAYS \
    --bt \
    --firth --approx \
    --pred ${OUTDIR}/regenie_step1_pred.list \
    --anno-file ${OUTDIR}/variant_gene_annotations.txt \
    --set-list ${OUTDIR}/variant_gene_annotations.txt \
    --aaf-file ${OUTDIR}/variant_gene_annotations.txt \
    --vc-tests skato,acato,burden \
    --aaf-bins 0.01,0.05 \
    --bsize 200 \
    --minMAC 1 \
    --out ${OUTDIR}/regenie_step2_burden \
    --threads 8

echo "======================================================"
echo "STEP 7: Filter and summarize results"
echo "======================================================"
# Genome-wide significance threshold: p < 5e-8 (single variant)
# Gene-based significance: p < 2.5e-6 (Bonferroni for ~20k genes)

echo "=== Top single-variant hits (p < 1e-5) ==="
awk 'NR==1 || $12 < 1e-5' ${OUTDIR}/regenie_step2_single_PHENO.regenie \
    | sort -k12,12g \
    > ${OUTDIR}/top_single_variants.txt

wc -l ${OUTDIR}/top_single_variants.txt
head -20 ${OUTDIR}/top_single_variants.txt

echo "=== Top gene-based hits (p < 0.001) ==="
for test in skato acato burden; do
    echo "--- ${test} ---"
    awk -v t=$test 'NR==1 || $12 < 0.001' \
        ${OUTDIR}/regenie_step2_burden_${test}_PHENO.regenie \
        | sort -k12,12g \
        | head -20
done

echo "======================================================"
echo "ALL DONE at $(date)"
echo "Key output files:"
echo "  Single-variant: ${OUTDIR}/regenie_step2_single_PHENO.regenie"
echo "  Gene-based:     ${OUTDIR}/regenie_step2_burden_*.regenie"
echo "  Top hits:       ${OUTDIR}/top_single_variants.txt"
echo "======================================================"
