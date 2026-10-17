# UBE2L3 Sepsis GWAS

**Whole-genome sequencing GWAS identifies Arab founder loss-of-function variant in *UBE2L3* associated with sepsis susceptibility**

[![medRxiv](https://img.shields.io/badge/medRxiv-2026.09.25.26363954-red)](https://doi.org/10.64898/2026.09.25.26363954)
[![GWAS Catalog](https://img.shields.io/badge/GWAS_Catalog-GCST91000424-blue)](https://www.ebi.ac.uk/gwas/studies/GCST91000424)

## Authors
Rostami MR, Hssain AA, Crystal RG, Rodriguez-Flores JL  
Weill Cornell Medicine | Precision Pharmacogenomics Programme

## Overview
Genome-wide association study of sepsis susceptibility using whole-genome
sequencing (WGS) in 210 Arab ICU patients (68 sepsis cases, 142 controls)
from Qatar (GRCh38). We identified a rare Arab founder start-loss variant
in *UBE2L3* (chr22:21,567,746 T>G; c.2T>G; p.Met1?; CADD=25.6) present
exclusively in sepsis cases, with independent gene-level replication in
26,939 bacterial infection cases in the NIH All of Us Research Program
(SKAT-O P=2.38×10⁻³; OR=23.10 [95% CI 4.2–125.8]).

## Key Finding
| Variant | Gene | Effect | P-value | OR | AF |
|---|---|---|---|---|---|
| chr22:21,567,746 T>G | *UBE2L3* | start_lost (p.Met1?) | 0.019 | 20.03 | 0.007 |

- **3/68 cases vs 0/142 controls** (Fisher's exact P=0.026)
- Absent from gnomAD (251,000 individuals) and All of Us (700,000+)
- All 3 carriers: male, sepsis cases, SOFA ≥2
- Independent replication: All of Us SKAT-O P=2.38×10⁻³ (n=26,939)
- GWAS Catalog accession: **GCST91000424**

## Repository Contents
scripts/
├── regenie_sepsis_pipeline.sh   # Bash: WGS QC → phasing → annotation → REGENIE GWAS + burden
└── UBE2L3_sepsis_analysis.R    # R: Manhattan, QQ, LD plot, fgsea, gprofiler2, carrier analysis

data/
└── AllofUs_SKATO_replication.csv  # All of Us burden test summary (no individual-level data)

results/
└── README.md  → GWAS Catalog GCST91000424

## Methods Summary
- **WGS:** 30× GRCh38, Illumina
- **Phasing:** BEAGLE v5.1 (1000 Genomes reference)
- **Annotation:** SnpEff + dbNSFP v5.3.1a (CADD, SIFT, AlphaMissense, gnomAD v4.1)
- **GWAS:** REGENIE v3.2.9, Firth logistic regression, 16,598,149 variants
- **Burden test:** SKAT-O, 4 functional masks × 4 AAF bins
- **Replication:** All of Us Researcher Workbench (Registered Tier)

## Data Availability
GWAS summary statistics are deposited in the NHGRI-EBI GWAS Catalog
(accession **GCST91000424**) under embargo until publication.
Individual-level WGS data cannot be shared due to participant privacy.
All analysis code is provided here under the MIT License.
## Citation
Mahboubeh R. Rostami, Ali Ait Hssain, Amal Robay, Edward J. Schenck,
Jason J. Mezey, Ronald G. Crystal, Juan Rodriguez-Flores.
Genome-Wide Association Study of Sepsis Susceptibility in ICU Patients
Identifies Arab Founder Loss-of-Function Variant in UBE2L3.
medRxiv 2026. https://doi.org/10.64898/2026.09.25.26363954
## License
MIT License
