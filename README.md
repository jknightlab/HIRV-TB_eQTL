# HIRV-TB_eQTL

## Overview

This repository contains scripts for the paper:

“**Common genetic determinants of dynamic human *in vivo* immune response to Mycobacterium tuberculosis.**”

Ping Zhang, Carolin T Turner, Aneesh Chandran, Joshua Rosenheim, Jana Jiang, Lucy K Bell, Santino Capocci, Marc Lipman, Heinke Kunst, Stefan Lozewicz, Gillian S Tomlinson, Julian C Knight and Mahdad Noursadeghi

Abstract

We are yet to unlock sufficient understanding of the immunological determinants of protection and pathogenesis in tuberculosis (TB) to stratify disease-risk following infection or develop an effective vaccine. Given evidence of common genetic susceptibility to TB, we sought to identify single nucleotide polymorphisms (SNPs) affecting transcriptional regulation of immune responses to Mycobacterium tuberculosis (Mtb) in order to relate immune response variation to disease-risk. We undertook a genome-wide expression quantitative trait loci (eQTL) study of human in vivo immune responses at transcriptional level to a standardised immunogenic challenge using the tuberculin skin test (TST), which recapitulates the molecular repertoire of cell mediated immunity at the site of TB disease. 267 individuals with latent or active TB underwent genotyping paired with 415 RNAseq of TST biopsies at day 2 to capture maximal inflammatory responses, and/or day 7 to capture proliferation/enrichment of antigen-specific T cells. We identified common cis-acting SNPs for 1719 of 5526 TST response genes, most enriched for antigen presentation and T cell activation pathways due to SNPs associated with HLA class II expression levels. We found associations between non-coding QTLs for HLA-DR expression with expansion of HLA-DR restricted Mtb-reactive T cells in day 7 TSTs, suggesting a potential interaction with coding variants to modulate antigen-specific responses. Interestingly, the greatest number of trans-acting QTLs were identified in day 7 TST data and enriched for cell-cycle/proliferation pathways, via a cis-effect on NCAPD3 (involved in mitotic chromosome assembly and segregation) that also associated with Mtb-reactive T cell responses. Co-localisation analysis with GWAS data revealed potential new disease susceptibility loci. We have established a catalogue of common genetic polymorphisms that influence variation of human in vivo immune responses, to enable investigation of novel mechanistic insights into the immunological determinants of TB. Variation in antigen-specific T cell proliferation may be regulated by coding/non-coding SNPs in HLA class II genes, and by non-coding SNPs that regulate cell proliferation machinery independent of antigen.


## Pipeline

| Step    | Script                                            | Description                                                                         |
| ------- | ------------------------------------------------- | ----------------------------------------------------------------------------------- |
| **1.1** | `vcf_for_imputation.sh`                           | Prepare VCF for imputation by filtering and splitting per chromosome                |
| **1.2** | `for.PCA.sh`                                      | Filter imputed VCFs and convert to PLINK format for PCA                             |
| **1.3** | `pca.sh`                                          | Perform PCA on merged genotype data with 1000 Genomes reference data                |
| **2.1** | `hisat2_mapping.sh`                               | RNA-seq read mapping using HISAT2                                                   |
| **2.2** | `featurecount.sh`                                 | Generate read counts using featureCounts                                            |
| **3.1** | `dedup_rg_reorder.sh`                             | Remove duplicates, add read groups, and reorder BAM files                           |
| **3.2** | `crosscheck.sh`                                   | Perform cross-checking of DNA and RNA samples using Picard’s CrosscheckFingerprints |
| **4.1** | `MatrixeQTL_cis.R`                                | Run MatrixEQTL using different expression PCs                                       |
| **4.2** | `trans-day7.PC25_subset.SNPs.R`                   | trans-eQTL mapping using day7 lead cis-eQTLs                                        |
| **5.1**   | `run_fastqtl.sh`                                  | Run FastQTL for permutations and nominal passes                                     |
| **6.1** | `lmm.R`                                           | LMM model for cis-eQTL mapping                                                      |
| **6.2** | `eigenMT-input.r`                                 | Prepare eigenMT input                                                               |
| **6.3** | `eigenMT_python3.py`                              | Python eigenMT script                                                               |
| **6.4** | `run_eigenMT_array.sh`                            | SLURM array script to run eigenMT per chromosome                                    |
| **6.5** | `eigenMT-qvalue-significance.threshold.r`         | Compute q-values for eigenMT results and define significance thresholds             |
| **6.6** | `conditional.analysis_lmerTest-array-ALL.pairs.r` | Stepwise forward/backward conditional mapping using lmerTest                        |
| **7.1** | `Coloc_run_arrayR-eQTL.r`                         | Colocalisation analysis for eQTL Catalogue vs TB summary statistics                 |
| **7.2** | `Coloc_SGE_array_LMM.resul_with.GWAS_opti.r`      | Colocalisation between TB eQTLs and GWAS summary statistics                         |
| **8.1** | `LD Score (ldsc) regression analysis`             | Partitioned heritability estimation                                                 |
| **9.1** | `interaction.analysis-Day7vsDay2.r`               | SNP × TB_Status (Day7 vs Day2) interaction analysis                                 |
| **9.2** | `interaction.analysis_cell.module.r`              | SNP x cell module interaction analysis                                              |
---

** directory layout:**

```
.
├── README.md
├── notebooks
│   ├── 1.TB_DESeq2_424.TB.Rmd
│   ├── 1.TB_DESeq2_424.TB.html
│   ├── 10.DeepSEA-Sei_LMM_results.Rmd
│   ├── 10.DeepSEA-Sei_LMM_results.html
│   ├── 2.Crosscheck_ALL.Rmd
│   ├── 2.Crosscheck_ALL.html
│   ├── 3.1KG_PCAs_267-freeze.TB.html
│   ├── 3.1KG_PCAs_267-freeze.TBs.Rmd
│   ├── 4.individual_eQTLs.plot_lmm_interaction.Rmd
│   ├── 4.individual_eQTLs.plot_lmm_interaction.html
│   ├── 5.Manhattan_LMM.Rmd
│   ├── 5.Manhattan_LMM.html
│   ├── 6.eGene.venn.and.XGR.pathway.Rmd
│   ├── 6.eGene.venn.and.XGR.pathway.html
│   ├── 7.coloc_TB.results.Rmd
│   ├── 7.coloc_TB.results.html
│   ├── 8.circular.plot.Rmd
│   ├── 8.circular.plot.html
│   ├── 9.ldsc.results.plot-opti.Rmd
│   └── 9.ldsc.results.plot-opti.html
└── scripts
    ├── 1.1.split.vcf.for.imputation.sh
    ├── 1.2.filter.and.vcf.to.plink.format_MIS_imputed.sh
    ├── 1.3.pca.with.1KG.data.sh
    ├── 2.1.hisat2.mapping.sh
    ├── 2.2.featureCounts.sh
    ├── 3.1.Crosscheck.RG.and.sort.sh
    ├── 3.2.Crosscheck.DNA.RNA.sh
    ├── 4.1.MatrixeQTL_cis.r
    ├── 4.2.trans-day7.PC25_subset.SNPs.r
    ├── 5.1.run.fastQTL.sh
    ├── 6.1.lmm.r
    ├── 6.2.eigenMT-input.r
    ├── 6.3.eigenMT_python3.py
    ├── 6.4.run_eigenMT_array.sh
    ├── 6.5.eigenMT-qvalue-significance.threshold.r
    ├── 6.6.conditional.analysis_lmerTest-array-ALL.pairs.r
    ├── 7.1.Coloc_SGE_array_LMM.resul_eQTL.R
    ├── 7.2.Coloc_SGE_array_LMM.resul_with.GWAS_opti.R
    ├── 8.1.LD Score (ldsc) regression analysis.txt
    ├── 9.1.interaction.analysis-Day7vsDay2.r
    └── 9.2.interaction.analysis_cell.module.r
3 directories, 42 files

```


---

## Environment & dependencies

### R (≥ 4.2.1)

| Package | Purpose |
|----------|----------|
| **argparse** | Command-line argument parsing for R scripts |
| **data.table** | Fast file I/O and manipulation |
| **tidyverse** | Core data wrangling (`dplyr`, `tidyr`, `ggplot2`, `readr`) |
| **reshape** | Parse variant IDs (`colsplit`) |
| **MatrixEQTL** | cis/trans-eQTL mapping |
| **lme4**, **lmerTest** | Linear mixed models and p-values |
| **foreach**, **doMC** | Parallel loops for SNP/gene scans |
| **coloc** | Bayesian colocalisation |
| **seqminer** | Tabix-based read of summary statistics |
| **GenomicRanges** | Interval operations for eQTLs |
| **ggplot2**, **circlize**, **scales** | Visualisation |

---

### Python & command-line tools

| Tool | Version (example) | Purpose |
|------|-------------------|---------|
| **FastQTL** | ≥ 2.184 | cis/trans eQTL discovery |
| **HISAT2** | ≥ 2.1.0 | RNA-seq alignment |
| **featureCounts** | ≥ 1.6.4 | Gene quantification |
| **samtools** | ≥ 1.9 | BAM processing |
| **bcftools** | ≥ 1.10 | VCF filtering / annotation |
| **htslib** (bgzip, tabix) | ≥ 1.8 | Compression / indexing |
| **Picard Tools** | ≥ 2.21.1 | BAM deduplication, read-group assignment |
| **eigenMT** | Python 3.12.4 | Conditional cis-eQTL fine-mapping |
| **ldsc** | Python 2.7.15 | Estimating heritability and genetic correlation from GWAS summary statistics |
| **conda** | v23.7.2 | Package & environment management |

---

### HPC environment

| Tool | Purpose |
|------|----------|
| **SBATCH (SLURM)** | Array & batch job submission |
| **module load / purge** | HPC environment setup |

---

## ✉️ Contact
For questions, please contact:  
**Ping Zhang**  (<ping.zhang@well.ox.ac.uk>)
