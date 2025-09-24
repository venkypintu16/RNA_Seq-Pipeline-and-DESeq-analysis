# RNA_Seq-Pipeline-and-DESeq-analysis
# A Simple, Cost-Effective, and Robust Method for rRNA Depletion in RNA-Sequencing Studies

This repository documents a two-phase project aimed at developing a **DIY rRNA depletion method** as a cost-effective and efficient alternative to the discontinued **Ribo-Zero kit (Illumina)**.  
The work spans from **pipeline setup and validation (Phase 1)** to **differential expression analysis and benchmarking (Phase 2)**.

---

## Background
- RNA sequencing (RNA-seq) enables global gene expression profiling.  
- rRNA constitutes ~90% of total RNA, making depletion essential.  
- Commercial kits (e.g., **Ribo-Zero**) are costly and discontinued, necessitating alternative approaches.  
- This project develops a **DIY kit** using biotinylated oligonucleotides against 23S, 16S, and 5S rRNAs with **magnetic streptavidin-coated beads** for depletion.  
- Cost comparison:  
  - **DIY Kit**: ~$10/reaction  
  - **Ribo-Zero Kit**: ~$80/reaction  

---

#  Phase 1: RNA_Seq Pipeline 
###  Objectives
- Develop a DIY rRNA depletion method.  
- Benchmark performance against Ribo-Zero across three bacterial species:  
  - **E. coli**  
  - **B. subtilis**  
  - **C. crescentus**

###  Materials & Tools
- **Datasets:** FASTQ files (GEO/EBI), FASTA, GFF/GTF from NCBI.  
- **Tools:**  
  - Bowtie2  
  - FastQC  
  - Samtools  
  - Subread (featureCounts)  
  - IGV Viewer  
  - R / RStudio (RPKM normalization, scatterplots, boxplots, histograms)  

###  Methods
1. Quality control of raw reads using **FastQC**.  
2. Alignment to reference genomes using **Bowtie2**.  
3. Conversion to BAM files with **Samtools**.  
4. Read quantification using **featureCounts (Subread)**.  
5. Normalization to **RPKM** and visualization in **R**.  
6. Comparison of DIY Kit vs Ribo-Zero efficiency.  

###  Results
- Reads mapped with high per-base quality.  
- DIY depletion yielded **75–80% mRNA reads**, comparable to Ribo-Zero.  
- Scatterplots and IGV confirmed effective depletion.  
- **Cost savings**: DIY kit significantly reduced per-sample cost.  

---

#  Phase 2: Analysis Part with R using DESeq package
###  Objectives
- Extend Phase 1 by incorporating **differential expression analysis**.  
- Confirm depletion efficiency of the DIY kit under different conditions.  
- Compare **DIY Kit vs Ribo-Zero** in control and antibiotic-treated **E. coli** samples.  

###  Materials & Tools
- **Data Sources:**  
  - RPKM and feature counts (Phase 1 outputs).  
  - E. coli datasets (control, rifampicin, chloramphenicol).  
- **Additional Tools:**  
  - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)  

###  Methods
1. Generate RPKM plots for depleted vs total RNA.  
2. Compare **DIY vs Ribo-Zero** in antibiotic-treated E. coli samples.  
3. Perform **DESeq2 differential expression analysis**:  
   - Fit negative binomial GLMs.  
   - Estimate dispersion and fold changes.  
   - Generate scatterplots, boxplots, heatmaps.  

###  Results
- **Phase 2 confirmed Phase 1 findings**: DIY kit effectively depleted rRNA.  
- Ribo-Zero showed higher read counts, but DIY kit exhibited **enhanced gene expression signal**.  
- **DESeq2** revealed treatment-specific differential expression.  
- Antibiotic comparisons: control vs rifampicin vs chloramphenicol consistent across both kits.  

###  Future Directions
- Optimize DIY chemistry to improve efficiency.  
- Compare DIY kit with **new market kits** (e.g., Zymo-Seq RiboFree).  
- Extend workflow to **metatranscriptomics** and **eukaryotic RNA-seq**.  

 
---

##  Reference
Culviner, P. H., Guegler, C. K., & Laub, M. T. (2020).  
*A Simple, Cost-Effective, and Robust Method for rRNA Depletion in RNA-Sequencing Studies.*  
mBio, 11(2), e00010-20. [https://doi.org/10.1128/mBio.00010-20](https://doi.org/10.1128/mBio.00010-20)

---

##  Example Workflows

### RNA-Seq Alignment and RPKM (Phase 1)
```bash
# Quality check
fastqc sample.fastq

# Alignment
bowtie2 -x reference -U sample.fastq -S sample.sam

# Convert to BAM
samtools view -bS sample.sam > sample.bam

# Feature counts
featureCounts -a annotation.gtf -o counts.txt sample.bam

