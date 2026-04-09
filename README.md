# Metagenomics Analysis Notebook

---

## 1. Abundance Table
- **Definition**: Matrix of counts (OTUs/ASVs × samples).
- **Statistical tests**:
  - Parametric: t-test, ANOVA
  - Non-parametric: Wilcoxon, Kruskal–Wallis
  - Multiple testing correction: FDR, Bonferroni
- **Differential abundance methods**:
  - DESeq2 (negative binomial model)
  - ANCOM (accounts for compositionality)
  - ALDEx2 (CLR-based Bayesian approach)

Reference: Love et al. (2014), *Genome Biology*; Mandal et al. (2015), *Microbiome*

---

## 2. Core Challenges
1. **Sparsity**  
   - Many zero counts due to rare taxa.  
   - Requires zero-inflated models or filtering.

2. **Dimensionality**  
   - Thousands of features vs. limited samples.  
   - Use dimensionality reduction (PCA, PCoA, NMDS).

3. **Compositionality**  
   - Data are relative, not absolute.  
   - Use log-ratio transformations (CLR, ALR).

Reference: Gloor et al. (2017), *Frontiers in Microbiology*

---

## 3. Batch Effects
- **PCR-based studies**: primer bias, annealing temperature, PCR cycles, post-PCR processes.
- **Shotgun sequencing**: less PCR bias, but still subject to library prep and sequencing run effects.
- **Correction**: ComBat (empirical Bayes).

| | Before correction | After correction |
|---|---|---|
| A | 50 | 50 |
| B | 35 | 50 |

Reference: Johnson et al. (2007), *Biostatistics*

---

## 4. Data Transformation
- Raw counts
- Relative abundance
- Log(x+1)
- CLR: `log(x / gm(x))`
- Rarefaction (controversial)

Reference: McMurdie & Holmes (2014), *PLoS Computational Biology*

---

## 5. Rarefaction
- **Pros**: balances sequencing depth, useful for alpha diversity.
- **Cons**: discards reads, reduces statistical power, biases rare taxa detection.
- **Alternatives**: normalization, variance-stabilizing transformations.

---

## 6. Alpha Diversity (within-sample)
- **Richness**: observed species, Chao1 estimator.
- **Evenness**: Pielou’s evenness, Shannon index (H’).
- **Phylogenetic diversity**: Faith’s PD.
- **Tests**: parametric (t-test) or non-parametric depending on distribution.

| | Before correction | After correction |
|---|---|---|
| A | Low | Moderate | High |
| B | Low | Moderate | High |

Reference: Chao (1984), *Scand J Stat*; Faith (1992), *Biological Conservation*

---

## 7. Beta Diversity (between-sample)
- **Metrics**: Bray–Curtis, Jaccard, UniFrac.
- **Ordination**: PCoA, NMDS.
- **Constrained ordination**:
  - RDA (linear, CLR-transformed data).
  - CCA (unimodal).
- **PERMANOVA**: tests group differences in community composition.

Reference: Anderson (2001), *Austral Ecology*

---

## 8. Differential Abundance
- **Challenges**: sparsity, compositionality, dimensionality, multiple testing.
- **Methods**: DESeq2, ANCOM, ALDEx2.
- **Goal**: identify taxa significantly associated with conditions.

---

## 9. Summary Formula
```
Observed data = Biological signal + Batch effect + Error
Healthy group = Healthy effect + Noise
Disease group = Disease effect + Noise
```

---

## 10. Workflow Diagram (Mermaid)

```mermaid
flowchart TD
    A[Raw Counts] --> B[Data Transformation]
    B --> C[Alpha Diversity]
    B --> D[Beta Diversity]
    B --> E[Differential Abundance]
    B --> F[Batch Effect Correction]
    F --> D
    D --> G[Ordination (PCoA, PCA, NMDS)]
    E --> H[Statistical Testing]
```
## 11. Practical Notes

Always check sequencing depth before analysis.

Normalize or transform data before ordination.

Use PERMANOVA for beta diversity hypothesis testing.

Apply multiple testing correction for differential abundance.

Document batch correction steps (ComBat, limma).
---
## 12. References
```mermaid
    - Love et al. (2014), Genome Biology
    - Mandal et al. (2015), Microbiome
    Gloor et al. (2017), Frontiers in Microbiology
    Johnson et al. (2007), Biostatistics
    McMurdie & Holmes (2014), PLoS Computational Biology
    Chao (1984), Scand J Stat
    Faith (1992), Biological Conservation
    Anderson (2001), Austral Ecology
```
---

# Quick Summary of Steps
```
Step 0: Load/install libraries.
Step 1: Import OTU, taxonomy, metadata, tree → build phyloseq object.
Step 2: Compute alpha diversity indices, visualize, test group differences.
Step 3: Compute beta diversity distances, ordination, PERMANOVA.
Step 4: Run constrained ordination (CCA) to link metadata with OTU distribution.
Step 5: Perform differential abundance analysis with ALDEx2, visualize volcano plot.
Step 6: Build and analyze microbial co-occurrence networks.
```

