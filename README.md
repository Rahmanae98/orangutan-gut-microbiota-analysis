# Microbiome Analysis of Wild Orangutans (Pongo abelii)
This repository presents a comprehensive R-based workflow for analyzing the gut microbiota of wild Sumatran orangutans (Pongo abelii) from the Suaq Balimbing Research Station, Gunung Leuser National Park, Indonesia. The project aims to understand how ecological and host-related factors influence microbial diversity in natural habitats.
## Background
Understanding the relationship between host ecology and gut microbiota is crucial in wildlife conservation. In this study, we explore how factors such as:
- Food availability (phenology and Fruit Availability Index)
- Climate variables (rainfall/temperature)
- Host traits (sex, age, class/social status)
affect the composition and diversity of orangutan gut microbiota. The dataset consists of gut microbiota profiles based on 16S rRNA amplicon sequencing, collected alongside ecological and behavioral observations.
## Analyses Overview
### Alpha Diversity
- Metrics: Shannon Index and ASV Richness
- Tests: Kruskal-Wallis to assess diversity differences by sex, class, and month
- Goal: Determine richness and evenness patterns among host groups
### Beta Diversity
- Dissimilarity: Bray-Curtis index
- Ordination: PCoA visualization
- Statistical Test: PERMANOVA (via `adonis2`) to test differences between groups
### Visualization
- Barplots, boxplots, and ordination plots created with `ggplot2`, `vegan`, `phyloseq`, and `microbiome` R packages.


