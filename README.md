---
editor_options: 
  markdown: 
    wrap: 72
---

My analysis aims to apply Cox regression models to microbiome data with
different transformations to identify significant associations with
survival outcomes.

The analysis is based on a **`phyloseq` object**, which is a structured
R object used for microbiome data integration and analysis. It contains
multiple types of information commonly used in microbiome studies,
including:

1\. Install Required Packages

For handling microbiome data (OTU tables, taxonomy, sample data):

```{r}
install.packages("BiocManager")
BiocManager::install("phyloseq")
```

For compositional data analysis and penalized Cox regression models:

```{r}
install.packages("coda4microbiome")
```

For survival analysis functions

```{r}
install.packages("survival")
```

2-Source Custom Function:\
The script includes a custom function abund_coxnet2, which implements a
penalized Cox regression model on different data transformations (e.g.,
ILR, ALR, CLR, and relative abundance). This function is defined in a
separate script and needs to be sourced before use:

```{r}
source("abund_coxnet2.R")
```
