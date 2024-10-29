# RGlycopacity tool

```{r}
devtools::install_github('CopenhagenCenterForGlycomics/RGlycopacity')

# OR

remotes::install_github('CopenhagenCenterForGlycomics/RGlycopacity')
```

```{r}
library(Seurat)
pre_processed_rds = readRDS('my_seurat.Rds')
```

```{r}
library(RGlycopacity)
data_pseudobulk = calculatePseudobulks(pre_processed_rds)
expression_range_table = calculateExpressedRange(data_pseudobulks)
```
