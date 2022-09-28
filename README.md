## Dozer
**Dozer** is a lightweight R package for personalized co-expression network analysis from population-scale single cell RNA-sequencing datasets. Given a single cell RNA-seq dataset of cells from multiple subjects and phenotypic information of each subject, Dozer is able to construct subject specific co-expression networks, extract network traits, and conduct association studies between network traits and subject specific information.

**Dozer** can be downloaded and installed in R by: 

```r
## install.packages("devtools")
devtools::install_github("shanlu01/Dozer", build_vignettes = TRUE)
```
As a demonstration of the analysis pipeline, we included an example dataset and a vignette in our package. The featured analyses include identifying genes showing differential centrality and gene modules having difference in connectivity between phenotypic groups.  The dataset and vignette can be access through:
```r
load(system.file("extdata", "Jerber_demo.rda", package = "Dozer"))
browseVignettes("Dozer")
```
## Reference
S. Lu and S. Keles, "Dozer: Debiased personalized gene co-expression networks for population-scale scRNA-seq data".
