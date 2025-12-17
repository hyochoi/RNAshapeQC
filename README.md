# RNAshapeQC

**RNA** Coverage-**Shape**-Based **Quality Control** Metrics


## Description

_RNAshapeQC_ provides coverage-shape-based quality control (QC) metrics for mRNA-seq and total RNA-seq data. It supports per-gene pileup construction from BAM files as well as toy datasets for quick-start examples. The package implements protocol-specific metrics, including decay rate (DR), degradation score (DS), mean coverage depth (MCD), window coefficient of variation (wCV), area under the curve (AUC), and shape-based sample-level indices. 

_RNAshapeQC_ also includes HPC-friendly functions for per-gene batch processing and cross-study pileup generation. This package enables interpretable, protocol-specific QC assessments for diverse RNA-seq workflows.


## Documentation

Comprehensive documentation is available in the [_Tutorial_](https://miyeonyeon.github.io/bioc-vignettes/RNAshapeQC_intro.html) and the [_User Manual_](https://miyeonyeon.github.io/bioc-vignettes/RNAshapeQCdocs_book/).


## Installation

### From GitHub
```r
if (!requireNamespace("devtools", quietly=TRUE)) {
  install.packages("devtools")
}
devtools::install_github("hyochoi/RNAshapeQC", dependencies=TRUE)
```

### From Bioconductor
```r
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("RNAshapeQC", dependencies=TRUE)
```


## Citation

If you use _RNAshapeQC_ in your research, please cite the package using:

```r
citation("RNAshapeQC")
```

A manuscript describing the methodology and applications of _RNAshapeQC_ is currently in preparation. Citation information will be updated upon publication.
