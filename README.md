# shiftR: Fast Enrichment Analysis via Circular Permutations

Fast enrichment analysis for locally correlated statistics
via circular permutations.
The analysis can be performed at multiple significance thresholds
for both primary and auxiliary data sets with 
with efficient correction for multiple testing.

## Installation

### Install CRAN Version

To install the
[CRAN version](https://CRAN.R-project.org/package=shiftR)
of `shiftR`, run

```
install.packages("shiftR")
```

### Install GitHub Version

To install `shiftR` directly from GitHub, run

```
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("andreyshabalin/shiftR")
```
