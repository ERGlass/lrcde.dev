# LRCDE ReadMe

## Getting started guide to using the LRCDE package

### Install the LRCDE package from GitHub.

I recommend installing the "devtools" package by Hadley Wickham.

The devtools package contains a function which allows painless installation of R packages from GitHub:

Just issue the following command:
```{r}
install_github("ERGlass/lrcde.dev", build_vignettes=TRUE)
```

Note that LRCDE depends upon the dplyer and pROC packages.

### Test the LRCDE installation.

The minimum required input to the LRCDE function is:
1. A heterogeneous matrix of genomic measures - a sample size by genomic site (rows by columns): het.matrix
2. A matrix of measured relative cell proportions - sample size by cell type: cell.props
3. A group membership vector (1's and 2's) indicating controls and cases: groups

The default output file will be "LRCDE_power_analysis.csv" which will be written to the current working directory.

Remaining input parameters to the lrcde function are already set to the recommended defaults.

issue the command:
```{r}
results.frame = lrcde( het.matrix, cell.props, groups )
```

The results of lrcde will be in the results.frame data frame and in the LRCEE_power_analysis.csv file in the current working directory.

## Quick Start

type:
```{r}
 browseVignettes('lrcde')
```

Read the vignette and use the code examples to run your own quick simulation.

