# LRCDE ReadMe

## Getting started guide to using the LRCDE package

### Install the LRCDE package from GitHub.

I recommend installing the "devtools" package by Hadley Wickham.

The devtools package contains a function which allows painless installation of R packages from GitHub:

I recommend installing the package "knitr" in order to build the vignettes contained in the package. The vignettes provide getting started instructions for the LRCDE package and will help demystify the operation and features that LRCDE contains.

Just issue the following command:

```{r}
devtools::install_github("ERGlass/lrcde.dev", build_vignettes=TRUE)
```

Note that LRCDE depends upon the `dplyr` and `pROC` packages.

### Test the LRCDE installation.

The minimum required input to the LRCDE function is:
1. `het.matrix` - A heterogeneous matrix of genomic measures, sample by gene (rows by columns).
2. `cell.props` - A matrix of measured relative cell proportions, sample by cell type.
3. `group` - A group membership vector (1's and 2's) indicating controls and cases.

See `?lrcde` for examples.

After you have setup your data matrices `het.matrix` and `cell.props` and the `group` membership vector, issue the command:

```{r}
results.frame = lrcde( het.matrix, cell.props, groups )
```

`lrcde` returns a `data.frame`, and saves it into a file. The default output file will be `LRCDE_power_analysis.csv`, unless a different `output.file` is specified.

Remaining input parameters to the lrcde function are already set to the recommended defaults.

## Quick Start

The "Using LRCDE" vignette will provide a quick tutorial and help simulate a small data set to test the LRCDE package installation.

Type:

```{r}
 browseVignettes('lrcde')
```

Read the "[Using LRCDE](vignettes/using_lrcde.Rmd)" vignette and use the code examples to run your own quick simulation.

The "LRCDE Permutations for AUC Variability" is a more elaborate demonstration.  Caution: The permutations will take about 20 minutes or so to run. This is the code which demonstrates AUROC variability when condition number (kappa) of the cell proportions predictor matrix is very high. The plotting facility in the permutations vignette depends upon the package `Hmisc`.
