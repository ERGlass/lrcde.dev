## Getting started with the LRCDE package

Install the LRCDE package from GitHub.

```{r}
devtools::install_github("ERGlass/lrcde.dev", build_vignettes=TRUE)
```

LRCDE package performs cell type-specific differential expression analysis provided three components:

1. `het.matrix` - A heterogeneous matrix of genomic measures, sample by gene (rows by columns).
2. `cell.props` - A matrix of measured relative cell proportions, sample by cell type. It can be estimated using cell signatures using the `decon.cell.props` function
3. `group` - A group membership vector (1's and 2's) indicating controls and cases.

See `?lrcde` for real life examples. Run the simulation example in the "[Using LRCDE](vignettes/using_lrcde.Rmd)" vignette.

```{r}
results.frame = lrcde( het.matrix, cell.props, groups )
```

`lrcde` returns a list with the first element containing a `data.frame` of results. The results are saved in `LRCDE_power_analysis.csv` file, unless a different `output.file` is specified.

## Getting help

```{r}
 browseVignettes('lrcde')
```

- "[Using LRCDE](vignettes/using_lrcde.Rmd)" - how to run a quick simulation

The "LRCDE Permutations for AUC Variability" is a more elaborate demonstration.  Caution: The permutations will take about 20 minutes or so to run. This is the code which demonstrates AUROC variability when condition number (kappa) of the cell proportions predictor matrix is very high. The plotting facility in the permutations vignette depends upon the package `Hmisc`.

## Functions

- `lrcde` - main hunction to perform cell type-specific differential expression analysis.