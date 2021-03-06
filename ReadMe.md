LRCDE package performs cell type-specific differential expression analysis provided three components:

1. `het.matrix` - A heterogeneous matrix of genomic measures, sample by gene (rows by columns).
2. `cell.props` - A matrix of measured relative cell proportions, sample by cell type. It can be estimated using cell signatures using the `decon.cell.props` function
3. `group` - A group membership vector (1's and 2's) indicating controls and cases.

LRCDE functionality applies to any type of high-throughput data other than gene expression. E.g., methylation profiles (M-values recommended), RNA-seq data (FPKM or TPM values recommended).

## Getting started with the LRCDE package

Install the LRCDE package from GitHub.

```{r}
devtools::install_github("ERGlass/lrcde.dev", build_vignettes=TRUE)
```
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
- "[GEO65391_deconvolve_proportions](vignettes/GEO65391_deconvolve_proportions.Rmd)" - how to estimate cell proportions from cell signatures
- "[random_data_simulation_lrcde](vignettes/random_data_simulation_lrcde.Rmd)" - running `lrsde` on completely random data
- "[replicating_cssam_results](vignettes/replicating_cssam_results.Rmd)" - replicating csSAM results published in [https://www.ncbi.nlm.nih.gov/pubmed/20208531](https://www.ncbi.nlm.nih.gov/pubmed/20208531)

The "[LRCDE Permutations for AUC Variability](vignettes/lrcde.permutations_for_auc_variability.Rmd)" is a more elaborate demonstration.  Caution: The permutations will take about 20 minutes or so to run. This is the code which demonstrates AUROC variability when condition number (kappa) of the cell proportions predictor matrix is very high. The plotting facility in the permutations vignette depends upon the package `Hmisc`.

## Functions

- `lrcde` - main function to perform cell type-specific differential expression analysis.
- `decon.cell.props` - a function to estimate cell proportions given heterogeneous gene expression matrix and cell type-specific gene expression signatures.

## Internal functions

- `do.dual.decon` - a function to actually perform two-group cell type-specific differential expression analysis. Used by the main `lrcde` function
- `cell.props.target` - a function to simulate cell proportion matrix with pre-defined cell proportion standard deviation for the target cell type, and the condition number for the whole matrix
- `custom.resid.synthetic` - a function to create a matrix of simulated residuals for each sample and cell type

## Misc 

- Frequencies of Cell Types in Human Peripheral Blood, [https://www.stemcell.com/media/files/wallchart/WA10006-Frequencies_Cell_Types_Human_Peripheral_Blood.pdf](https://www.stemcell.com/media/files/wallchart/WA10006-Frequencies_Cell_Types_Human_Peripheral_Blood.pdf)
