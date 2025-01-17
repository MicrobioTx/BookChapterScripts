MelonnPan - Model-based Genomically Informed High-dimensional Predictor of Microbial Community Metabolic Profiles
================
Himel Mallick
2023-02-14

-   [Introduction](#introduction)
-   [How to Install](#how-to-install)
    -   [From R](#from-r)
    -   [From the command line](#from-the-command-line)
-   [Usage](#usage)
    -   [Input](#input)
    -   [Output](#output)
-   [References](#references)
-   [Citation](#citation)

Introduction
------------

MelonnPan is a computational method for predicting metabolite compositions from microbiome sequencing data.

![Overview of MelonnPan](./Figure_0.jpg)

MelonnPan is composed of two high-level workflows: `MelonnPan-Predict` and `MelonnPan-Train`.

The `MelonnPan-Predict` workflow takes a table of microbial sequence features (i.e., taxonomic or functional abundances on a per sample basis) as input, and outputs a predicted metabolomic table (i.e., relative abundances of metabolite compounds across samples).

The `MelonnPan-Train` workflow creates an weight matrix that links an optimal set of sequence features to a subset of predictable metabolites following rigorous internal validation, which is then used to generate a table of predicted metabolite compounds (i.e., relative abundances of metabolite compounds per sample). When sufficiently accurate, these predicted metabolite relative abundances can be used for downstream statistical analysis and end-to-end biomarker discovery.

How to Install
--------------

There are two options for installing `MelonnPan`:

### From R

In R, you can install `MelonnPan` using the `devtools` package as follows (execute from within a fresh R session):

``` r
install.packages('devtools') # Install devtools if not installed already
library(devtools) # Load devtools
devtools::install_github("biobakery/melonnpan") # Install MelonnPan
```

### From the command line

Clone the repository using `git clone`, which downloads the package as its own directory called `melonnpan`.

``` bash
git clone https://github.com/biobakery/melonnpan.git
```

Then, install MelonnPan using `R CMD INSTALL`.

``` bash
R CMD INSTALL melonnpan
```

Usage
-----

MelonnPan can be run from the command line or from within R. Both methods require the same arguments, have the same options, and use the same default settings. Check out the [MelonnPan tutorial](https://github.com/biobakery/biobakery/wiki/melonnpan) for an example application.

-   The default `MelonnPan-Predict` function can be run by executing the script `predict_metabolites.R` from the command line or within R using the function `melonnpan.predict()`. Currently it uses a [pre-trained model](https://github.com/biobakery/melonnpan/blob/master/data/melonnpan.trained.model.txt) from the human gut based on UniRef90 gene families (functionally profiled by [HUMAnN2](http://huttenhower.sph.harvard.edu/humann2)), as described in Franzosa et al. (2019) and the original MelonnPan paper (Mallick et al., 2019), which is included in the package and can also be downloaded from the [`data/`](https://github.com/biobakery/melonnpan/blob/master/data) sub-directory (**melonnpan.trained.model.txt**).

-   If you have paired metabolite and microbial sequencing data (possibly measured from the same biospecimen), you can also train a MelonnPan model by running the script `train_metabolites.R` from the command line or within R using the function `melonnpan.train()`.

-   MelonnPan currently requires input data that is specified using UniRef90 gene families (functionally profiled by [HUMAnN2](http://huttenhower.sph.harvard.edu/humann2)). If you do not have functionally profiled UniRef90 gene families from the human gut or other environments, you may need to first train a MelonnPan model using the `MelonnPan-Train` workflow and supply the resulting weights to the `MelonnPan-Predict` module to get the relevant predictions.

### Input

-   `MelonnPan-Predict` workflow requires the following input:
    -   a table of microbial sequence features' relative abundances (samples in rows)
-   `MelonnPan-Train` workflow requires the following inputs:
    -   a table of metabolite relative abundances (samples in rows)
    -   a table of microbial sequence features' relative abundances (samples in rows)
-   For a complete description of the possible parameters for specific `MelonnPan` functions and their default values and output, run the help within R with the `?` operator.

### Output

-   The `MelonnPan-Predict` workflow outputs the following:
    -   **MelonnPan\_Predicted\_Metabolites.txt**: Predicted relative abundances of metabolites as determined by `MelonnPan-Predict`.
    -   **MelonnPan\_RTSI.txt**: Table summarizing RTSI scores per sample.
-   Similarly, the `MelonnPan-Train` workflow outputs the following:
    -   **MelonnPan\_Training\_Summary.txt**: Significant compounds list with per-compound prediction accuracy (correlation coefficient) and the associated p-value and q-value.
    -   **MelonnPan\_Trained\_Metabolites.txt**: Predicted relative abundances of statisticially significant metabolites as determined by `MelonnPan-Train`.
    -   **MelonnPan\_Trained\_Weights.txt**: Table summarizing coefficient estimates (weights) per compound.

References
----------

Zou H, Hastie T (2005). Regularization and Variable Selection via the Elastic Net. *Journal of the Royal Statistical Society. Series B (Methodological)* 67(2):301–320.

Franzosa EA et al. (2019). [Gut microbiome structure and metabolic activity in inflammatory bowel disease](https://www.ncbi.nlm.nih.gov/pubmed/30531976). *Nature Microbiology* 4(2):293–305.

Citation
--------

Mallick H, Franzosa EA, McIver LJ, Banerjee S, Sirota-Madi A, Kostic AD, Clish CB, Vlamakis H, Xavier R, Huttenhower C (2019). [Predictive metabolomic profiling of microbial communities using amplicon or metagenomic sequences](https://www.ncbi.nlm.nih.gov/pubmed/31316056). *Nature Communications* 10(1):3136-3146.
