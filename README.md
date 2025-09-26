

<!-- README.md is generated from README.qmd. Please edit that file -->

# dina

<!-- badges: start -->

[![R-CMD-check](https://github.com/tmsalab/dina/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tmsalab/dina/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Estimate the Deterministic Input, Noisy And Gate (DINA) cognitive
diagnostic model parameters using the Gibbs sampler described by
Culpepper (2015) \<doi: 10.3102/1076998615595403\>.

## Installation

You can install `dina` from CRAN using:

``` r
install.packages("dina")
```

Or, you can be on the cutting-edge development version on GitHub using:

``` r
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("tmsalab/dina")
```

## Usage

To use the `dina` package, load it into *R* using:

``` r
library("dina")
```

From there, the DINA CDM can be estimated using:

``` r
dina_model = dina(<data>, <q>, chain_length = 10000)
```

To simulate item data under DINA, use:

``` r
# Set a seed for reproducibility
set.seed(888)

# Setup Parameters
N = 15   # Number of Examinees / Subjects
J = 10   # Number of Items
K = 2    # Number of Skills / Attributes

# Assign slipping and guessing values for each item
ss = gs = rep(.2, J)

# Simulate identifiable Q matrix
Q = sim_q_matrix(J, K)

# Simulate subject attributes
subject_alphas = sim_subject_attributes(N, K)

# Simulate Item Data
items_dina = sim_dina_items(subject_alphas, Q, ss, gs)
```

## Authors

Steven Andrew Culpepper and James Joseph Balamuta

## Citing the `dina` package

To ensure future development of the package, please cite `dina` package
if used during an analysis or simulation studies. Citation information
for the package may be acquired by using in *R*:

``` r
citation("dina")
```

## License

GPL (\>= 2)
