
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/tmsalab/dina.svg)](https://travis-ci.org/tmsalab/dina)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=2\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN Version
Badge](http://www.r-pkg.org/badges/version/dina)](https://cran.r-project.org/package=dina)
[![CRAN
Status](https://cranchecks.info/badges/worst/dina)](https://cran.r-project.org/web/checks/check_results_dina.html)
[![RStudio CRAN Mirror’s Monthly
Downloads](http://cranlogs.r-pkg.org/badges/dina?color=brightgreen)](http://www.r-pkg.org/pkg/dina)
[![RStudio CRAN Mirror’s Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/dina?color=brightgreen)](http://www.r-pkg.org/pkg/dina)
[![Coverage
status](https://codecov.io/gh/tmsalab/dina/branch/master/graph/badge.svg)](https://codecov.io/github/tmsalab/dina?branch=master)

# `dina` R package

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

# Item data
items_dina = sim_dina_items(subject_alphas, Q, ss, gs)
```

## Authors

Steven Andrew Culpepper and James Joseph Balamuta

## Citing the `dina` package

To ensure future development of the package, please cite `dina` package
if used during an analysis or simulations. Citation information for the
package may be acquired by using in *R*:

``` r
citation("dina")
```

## License

GPL (\>= 2)
