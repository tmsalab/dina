# dina 2.0.0

## API Breakage

- Deprecated `DINA_Gibbs()` in favor of `dina()`, which generates the correct
  alpha matrix (`Amat`) inside of the function instead of relying on the user
  to set it up.
- The call to estimate with the gibbs sampling technique is now: `dina(Y, Q, chain_length)`

## Changes

- Switched internal portions of the package to use the `simcdm` _C++_ routines
  and imported _R_ level-routines.
- Switched from `src/init.c` to autogeneration via Rcpp 0.12.15
- Removed miscellaneous RNG seed. 

## Documentation

- Enabled markdown for inline documentiona with roxygen2.
- Improved documentation flow

## Deployment

- Added TMSA Lab's Travis-CI configuration for testing across R versions.
- Added Unit Tests for model reproducibility.
- Added code coverage results.

# dina 1.0.2

- Addressed R 3.4 registration requirements.
- Added URL to GitHub repository.

# dina 1.0.1

- Fixed notation in a few examples and computation code.

# dina 1.0.0

- Initial release of the `dina` package.
