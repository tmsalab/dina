#' @useDynLib "dina", .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @aliases dina-pkg
"_PACKAGE"

# -------- Model Simulation

#' @importFrom simcdm sim_dina
#' @export
simcdm::sim_dina