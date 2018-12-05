#' @useDynLib "dina", .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @aliases dina-pkg
"_PACKAGE"

# -------- Model Simulation

#' @importFrom simcdm sim_dina_items
#' @export
simcdm::sim_dina_items

#' @importFrom simcdm sim_alpha_matrix
#' @export
simcdm::sim_alpha_matrix