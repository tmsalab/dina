#include <RcppArmadillo.h>
#include <rgen.h>
#include <simcdm.h>

//' Update attributes and latent class probabilities
//'
//' Update attributes and latent class probabilities by sampling from full
//' conditional distribution.
//'
//' @param Amat   A \eqn{C \times K}{C x K} `matrix` of latent classes.
//' @param Q      A \eqn{N \times K}{N x K} `matrix` indicating which
//'               skills are required for which items.
//' @param ss     A \eqn{J} `vector` of item slipping parameters.
//' @param gs     A \eqn{J} `vector` of item guessing parameters.
//' @param Y      A \eqn{N \times J}{N x J} `matrix` of observed responses.
//' @param PIs    A \eqn{C} `vector` of latent class probabilities.
//' @param ALPHAS A \eqn{N \times K}{N x K} `matrix` of latent attributes.
//' @param delta0 A \eqn{J} `vector` of Dirichlet prior parameters.
//'
//' @return
//' A \eqn{N \times K}{N x K} `matrix` of attributes and a C `vector` of
//' class probabilities.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List update_alpha(const arma::mat &Amat, const arma::mat &Q,
                        const arma::vec &ss, const arma::vec &gs,
                        const arma::mat &Y, const arma::vec &PIs,
                        arma::mat &ALPHAS, const arma::vec &delta0)
{

    unsigned int N = Y.n_rows;
    unsigned int J = Q.n_rows;
    unsigned int C = Amat.n_rows;
    unsigned int K = Q.n_cols;
    double ci;

    arma::mat alphas_new = arma::zeros<arma::mat>(N, K);
    arma::vec PYCS(C);
    arma::vec CLASSES(N);
    arma::vec PS;
    arma::vec Ncs = arma::zeros<arma::vec>(C);
    double etaij;

    for (unsigned int i = 0; i < N; i++) {
        PYCS = arma::ones<arma::vec>(C);

        for (unsigned int c = 0; c < C; c++) {

            for (unsigned int j = 0; j < J; j++) {
                etaij = 1.0;
                if (arma::dot(Amat.row(c), Q.row(j)) <
                    arma::dot(Q.row(j), Q.row(j))) {
                    etaij = 0.0;
                }
                if (etaij == 1.0 and Y(i, j) == 1.0) {
                    PYCS(c) = (1.0 - ss(j)) * PYCS(c);
                }
                if (etaij == 0.0 and Y(i, j) == 1.0) {
                    PYCS(c) = gs(j) * PYCS(c);
                }
                if (etaij == 1.0 and Y(i, j) == 0.0) {
                    PYCS(c) = ss(j) * PYCS(c);
                } else {
                    PYCS(c) = (1.0 - gs(j)) * PYCS(c);
                }
            }
        }
        PS = PYCS % PIs / (arma::as_scalar(PYCS.t() * PIs));
        ci = rgen::rmultinomial(PS);
        ALPHAS.row(i) = Amat.row(ci);
        Ncs(ci) = 1.0 + Ncs(ci);
        CLASSES(i) = ci;
    }
    return Rcpp::List::create(
        Rcpp::Named("PYCS") = PYCS, Rcpp::Named("PS") = PS,
        Rcpp::Named("PIs_new") = rgen::rdirichlet(Ncs + delta0),
        Rcpp::Named("CLASSES") = CLASSES);
}

//' Update item parameters
//'
//' Update guessing and slipping parameters from full conditional distribution.
//'
//' @param Y      A N by J `matrix` of observed responses.
//' @param Q      A N by K `matrix` indicating which skills are required 
//'               for which items. 
//' @param ALPHAS A N by K `matrix` of latent attributes. 
//' @param ss_old A J `vector` of item slipping parameters from prior iteration.
//' @param as0    Slipping prior alpha parameter for Beta distribution.
//' @param bs0    Slipping prior beta parameter for Beta distribution. 
//' @param ag0    Guessing prior alpha parameter for Beta distribution. 
//' @param bg0    Guessing prior beta parameter for Beta distribution.
//'
//' @return
//' A `list` with two J `vectors` of guessing and slipping parameters.
//'
//' @author Steven Andrew Culpepper
//' @noRd
// [[Rcpp::export]]
Rcpp::List update_sg(const arma::mat &Y, const arma::mat &Q,
                     const arma::mat &ALPHAS, const arma::vec &ss_old,
                     double as0, double bs0, double ag0, double bg0)
{

    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;

    arma::vec ETA;
    arma::vec ss_new(J);
    arma::vec gs_new(J);
    arma::mat AQ = ALPHAS * Q.t();
    double T, S, G, y_dot_eta, qq, ps, pg;
    double ug, us;

    for (unsigned int j = 0; j < J; j++) {
        us = R::runif(0, 1);
        ug = R::runif(0, 1);
        ETA = arma::zeros<arma::vec>(N);
        qq = arma::as_scalar(Q.row(j) * (Q.row(j)).t());
        ETA.elem(arma::find(AQ.col(j) == qq)).fill(1.0);

        y_dot_eta = arma::as_scalar((Y.col(j)).t() * ETA);
        T = sum(ETA);
        S = T - y_dot_eta;
        G = sum(Y.col(j)) - y_dot_eta;

        // sample s and g as linearly truncated bivariate beta

        // draw g conditoned upon s_t-1
        pg = R::pbeta(1.0 - ss_old(j), G + ag0, N - T - G + bg0, 1, 0);
        gs_new(j) = R::qbeta(ug * pg, G + ag0, N - T - G + bg0, 1, 0);
        // draw s conditoned upon g
        ps = R::pbeta(1.0 - gs_new(j), S + as0, T - S + bs0, 1, 0);
        ss_new(j) = R::qbeta(us * ps, S + as0, T - S + bs0, 1, 0);
    }
    return Rcpp::List::create(Rcpp::Named("ss_new") = ss_new,
                              Rcpp::Named("gs_new") = gs_new);
}

//' Generate Posterior Distribution with Gibbs sampler
//'
//' Function for sampling parameters from full conditional distributions.
//' The function returns a list of arrays or matrices with parameter posterior
//' samples. Note that the output includes the posterior samples in objects.
//'
//' @param Y            A \eqn{N \times J}{N x J} `matrix` of observed
//'                     responses. 
//' @param Amat         A \eqn{C \times K}{C x K} `matrix` of latent
//'                     classes. 
//' @param Q            A \eqn{N \times K}{N x K} `matrix` indicating
//'                     which skills are required for which items. 
//' @param chain_length Number of MCMC iterations.
//'
//' @return
//' A `list` with samples from the posterior distribution with each
//' entry named:
//'
//' - `CLASSES` = individual attribute profiles,
//' - `PIs` = latent class proportions,
//' - `SigS` = item slipping parameters, and
//' - `GamS` = item guessing parameters.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [simcdm::sim_dina_items()] and [simcdm::attribute_classes()]
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List DINA_Gibbs_cpp(const arma::mat &Y, const arma::mat &Q,
                          unsigned int chain_length = 10000)
{

    // Number of Observations
    unsigned int N = Y.n_rows;

    // Number of Items
    unsigned int J = Y.n_cols;

    // Number of Attributes
    unsigned int K = Q.n_cols;

    // Number of Latent Classes (2^k)
    unsigned int C = static_cast<unsigned int>(pow(2.0, static_cast<double>(K)));

    // Generate the latent class alpha matrix
    arma::mat Amat = simcdm::attribute_classes(K);

    // Prior values for betas and Dirichlet distribution
    arma::vec delta0 = arma::ones<arma::vec>(C);
    double as0 = 1.0;
    double bs0 = 1.0;
    double ag0 = 1.0;
    double bg0 = 1.0;

    arma::vec pil0 = arma::ones<arma::vec>(C) / double(C); // prior probability

    // Saving output
    arma::mat SigS(J, chain_length);
    arma::mat GamS(J, chain_length);
    arma::mat US(J, chain_length);
    arma::mat PIs(C, chain_length);
    arma::mat CLASSES(N, chain_length);
    arma::cube QS(J, K, chain_length);

    // Need to initialize, alphas, ss, gs, and pis
    //  arma::mat alphas = arma::zeros<arma::mat>(N,K); //K>1 is assumed
    arma::mat alphas = arma::randu<arma::mat>(N, K); // K>1 is assumed
    alphas.elem(find(alphas > 0.5)).ones();
    alphas.elem(find(alphas <= 0.5)).zeros();

    arma::vec ss = arma::randu<arma::vec>(J);
    arma::vec gs = arma::randu<arma::vec>(J);
    arma::vec pis = arma::randu<arma::vec>(C);

    // Start Markov chain
    for (unsigned int t = 0; t < chain_length; t++) {

        // updata alpha and pi
        Rcpp::List step1a =
            update_alpha(Amat, Q, ss, gs, Y, pis, alphas, delta0);

        // update value for pis. alphas are updated via pointer. save classes
        // and PIs
        pis = Rcpp::as<arma::vec>(step1a[2]);
        CLASSES.col(t) = Rcpp::as<arma::vec>(step1a[3]);
        PIs.col(t) = pis;

        // update s and g
        Rcpp::List step1b = update_sg(Y, Q, alphas, ss, as0, bs0, ag0, bg0);

        // update value for ss and gs.
        ss = Rcpp::as<arma::vec>(step1b[0]);
        gs = Rcpp::as<arma::vec>(step1b[1]);
        SigS.col(t) = ss;
        GamS.col(t) = gs;
    }

    return Rcpp::List::create(
        Rcpp::Named("CLASSES", CLASSES), Rcpp::Named("PIs", PIs),
        Rcpp::Named("SigS", SigS), Rcpp::Named("GamS", GamS));
}
