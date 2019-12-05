
#include <RcppArmadillo.h>
#include <simcdm.h>
#include <rgen.h>

class DINA
{
private:
    const arma::mat& Y, Q;
    arma::mat Amat, alphas, SigS, GamS, US, PIs, CLASSES;
    arma::cube QS;
    arma::vec ss, gs, pis, delta0, pil0, CLASS;
    unsigned int chain_length = 10000;
    unsigned int N, J, K, C;
    
    double as0 = 1.0, bs0 = 1.0, ag0 = as0, bg0 = bs0;
    
public:
    DINA(const arma::mat &in_Y, 
         const arma::mat &in_Q, 
         unsigned int in_chain_length): 
       Y(in_Y), Q(in_Q), chain_length(in_chain_length)
    {
        // Number of Observations
        N = Y.n_rows;
        
        // Number of Items
        J = Y.n_cols;
        
        // Number of Attributes
        K = Q.n_cols;
        
        // Number of Latent Classes (2^k)
        C = static_cast<unsigned int>(pow(2.0, static_cast<double>(K)));
        
        // Generate the latent class alpha matrix
        Amat = simcdm::attribute_classes(K);
        
        // Prior values for betas and Dirichlet distribution
        delta0 = arma::ones<arma::vec>(C);
        as0 = 1.0;
        bs0 = 1.0;
        ag0 = 1.0;
        bg0 = 1.0;
        
        pil0 = arma::ones<arma::vec>(C) / double(C); // prior probability
        
        CLASS = arma::vec(N);
        
        // Saving output
        SigS = arma::mat(J, chain_length);
        GamS = arma::mat(J, chain_length);
        US = arma::mat(J, chain_length);
        PIs = arma::mat(C, chain_length);
        CLASSES = arma::mat(N, chain_length);
        QS = arma::cube(J, K, chain_length);
        
        // Need to initialize, alphas, ss, gs, and pis
        alphas = arma::randu<arma::mat>(N, K); // K>1 is assumed
        alphas.elem(find(alphas > 0.5)).ones();
        alphas.elem(find(alphas <= 0.5)).zeros();
        
        ss = arma::randu<arma::vec>(J);
        gs = arma::randu<arma::vec>(J);
        pis = arma::randu<arma::vec>(C);
        
    }
    
    void update_alpha () { 
        double ci;
        
        arma::vec PYCS(C);
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
            PS = PYCS % pis / (arma::conv_to<double>::from(PYCS.t() * pis));
            ci = rgen::rmultinomial(PS);
            alphas.row(i) = Amat.row(ci);
            Ncs(ci) = 1.0 + Ncs(ci);
            CLASS(i) = ci;
        }
        
        // Update pi
        pis = rgen::rdirichlet(Ncs + delta0);
    }
    
    void update_guessing_and_slipping () {
        arma::vec ETA;
        arma::mat AQ = alphas * Q.t();
        double T, S, G, y_dot_eta, qq, ps, pg;
        double ug, us;
        
        for (unsigned int j = 0; j < J; j++) {
            us = R::runif(0, 1);
            ug = R::runif(0, 1);
            
            ETA = arma::zeros<arma::vec>(N);
            qq = arma::conv_to<double>::from(Q.row(j) * (Q.row(j)).t());
            ETA.elem(arma::find(AQ.col(j) == qq)).fill(1.0);
            
            y_dot_eta = arma::conv_to<double>::from((Y.col(j)).t() * ETA);
            T = sum(ETA);
            S = T - y_dot_eta;
            G = sum(Y.col(j)) - y_dot_eta;
            
            // sample s and g as linearly truncated bivariate beta
            
            // draw g conditoned upon s_t-1
            pg = R::pbeta(1.0 - ss(j), G + ag0, N - T - G + bg0, 1, 0);
            gs(j) = R::qbeta(ug * pg, G + ag0, N - T - G + bg0, 1, 0);
            // draw s conditoned upon g
            ps = R::pbeta(1.0 - gs(j), S + as0, T - S + bs0, 1, 0);
            ss(j) = R::qbeta(us * ps, S + as0, T - S + bs0, 1, 0);
        }
    }
    
    void run_estimation() {
        // Start Markov chain
        for (unsigned int t = 0; t < chain_length; t++) {

            // Step1a: updata alpha and pi
            update_alpha();
            
            // Update values for pis and alphas are updated via pointer. 
            // Save classes and PIs
            CLASSES.col(t) = CLASS;
            PIs.col(t) = pis;

            // Step1b: update s and g
            update_guessing_and_slipping();
            
            // Save values for ss and gs.
            SigS.col(t) = ss;
            GamS.col(t) = gs;
        }
    }
    
    
    
    unsigned int get_chain_length() const { return chain_length; }
    unsigned int get_observations() const { return N; }
    unsigned int get_items() const { return J; }
    unsigned int get_attributes() const { return K; }
    
    arma::mat get_classes() const { return CLASSES; }
    arma::mat get_pis() const { return PIs; }
    arma::mat get_slipping() const { return SigS; }
    arma::mat get_guessing() const { return GamS; }
    
};


//' @export
// [[Rcpp::export]]
Rcpp::List dina_cpp_class(const arma::mat& Y, const arma::mat& Q, 
                          unsigned int chain_length) {
    DINA dina_model(Y, Q, chain_length);
    
    dina_model.run_estimation();
    
    return Rcpp::List::create(
        Rcpp::Named("CLASSES", dina_model.get_classes()),
        Rcpp::Named("PIs", dina_model.get_pis()),
        Rcpp::Named("SigS", dina_model.get_slipping()), 
        Rcpp::Named("GamS", dina_model.get_guessing()));
    
    }