#include <RcppArmadillo.h>

//' Generate Multinomial Random Variable
//' 
//' Sample a multinomial random variable for given probabilities. 
//' @param ps A \code{vector} for the probability of each category.
//' @return A \code{vector} from a multinomial with probability ps.
//' @author Steven Andrew Culpepper
//' @export
// [[Rcpp::export]]
double rmultinomial(const arma::vec& ps){
  unsigned int C = ps.n_elem;
  double u = R::runif(0,1);
  arma::vec cps = cumsum(ps);
  arma::vec Ips = arma::zeros<arma::vec>(C);
  
  Ips.elem(arma::find(cps < u) ).fill(1.0);

  return sum(Ips);
}

//' Generate Dirichlet Random Variable
//' 
//' Sample a Dirichlet random variable. 
//' @param deltas A \code{vector} of Dirichlet parameters.
//' @return A \code{vector} from a Dirichlet.
//' @author Steven Andrew Culpepper
//' @export
// [[Rcpp::export]]
arma::vec rDirichlet(const arma::vec& deltas){
  unsigned int C = deltas.n_elem;
  arma::vec Xgamma(C);
  
  //generating gamma(deltac,1)
  for(unsigned int c=0;c<C;c++){
    Xgamma(c) = R::rgamma(deltas(c),1.0);
  }
  return Xgamma/sum(Xgamma);
}


//' Simulation Responses from the DINA model
//' 
//' Sample responses from the DINA model for given attribute profiles, Q matrix,
//' and item parmeters. Returns a \code{matrix} of dichotomous responses
//' generated under DINA model.
//' @param alphas A N by K \code{matrix} of latent attributes.
//' @param Q A N by K \code{matrix} indicating which skills are required for which items.
//' @param ss A J \code{vector} of item slipping parameters.
//' @param gs A J \code{vector} of item guessing parameters.
//' @return A N by J \code{matrix} of responses from the DINA model.
//' @author Steven Andrew Culpepper
//' @export
//' @seealso \code{\link{DINA_Gibbs} }
//' @examples
//' ###########################################
//' #de la Torre (2009) Simulation Replication
//' ###########################################
//' N = 200
//' K = 5
//' J = 30
//' delta0 = rep(1,2^K)
//'     
//' # Creating Q matrix
//' Q = matrix(rep(diag(K),2),2*K,K,byrow=TRUE)
//' for(mm in 2:K){
//'     temp = combn(1:K,m=mm)
//'     tempmat = matrix(0,ncol(temp),K)
//'     for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
//'     Q = rbind(Q,tempmat)
//' }
//' Q = Q[1:J,]
//'         
//' # Setting item parameters and generating attribute profiles
//' ss = gs = rep(.2,J)
//' PIs = rep(1/(2^K),2^K)
//' CLs = c((1:(2^K))\%*\%rmultinom(n=N,size=1,prob=PIs) )
//'             
//' # Defining matrix of possible attribute profiles
//' As = rep(0,K)
//' for(j in 1:K){
//'     temp = combn(1:K,m=j)
//'     tempmat = matrix(0,ncol(temp),K)
//'     for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
//'     As = rbind(As,tempmat)
//' }
//' As = as.matrix(As)
//'                     
//' # Sample true attribute profiles
//' Alphas = As[CLs,]
//'                     
//' # Simulate data under DINA model 
//' gen = DINAsim(Alphas,Q,ss,gs)
//' Y_sim = gen$Y
//'                         
//' # Execute MCMC
//' # NOTE: small chain length used to reduce computation time for pedagogical example.
//' chainLength = 200
//' burnin = 100
//' outchain <- DINA_Gibbs(Y_sim, Amat = As, Q, chain_length = chainLength)
//'                             
//' # Summarize posterior samples for g and 1-s
//' mGs = apply(outchain$GamS[,burnin:chainLength],1,mean)
//' sGs = apply(outchain$GamS[,burnin:chainLength],1,sd)
//' m1mSS = 1-apply(outchain$SigS[,burnin:chainLength],1,mean)
//' s1mSS = apply(outchain$SigS[,burnin:chainLength],1,sd)
//' output=cbind(mGs,sGs,m1mSS,s1mSS)
//' colnames(output) = c('g Est','g SE','1-s Est','1-s SE')
//' rownames(output) = paste0('Item ',1:J)
//' print(output,digits=3)
//'                                 
//' # Summarize marginal skill distribution using posterior samples for latent class proportions
//' PIoutput = cbind(apply(outchain$PIs,1,mean),apply(outchain$PIs,1,sd))
//' colnames(PIoutput) = c('EST','SE')
//' rownames(PIoutput) = apply(As,1,paste0,collapse='')
//' print(PIoutput,digits=3)
// [[Rcpp::export]]
Rcpp::List DINAsim(const arma::mat& alphas, const arma::mat& Q,
                   const arma::vec& ss, const arma::vec& gs){
  unsigned int N = alphas.n_rows;
  unsigned int J = Q.n_rows;
  
  arma::mat Y = arma::zeros<arma::mat>(N,J);
  arma::mat ETA(N,J);
  double uij;
  
  for(unsigned int j=0;j<J;j++){
    for(unsigned int i=0;i<N;i++){
      ETA(i,j)=1.0;
      uij = R::runif(0.,1.);
      
      if(arma::dot(alphas.row(i),Q.row(j)) <arma::dot(Q.row(j),Q.row(j))  ){
        ETA(i,j) = 0.0;
      }
      
      if(pow(1.0-ss(j),ETA(i,j)) * pow(gs(j),1.0-ETA(i,j))>uij){
        Y(i,j) = 1;
      }
      else{
        Y(i,j) = 0;
      }
    }
  }
    return Rcpp::List::create(Rcpp::Named("Y")=Y,
                              Rcpp::Named("ETA")=ETA
                              );
}

//' Update attributes and latent class probabilities
//' 
//' Update attributes and latent class probabilities by sampling from full
//' conditional distribution. 
//' @param Amat A C by K \code{matrix} of latent classes.
//' @param Q A N by K \code{matrix} indicating which skills are required for which items.
//' @param ss A J \code{vector} of item slipping parameters.
//' @param gs A J \code{vector} of item guessing parameters.
//' @param Y A N by J \code{matrix} of observed responses.
//' @param PIs A C \code{vector} of latent class probabilities.
//' @param ALPHAS A N by K \code{matrix} of latent attributes.
//' @param delta0 A J \code{vector} of Dirichlet prior parameters.
//' @return A N by K \code{matrix} of attributes and a C \code{vector} of class probabilities.
//' @author Steven Andrew Culpepper
//' @export
// [[Rcpp::export]]
Rcpp::List update_alpha(const arma::mat& Amat, const arma::mat& Q,
                        const arma::vec& ss, const arma::vec& gs,
                        const arma::mat& Y, const arma::vec& PIs,
                        arma::mat& ALPHAS, const arma::vec& delta0) {
    
  unsigned int N = Y.n_rows;
  unsigned int J = Q.n_rows;
  unsigned int C = Amat.n_rows;
  unsigned int K = Q.n_cols;
  double ci;
  
  arma::mat alphas_new = arma::zeros<arma::mat>(N,K);
  arma::vec PYCS(C);
  arma::vec CLASSES(N);
  arma::vec PS;
  arma::vec Ncs = arma::zeros<arma::vec>(C);
  double etaij;

  for(unsigned int i=0;i<N;i++){
    PYCS = arma::ones<arma::vec>(C);
    
    for(unsigned int c=0;c<C;c++){
      
      for(unsigned int j=0;j<J;j++){
        etaij=1.0;
        if(arma::dot(Amat.row(c),Q.row(j)) < arma::dot(Q.row(j),Q.row(j))  ){
          etaij = 0.0;
        }
        if(etaij == 1.0 and Y(i,j) == 1.0){
          PYCS(c) = (1.0-ss(j))*PYCS(c);
        }
        if(etaij == 0.0 and Y(i,j) == 1.0){
          PYCS(c) = gs(j)*PYCS(c);
        }
        if(etaij == 1.0 and Y(i,j) == 0.0){
          PYCS(c) = ss(j)*PYCS(c);
        }
        else{
          PYCS(c) = (1.0-gs(j))*PYCS(c);
        }
      }
    }
    PS = PYCS % PIs/(arma::conv_to< double >::from(PYCS.t()*PIs));
    ci = rmultinomial(PS);
    ALPHAS.row(i) = Amat.row(ci);
    Ncs(ci) = 1.0 + Ncs(ci);
    CLASSES(i) = ci;
  }  
    return Rcpp::List::create(Rcpp::Named("PYCS")=PYCS,
                              Rcpp::Named("PS")=PS,
                              Rcpp::Named("PIs_new")=rDirichlet(Ncs+delta0),
                              Rcpp::Named("CLASSES")=CLASSES 
                              );
}

//' Update item parameters
//' 
//' Update guessing and slipping parameters from full conditional distribution. 
//' @param Y A N by J \code{matrix} of observed responses.
//' @param Q A N by K \code{matrix} indicating which skills are required for which items.
//' @param ALPHAS A N by K \code{matrix} of latent attributes.
//' @param ss_old A J \code{vector} of item slipping parameters from prior iteration.
//' @param as0 Slipping prior alpha parameter for Beta distribution.
//' @param bs0 Slipping prior beta parameter for Beta distribution.
//' @param ag0 Guessing prior alpha parameter for Beta distribution.
//' @param bg0 Guessing prior beta parameter for Beta distribution.
//' @return A list with two J \code{vectors} of guessing and slipping parameters.
//' @author Steven Andrew Culpepper
//' @export
// [[Rcpp::export]]
Rcpp::List update_sg(const arma::mat& Y, const arma::mat& Q, 
                     const arma::mat& ALPHAS, const arma::vec& ss_old,
                     double as0, double bs0, double ag0, double bg0) {

  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;

  arma::vec ETA;
  arma::vec ss_new(J);
  arma::vec gs_new(J);
  arma::mat AQ = ALPHAS * Q.t();
  double T,S,G,y_dot_eta,qq,ps,pg;
  double ug,us;
  
  for(unsigned int j=0;j<J;j++){
    us=R::runif(0,1);
    ug=R::runif(0,1);
    ETA=arma::zeros<arma::vec>(N);
    qq = arma::conv_to< double >::from( Q.row(j) * (Q.row(j)).t() );
    ETA.elem(arma::find(AQ.col(j) == qq) ).fill(1.0);
    
    y_dot_eta = arma::conv_to< double >::from( (Y.col(j)).t() * ETA );
    T = sum(ETA);
    S = T - y_dot_eta;
    G = sum(Y.col(j)) - y_dot_eta;
    
  //sample s and g as linearly truncated bivariate beta

   //draw g conditoned upon s_t-1
      pg = R::pbeta(1.0-ss_old(j),G+ag0,N-T-G+bg0,1,0);
      gs_new(j) = R::qbeta(ug*pg,G+ag0,N-T-G+bg0,1,0);
    //draw s conditoned upon g
      ps = R::pbeta(1.0-gs_new(j),S+as0,T-S+bs0,1,0);
      ss_new(j) = R::qbeta(us*ps,S+as0,T-S+bs0,1,0);
  }
    return Rcpp::List::create(Rcpp::Named("ss_new")=ss_new,
                              Rcpp::Named("gs_new")=gs_new
                              );
}


//' Generate Posterior Distribution with Gibbs sampler
//' 
//' Function for sampling parameters from full conditional distributions.
//' The function returns a list of arrays or matrices with parameter posterior
//' samples. Note that the output includes the posterior samples in objects 
//' named: \code{CLASSES} = individual attribute profiles,
//' \code{PIs} = latent class proportions,
//' \code{SigS} = item slipping parameters, and 
//' \code{GamS} = item guessing parameters.
//' @param Y            A N by J \code{matrix} of observed responses.
//' @param Amat         A C by K \code{matrix} of latent classes.
//' @param Q            A N by K \code{matrix} indicating which skills are
//'                     required for which items.
//' @param chain_length Number of MCMC iterations.
//' @return A list with samples from the posterior distribution.
//' @author Steven Andrew Culpepper
//' @export
//' @seealso \code{\link{DINAsim}} 
//' @examples
//' \dontrun{
//' ###################################
//' #Tatsuoka Fraction Subtraction Data
//' ###################################
//' require(CDM)
//' data(fraction.subtraction.data)
//' Y_1984 = as.matrix(fraction.subtraction.data)
//' Q_1984 = as.matrix(fraction.subtraction.qmatrix)
//' K_1984 = ncol(fraction.subtraction.qmatrix)
//' J_1984 = ncol(Y_1984)
//'     
//' # Creating matrix of possible attribute profiles
//' As_1984 = rep(0,K_1984)
//' for(j in 1:K_1984){
//'     temp = combn(1:K_1984,m=j)
//'     tempmat = matrix(0,ncol(temp),K_1984)
//'     for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
//'     As_1984 = rbind(As_1984,tempmat)
//' }
//' As_1984 = as.matrix(As_1984)
//'             
//' # Generate samples from posterior distribution
//' # May take 8 minutes
//' chainLength <- 5000
//' burnin <- 1000
//' outchain_1984 <- DINA_Gibbs(Y = Y_1984, Amat = As_1984,
//'                             Q_1984, chain_length = chainLength)
//'                 
//' # Summarize posterior samples for g and 1-s
//' mgs_1984 = apply(outchain_1984$GamS[,burnin:chainLength],1,mean)
//' sgs_1984 = apply(outchain_1984$GamS[,burnin:chainLength],1,sd)
//' mss_1984 = 1-apply(outchain_1984$SigS[,burnin:chainLength],1,mean)
//' sss_1984 = apply(outchain_1984$SigS[,burnin:chainLength],1,sd)
//' output_1984=cbind(mgs_1984,sgs_1984,mss_1984,sss_1984)
//' colnames(output_1984) = c('g Est','g SE','1-s Est','1-s SE')
//' rownames(output_1984) = colnames(Y_1984)
//' print(output_1984,digits=3)
//'                     
//' # Summarize marginal skill distribution using posterior samples for latent class proportions
//' marg_PIs = t(As_1984)\%*\%outchain_1984$PIs
//' PI_Est = apply(marg_PIs[,burnin:chainLength],1,mean)
//' PI_Sd = apply(marg_PIs[,burnin:chainLength],1,sd)
//' PIoutput = cbind(PI_Est,PI_Sd)
//' colnames(PIoutput) = c('EST','SE')
//' rownames(PIoutput) = paste0('Skill ',1:K_1984)
//' print(PIoutput,digits=3)
//' 
//' #####################################################
//' #de la Torre (2009) Simulation Replication w/ N = 200
//' #####################################################
//' N = 200
//' K = 5
//' J = 30
//' delta0 = rep(1,2^K)
//' 
//' #Creating Q matrix
//' Q = matrix(rep(diag(K),2),2*K,K,byrow=TRUE)
//' for(mm in 2:K){
//'     temp = combn(1:K,m=mm)
//'     tempmat = matrix(0,ncol(temp),K)
//'     for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
//'     Q = rbind(Q,tempmat)
//' }
//' Q = Q[1:J,]
//'     
//' # Setting item parameters and generating attribute profiles
//' ss = gs = rep(.2,J)
//' PIs = rep(1/(2^K),2^K)
//' CLs = c((1:(2^K))\%*\%rmultinom(n=N,size=1,prob=PIs) )
//'         
//' # Defining matrix of possible attribute profiles
//' As = rep(0,K)
//' for(j in 1:K){
//'     temp = combn(1:K,m=j)
//'     tempmat = matrix(0,ncol(temp),K)
//'     for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
//'     As = rbind(As,tempmat)
//' }
//' As = as.matrix(As)
//'                 
//' # Sample true attribute profiles
//' Alphas = As[CLs,]
//'                 
//' # Simulate data under DINA model 
//' gen = DINAsim(Alphas,Q,ss,gs)
//' Y_sim = gen$Y
//'                     
//' # Execute MCMC
//' # NOTE small chain length used to reduce computation time for pedagogical example.
//' chainLength = 200
//' burnin = 100
//'                     
//' outchain <- DINA_Gibbs(Y_sim,Amat=As,Q,chain_length=chainLength)
//'                         
//' # Summarize posterior samples for g and 1-s  
//' mGs = apply(outchain$GamS[,burnin:chainLength],1,mean)
//' sGs = apply(outchain$GamS[,burnin:chainLength],1,sd)
//' m1mSS = 1 - apply(outchain$SigS[,burnin:chainLength],1,mean)
//' s1mSS = apply(outchain$SigS[,burnin:chainLength],1,sd)
//' output = cbind(mGs,sGs,m1mSS,s1mSS)
//' colnames(output) = c('g Est','g SE','1-s Est','1-s SE')
//' rownames(output) = paste0('Item ',1:J)
//' print(output, digits=3)
//'                             
//' # Summarize marginal skill distribution using posterior samples for latent class proportions
//' PIoutput = cbind(apply(outchain$PIs,1,mean),apply(outchain$PIs,1,sd))
//' colnames(PIoutput) = c('EST','SE')
//' rownames(PIoutput) = apply(As,1,paste0,collapse='')
//' print(PIoutput,digits=3)
//' }
// [[Rcpp::export]]
Rcpp::List DINA_Gibbs(const arma::mat& Y, 
                      const arma::mat& Amat,
                      const arma::mat& Q,
                      unsigned int chain_length = 10000){
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  unsigned int K = Amat.n_cols;
  unsigned int C = Amat.n_rows;
  
  //Prior values for betas and Dirichlet distribution
  arma::vec delta0 = arma::ones<arma::vec>(C);
  double as0=1.0;
  double bs0=1.0;
  double ag0=1.0;
  double bg0=1.0;
  arma::vec pil0=arma::ones<arma::vec>(C)/double(C);//prior probability

  //Savinging output
  arma::mat SigS(J,chain_length);
  arma::mat GamS(J,chain_length);
  arma::mat US(J,chain_length);
  arma::mat PIs(C,chain_length);
  arma::mat CLASSES(N,chain_length);
  arma::cube QS(J,K,chain_length);

  //need to initialize, alphas, ss, gs,pis 
//  arma::mat alphas = arma::zeros<arma::mat>(N,K); //K>1 is assumed
  arma::mat alphas = arma::randu<arma::mat>(N,K); //K>1 is assumed   
  alphas.elem( find(alphas > 0.5) ).ones();
  alphas.elem( find(alphas <= 0.5) ).zeros();
  
  arma::vec ss = arma::randu<arma::vec>(J);
  arma::vec gs = arma::randu<arma::vec>(J);
  arma::vec pis = arma::randu<arma::vec>(C);

  //Start Markov chain
  for(unsigned int t = 0; t < chain_length; t++){

    //updata alpha and pi
    Rcpp::List step1a = update_alpha(Amat,Q,ss,gs,Y,pis,alphas,delta0);

      //update value for pis. alphas are updated via pointer. save classes and PIs
      pis             = Rcpp::as<arma::vec>(step1a[2]);
      CLASSES.col(t)  = Rcpp::as<arma::vec>(step1a[3]);
      PIs.col(t)      = pis;

    //update s and g
    Rcpp::List step1b = update_sg(Y,Q,alphas,ss,as0,bs0,ag0,bg0);
    
      //update value for ss and gs. 
      ss              = Rcpp::as<arma::vec>(step1b[0]);
      gs              = Rcpp::as<arma::vec>(step1b[1]);
      SigS.col(t)       = ss;
      GamS.col(t)       = gs;
  }

 return Rcpp::List::create(Rcpp::Named("CLASSES",CLASSES),
                          Rcpp::Named("PIs",PIs),
                          Rcpp::Named("SigS",SigS),
                          Rcpp::Named("GamS",GamS)
                          );
}
  
