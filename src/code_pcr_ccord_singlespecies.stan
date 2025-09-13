data {
  // number of sites
  int<lower = 1> n;
  int<lower = 0> ncov_psi;
  matrix[n, ncov_psi] X_psi;
  int<lower = 0> ncov_ord;
  matrix[n, ncov_ord] X_ord;
  
  // sample-level detection covariates
  int<lower = 1> N;
  int<lower = 1> M[n];
  int<lower = 0> sumM[n];
  int<lower = 1> ncov_theta;
  matrix[N, ncov_theta] X_theta;
  
  // survey level information  
  int<lower = 0> y[N];
  int<lower = 1> k_samples[N];
  
  // priors
  real prior_beta_psi;
  real<lower = 0> prior_beta_psi_sd;
  real prior_beta_theta;
  real<lower = 0> prior_beta_theta_sd;
  real<lower = 0> prior_atheta0;
  real<lower = 0> prior_btheta0;
  real<lower = 0> prior_ap;
  real<lower = 0> prior_bp;
  real<lower = 0> prior_aq;
  real<lower = 0> prior_bq;
  
}
parameters {
  vector[ncov_psi] beta_psi;
  
  vector[ncov_theta] beta_theta;
  real<lower = 0, upper = 1> theta0;
  real<lower = 0, upper = 1> p;
  real<lower = 0, upper = 1> q;
}
transformed parameters {
  // matrix[d, S] L;
  // for(k in 1:d){
    //   L[k,1] = 0;
    //   for(s in 1:(S-1)){
      //     L[k,s + 1] = L_params[k, s];
      //   }
    // }
 
  vector[n] logit_psi = X_psi * beta_psi;
  vector[N] logit_theta = X_theta * beta_theta;
  
  vector[n] log_psi = log_inv_logit(logit_psi);
  vector[n] log1m_psi = log1m_inv_logit(logit_psi);
  
  vector[N] log_theta = log_inv_logit(logit_theta);
  vector[N] log1m_theta = log1m_inv_logit(logit_theta);
}
model {
  
  real log_p_yz1;
  real log_p_yz0;
  
  for(cov in 1:ncov_psi){
    beta_psi[cov] ~ normal(prior_beta_psi, prior_beta_psi_sd);
  }
  
  for(cov in 1:ncov_theta){
    beta_theta[cov] ~ normal(prior_beta_theta, prior_beta_theta_sd);
  }
  
  theta0 ~ beta(prior_atheta0, prior_btheta0);
  p ~ beta(prior_ap, prior_bp);
  q ~ beta(prior_aq, prior_bq);
  
  
  for (i in 1:n) {
      
      log_p_yz1 = 0;
      log_p_yz0 = 0;
      
      for(m in 1:M[i]){
        
         // z = 1
        log_p_yz1 = log_p_yz1 + 
        log_sum_exp(
           // delta = 1
          log_theta[sumM[i] + m] +  
          binomial_lpmf(
            y[sumM[i] + m] | 
            k_samples[sumM[i] + m], 
            p),
             // delta = 0
            log1m_theta[sumM[i] + m] +  
            binomial_lpmf(
              y[sumM[i] + m] | 
              k_samples[sumM[i] + m], 
              q)
              );
        
         // z = 0
        log_p_yz0 = log_p_yz0 + 
        log_sum_exp(
           // delta = 1
          log(theta0) +  
          binomial_lpmf(
          y[sumM[i] + m] | 
          k_samples[sumM[i] + m], 
          p),
           // delta = 0
          log(1 - theta0) +  
          binomial_lpmf(
          y[sumM[i] + m] |
          k_samples[sumM[i] + m], 
          q)
          );
      }
      
      target += log_sum_exp(
          log_psi[i] +  log_p_yz1,
          log1m_psi[i] + log_p_yz0
          );
      
    }
  
}


generated quantities {
  
  vector[n] log_lik;
  
  real log_p_yz1;
  real log_p_yz0;
  
  for (i in 1:n) {
    
    log_p_yz1 = 0;
    log_p_yz0 = 0;
    
    for(m in 1:M[i]){
      
      // z = 1
      log_p_yz1 = log_p_yz1 + 
      log_sum_exp(
        // delta = 1
        log_theta[sumM[i] + m] +  
        binomial_lpmf(
          y[sumM[i] + m] | 
          k_samples[sumM[i] + m], 
          p),
          // delta = 0
          log1m_theta[sumM[i] + m] +  
          binomial_lpmf(
            y[sumM[i] + m] | 
            k_samples[sumM[i] + m], 
            q)
            );
            
            // z = 0
            log_p_yz0 = log_p_yz0 + 
            log_sum_exp(
              // delta = 1
              log(theta0) +  
              binomial_lpmf(
                y[sumM[i] + m] | 
                k_samples[sumM[i] + m], 
                p),
                // delta = 0
                log(1 - theta0) +  
                binomial_lpmf(
                  y[sumM[i] + m] |
                  k_samples[sumM[i] + m], 
                  q)
                  );
    }
    
    log_lik[i] = log_sum_exp(
      log_psi[i] +  log_p_yz1,
      log1m_psi[i] + log_p_yz0
      );
      
  }
  
  
}
