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
  
  // number of species
  int<lower = 1> S;
  
  // number of factors
  int<lower = 1> d;
  
  // survey level information  
  int<lower = 0> y[N, S];
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
  matrix[ncov_psi, S] beta_psi;
  matrix[ncov_ord, d] beta_ord;
  matrix[n, d] E;
  matrix[d, S] L_params;
  // vector[n] U;
  // vector[n] E;
  // vector[S-1] L_params;
  
  matrix[ncov_theta, S] beta_theta;
  real<lower = 0, upper = 1> theta0[S];
  real<lower = 0, upper = 1> p[S];
  real<lower = 0, upper = 1> q[S];
  
}

transformed parameters {
  
  matrix[d, S] L;
  for(k in 1:d){
    if(k > 1){
      for(s in 1:(k-1)){
        L[k,s] = 0;
      }
    } 
    L[k,k] = 1;
    for(s in (k+1):S){
      L[k,s] = L_params[k, s];
    }
  }
  
  // vector[S] L;
  // L[1] = 1;
  // for(s in 2:S){
  //   L[s] = L_params[s-1];
  // }
  
  // matrix[n, d] U = E;
  // matrix[n, S] UL = U * L;
  matrix[n, d] U = X_ord * beta_ord + E;
  // matrix[n, S] logit_psi = X_psi * beta_psi;
  matrix[n, S] logit_psi = X_psi * beta_psi + U * L;
  matrix[N, S] logit_theta = X_theta * beta_theta;
}

model {
  matrix[n, S] log_psi = log_inv_logit(logit_psi);
  matrix[n, S] log1m_psi = log1m_inv_logit(logit_psi);
  
  matrix[N, S] log_theta = log_inv_logit(logit_theta);
  matrix[N, S] log1m_theta = log1m_inv_logit(logit_theta);
  
  real log_p_yz1;
  real log_p_yz0;
  
  for(cov in 1:ncov_ord){
    for(k in 1:d){
      
    beta_ord[cov, k] ~ normal(prior_beta_psi, prior_beta_psi_sd);
    }
  }
  
  for(s in 1:S){
    
    for(cov in 1:ncov_psi){
      beta_psi[cov, s] ~ normal(prior_beta_psi, prior_beta_psi_sd);
    }
    
    for(cov in 1:ncov_theta){
      beta_theta[cov, s] ~ normal(prior_beta_theta, prior_beta_theta_sd);
    }
    
    theta0[s] ~ beta(prior_atheta0, prior_btheta0);
    p[s] ~ beta(prior_ap, prior_bp);
    q[s] ~ beta(prior_aq, prior_bq);
    
    // for(k in 1:d){
    //   L[k,s] ~ normal(0, 1);
    // }  
    
  }
  
  for(k in 1:d){
    for(s in (k+1):S){
      L[k,s] ~ normal(0, 1);
    }
  }
  
  for(i in 1:n){
    for(k in 1:d){
      E[i,k] ~ normal(0, 1);
    }
  }
  
  for(s in 1:S){
    
    for (i in 1:n) {
      
      log_p_yz1 = 0;
      log_p_yz0 = 0;
      
      for(m in 1:M[i]){
        
         // z = 1
        log_p_yz1 = log_p_yz1 + 
        log_sum_exp(
           // delta = 1
          log_theta[sumM[i] + m,s] +  
          binomial_lpmf(
            y[sumM[i] + m,s] | 
            k_samples[sumM[i] + m], 
            p[s]),
             // delta = 0
            log1m_theta[sumM[i] + m,s] +  
            binomial_lpmf(
              y[sumM[i] + m,s] | 
              k_samples[sumM[i] + m], 
              q[s])
              );
        
         // z = 0
        log_p_yz0 = log_p_yz0 + 
        log_sum_exp(
           // delta = 1
          log(theta0[s]) +  
          binomial_lpmf(
          y[sumM[i] + m,s] | 
          k_samples[sumM[i] + m], 
          p[s]),
           // delta = 0
          log(1 - theta0[s]) +  
          binomial_lpmf(
          y[sumM[i] + m,s] |
          k_samples[sumM[i] + m], 
          q[s])
          );
      }
      
      target += log_sum_exp(
          log_psi[i, s] +  log_p_yz1,
          log1m_psi[i, s] + log_p_yz0
          );
      
    }
  }
  
}
