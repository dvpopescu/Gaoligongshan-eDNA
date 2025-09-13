
logit <- function(x){
  log(x / (1 - x))
}

logistic <- function(x) 1 / (1 + exp(-x))

simulateData <- function(n, M, K, S, d, 
                         ncov_psi, ncov_theta, ncov_ord,
                         meanPsi, sdPsi, 
                         meanTheta, sdTheta,
                         a_theta0, b_theta0,
                         a_p, b_p, a_q, b_q){
  
  N <- sum(M)
  
  p_true <- rbeta(S, a_p, b_p)
  q_true <- rbeta(S, a_q, b_q)
  theta0_true <- rbeta(S, a_theta0, b_theta0)
    
  beta_psi_true <- matrix(sample(c(-1,1), ncov_psi * S, replace = T), ncov_psi, S)
  beta_theta0_true <- rnorm(S, logit(meanTheta), sdTheta)
  beta_theta_true <- matrix(sample(c(-1,1), ncov_theta * S, replace = T), ncov_theta, S)
  beta_theta_true <- rbind(beta_theta0_true, beta_theta_true)
  beta_ord_true <- matrix(sample(c(-1,1), ncov_ord * d, replace = T), ncov_ord, d)

  X_psi <- matrix(rnorm(n * ncov_psi), n, ncov_psi)
  X_theta <- matrix(rnorm(N * ncov_theta), N, ncov_theta)
  X_theta <- cbind(1, X_theta)
  X_ord <- matrix(rnorm(n * ncov_ord), n, ncov_ord)
  
  E_true <- matrix(rnorm(n * d, sd = .5), n, d)
  U_true <- X_ord %*% beta_ord_true + E_true
  
  L_true <- matrix(rnorm(d * S), d, S)
  L_true[lower.tri(L_true)] <- 0
  diag(L_true) <- 1
  
  psi <- logistic(X_psi %*% beta_psi_true + U_true %*% L_true)
  theta <- logistic(X_theta %*% beta_theta_true)
  
  z <- t(sapply(1:n, function(i){
    sapply(1:S, function(j){
      rbinom(1, 1, psi[i,j])
    })
  }))
  
  delta <- matrix(0, N, S)
  for (j in 1:S) {
    for (i in 1:n) {
      for (m in 1:M[i]) {
        if(z[i,j] == 1){
          delta[sum(M[seq_len(i-1)]) + m,j] <- 
            rbinom(1, 1, theta[sum(M[seq_len(i-1)]) + m,j])
        } else {
          delta[sum(M[seq_len(i-1)]) + m,j] <- 
            rbinom(1, 1, theta0_true[j])
        }  
      }
    }
  }
  
  y <- matrix(0, N, S)
  for (j in 1:S) {
    for (i in 1:N) {
      y[i,j] <- rbinom(1, K, p_true[j] * delta[i,j] + q_true[j] * (1 - delta[i,j]))
    }
  }
  
  data <- list("X_psi" = X_psi,
               "X_ord" = X_ord,
               "X_theta" = X_theta,
               "y" = y)
  
  trueParams <- list(
    "beta_psi" = beta_psi_true,
    "beta_theta" = beta_theta_true,
    "beta_ord" = beta_ord_true,
    "L" = L_true,
    "U" = U_true,
    "p" = p_true,
    "q" = q_true,
    "theta0" = theta0_true)  
  
  list("trueParams" = trueParams,
       "data" = data)
}

simulateDataLogistic <- function(n, K, S, d, 
                         ncov_psi, ncov_ord,
                         meanPsi, sdPsi){
  
  beta_psi_true <- matrix(sample(c(-1,1), ncov_psi * S, replace = T), ncov_psi, S)
  beta_ord_true <- matrix(sample(c(-1,1), ncov_ord * d, replace = T), ncov_ord, d)

  X_psi <- matrix(rnorm(n * ncov_psi), n, ncov_psi)
  X_ord <- matrix(rnorm(n * ncov_ord), n, ncov_ord)
  
  E_true <- matrix(rnorm(n * d, sd = .5), n, d)
  U_true <- X_ord %*% beta_ord_true + E_true
  
  L_true <- matrix(rnorm(d * S), d, S)
  L_true[lower.tri(L_true)] <- 0
  diag(L_true) <- 1
  
  psi <- logistic(X_psi %*% beta_psi_true + U_true %*% L_true)
  
  y <- t(sapply(1:n, function(i){
    sapply(1:S, function(j){
      rbinom(1, K, psi[i,j])
    })
  }))
 
  data <- list("X_psi" = X_psi,
               "X_ord" = X_ord,
               "y" = y)
  
  trueParams <- list(
    "beta_psi" = beta_psi_true,
    "beta_ord" = beta_ord_true,
    "L" = L_true,
    "U" = U_true)  
  
  list("trueParams" = trueParams,
       "data" = data)
}

# data is a list with elements:
# info: dataframe with Site variables and covaraites
# OTU
runOccPlus <- function(data,
                       d,
                       occCovariates = c(),
                       ordCovariates = c(),
                       detCovariates = c()){  
  # clean data for stan
  
  data_info <- as.data.frame(data$info)
  OTU <- data$OTU
  K <- data$K

  {
    n <- length(unique(data_info$Site))
    M <- data_info %>% 
      group_by(Site) %>%
      summarise(M = n()) %>% 
      select(M) 
    M <- M$M
    N <- sum(M)
    if(is.null(dim(OTU))){
      S <- 1
      OTU <- unlist(OTU)
    } else {
      S <- ncol(OTU)
    }
    
    if(length(ordCovariates) > 0){
      
      X_ord <- data_info %>% 
        group_by(Site) %>% 
        summarise(across(all_of(ordCovariates), 
                         function(x) {x[1]}))
      
      sitesNames <- X_ord$Site
      
      X_ord <- X_ord %>% 
        select(-Site) %>% 
        mutate_if(is.numeric, scale) %>% 
        mutate(across(where(~ !is.numeric(.x)), as.factor)) %>% 
        mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>% 
        model.matrix(~., .) 
      
      X_ord <- X_ord[,-1]
      
      rownames(X_ord) <- sitesNames
      
    } else {
      X_ord <- matrix(0, n, 0)
    }
    
    if(length(occCovariates) > 0){
      
      X_psi <- data_info %>% 
        group_by(Site) %>% 
        summarise(across(all_of(occCovariates), 
                         function(x) {x[1]})) %>% 
        select(-Site) %>% 
        mutate_if(is.numeric, scale) %>% 
        mutate(across(where(~ !is.numeric(.x)), as.factor)) %>% 
        model.matrix(~., .) 
      
      X_psi <- X_psi[,-1]
    
    } else {
      
      X_psi <- matrix(0, n, 0)
      
    }
    
    if(length(detCovariates) > 0){
      
      X_theta <- data_info %>% 
        select(detCovariates) %>% 
        mutate_if(is.numeric, scale) %>% 
        mutate(across(where(~ !is.numeric(.x)), as.factor)) %>% 
        model.matrix(~., .) 
      
      X_theta <- X_theta[,-1,drop=F]
      
    } else {
      
      X_theta <- matrix(0, N, 0)
      
    }
    
    sumM <- c(0, cumsum(M)[-n])
    k_samples <- rep(K, N)
    
    ncov_psi <- ncol(X_psi)
    ncov_theta <- ncol(X_theta)
    ncov_ord <- ncol(X_ord)
    X_theta <- cbind(1, X_theta)
    
    sites <- unique(data_info$Site)
    
    start_idx <- sumM + 1
    end_idx <- c(sumM[-1], N)
    
  }
  
  # priors 
  {
    prior_beta_psi <- 0
    prior_beta_psi_sd <- 1
    prior_beta_theta <- 0
    prior_beta_theta_sd <- 1
    prior_atheta0 <- 1
    prior_btheta0 <- 10
    prior_ap <- 10
    prior_bp <- 1
    prior_aq <- 1
    prior_bq <- 10
  }
  
  occ_dat <- list(
    y = OTU,
    n = n,
    N = N,
    S = S,
    M = M,
    K = K,
    d = d,
    sumM = sumM,
    X_psi = X_psi,
    ncov_psi = ncov_psi,
    X_ord = X_ord,
    ncov_ord = ncov_ord,
    X_theta = X_theta,
    ncov_theta = ncov_theta + 1,
    k_samples = k_samples,
    prior_beta_psi = prior_beta_psi,
    prior_beta_psi_sd = prior_beta_psi_sd,
    prior_beta_theta = prior_beta_theta,
    prior_beta_theta_sd = prior_beta_theta_sd,
    prior_atheta0 = prior_atheta0,
    prior_btheta0 = prior_btheta0,
    prior_ap = prior_ap,
    prior_bp = prior_bp,
    prior_aq = prior_aq,
    prior_bq = prior_bq
  )
  
  init_fun <- function(...) list(
    theta0 = rep(0.05, S),
    p = rep(.95, S), 
    q = rep(.05, S)
  )
  
  if(S == 1){
    model <- stan_model(file = "code_pcr_ccord_singlespecies.stan", verbose = T)  
  } else {
    model <- stan_model(file = "code_pcr_ccord2.stan", verbose = T)
  }
  
  # model <- stan_model(file = here('Stan','code_pcr_ccord2.stan'), verbose = T)
  
  params <- c("beta_psi","beta_theta",
              "p","q","theta0",
              "log_lik")
  if(S > 1){
    params <- c(params, c("beta_ord","U","L","E"))
  }
  
  vb_fit <- vb(model, data = occ_dat,
               algorithm = "meanfield",
               pars = params,
               elbo_samples = 500,
               init = init_fun,
               tol_rel_obj = 0.00001,
               output_samples = 10000)
  
  matrix_of_draws <- as.matrix(vb_fit)
  
  list(
    "vb_fit" = vb_fit,
    "matrix_of_draws" = matrix_of_draws,
       "d" = d,
       "X_ord" = X_ord,
       "X_theta" = X_theta,
       "X_psi" = X_psi) 

}
