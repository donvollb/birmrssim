data {
  int<lower=1> N;                     // Number of persons
  int<lower=1> theta_n;               // Number of theta dimensions
  array[theta_n] int<lower=1> T;      // Items per theta
  vector[sum(T)] x;                   // Reverse-scored item indicators; not relevant here
  matrix[N, sum(T)] r;                // Item responses
  real ars_prior;                     // SD-Prior for ARS; not relevant here
}

transformed data {
  int<lower=1> K = sum(T);            // Total number of items
}

parameters {
  array[N] vector[theta_n + 1] theta; // Latent traits (last column = ERS)
  vector[K] delta;                    // Item difficulties
  vector<lower=0.1>[K] tau;           // Slightly relaxed lower bound
  vector<lower=0>[theta_n + 1] sigma; // Standard deviations
  corr_matrix[theta_n + 1] Sigma_corr;// Correlation matrix
}

transformed parameters {
  matrix[N, K] m;                     // first shape parameter for the Beta-distribution
  matrix[N, K] n;                     // second shape parameter for the Beta-distribution
  matrix[theta_n + 1, theta_n + 1] L_Sigma = diag_pre_multiply(sigma, cholesky_decompose(Sigma_corr));    // build the covariance matrix for theta(s) and ERS
  vector[N] ers_exp;                  // create a vector for exp(ers)

  for (j in 1:N)
    ers_exp[j] = exp(theta[j, theta_n + 1]);  // fill the vector with exp(ers)
  
  {

    // In the following lines, we constrain the exponents for the beta distribution shape parameters 
    // to lie within the range of exp(-7) and exp(7). This constraint is purely computational: 
    // values outside this range have negligible practical influence on the beta distribution but significantly increase computational complexity.
    // Importantly, this restriction is intentionally chosen to be non-informative and may, in fact, slightly degrade parameter estimates. 
    real MAX_ARG = 7;                 
    real MIN_ARG = -7;


    int k_lower = 1;
    int k_upper = T[1];
    
    for (i in 1:theta_n) {
      for (k in k_lower:k_upper) {
        vector[N] core; // just a placeholder variable for the core part of the equation


        for (j in 1:N)        
        core[j] = ers_exp[j] * (theta[j, i] - delta[k]); // compute the core
      
        m[, k] = exp(fmin(fmax((core + tau[k]) / 2, MIN_ARG), MAX_ARG)); // compute the shape parameters
        n[, k] = exp(fmin(fmax((-core + tau[k]) / 2, MIN_ARG), MAX_ARG)); // compute the shape parameters

      }
      if (i != theta_n) {
        k_lower += T[i];
        k_upper += T[i + 1];
      }
    }
  }
}

model {
  // Priors
  Sigma_corr ~ lkj_corr(2.0);        // Mildly informative prior
  delta ~ normal(0, 3);              // wide prior
  tau ~ normal(1, 2);                // wide prior
  sigma ~ normal(0, 1);              // standard normal prior for SD
  
  // Vectorized and Cholesky factorized latent traits
  theta ~ multi_normal_cholesky(rep_vector(0, theta_n + 1), L_Sigma);  // multivariate prior for theta(s) and ERS

  // Likelihood 
  to_vector(r) ~ beta(to_vector(m), to_vector(n));
}

generated quantities {

  // log-likelihood
  matrix[N, K] log_lik;

  for (j in 1:N) {
    for (k in 1:K) {
      log_lik[j, k] = beta_lpdf(r[j, k] | m[j, k], n[j, k]);
    }
  }
}
