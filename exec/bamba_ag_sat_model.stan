// BAMBA model
// 2-component grouped mixture model for BAMA or Fc array data
// 1 component is centered at mu0 (unique per timepoint)
// The other component is centered at a group mean
// which is mu0 + an antigen (ag) effect (positive)
// Each observation comes either from the group centered at mu0
// or from the particular group it is a part of
// given by ag and timepoint
// Prior probability of is a product of ag and timepoint probabilities

// y_i = [y_{i1}, ..., y_{iT}]^T
// f(y_{it} | \theta) = N(\mu_{0t}, \sigma_t)^{1-z_{it}} N(\mu_{0t} + mu_{at}, \sigma_t)^{z_{it}}
// z_{it} ~ Bernoulli(\omega_{ag} * \omega_t)
// \mu_{0t} ~ N(0, 0.25) (variance is 0.25)
// \omega_{ag} ~ Beta(a_{ag}, b_{ag})
// a_{ag} = \lambda_{ag} \phi_{ag}
// b_{ag} = \lambda_{ag} (1 - \phi_{ag})
// \lambda_{ag} ~ Pareto(0.1, 1.5)
// \phi_{ag} ~ Beta(1, 1)
// \omega_t ~ Beta(a_t, b_t)
// \mu_{\at} ~ Gamma(\alpha_{ag}, \beta_{ag}
// \alpha_{ag} ~ Uniform(0, 1000)
// \beta_{ag} ~ Uniform(0, 1000)
// \sigma_t ~ Cauchy+(0, 10)

functions {
  real normal_mix_lpdf(real y, real theta, real mu1, real mu2, real sigma) { 
    return log_mix(theta, normal_lpdf(y | mu1, sigma), 
                          normal_lpdf(y | mu2, sigma)); 
  } 
  real normal_mix_lccdf(real y, real theta, real mu1, real mu2, real sigma) { 
    return log_mix(theta, normal_lccdf(y | mu1, sigma), 
                          normal_lccdf(y | mu2, sigma)); 
  } 
}

data {
  int<lower=1> N; // Number of observations
  int<lower=1> T; // Number of timepoints (dimension of observations)
  int<lower=1> obs_to_t[N];
  int<lower=1> N_ag; // grouping for antigen
  int<lower=1> obs_to_ag[N];
  int<lower=0> N_sat; // Number of saturated observations
  int<lower=1> sat_ind[N_sat]; // indexes of saturated observations
  vector[N] y;
}

parameters {
  // Overall probability of an observation being in a non-zero group
  real<lower=0, upper=1> omega_ag[N_ag];
  real<lower=0, upper=1> omega_t[T];
  
  vector<lower=0>[T] sigma; // Group variance (stddev)
  real<lower=0> alpha_ag; // mean parameter for mu_ag gamma prior
  real<lower=0> beta_ag; // variance parameter for mu_ag gamma prior
  matrix<lower=0>[N_ag, T] mu_ag; // antigen response mean
  real<lower=0, upper=1> phi_ag; // hyperprior for omega (as found in Stan manual)
  real<lower=0.1> lambda_ag; // hyperprior for omega (as found in Stan manual)
  real<lower=0, upper=1> phi_t; // hyperprior for omega (as found in Stan manual)
  real<lower=0.1> lambda_t; // hyperprior for omega (as found in Stan manual)
  vector[T] mu0;
  
  vector<lower=0>[N_sat] delta_sat; // extra to add onto saturated observations
}

transformed parameters {
  real<lower=0> a_ag;
  real<lower=0> b_ag;
  real<lower=0> a_t;
  real<lower=0> b_t;
  vector<upper=0>[2] soft_z[N]; // used to track the group membership

  vector[N] ystar;

  a_ag = lambda_ag * phi_ag;
  b_ag = lambda_ag * (1 - phi_ag);
  a_t = lambda_t * phi_t;
  b_t = lambda_t * (1 - phi_t);

  ystar = y;
  if (N_sat > 0) {
    for (n_s in 1:N_sat) {
      ystar[sat_ind[n_s]] = y[sat_ind[n_s]] + delta_sat[n_s];
    }
  }

  // Calculate log probabilities for each observation in the zero and non-zero groups
  for (n in 1:N) {
      soft_z[n][1] = log(1-omega_ag[obs_to_ag[n]] * omega_t[obs_to_t[n]]) +
			normal_lpdf(ystar[n] | mu0[obs_to_t[n]], sigma[obs_to_t[n]]);
      soft_z[n][2] = log(omega_ag[obs_to_ag[n]] * omega_t[obs_to_t[n]]) +
                    normal_lpdf(ystar[n] | mu0[obs_to_t[n]] + mu_ag[obs_to_ag[n]][obs_to_t[n]],
				       sigma[obs_to_t[n]]);
  }
}

model {
  sigma ~ cauchy(0, 10);
  phi_ag ~ beta(1, 1);          // as in Stan manual
  lambda_ag ~ pareto(0.1, 1.5); // as in Stan manual
  phi_t ~ beta(1, 1);          // as in Stan manual
  lambda_t ~ pareto(0.1, 1.5); // as in Stan manual
  omega_ag ~ beta(a_ag, b_ag);
  omega_t ~ beta(a_t, b_t);
  alpha_ag ~ uniform(0, 1000);
  beta_ag ~ uniform(0, 1000);

  mu0 ~ normal(0, 0.5);

  for (n_a in 1:N_ag) {
    mu_ag[n_a,] ~ gamma(alpha_ag, beta_ag);
  }

  if (N_sat > 0) {
    for (n_s in 1:N_sat) {
      int n = sat_ind[n_s];
      delta_sat[n_s] ~ normal_mix(omega_ag[obs_to_ag[n]] * omega_t[obs_to_t[n]],
                                  mu0[obs_to_t[n]] + mu_ag[obs_to_ag[n]][obs_to_t[n]] - y[n],
                                  mu0[obs_to_t[n]] - y[n],
                                  sigma[obs_to_t[n]]);
      if (delta_sat[n_s] < 0) {
        target += negative_infinity();
      }
      else {
        target += -normal_mix_lccdf(0 | omega_ag[obs_to_ag[n]] * omega_t[obs_to_t[n]],
                                    mu0[obs_to_t[n]] + mu_ag[obs_to_ag[n]][obs_to_t[n]] - y[n],
                                    mu0[obs_to_t[n]] - y[n],
                                    sigma[obs_to_t[n]]);
      }

    }
  }

  for (n in 1:N) {
      target += log_sum_exp(soft_z[n]);
  }
}

// compute group membership probabilities
// z[i,1] is probability that the observation is from the zero-centered group
// z[i,2] is probability that the observation is from the non-zero group
generated quantities {
  vector<lower=0, upper=1>[2] z[N];
  for (n in 1:N) {
    z[n] = softmax(soft_z[n]);
  }
}
