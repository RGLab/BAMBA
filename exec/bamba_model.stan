// BAMBA model
// 2-component grouped mixture model for BAMA or Fc array data
// 1 component is centered at mu0 (unique per timepoint)
// The other component is centered at a group mean
// which is mu0 + an antigen (ag) effect + an Fc variable (re) effect (both positive)
// Each observation comes either from the group centered at mu0
// or from the particular group it is a part of
// given by ag and re
// Prior probability of is a product of group, ag, re, and timepoint probabilities

// y_i = [y_{i1}, ..., y_{iT}]^T
// f(y_{it} | \theta) = N(\mu_{0t}, \sigma_t)^{1-z_{it}} N(\mu_{0t} + mu_{art}, \sigma_t)^{z_{it}}
// z_{it} ~ Bernoulli(\omega_{grp} * \omega_{ag} * \omega_{re} * \omega_t)
// \mu_{0t} ~ N(0, 0.25) (variance is 0.25)
// \omega_{grp} ~ Beta(a_{grp}, b_{grp})
// a_{grp} = \lambda_{grp} \phi_{grp}
// b_{grp} = \lambda_{grp} (1 - \phi_{grp})
// \lambda_{grp} ~ Pareto(0.1, 1.5)
// \phi_{grp} ~ Beta(1, 1)
// \omega_{ag} ~ Beta(a_{ag}, b_{ag})
// a_{ag} = \lambda_{ag} \phi_{ag}
// b_{ag} = \lambda_{ag} (1 - \phi_{ag})
// \lambda_{ag} ~ Pareto(0.1, 1.5)
// \phi_{ag} ~ Beta(1, 1)
// \omega_{re} ~ Beta(a_{re}, b_{re})
// a_{re} = \lambda_{re} \phi_{re}
// b_{re} = \lambda_{re} (1 - \phi_{re})
// \lambda_{re} ~ Pareto(0.1, 1.5)
// \phi_{re} ~ Beta(1, 1)
// \omega_t ~ Beta(a_t, b_t)
// \mu_{\art} = \mu_{ag t} + \mu_{re t}
// \mu_{\art} ~ Gamma(\alpha_{ar}, \beta_{ar}
// \alpha_{ar} ~ Uniform(0, 1000)
// \beta_{ar} ~ Uniform(0, 1000)
// \sigma_t ~ Cauchy+(0, 10)

data {
  int<lower=1> N; // Number of observations
  int<lower=1> T; // Number of timepoints (dimension of observations)
  int<lower=1> obs_to_t[N];
  int<lower=1> N_grp; // grouping for mixture probability
  int<lower=1> obs_to_grp[N];
  int<lower=1> N_ag; // grouping for antigen
  int<lower=1> obs_to_ag[N];
  int<lower=1> N_re; // grouping for Fc variable
  int<lower=1> obs_to_re[N];
  int<lower=0> N_sat; // Number of saturated observations
  int<lower=1> sat_ind[N_sat]; // indexes of saturated observations
  vector[N] y;
}

parameters {
  // Overall probability of an observation being in a non-zero group
  real<lower=0, upper=1> omega_grp[N_grp];
  real<lower=0, upper=1> omega_ag[N_ag];
  real<lower=0, upper=1> omega_re[N_re];
  real<lower=0, upper=1> omega_t[T];
  
  vector<lower=0>[T] sigma; // Group variance (stddev)
  real<lower=0> alpha_ag; // mean parameter for mu_ag gamma prior
  real<lower=0> beta_ag; // variance parameter for mu_ag gamma prior
  real<lower=0> alpha_re; // mean parameter for mu_re gamma prior
  real<lower=0> beta_re; // variance parameter for mu_re gamma prior
  matrix<lower=0>[N_ag, T] mu_ag; // antigen response mean
  matrix<lower=0>[N_re, T] mu_re; // Fc variable response mean
  real<lower=0, upper=1> phi_grp; // hyperprior for omega (as found in Stan manual)
  real<lower=0.1> lambda_grp; // hyperprior for omega (as found in Stan manual)
  real<lower=0, upper=1> phi_ag; // hyperprior for omega (as found in Stan manual)
  real<lower=0.1> lambda_ag; // hyperprior for omega (as found in Stan manual)
  real<lower=0, upper=1> phi_re; // hyperprior for omega (as found in Stan manual)
  real<lower=0.1> lambda_re; // hyperprior for omega (as found in Stan manual)
  real<lower=0, upper=1> phi_t; // hyperprior for omega (as found in Stan manual)
  real<lower=0.1> lambda_t; // hyperprior for omega (as found in Stan manual)
  vector[T] mu0;
  
  vector<lower=0>[N_sat] delta_sat; // extra to add onto saturated observations
  real<lower=0> alpha_sat; // mean parameter for delta_sat gamma prior
  real<lower=0> beta_sat; // variance parameter for delta_sat gamma prior
}

transformed parameters {
  real<lower=0> a_grp;
  real<lower=0> b_grp;
  real<lower=0> a_ag;
  real<lower=0> b_ag;
  real<lower=0> a_re;
  real<lower=0> b_re;
  real<lower=0> a_t;
  real<lower=0> b_t;
  vector<upper=0>[2] soft_z[N]; // used to track the group membership
  matrix<lower=0>[N_ag, N_re] mu_ar[T];
  
  vector[N] ystar;

  a_grp = lambda_grp * phi_grp;
  b_grp = lambda_grp * (1 - phi_grp);
  a_ag = lambda_ag * phi_ag;
  b_ag = lambda_ag * (1 - phi_ag);
  a_re = lambda_re * phi_re;
  b_re = lambda_re * (1 - phi_re);
  a_t = lambda_t * phi_t;
  b_t = lambda_t * (1 - phi_t);

  for (n_a in 1:N_ag) {
    for (n_r in 1:N_re) {
      for (t in 1:T) {
        mu_ar[t][n_a, n_r] = mu_ag[n_a, t] + mu_re[n_r, t];
      }
    }
  }
  
  ystar = y;
  if (N_sat > 0) {
    for (n_s in 1:N_sat) {
      ystar[sat_ind[n_s]] = ystar[sat_ind[n_s]] + delta_sat[n_s];
    }
  }

  // Calculate log probabilities for each observation in the zero and non-zero groups
  for (n in 1:N) {
      soft_z[n][1] = log(1-omega_grp[obs_to_grp[n]] * omega_ag[obs_to_ag[n]] *
                             omega_re[obs_to_re[n]] * omega_t[obs_to_t[n]]) +
			normal_lpdf(ystar[n] | mu0[obs_to_t[n]], sigma[obs_to_t[n]]);
      soft_z[n][2] = log(omega_grp[obs_to_grp[n]] * omega_ag[obs_to_ag[n]] *
                           omega_re[obs_to_re[n]] * omega_t[obs_to_t[n]]) +
                    normal_lpdf(ystar[n] | mu0[obs_to_t[n]] +
		                       mu_ar[obs_to_t[n]][obs_to_ag[n], obs_to_re[n]],
				       sigma[obs_to_t[n]]);
  }
}

model {
  sigma ~ cauchy(0, 10);
  phi_grp ~ beta(1, 1);          // as in Stan manual
  lambda_grp ~ pareto(0.1, 1.5); // as in Stan manual
  phi_ag ~ beta(1, 1);          // as in Stan manual
  lambda_ag ~ pareto(0.1, 1.5); // as in Stan manual
  phi_re ~ beta(1, 1);          // as in Stan manual
  lambda_re ~ pareto(0.1, 1.5); // as in Stan manual
  phi_t ~ beta(1, 1);          // as in Stan manual
  lambda_t ~ pareto(0.1, 1.5); // as in Stan manual
  omega_grp ~ beta(a_grp, b_grp);
  omega_ag ~ beta(a_ag, b_ag);
  omega_re ~ beta(a_re, b_re);
  omega_t ~ beta(a_t, b_t);
  alpha_ag ~ uniform(0, 1000);
  beta_ag ~ uniform(0, 1000);
  alpha_re ~ uniform(0, 1000);
  beta_re ~ uniform(0, 1000);
  alpha_sat ~ uniform(0, 1000);
  beta_sat ~ uniform(0, 1000);

  mu0 ~ normal(0, 0.5);

  for (n_a in 1:N_ag) {
    mu_ag[n_a,] ~ gamma(alpha_ag, beta_ag);
  }

  for (n_r in 1:N_re) {
    mu_re[n_r,] ~ gamma(alpha_re, beta_re);
  }
  
  delta_sat ~ gamma(alpha_sat, beta_sat);

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
