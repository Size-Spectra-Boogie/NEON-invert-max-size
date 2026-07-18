functions {
  real positive_trunc_normal_quantile(real p, real mu, real sigma) {
    real p_zero = normal_cdf(0 | mu, sigma);
    real p_normal = p_zero + p * (1 - p_zero);

    p_normal = fmin(1 - 1e-12, fmax(1e-12, p_normal));
    return mu + sigma * inv_Phi(p_normal);
  }

  real nb_max_quantile_conditional(
      real q,
      real mu,
      real sigma,
      real log_lambda,
      real phi,
      int k_ref
  ) {
    real total_shape = k_ref * phi;
    real log_p_zero =
      -total_shape * log1p_exp(log_lambda - log(phi));
    real log_target_cdf = log_sum_exp(
      log_p_zero,
      log(q) + log1m_exp(log_p_zero)
    );
    real one_minus_individual_cdf =
      exp(log(phi) - log_lambda) *
      expm1(-log_target_cdf / total_shape);
    real individual_cdf = 1 - one_minus_individual_cdf;

    individual_cdf =
      fmin(1 - 1e-12, fmax(1e-12, individual_cdf));

    return positive_trunc_normal_quantile(
      individual_cdf, mu, sigma
    );
  }
}

data {
  int<lower=1> S;
  int<lower=1> K;
  int<lower=1> n_obs;

  vector<lower=0>[n_obs] x;
  array[K] int<lower=0> n_per_sample;
  array[K] int<lower=1, upper=n_obs + 1> start_idx;
  array[K] int<lower=1, upper=S> site_id;

  int<lower=1> k_ref;

  /*
   * Generous, scientifically plausible upper bounds for site-level
   * truncated-normal parameters. These prevent exp(log_parameter)
   * from overflowing during warmup.
   */
  real<lower=1e-6> mu_upper;
  real<lower=1e-6> sigma_upper;

/*
 * create a switch to sample from the prior for checks
 */
  int<lower=0, upper=1> prior_only;
}

transformed data {
  int next_idx = 1;

  if (sum(n_per_sample) != n_obs) {
    reject(
      "sum(n_per_sample) must equal n_obs; received ",
      sum(n_per_sample), " and ", n_obs
    );
  }

  for (j in 1:K) {
    if (start_idx[j] != next_idx) {
      reject(
        "start_idx is inconsistent at event ", j,
        "; expected ", next_idx,
        " but received ", start_idx[j]
      );
    }
    next_idx += n_per_sample[j];
  }

  if (next_idx != n_obs + 1) {
    reject("Event indexing does not exhaust x.");
  }
}

parameters {
  /* Across-site centers on the log-size scale. */
  real alpha_log_mu;
  real alpha_log_sigma;

  /*
   * Site parameters are sampled directly on bounded log scales.
   * This centered hierarchy is appropriate because each site usually
   * contributes many individual size observations.
   */
  vector<lower=log(1e-8), upper=log(mu_upper)>[S] log_mu_site;
  vector<lower=log(1e-8), upper=log(sigma_upper)>[S] log_sigma_site;

  real<lower=1e-4, upper=2.5> tau_log_mu;
  real<lower=1e-4, upper=2.5> tau_log_sigma;

  /* Count hierarchy remains non-centered. */
  real alpha_log_lambda;
  real<lower=0, upper=3.5> tau_log_lambda;
  vector[S] z_lambda;

  real<lower=-4, upper=12> log_phi;
}

transformed parameters {
  vector<lower=0, upper=mu_upper>[S] mu = exp(log_mu_site);
  vector<lower=0, upper=sigma_upper>[S] sigma = exp(log_sigma_site);
  vector[S] log_lambda =
    alpha_log_lambda + tau_log_lambda * z_lambda;
  real<lower=0> phi = exp(log_phi);
}

model {
  /* Hyperpriors for the across-site body-size distribution. */
  alpha_log_mu ~ student_t(3, log(3), 0.5);
  alpha_log_sigma ~ normal(log(10), 1.00);

  tau_log_mu ~ normal(0, 0.35);
  tau_log_sigma ~ normal(0, 0.75);

  /*
   * Correctly normalized bounded hierarchical priors on log(mu_s)
   * and log(sigma_s). The normalization depends on the hyperparameters
   * and therefore must be included.
   */
  for (s in 1:S) {
    target += normal_lpdf(
      log_mu_site[s] | alpha_log_mu, tau_log_mu
    ) - log_diff_exp(
      normal_lcdf(log(mu_upper) | alpha_log_mu, tau_log_mu),
      normal_lcdf(log(1e-8) | alpha_log_mu, tau_log_mu)
    );

    target += normal_lpdf(
      log_sigma_site[s] | alpha_log_sigma, tau_log_sigma
    ) - log_diff_exp(
      normal_lcdf(log(sigma_upper) | alpha_log_sigma, tau_log_sigma),
      normal_lcdf(log(1e-8) | alpha_log_sigma, tau_log_sigma)
    );
  }

  /* Broad count hierarchy for large differences among sites. */
  alpha_log_lambda ~ normal(log(500), 2.00);
  tau_log_lambda ~ normal(0, 1.50);
  z_lambda ~ std_normal();

  log_phi ~ normal(log(20), 1.00);

  for (j in 1:K) {
    int s = site_id[j];

    target += neg_binomial_2_log_lpmf(
      n_per_sample[j] | log_lambda[s], phi
    );

    if (n_per_sample[j] > 0) {
      target += normal_lpdf(
        segment(x, start_idx[j], n_per_sample[j]) |
        mu[s], sigma[s]
      );

      target += -n_per_sample[j] *
        normal_lccdf(0 | mu[s], sigma[s]);
    }
  }
}

generated quantities {
  vector[K] log_lik;
  vector[K] log_lik_count;
  vector[K] log_lik_size;

  vector[S] lambda;
  vector[S] expected_n_ref;
  vector[S] prob_zero_ref;

  vector[S] size_mean;
  vector[S] size_median;

  vector[S] max_ref_q025;
  vector[S] max_ref_q50;
  vector[S] max_ref_q975;
  vector[S] max_ref_rep;

  for (j in 1:K) {
    int s = site_id[j];

    log_lik_count[j] = neg_binomial_2_log_lpmf(
      n_per_sample[j] | log_lambda[s], phi
    );

    if (n_per_sample[j] > 0) {
      log_lik_size[j] = normal_lpdf(
        segment(x, start_idx[j], n_per_sample[j]) |
        mu[s], sigma[s]
      ) - n_per_sample[j] *
        normal_lccdf(0 | mu[s], sigma[s]);
    } else {
      log_lik_size[j] = 0;
    }

    log_lik[j] = log_lik_count[j] + log_lik_size[j];
  }

  for (s in 1:S) {
    real a = -mu[s] / sigma[s];
    real inverse_mills = exp(
      std_normal_lpdf(a) - std_normal_lccdf(a)
    );
    real log_prob_zero_count =
      -(k_ref * phi) * log1p_exp(log_lambda[s] - log(phi));

    lambda[s] = exp(log_lambda[s]);
    expected_n_ref[s] = k_ref * lambda[s];
    prob_zero_ref[s] = exp(log_prob_zero_count);

    size_mean[s] = mu[s] + sigma[s] * inverse_mills;
    size_median[s] = positive_trunc_normal_quantile(
      0.5, mu[s], sigma[s]
    );

    max_ref_q025[s] = nb_max_quantile_conditional(
      0.025, mu[s], sigma[s], log_lambda[s], phi, k_ref
    );
    max_ref_q50[s] = nb_max_quantile_conditional(
      0.5, mu[s], sigma[s], log_lambda[s], phi, k_ref
    );
    max_ref_q975[s] = nb_max_quantile_conditional(
      0.975, mu[s], sigma[s], log_lambda[s], phi, k_ref
    );
    max_ref_rep[s] = nb_max_quantile_conditional(
      uniform_rng(1e-12, 1 - 1e-12),
      mu[s], sigma[s], log_lambda[s], phi, k_ref
    );
  }
}
