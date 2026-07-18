functions {
  // Quantile of a Normal(mu, sigma) distribution truncated below at zero.
  real positive_trunc_normal_quantile(real p, real mu, real sigma) {
    real p_zero = normal_cdf(0 | mu, sigma);
    real p_normal = p_zero + p * (1 - p_zero);

    // Protect inv_Phi() from probabilities numerically equal to 0 or 1.
    p_normal = fmin(1 - 1e-12, fmax(1e-12, p_normal));

    return mu + sigma * inv_Phi(p_normal);
  }

  /*
   * Quantile of the maximum body size across k_ref independent events,
   * conditional on at least one individual being collected.
   *
   * Event counts follow NegBinomial2(mean = exp(log_lambda), shape = phi).
   * The sum across k_ref events is NegBinomial2 with
   * mean = k_ref * exp(log_lambda) and shape = k_ref * phi.
   */
  real nb_max_quantile_conditional(
      real q,
      real mu,
      real sigma,
      real log_lambda,
      real phi,
      int k_ref
  ) {
    real total_shape = k_ref * phi;

    // Probability of zero individuals across all k_ref events.
    real log_p_zero =
      -total_shape * log1p_exp(log_lambda - log(phi));

    // Unconditional maximum CDF value corresponding to conditional quantile q.
    real log_target_cdf = log_sum_exp(
      log_p_zero,
      log(q) + log1m_exp(log_p_zero)
    );

    // Solve the negative-binomial probability-generating function for F(x).
    real one_minus_individual_cdf =
      exp(log(phi) - log_lambda) *
      expm1(-log_target_cdf / total_shape);

    real individual_cdf = 1 - one_minus_individual_cdf;

    individual_cdf =
      fmin(1 - 1e-12, fmax(1e-12, individual_cdf));

    return positive_trunc_normal_quantile(
      individual_cdf,
      mu,
      sigma
    );
  }
}

data {
  int<lower=1> S;                 // Number of sites
  int<lower=1> K;                 // Total number of sampling events
  int<lower=1> n_obs;             // Total number of measured individuals

  // Individual body sizes, contiguous within sampling events.
  vector<lower=0>[n_obs] x;

  // Number collected in each event; zero-count events are allowed.
  array[K] int<lower=0> n_per_sample;

  /*
   * Starting position in x for each event. Zero-count events use the
   * next available index and therefore may share an index with the next event.
   */
  array[K] int<lower=1, upper=n_obs + 1> start_idx;

  // Site associated with each event.
  array[K] int<lower=1, upper=S> site_id;

  // Number of equal-effort events used for standardized maximum predictions.
  int<lower=1> k_ref;
}

transformed data {
  int next_idx = 1;

  if (sum(n_per_sample) != n_obs) {
    reject(
      "sum(n_per_sample) must equal n_obs; received ",
      sum(n_per_sample), " and ", n_obs
    );
  }

  // Require x to be partitioned into contiguous, non-overlapping event blocks.
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
  // Across-site log-scale centers.
  real alpha_log_mu;
  real alpha_log_sigma;
  real alpha_log_lambda;

  // Among-site standard deviations on the log scale.
  real<lower=0> tau_log_mu;
  real<lower=0> tau_log_sigma;
  real<lower=0> tau_log_lambda;

  // Non-centered site deviations.
  vector[S] z_mu;
  vector[S] z_sigma;
  vector[S] z_lambda;

  /*
   * Shared negative-binomial shape/precision on the log scale.
   * The bounds are extremely broad but prevent numerical overflow.
   */
  real<lower=-4, upper=12> log_phi;
}

transformed parameters {
  vector<lower=0>[S] mu;
  vector<lower=0>[S] sigma;
  vector[S] log_lambda;
  real<lower=0> phi;

  mu = exp(alpha_log_mu + tau_log_mu * z_mu);
  sigma = exp(alpha_log_sigma + tau_log_sigma * z_sigma);

  // Keep expected event counts on the log scale inside the likelihood.
  log_lambda = alpha_log_lambda + tau_log_lambda * z_lambda;

  phi = exp(log_phi);
}

model {
  // Across-site size-distribution priors.
  alpha_log_mu ~ normal(log(20), 0.75);
  alpha_log_sigma ~ normal(log(10), 1.00);

  // Allows typical event counts spanning roughly tens to tens of thousands.
  alpha_log_lambda ~ normal(log(500), 2.00);

  // Among-site heterogeneity.
  tau_log_mu ~ normal(0, 0.50);
  tau_log_sigma ~ normal(0, 0.75);
  tau_log_lambda ~ normal(0, 1.50);

  // Non-centered site effects.
  z_mu ~ std_normal();
  z_sigma ~ std_normal();
  z_lambda ~ std_normal();

  // Shared within-site count overdispersion.
  log_phi ~ normal(log(20), 1.00);

  // Joint event-count and individual-size likelihood.
  for (j in 1:K) {
    int s = site_id[j];

    // E[n_j] = exp(log_lambda[s]); Var[n_j] = mean + mean^2 / phi.
    target += neg_binomial_2_log_lpmf(
      n_per_sample[j] | log_lambda[s], phi
    );

    // Every measured size informs the positive-truncated Normal distribution.
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
  // Event-level log likelihoods for diagnostics or LOO.
  vector[K] log_lik;
  vector[K] log_lik_count;
  vector[K] log_lik_size;

  // Site-level quantities on interpretable scales.
  vector[S] lambda;
  vector[S] expected_n_ref;
  vector[S] prob_zero_ref;

  // Mean and median of the fitted positive-truncated size distribution.
  vector[S] size_mean;
  vector[S] size_median;

  /*
   * Conditional quantiles and one posterior-predictive draw of the maximum
   * across k_ref events. These integrate over negative-binomial count variation
   * and are conditional on at least one individual being collected.
   */
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
      0.5,
      mu[s],
      sigma[s]
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
