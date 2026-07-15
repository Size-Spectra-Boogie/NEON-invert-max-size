data {
  int<lower=1> S;               // Number of sites
  int<lower=1> K;               // Total number of sampling events
  int<lower=1> n_obs;           // Total number of measured individuals

  // Individual sizes, arranged contiguously within sampling events
  vector<lower=0>[n_obs] x;

  // Sampling-event information
  array[K] int<lower=1> n_per_sample;
  array[K] int<lower=1, upper=n_obs> start_idx;

  // Site associated with each sampling event
  array[K] int<lower=1, upper=S> site_id;
}

transformed data {
  vector[K] smallest_in_sample;

  /*
   * Calculate the smallest recorded individual in each sampling event.
   * The observations belonging to event j must occupy:
   *
   * start_idx[j] :
   * (start_idx[j] + n_per_sample[j] - 1)
   */
  for (j in 1:K) {
    int end_idx =
      start_idx[j] + n_per_sample[j] - 1;

    if (end_idx > n_obs) {
      reject(
        "Sampling event ", j,
        " extends beyond the end of x."
      );
    }

    smallest_in_sample[j] =
      min(
        segment(
          x,
          start_idx[j],
          n_per_sample[j]
        )
      );
  }
}

parameters {
  /*
   * Population-level arithmetic means across sites.
   * These are sampled on the log scale.
   */
  real log_mu_pop;
  real log_sigma_pop;
  real log_lambda_pop;

  /*
   * Among-site standard deviations on the log scale.
   */
  real<lower=0> tau_log_mu;
  real<lower=0> tau_log_sigma;
  real<lower=0> tau_log_lambda;

  /*
   * Non-centered standardized site effects.
   */
  vector[S] z_mu;
  vector[S] z_sigma;
  vector[S] z_lambda;
}

transformed parameters {
  /*
   * Population-level arithmetic means.
   */
  real<lower=0> mu_pop;
  real<lower=0> sigma_pop;
  real<lower=0> lambda_pop;

  /*
   * Site-specific parameters.
   */
  vector<lower=0>[S] mu;
  vector<lower=0>[S] sigma;
  vector<lower=0>[S] lambda;

  mu_pop = exp(log_mu_pop);
  sigma_pop = exp(log_sigma_pop);
  lambda_pop = exp(log_lambda_pop);

  /*
   * The -0.5 * tau^2 correction means:
   *
   * E(mu[s] | mu_pop, tau_log_mu) = mu_pop
   *
   * with equivalent interpretations for sigma and lambda.
   */
  mu =
    exp(
      rep_vector(
        log_mu_pop - 0.5 * square(tau_log_mu),
        S
      ) +
      tau_log_mu * z_mu
    );

  sigma =
    exp(
      rep_vector(
        log_sigma_pop - 0.5 * square(tau_log_sigma),
        S
      ) +
      tau_log_sigma * z_sigma
    );

  lambda =
    exp(
      rep_vector(
        log_lambda_pop - 0.5 * square(tau_log_lambda),
        S
      ) +
      tau_log_lambda * z_lambda
    );
}

model {
  /*
   * Population-level priors.
   *
   * The prior for the across-site mean of mu is centered at 20,
   * but the lognormal form allows substantially larger values.
   */
  log_mu_pop ~ normal(log(20), 0.60);

  /*
   * These retain approximately the centers of the original priors.
   * They should be adjusted if sigma or lambda differ substantially
   * among taxa or sampling methods.
   */
  log_sigma_pop ~ normal(log(20), 0.50);
  log_lambda_pop ~ normal(log(10000), 0.50);

  /*
   * Half-normal priors on among-site heterogeneity.
   */
  tau_log_mu ~ normal(0, 0.50);
  tau_log_sigma ~ normal(0, 0.40);
  tau_log_lambda ~ normal(0, 0.50);

  /*
   * Non-centered site effects.
   */
  z_mu ~ std_normal();
  z_sigma ~ std_normal();
  z_lambda ~ std_normal();

  /*
   * Exact finite-sampling likelihood.
   */
  for (j in 1:K) {
    int s = site_id[j];

    /*
     * Probability that an untruncated normal observation exceeds zero.
     *
     * This is the normalizing constant for the positive-truncated
     * normal distribution.
     */
    real log_norm_const =
      normal_lccdf(0 | mu[s], sigma[s]);

    /*
     * Numerator of the positive-truncated normal CDF evaluated at
     * the smallest observed size:
     *
     * log[P(0 < X <= smallest)]
     */
    real log_cdf_numerator =
      log_diff_exp(
        normal_lcdf(
          smallest_in_sample[j] |
          mu[s],
          sigma[s]
        ),
        normal_lcdf(
          0 |
          mu[s],
          sigma[s]
        )
      );

    /*
     * Conditional truncated-normal CDF:
     *
     * P(X <= smallest | X > 0)
     */
    real trunc_cdf =
      exp(
        log_cdf_numerator -
        log_norm_const
      );

    /*
     * Original likelihood contribution:
     *
     * -lambda + n_j log(lambda)
     * + lambda F_T(x_min)
     * + sum_i log f_T(x_i)
     */
    target += -lambda[s];

    target +=
      n_per_sample[j] *
      log(lambda[s]);

    target +=
      lambda[s] *
      trunc_cdf;

    target +=
      normal_lpdf(
        segment(
          x,
          start_idx[j],
          n_per_sample[j]
        ) |
        mu[s],
        sigma[s]
      );

    target +=
      -n_per_sample[j] *
      log_norm_const;
  }
}

generated quantities {
  /*
   * Event-level log likelihood, useful for diagnostics and LOO.
   */
  vector[K] log_lik;

  for (j in 1:K) {
    int s = site_id[j];

    real log_norm_const =
      normal_lccdf(0 | mu[s], sigma[s]);

    real log_cdf_numerator =
      log_diff_exp(
        normal_lcdf(
          smallest_in_sample[j] |
          mu[s],
          sigma[s]
        ),
        normal_lcdf(
          0 |
          mu[s],
          sigma[s]
        )
      );

    real trunc_cdf =
      exp(
        log_cdf_numerator -
        log_norm_const
      );

    log_lik[j] =
      -lambda[s] +
      n_per_sample[j] * log(lambda[s]) +
      lambda[s] * trunc_cdf +
      normal_lpdf(
        segment(
          x,
          start_idx[j],
          n_per_sample[j]
        ) |
        mu[s],
        sigma[s]
      ) -
      n_per_sample[j] * log_norm_const;
  }
}
