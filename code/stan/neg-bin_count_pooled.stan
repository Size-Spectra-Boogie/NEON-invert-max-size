functions {
  /*
   * Quantile function for a normal distribution truncated below at zero.
   * p is a probability on the truncated-normal scale.
   */
  real positive_trunc_normal_quantile(
      real p,
      real mu,
      real sigma
  ) {
    real p_zero;
    real p_normal;

    p_zero = normal_cdf(0 | mu, sigma);

    // Convert truncated-normal probability to ordinary-normal probability
    p_normal = p_zero + p * (1 - p_zero);

    // Protect inv_Phi() against numerical values exactly equal to 0 or 1
    p_normal = fmin(
      1 - 1e-12,
      fmax(1e-12, p_normal)
    );

    return mu + sigma * inv_Phi(p_normal);
  }
}

data {
  int<lower=1> S;                 // Number of sites
  int<lower=1> K;                 // Total number of sampling events
  int<lower=1> n_obs;             // Total number of measured individuals

  // Individual body sizes, contiguous within sampling events
  vector<lower=0>[n_obs] x;

  // Number collected during each event; zero counts are permitted
  array[K] int<lower=0> n_per_sample;

  /*
   * Starting index in x for each event.
   * For an event with zero observations, this index is not used.
   */
  array[K] int<lower=1> start_idx;

  // Site associated with each sampling event
  array[K] int<lower=1, upper=S> site_id;

  // Standardized number of events used for maximum-size prediction
  int<lower=1> k_ref;
}

transformed data {
  /*
   * Every observed individual must belong to exactly one event.
   */
  if (sum(n_per_sample) != n_obs) {
    reject(
      "sum(n_per_sample) must equal n_obs. ",
      "sum(n_per_sample) = ", sum(n_per_sample),
      "; n_obs = ", n_obs
    );
  }

  /*
   * Check event indexing.
   */
  for (j in 1:K) {
    if (n_per_sample[j] > 0) {
      int end_idx =
        start_idx[j] + n_per_sample[j] - 1;

      if (end_idx > n_obs) {
        reject(
          "Sampling event ", j,
          " extends beyond the end of x."
        );
      }
    }
  }
}

parameters {
  /*
   * Across-site log-scale centers.
   *
   * exp(alpha_log_mu) is the median mu across sites,
   * conditional on the hyperparameters.
   */
  real alpha_log_mu;
  real alpha_log_sigma;
  real alpha_log_lambda;

  // Among-site heterogeneity on the log scale
  real<lower=0> tau_log_mu;
  real<lower=0> tau_log_sigma;
  real<lower=0> tau_log_lambda;

  // Non-centered site deviations
  vector[S] z_mu;
  vector[S] z_sigma;
  vector[S] z_lambda;

  // Negative-binomial dispersion on the log scale
  // this is not estimated for each site but is shared
  // across sites
  real log_phi;
}

transformed parameters {
  vector<lower=0>[S] mu;
  vector<lower=0>[S] sigma;
  vector<lower=0>[S] lambda;

  real<lower=0> phi;

  mu = exp(
    alpha_log_mu +
    tau_log_mu * z_mu
  );

  sigma = exp(
    alpha_log_sigma +
    tau_log_sigma * z_sigma
  );

  lambda = exp(
    alpha_log_lambda +
    tau_log_lambda * z_lambda
  );

  phi = exp(log_phi);
}

model {
  /*
   * Hyperpriors
   *
   * Across-site median mu is centered at 20.
   * A log SD of 0.75 gives a long upper tail while retaining
   * meaningful prior concentration near 20.
   */
  alpha_log_mu ~ normal(log(20), 0.75);

  /*
   * Broad prior for the across-site center of sigma.
   * Change the center if there is stronger prior information.
   */
  alpha_log_sigma ~ normal(log(10), 0.80);

  /*
   * lambda is now the expected number collected per standardized event,
   * not an unobserved total population size.
   *
   * This prior is centered at 200 individuals per event but is broad.
   */
  alpha_log_lambda ~ normal(log(200), 1.00);

  // Among-site variation
  tau_log_mu ~ normal(0, 0.50);
  tau_log_sigma ~ normal(0, 0.50);
  tau_log_lambda ~ normal(0, 0.75);

  // Non-centered site effects
  z_mu ~ std_normal();
  z_sigma ~ std_normal();
  z_lambda ~ std_normal();

  /*
   * Shared negative-binomial dispersion.
   *
   * Median phi = 20, with a broad prior. Smaller phi means more
   * event-to-event overdispersion; large phi approaches Poisson.
   */
  log_phi ~ normal(log(20), 1.00);

  /*
   * Joint count and body-size likelihood.
   */
  for (j in 1:K) {
    int s = site_id[j];

    /*
     * Event count model:
     *
     * E[n_j]   = lambda[s]
     * Var[n_j] = lambda[s] + lambda[s]^2 / phi
     */
    n_per_sample[j] ~
      neg_binomial_2(lambda[s], phi);

    /*
     * All individual body sizes inform the site-specific size
     * distribution. The likelihood is conditioned on x > 0.
     */
    if (n_per_sample[j] > 0) {
      target += normal_lpdf(
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
        normal_lccdf(
          0 |
          mu[s],
          sigma[s]
        );
    }
  }
}

generated quantities {
  /*
   * Event-level log likelihoods.
   * The combined likelihood uses the sampling event as the
   * pointwise unit for LOO or related diagnostics.
   */
  vector[K] log_lik;
  vector[K] log_lik_count;
  vector[K] log_lik_size;

  // Posterior-predictive count for each observed event
  array[K] int n_rep;

  /*
   * Standardized predictions over k_ref events.
   */
  vector[S] expected_n_ref;
  array[S] int n_ref_rep;
  vector[S] max_ref_rep;

  for (j in 1:K) {
    int s = site_id[j];

    log_lik_count[j] =
      neg_binomial_2_lpmf(
        n_per_sample[j] |
        lambda[s],
        phi
      );

    if (n_per_sample[j] > 0) {
      log_lik_size[j] =
        normal_lpdf(
          segment(
            x,
            start_idx[j],
            n_per_sample[j]
          ) |
          mu[s],
          sigma[s]
        ) -
        n_per_sample[j] *
        normal_lccdf(
          0 |
          mu[s],
          sigma[s]
        );
    } else {
      log_lik_size[j] = 0;
    }

    log_lik[j] =
      log_lik_count[j] +
      log_lik_size[j];

    n_rep[j] =
      neg_binomial_2_rng(
        lambda[s],
        phi
      );
  }

  /*
   * Draw a total count and realized maximum for k_ref new
   * equal-effort sampling events at each site.
   */
  for (s in 1:S) {
    expected_n_ref[s] =
      k_ref * lambda[s];

    n_ref_rep[s] = 0;

    for (r in 1:k_ref) {
      n_ref_rep[s] +=
        neg_binomial_2_rng(
          lambda[s],
          phi
        );
    }

    if (n_ref_rep[s] > 0) {
      real u;
      real p_max;

      /*
       * Conditional on N, the maximum CDF is F(x)^N.
       * Therefore, F(M) = U^(1/N).
       */
      u = uniform_rng(
        1e-12,
        1 - 1e-12
      );

      p_max =
        exp(
          log(u) /
          n_ref_rep[s]
        );

      max_ref_rep[s] =
        positive_trunc_normal_quantile(
          p_max,
          mu[s],
          sigma[s]
        );
    } else {
      /*
       * There is no defined maximum if all k_ref events contain
       * zero individuals. Zero is used as a sentinel value.
       */
      max_ref_rep[s] = 0;
    }
  }
}
