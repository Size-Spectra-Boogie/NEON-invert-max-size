functions {
  /*
   * Quantile function for a normal distribution truncated below at zero.
   *
   * p is a probability on the positive-truncated-normal scale.
   */
  real positive_trunc_normal_quantile(
      real p,
      real mu,
      real sigma
  ) {
    real p_zero =
      normal_cdf(0 | mu, sigma);

    real p_normal =
      p_zero +
      p * (1 - p_zero);

    /*
     * Protect inv_Phi() against probabilities that round
     * exactly to zero or one.
     */
    p_normal =
      fmin(
        1 - 1e-12,
        fmax(1e-12, p_normal)
      );

    return mu + sigma * inv_Phi(p_normal);
  }

  /*
   * Quantile of the maximum across k_ref negative-binomial
   * sampling events, conditional on at least one individual
   * being observed.
   *
   * Count model for one event:
   *
   * N_r ~ NegBinomial2(lambda, phi)
   *
   * with log_lambda = log(lambda).
   */
  real nb_max_quantile_conditional(
      real q,
      real mu,
      real sigma,
      real log_lambda,
      real phi,
      int k_ref
  ) {
    real total_shape =
      k_ref * phi;

    /*
     * Probability that all k_ref events contain zero
     * individuals.
     */
    real log_p_zero =
      -total_shape *
      log1p_exp(
        log_lambda -
        log(phi)
      );

    /*
     * Convert the conditional maximum probability q into
     * the corresponding unconditional maximum CDF.
     */
    real log_target_cdf =
      log_sum_exp(
        log_p_zero,
        log(q) +
        log1m_exp(log_p_zero)
      );

    /*
     * Recover the individual-size CDF associated with the
     * requested maximum quantile.
     */
    real one_minus_individual_cdf =
      exp(
        log(phi) -
        log_lambda
      ) *
      expm1(
        -log_target_cdf /
        total_shape
      );

    real individual_cdf =
      1 -
      one_minus_individual_cdf;

    individual_cdf =
      fmin(
        1 - 1e-12,
        fmax(
          1e-12,
          individual_cdf
        )
      );

    return positive_trunc_normal_quantile(
      individual_cdf,
      mu,
      sigma
    );
  }
}

data {
  /*
   * This model must receive exactly one site.
   *
   * Retaining S and site_id allows the same data interface
   * and output naming as the hierarchical model.
   */
  int<lower=1, upper=1> S;

  int<lower=1> K;
  int<lower=1> n_obs;

  /*
   * Individual body sizes, contiguous within events.
   */
  vector<lower=0>[n_obs] x;

  /*
   * Number collected during each event.
   */
  array[K] int<lower=0> n_per_sample;

  /*
   * Starting location of each event in x.
   *
   * For a zero-count event, start_idx points to the next
   * available position but is not used in the likelihood.
   */
  array[K] int<
    lower=1,
    upper=n_obs + 1
  > start_idx;

  /*
   * All values must equal one because S = 1.
   */
  array[K] int<
    lower=1,
    upper=S
  > site_id;

  /*
   * Standardized number of sampling events used for
   * maximum-size prediction.
   */
  int<lower=1> k_ref;

  /*
   * Upper safeguards for the latent normal location and
   * scale parameters.
   */
  real<lower=1e-6> mu_upper;
  real<lower=1e-6> sigma_upper;

  /*
   * prior_only = 0: posterior sampling
   * prior_only = 1: prior-only sampling
   */
  int<lower=0, upper=1> prior_only;
}

transformed data {
  int next_idx = 1;

  if (sum(n_per_sample) != n_obs) {
    reject(
      "sum(n_per_sample) must equal n_obs; received ",
      sum(n_per_sample),
      " and ",
      n_obs
    );
  }

  /*
   * Confirm that x is stored contiguously by event.
   */
  for (j in 1:K) {
    if (start_idx[j] != next_idx) {
      reject(
        "start_idx is inconsistent at event ",
        j,
        "; expected ",
        next_idx,
        " but received ",
        start_idx[j]
      );
    }

    next_idx +=
      n_per_sample[j];
  }

  if (next_idx != n_obs + 1) {
    reject(
      "Event indexing does not exhaust x."
    );
  }
}

parameters {
  /*
   * Direct single-site parameters.
   *
   * These retain the same variable names and indexing as the
   * site-level parameters in the hierarchical model.
   */
  vector<
    lower=log(1e-8),
    upper=log(mu_upper)
  >[S] log_mu_site;

  vector<
    lower=log(1e-8),
    upper=log(sigma_upper)
  >[S] log_sigma_site;

  /*
   * Expected count per sampling event on the log scale.
   */
  vector[S] log_lambda;

  /*
   * Negative-binomial inverse-overdispersion parameter on
   * the log scale.
   */
  real<
    lower=-4,
    upper=12
  > log_phi;
}

transformed parameters {
  vector<lower=0>[S] mu =
    exp(log_mu_site);

  vector<lower=0>[S] sigma =
    exp(log_sigma_site);

  real<lower=0> phi =
    exp(log_phi);
}

model {
  /*
   * Direct single-site prior on the latent normal location.
   *
   * This matches the prior used for alpha_log_mu in the
   * hierarchical model.
   */
  log_mu_site[1] ~
    student_t(
      3,
      log(3),
      0.5
    );

  /*
   * Direct single-site prior on the latent normal scale.
   *
   * This matches the prior used for alpha_log_sigma.
   */
  log_sigma_site[1] ~
    normal(
      log(1),
      0.8
    );

  /*
   * Direct single-site prior on expected event abundance.
   *
   * This matches the prior used for alpha_log_lambda.
   */
  log_lambda[1] ~
    normal(
      log(500),
      2.0
    );

  /*
   * Negative-binomial dispersion prior.
   */
  log_phi ~
    normal(
      log(20),
      1.0
    );

  /*
   * Evaluate the observed-data likelihood only during
   * posterior sampling.
   */
  if (prior_only == 0) {
    for (j in 1:K) {
      int s =
        site_id[j];

      /*
       * Event count likelihood.
       *
       * E[N_j] =
       *   exp(log_lambda[s])
       *
       * Var[N_j] =
       *   exp(log_lambda[s])
       *   + exp(2 * log_lambda[s]) / phi
       */
      target +=
        neg_binomial_2_log_lpmf(
          n_per_sample[j] |
          log_lambda[s],
          phi
        );

      /*
       * Positive-truncated normal size likelihood.
       */
      if (n_per_sample[j] > 0) {
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
          normal_lccdf(
            0 |
            mu[s],
            sigma[s]
          );
      }
    }
  }
}

generated quantities {
  /*
   * Event-level log likelihoods.
   */
  vector[K] log_lik;
  vector[K] log_lik_count;
  vector[K] log_lik_size;

  /*
   * Expected count across k_ref events, retained on the
   * log scale.
   */
  vector[S] log_expected_n_ref;

  /*
   * Probability of zero total individuals across k_ref events.
   */
  vector[S] prob_zero_ref;

  /*
   * Population summaries of the positive-truncated normal.
   */
  vector[S] size_mean;
  vector[S] size_median;

  /*
   * Conditional quantiles of the maximum across k_ref events,
   * conditional on observing at least one individual.
   */
  vector[S] max_ref_q025;
  vector[S] max_ref_q50;
  vector[S] max_ref_q975;

  /*
   * One posterior- or prior-predictive realization of the
   * conditional k_ref-event maximum.
   */
  vector[S] max_ref_rep;

  /*
   * Initialize likelihood outputs so prior-only sampling
   * never evaluates the observed likelihood.
   */
  log_lik =
    rep_vector(0, K);

  log_lik_count =
    rep_vector(0, K);

  log_lik_size =
    rep_vector(0, K);

  /*
   * Calculate likelihood outputs only for posterior fits.
   */
  if (prior_only == 0) {
    for (j in 1:K) {
      int s =
        site_id[j];

      log_lik_count[j] =
        neg_binomial_2_log_lpmf(
          n_per_sample[j] |
          log_lambda[s],
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
      }

      log_lik[j] =
        log_lik_count[j] +
        log_lik_size[j];
    }
  }

  /*
   * Derived single-site quantities.
   *
   * The loop is retained so output names match the hierarchical
   * model exactly, for example max_ref_rep[1].
   */
  for (s in 1:S) {
    real a =
      -mu[s] /
      sigma[s];

    real inverse_mills =
      exp(
        std_normal_lpdf(a) -
        std_normal_lccdf(a)
      );

    real log_prob_zero_count =
      -(k_ref * phi) *
      log1p_exp(
        log_lambda[s] -
        log(phi)
      );

    /*
     * Expected total count across k_ref events:
     *
     * exp(log_expected_n_ref[s])
     */
    log_expected_n_ref[s] =
      log(k_ref) +
      log_lambda[s];

    prob_zero_ref[s] =
      exp(
        log_prob_zero_count
      );

    /*
     * Mean of the positive-truncated normal.
     */
    size_mean[s] =
      mu[s] +
      sigma[s] *
      inverse_mills;

    /*
     * Median of the positive-truncated normal.
     */
    size_median[s] =
      positive_trunc_normal_quantile(
        0.5,
        mu[s],
        sigma[s]
      );

    /*
     * Conditional maximum quantiles.
     */
    max_ref_q025[s] =
      nb_max_quantile_conditional(
        0.025,
        mu[s],
        sigma[s],
        log_lambda[s],
        phi,
        k_ref
      );

    max_ref_q50[s] =
      nb_max_quantile_conditional(
        0.50,
        mu[s],
        sigma[s],
        log_lambda[s],
        phi,
        k_ref
      );

    max_ref_q975[s] =
      nb_max_quantile_conditional(
        0.975,
        mu[s],
        sigma[s],
        log_lambda[s],
        phi,
        k_ref
      );

    max_ref_rep[s] =
      nb_max_quantile_conditional(
        uniform_rng(
          1e-12,
          1 - 1e-12
        ),
        mu[s],
        sigma[s],
        log_lambda[s],
        phi,
        k_ref
      );
  }
}
