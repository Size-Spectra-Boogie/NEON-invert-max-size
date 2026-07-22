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
  real<lower=0, upper=4.5> tau_log_lambda;
  vector[S] z_lambda;

  real<lower=-4, upper=12> log_phi;
}

transformed parameters {
  vector<lower=0>[S] mu = exp(log_mu_site);
  vector<lower=0>[S] sigma = exp(log_sigma_site);
  vector[S] log_lambda =
    alpha_log_lambda + tau_log_lambda * z_lambda;
  real<lower=0> phi = exp(log_phi);
}

model {
  /*
   * Hyperprior for the typical site-level truncated-normal
   * location parameter.
   *
   * Median near 3, strong concentration near 2–4, but with a
   * Student-t upper tail allowing rare values near or above 50.
   */
  alpha_log_mu ~ student_t( 3, log(3), 0.5 );

  /*
   * Hyperprior for the typical site-level standard deviation.
   *
   * Median sigma = 1. A log-scale SD of 0.8 remains broad,
   * allowing substantial variation among taxa.
   */
  alpha_log_sigma ~ normal( log(1), 0.8 );

  /*
   * Among-site variation in mu and sigma on the log scale.
   * Because these parameters are constrained positive, these
   * statements define half-normal priors.
   */
  tau_log_mu ~ normal( 0, 0.35 );

  tau_log_sigma ~ normal( 0, 0.50 );

  /*
   * Correctly normalized hierarchical priors for the bounded
   * site-level log(mu) and log(sigma) parameters.
   *
   * The normalizing constants must be retained because they
   * depend on the estimated hyperparameters.
   */
  for (s in 1:S) {
    target += normal_lpdf( log_mu_site[s] |
      alpha_log_mu, tau_log_mu
    );

    target += -log_diff_exp(
      normal_lcdf( log(mu_upper) |
        alpha_log_mu, tau_log_mu
      ),
      normal_lcdf( log(1e-8) |
        alpha_log_mu, tau_log_mu
      )
    );

    target += normal_lpdf( log_sigma_site[s] |
      alpha_log_sigma, tau_log_sigma
    );

    target += -log_diff_exp(
      normal_lcdf( log(sigma_upper) |
        alpha_log_sigma, tau_log_sigma
      ),
      normal_lcdf( log(1e-8) |
        alpha_log_sigma, tau_log_sigma
      )
    );
  }

  /*
   * Count hierarchy.
   *
   * alpha_log_lambda is the log expected count for a typical
   * equal-effort event. This prior is intentionally broad because
   * event counts may differ by hundreds-fold among sites.
   */
  alpha_log_lambda ~ normal( log(500), 2.0 );

  /*
   * Among-site variation in expected event counts.
   */
  tau_log_lambda ~ normal( 0, 1.5 );

  /*
   * Non-centered site effects for expected counts.
   */
  z_lambda ~ std_normal();

  /*
   * Shared negative-binomial dispersion.
   *
   * Smaller phi indicates greater event-to-event overdispersion;
   * large phi approaches a Poisson count distribution.
   */
  log_phi ~ normal( log(20), 1.0 );

  /*
   * Observed-data likelihood.
   *
   * Set prior_only = 1 to omit this section and sample from the
   * hierarchical prior. Set prior_only = 0 for posterior sampling.
   */
  if (prior_only == 0) {
    for (j in 1:K) {
      int s = site_id[j];

      /*
       * Count likelihood:
       *
       * E[n_j] = exp(log_lambda[s])
       * Var[n_j] =
       *   exp(log_lambda[s])
       *   + exp(2 * log_lambda[s]) / phi
       *
       * The log-location parameterization avoids explicitly
       * exponentiating log_lambda in the likelihood.
       */
      target += neg_binomial_2_log_lpmf(
        n_per_sample[j] |
        log_lambda[s],
        phi
      );

      /*
       * Positive-truncated normal body-size likelihood.
       * Zero-count events contribute only to the count model.
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

        /*
         * Normalize the normal density conditional on X > 0.
         */
        target += -n_per_sample[j] * normal_lccdf(
            0 | mu[s], sigma[s] );
      }
    }
  }
}

generated quantities {
  /*
   * Event-level log likelihood.
   *
   * These are calculated only for posterior sampling.
   * For prior-only sampling they are set to zero.
   */
  vector[K] log_lik;
  vector[K] log_lik_count;
  vector[K] log_lik_size;

  /*
   * Standardized expected count is retained on the log scale
   * to avoid overflow under broad prior draws.
   *
   * expected count over k_ref events:
   * exp(log_expected_n_ref[s])
   */
  vector[S] log_expected_n_ref;

  /*
   * Probability of observing zero individuals across k_ref
   * equal-effort events.
   */
  vector[S] prob_zero_ref;

  /*
   * Summaries of the positive-truncated normal size distribution.
   */
  vector[S] size_mean;
  vector[S] size_median;

  /*
   * Conditional quantiles of the maximum across k_ref events,
   * given that at least one individual is observed.
   */
  vector[S] max_ref_q025;
  vector[S] max_ref_q50;
  vector[S] max_ref_q975;

  /*
   * One posterior- or prior-predictive maximum realization
   * per site and draw, conditional on at least one individual.
   */
  vector[S] max_ref_rep;

  /*
   * Initialize likelihood outputs to zero. This ensures that
   * prior-only sampling never evaluates the observed likelihood.
   */
  log_lik = rep_vector(0, K);
  log_lik_count = rep_vector(0, K);
  log_lik_size = rep_vector(0, K);

  /*
   * Calculate observed-data log likelihood only during
   * posterior sampling.
   */
  if (prior_only == 0) {
    for (j in 1:K) {
      int s = site_id[j];

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
   * Site-level derived quantities do not depend on the observed
   * likelihood, so they are generated under both prior-only and
   * posterior sampling.
   */
  for (s in 1:S) {
    real a =
      -mu[s] / sigma[s];

    real inverse_mills =
      exp(
        std_normal_lpdf(a) -
        std_normal_lccdf(a)
      );

    /*
     * For one event:
     *
     * N ~ NegBinomial2(lambda, phi)
     *
     * The sum over k_ref events has shape k_ref * phi.
     * This expression calculates P(N_total = 0) on the log scale.
     */
    real log_prob_zero_count =
      -(k_ref * phi) *
      log1p_exp(
        log_lambda[s] -
        log(phi)
      );

    /*
     * Expected total count over k_ref events:
     *
     * log(k_ref * lambda_s)
     */
    log_expected_n_ref[s] =
      log(k_ref) +
      log_lambda[s];

    prob_zero_ref[s] =
      exp(log_prob_zero_count);

    /*
     * Population mean and median of the positive-truncated
     * normal size distribution.
     */
    size_mean[s] =
      mu[s] +
      sigma[s] * inverse_mills;

    size_median[s] =
      positive_trunc_normal_quantile(
        0.5,
        mu[s],
        sigma[s]
      );

    /*
     * Maximum-size distribution across k_ref events,
     * conditional on at least one individual being observed.
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
