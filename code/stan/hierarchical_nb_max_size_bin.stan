functions {
  /*
   * Log probability that a normal random variable falls in
   * (x_lower, x_upper].
   *
   * The complement representation is branch-free with respect to
   * model parameters:
   *
   * P(x_lower < X <= x_upper)
   *   = 1 - P(X <= x_lower) - P(X > x_upper).
   */
  real normal_interval_log_prob(
      real x_lower,
      real x_upper,
      real mu,
      real sigma
  ) {
    real log_outside;

    if (x_upper <= x_lower) {
      reject(
        "normal_interval_log_prob requires x_upper > x_lower; received ",
        x_lower, " and ", x_upper
      );
    }

    log_outside =
      log_sum_exp(
        normal_lcdf(
          x_lower |
          mu,
          sigma
        ),
        normal_lccdf(
          x_upper |
          mu,
          sigma
        )
      );

    /*
     * Floating-point rounding can rarely return exactly zero for
     * log_outside. Keep the argument inside the domain of log1m_exp.
     */
    return log1m_exp(
      fmin(
        -1e-15,
        log_outside
      )
    );
  }

  /*
   * Quantile of a normal distribution truncated below at size_lower.
   */
  real lower_trunc_normal_quantile(
      real p,
      real mu,
      real sigma,
      real size_lower
  ) {
    real p_safe =
      fmin(
        1 - 1e-12,
        fmax(
          1e-12,
          p
        )
      );

    real p_at_lower =
      normal_cdf(
        size_lower |
        mu,
        sigma
      );

    real p_normal =
      p_at_lower +
      p_safe *
      (
        1 -
        p_at_lower
      );

    p_normal =
      fmin(
        1 - 1e-12,
        fmax(
          1e-12,
          p_normal
        )
      );

    return
      mu +
      sigma *
      inv_Phi(p_normal);
  }

  /*
   * Mean of a normal distribution truncated below at size_lower.
   */
  real lower_trunc_normal_mean(
      real mu,
      real sigma,
      real size_lower
  ) {
    real alpha =
      (
        size_lower -
        mu
      ) /
      sigma;

    real inverse_mills =
      exp(
        std_normal_lpdf(alpha) -
        std_normal_lccdf(alpha)
      );

    return
      mu +
      sigma *
      inverse_mills;
  }

  /*
   * Quantile of the maximum across k_ref negative-binomial sampling
   * events, conditional on at least one individual being observed.
   *
   * Each event has:
   *
   *   N_r ~ NegBinomial2(lambda, phi)
   *
   * and individual sizes follow the lower-truncated normal.
   */
  real nb_max_quantile_conditional(
      real q,
      real mu,
      real sigma,
      real log_lambda,
      real phi,
      int k_ref,
      real size_lower
  ) {
    real total_shape =
      k_ref *
      phi;

    real log_p_zero =
      -total_shape *
      log1p_exp(
        log_lambda -
        log(phi)
      );

    real log_target_cdf =
      log_sum_exp(
        log_p_zero,
        log(q) +
        log1m_exp(log_p_zero)
      );

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

    return lower_trunc_normal_quantile(
      individual_cdf,
      mu,
      sigma,
      size_lower
    );
  }
}

data {
  int<lower=1> S;
  int<lower=1> K;
  int<lower=1> B_obs;

  /*
   * Event totals and site membership.
   */
  array[K] int<lower=0> n_per_sample;
  array[K] int<lower=1, upper=S> site_id;

  /*
   * Compressed positive event-by-bin cells.
   */
  array[B_obs] int<lower=1, upper=K> bin_event;
  array[B_obs] int<lower=1> bin_count;
  vector<lower=0>[B_obs] bin_lower;
  vector<lower=0>[B_obs] bin_upper;

  /*
   * Fixed lower truncation of the body-size distribution.
   */
  real<lower=0> size_lower;

  /*
   * Numerical lower floor for the fitted size scale. A value of
   * 0.01 is suitable for 1-mm bins and prevents floating-point
   * underflow to an invalid zero normal scale.
   */
  real<lower=1e-8> sigma_floor;

  /*
   * Number of equal-effort sampling events used for the
   * effort-dependent reference-maximum summaries.
   */
  int<lower=1> k_ref;

  /*
   * prior_only = 0: posterior sampling
   * prior_only = 1: prior sampling conditional on the data structure
   */
  int<lower=0, upper=1> prior_only;
}

transformed data {
  array[K] int reconstructed_event_count =
    rep_array(0, K);

  vector[K] bin_multinomial_log_coefficient;

  for (j in 1:K) {
    bin_multinomial_log_coefficient[j] =
      lgamma(
        n_per_sample[j] +
        1.0
      );
  }

  for (b in 1:B_obs) {
    int event =
      bin_event[b];

    if (bin_upper[b] <= bin_lower[b]) {
      reject(
        "bin_upper must exceed bin_lower at bin row ",
        b,
        "; received ",
        bin_lower[b], " and ", bin_upper[b]
      );
    }

    if (bin_lower[b] < size_lower) {
      reject(
        "Occupied bin row ",
        b,
        " extends below size_lower. bin_lower = ",
        bin_lower[b],
        "; size_lower = ",
        size_lower
      );
    }

    reconstructed_event_count[event] +=
      bin_count[b];

    bin_multinomial_log_coefficient[event] -=
      lgamma(
        bin_count[b] +
        1.0
      );
  }

  for (j in 1:K) {
    if (
      reconstructed_event_count[j] !=
      n_per_sample[j]
    ) {
      reject(
        "Binned counts do not equal n_per_sample for event ",
        j,
        ". Reconstructed = ",
        reconstructed_event_count[j],
        "; n_per_sample = ",
        n_per_sample[j]
      );
    }
  }
}

parameters {
  /*
   * Noncentered partial pooling for latent normal location.
   */
  real alpha_log_mu;
  real<lower=1e-4, upper=2.5> tau_log_mu;
  vector[S] z_mu;

  /*
   * Noncentered partial pooling for the log excess of the latent
   * normal scale above sigma_floor.
   */
  real alpha_log_sigma;
  real<lower=1e-4, upper=2.5> tau_log_sigma;
  vector[S] z_sigma;

  /*
   * Direct site-level expected event counts. These are deliberately
   * not pooled across sites.
   */
  vector[S] log_lambda_site;

  /*
   * Shared negative-binomial inverse-overdispersion.
   */
  real<lower=-4, upper=12> log_phi;
}

transformed parameters {
  vector[S] log_mu_site =
    alpha_log_mu +
    tau_log_mu *
    z_mu;

  vector[S] log_sigma_site =
    alpha_log_sigma +
    tau_log_sigma *
    z_sigma;

  /*
   * Tiny additive floors prevent exp() underflow from producing
   * invalid exact-zero values during warmup. For sigma,
   * log_sigma_site is the log excess above sigma_floor.
   */
  vector<lower=0>[S] mu =
    rep_vector(1e-8, S) +
    exp(log_mu_site);

  vector<lower=0>[S] sigma =
    rep_vector(sigma_floor, S) +
    exp(log_sigma_site);

  real<lower=0> phi =
    exp(log_phi);
}

model {
  /*
   * Size-distribution hierarchy.
   */
  alpha_log_mu ~
    student_t(
      3,
      log(3),
      0.5
    );

  tau_log_mu ~
    normal(
      0,
      0.35
    );

  z_mu ~
    std_normal();

  alpha_log_sigma ~
    normal(
      log(1),
      0.8
    );

  tau_log_sigma ~
    normal(
      0,
      0.50
    );

  z_sigma ~
    std_normal();

  /*
   * Independent site-level count means and shared dispersion.
   */
  log_lambda_site ~
    normal(
      log(500),
      2.0
    );

  log_phi ~
    normal(
      log(20),
      1.0
    );

  if (prior_only == 0) {
    vector[S] log_size_normalizer;

    for (s in 1:S) {
      log_size_normalizer[s] =
        normal_lccdf(
          size_lower |
          mu[s],
          sigma[s]
        );
    }

    /*
     * Event totals.
     */
    for (j in 1:K) {
      target +=
        neg_binomial_2_log_lpmf(
          n_per_sample[j] |
          log_lambda_site[
            site_id[j]
          ],
          phi
        );

      /*
       * Multinomial coefficient for the conditional bin-count vector.
       * Unoccupied bins have count zero and contribute no additional
       * term.
       */
      target +=
        bin_multinomial_log_coefficient[j];
    }

    /*
     * Exact grouped-bin size likelihood conditional on X > size_lower.
     */
    for (b in 1:B_obs) {
      int event =
        bin_event[b];

      int site =
        site_id[event];

      target +=
        bin_count[b] *
        (
          normal_interval_log_prob(
            bin_lower[b],
            bin_upper[b],
            mu[site],
            sigma[site]
          ) -
          log_size_normalizer[site]
        );
    }
  }
}

generated quantities {
  /*
   * Event-level joint likelihood components.
   */
  vector[K] log_lik;
  vector[K] log_lik_count;
  vector[K] log_lik_size;

  /*
   * Site-level event-count summaries.
   */
  vector[S] lambda_site =
    exp(log_lambda_site);

  vector[S] log_expected_n_ref;
  vector[S] prob_zero_ref;

  /*
   * Effort-independent body-size summaries.
   */
  vector[S] size_mean;
  vector[S] size_q50;
  vector[S] size_q75;
  vector[S] size_q95;
  vector[S] size_q99;

  /*
   * Effort-dependent quantiles of the maximum observed across k_ref
   * equal-effort sampling events, conditional on at least one
   * individual being observed.
   */
  vector[S] max_ref_q50;
  vector[S] max_ref_q75;
  vector[S] max_ref_q95;
  vector[S] max_ref_q99;

  /*
   * One conditional posterior- or prior-predictive maximum draw per
   * site and posterior draw.
   */
  vector[S] max_ref_rep;

  log_lik =
    rep_vector(0, K);

  log_lik_count =
    rep_vector(0, K);

  log_lik_size =
    rep_vector(0, K);

  if (prior_only == 0) {
    vector[S] log_size_normalizer;

    for (s in 1:S) {
      log_size_normalizer[s] =
        normal_lccdf(
          size_lower |
          mu[s],
          sigma[s]
        );
    }

    for (j in 1:K) {
      log_lik_count[j] =
        neg_binomial_2_log_lpmf(
          n_per_sample[j] |
          log_lambda_site[
            site_id[j]
          ],
          phi
        );

      log_lik_size[j] =
        bin_multinomial_log_coefficient[j];
    }

    for (b in 1:B_obs) {
      int event =
        bin_event[b];

      int site =
        site_id[event];

      log_lik_size[event] +=
        bin_count[b] *
        (
          normal_interval_log_prob(
            bin_lower[b],
            bin_upper[b],
            mu[site],
            sigma[site]
          ) -
          log_size_normalizer[site]
        );
    }

    log_lik =
      log_lik_count +
      log_lik_size;
  }

  for (s in 1:S) {
    real log_prob_zero_count =
      -(k_ref * phi) *
      log1p_exp(
        log_lambda_site[s] -
        log(phi)
      );

    log_expected_n_ref[s] =
      log(k_ref) +
      log_lambda_site[s];

    prob_zero_ref[s] =
      exp(log_prob_zero_count);

    size_mean[s] =
      lower_trunc_normal_mean(
        mu[s],
        sigma[s],
        size_lower
      );

    size_q50[s] =
      lower_trunc_normal_quantile(
        0.50,
        mu[s],
        sigma[s],
        size_lower
      );

    size_q75[s] =
      lower_trunc_normal_quantile(
        0.75,
        mu[s],
        sigma[s],
        size_lower
      );

    size_q95[s] =
      lower_trunc_normal_quantile(
        0.95,
        mu[s],
        sigma[s],
        size_lower
      );

    size_q99[s] =
      lower_trunc_normal_quantile(
        0.99,
        mu[s],
        sigma[s],
        size_lower
      );

    max_ref_q50[s] =
      nb_max_quantile_conditional(
        0.50,
        mu[s],
        sigma[s],
        log_lambda_site[s],
        phi,
        k_ref,
        size_lower
      );

    max_ref_q75[s] =
      nb_max_quantile_conditional(
        0.75,
        mu[s],
        sigma[s],
        log_lambda_site[s],
        phi,
        k_ref,
        size_lower
      );

    max_ref_q95[s] =
      nb_max_quantile_conditional(
        0.95,
        mu[s],
        sigma[s],
        log_lambda_site[s],
        phi,
        k_ref,
        size_lower
      );

    max_ref_q99[s] =
      nb_max_quantile_conditional(
        0.99,
        mu[s],
        sigma[s],
        log_lambda_site[s],
        phi,
        k_ref,
        size_lower
      );

    max_ref_rep[s] =
      nb_max_quantile_conditional(
        uniform_rng(
          1e-12,
          1 - 1e-12
        ),
        mu[s],
        sigma[s],
        log_lambda_site[s],
        phi,
        k_ref,
        size_lower
      );
  }
}
