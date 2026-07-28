functions {
  /*
   * Numerically stable log probability that a normal random variable
   * falls in the interval (lower, upper].
   */
  real normal_interval_log_prob(
      real x_lower,
      real x_upper,
      real mu,
      real sigma
  ) {
    if (x_upper <= x_lower) {
      reject(
        "normal_interval_log_prob requires upper > lower; received ",
        x_lower, " and ", x_upper
      );
    }

    if (x_upper <= mu) {
      return log_diff_exp(
        normal_lcdf(x_upper | mu, sigma),
        normal_lcdf(x_lower | mu, sigma)
      );
    } else if (x_lower >= mu) {
      return log_diff_exp(
        normal_lccdf(x_lower | mu, sigma),
        normal_lccdf(x_upper | mu, sigma)
      );
    } else {
      return log1m_exp(
        log_sum_exp(
          normal_lcdf(x_lower | mu, sigma),
          normal_lccdf(x_upper | mu, sigma)
        )
      );
    }
  }

  /*
   * Quantile function for a normal distribution truncated to
   * [lower, upper].
   */
  real bounded_trunc_normal_quantile(
      real p,
      real mu,
      real sigma,
      real x_lower,
      real x_upper
  ) {
    real p_safe =
      fmin(1 - 1e-12, fmax(1e-12, p));

    real p_normal;

    if (x_upper <= mu) {
      real log_p_normal =
        log_sum_exp(
          log1m(p_safe) +
            normal_lcdf(x_lower | mu, sigma),
          log(p_safe) +
            normal_lcdf(x_upper | mu, sigma)
        );

      p_normal =
        exp(log_p_normal);
    } else if (x_lower >= mu) {
      real log_survival =
        log_sum_exp(
          log1m(p_safe) +
            normal_lccdf(x_lower | mu, sigma),
          log(p_safe) +
            normal_lccdf(x_upper | mu, sigma)
        );

      p_normal =
        -expm1(log_survival);
    } else {
      real p_lower =
        normal_cdf(x_lower | mu, sigma);

      real p_upper =
        normal_cdf(x_upper | mu, sigma);

      p_normal =
        p_lower +
        p_safe * (p_upper - p_lower);
    }

    p_normal =
      fmin(
        1 - 1e-12,
        fmax(1e-12, p_normal)
      );

    return fmin(
      x_upper,
      fmax(
        x_lower,
        mu + sigma * inv_Phi(p_normal)
      )
    );
  }

  /*
   * Mean of a normal distribution truncated to [lower, upper].
   */
  real bounded_trunc_normal_mean(
      real mu,
      real sigma,
      real x_lower,
      real x_upper
  ) {
    real alpha =
      (x_lower - mu) / sigma;

    real beta =
      (x_upper - mu) / sigma;

    real log_z =
      normal_interval_log_prob(
        x_lower,
        x_upper,
        mu,
        sigma
      );

    real log_pdf_alpha =
      std_normal_lpdf(alpha);

    real log_pdf_beta =
      std_normal_lpdf(beta);

    real correction;

    if (log_pdf_alpha > log_pdf_beta) {
      correction =
        exp(
          log_diff_exp(
            log_pdf_alpha,
            log_pdf_beta
          ) -
          log_z
        );
    } else if (log_pdf_beta > log_pdf_alpha) {
      correction =
        -exp(
          log_diff_exp(
            log_pdf_beta,
            log_pdf_alpha
          ) -
          log_z
        );
    } else {
      correction = 0;
    }

    return fmin(
      x_upper,
      fmax(
        x_lower,
        mu + sigma * correction
      )
    );
  }

  /*
   * Quantile of the maximum across k_ref negative-binomial sampling
   * events, conditional on at least one individual being observed.
   */
  real nb_max_quantile_conditional_bounded(
      real q,
      real mu,
      real sigma,
      real log_lambda,
      real phi,
      int k_ref,
      real size_lower,
      real size_upper
  ) {
    real total_shape =
      k_ref * phi;

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

    return bounded_trunc_normal_quantile(
      individual_cdf,
      mu,
      sigma,
      size_lower,
      size_upper
    );
  }
}

data {
  /*
   * Number of sampling events for this single site.
   */
  int<lower=1> K;

  /*
   * Number of occupied event-by-size-bin cells.
   */
  int<lower=1> B_obs;

  /*
   * Total count during each event.
   */
  array[K] int<lower=0> n_per_sample;

  /*
   * Compressed binned-size data.
   *
   * bin_event[b] identifies the sampling event.
   * bin_count[b] is the number of individuals in that bin.
   * bin_lower[b] and bin_upper[b] are the bin edges.
   */
  array[B_obs] int<lower=1, upper=K> bin_event;
  array[B_obs] int<lower=1> bin_count;
  vector<lower=0>[B_obs] bin_lower;
  vector<lower=0>[B_obs] bin_upper;

  /*
   * Lower support of the size distribution.
   */
  real<lower=0> size_lower;

  /*
   * The taxon-level upper endpoint is estimated between the upper
   * edge of the largest occupied bin and this multiplier times that
   * edge. Use 1.3 to allow a maximum 30% offset.
   */
  real<lower=1> upper_multiplier_max;

  /*
   * Standardized number of events used for predictive maxima.
   */
  int<lower=1> k_ref;

  /*
   * Upper safeguards for the latent normal location and scale.
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
  array[K] int reconstructed_event_count =
    rep_array(0, K);

  real max_bin_edge_taxon =
    size_lower;

  /*
   * Multinomial coefficient for each event's vector of bin counts.
   */
  vector[K] bin_multinomial_log_coefficient;

  if (upper_multiplier_max <= 1) {
    reject(
      "upper_multiplier_max must be greater than 1; received ",
      upper_multiplier_max
    );
  }

  if (mu_upper <= size_lower) {
    reject(
      "mu_upper must exceed size_lower; received ",
      mu_upper, " and ", size_lower
    );
  }

  for (j in 1:K) {
    bin_multinomial_log_coefficient[j] =
      lgamma(n_per_sample[j] + 1.0);
  }

  for (b in 1:B_obs) {
    int j =
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

    reconstructed_event_count[j] +=
      bin_count[b];

    bin_multinomial_log_coefficient[j] -=
      lgamma(bin_count[b] + 1.0);

    if (bin_upper[b] > max_bin_edge_taxon) {
      max_bin_edge_taxon =
        bin_upper[b];
    }
  }

  for (j in 1:K) {
    if (reconstructed_event_count[j] != n_per_sample[j]) {
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
   * Direct single-site latent normal parameters.
   */
  real<
    lower=log(1e-8),
    upper=log(mu_upper)
  > log_mu;

  real<
    lower=log(1e-8),
    upper=log(sigma_upper)
  > log_sigma;

  /*
   * Position of the taxon-level endpoint between the largest
   * occupied bin edge and its upper_multiplier_max multiple.
   */
  real<lower=0, upper=1> upper_fraction;

  /*
   * Expected count per event on the log scale.
   */
  real log_lambda;

  /*
   * Negative-binomial inverse-overdispersion on the log scale.
   */
  real<lower=-4, upper=12> log_phi;
}

transformed parameters {
  real<lower=0> mu =
    exp(log_mu);

  real<lower=0> sigma =
    exp(log_sigma);

  real<lower=0> phi =
    exp(log_phi);

  real upper_multiplier =
    1 +
    (
      upper_multiplier_max - 1
    ) *
    upper_fraction;

  real size_upper_taxon =
    max_bin_edge_taxon *
    upper_multiplier;
}

model {
  /*
   * Direct priors corresponding to the centers used in the
   * hierarchical model.
   */
  log_mu ~
    student_t(
      3,
      log(3),
      0.5
    );

  log_sigma ~
    normal(
      log(1),
      0.8
    );

  upper_fraction ~
    beta(2, 2);

  log_lambda ~
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
    /*
     * Negative-binomial event totals.
     */
    for (j in 1:K) {
      target +=
        neg_binomial_2_log_lpmf(
          n_per_sample[j] |
          log_lambda,
          phi
        );

      target +=
        bin_multinomial_log_coefficient[j];
    }

    /*
     * Exact grouped-bin likelihood under a doubly truncated normal.
     */
    for (b in 1:B_obs) {
      target +=
        bin_count[b] *
        (
          normal_interval_log_prob(
            bin_lower[b],
            bin_upper[b],
            mu,
            sigma
          ) -
          normal_interval_log_prob(
            size_lower,
            size_upper_taxon,
            mu,
            sigma
          )
        );
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
   * Standardized count summaries.
   */
  real log_expected_n_ref;
  real prob_zero_ref;

  /*
   * Summaries and quantiles of the doubly truncated size
   * distribution.
   */
  real size_mean;
  real size_median;
  real size_q05;
  real size_q25;
  real size_q50;
  real size_q75;
  real size_q95;
  real size_q975;
  real size_q99;

  /*
   * Conditional quantiles and one predictive realization of the
   * maximum across k_ref events.
   */
  real max_ref_q025;
  real max_ref_q50;
  real max_ref_q975;
  real max_ref_rep;

  /*
   * Endpoint diagnostic.
   */
  real observed_max_bin_edge_taxon =
    max_bin_edge_taxon;

  real log_prob_zero_count =
    -(k_ref * phi) *
    log1p_exp(
      log_lambda -
      log(phi)
    );

  log_lik =
    rep_vector(0, K);

  log_lik_count =
    rep_vector(0, K);

  log_lik_size =
    rep_vector(0, K);

  if (prior_only == 0) {
    for (j in 1:K) {
      log_lik_count[j] =
        neg_binomial_2_log_lpmf(
          n_per_sample[j] |
          log_lambda,
          phi
        );

      log_lik_size[j] =
        bin_multinomial_log_coefficient[j];
    }

    for (b in 1:B_obs) {
      int j =
        bin_event[b];

      log_lik_size[j] +=
        bin_count[b] *
        (
          normal_interval_log_prob(
            bin_lower[b],
            bin_upper[b],
            mu,
            sigma
          ) -
          normal_interval_log_prob(
            size_lower,
            size_upper_taxon,
            mu,
            sigma
          )
        );
    }

    log_lik =
      log_lik_count +
      log_lik_size;
  }

  log_expected_n_ref =
    log(k_ref) +
    log_lambda;

  prob_zero_ref =
    exp(log_prob_zero_count);

  size_mean =
    bounded_trunc_normal_mean(
      mu,
      sigma,
      size_lower,
      size_upper_taxon
    );

  size_median =
    bounded_trunc_normal_quantile(
      0.50,
      mu,
      sigma,
      size_lower,
      size_upper_taxon
    );

  size_q05 =
    bounded_trunc_normal_quantile(
      0.05,
      mu,
      sigma,
      size_lower,
      size_upper_taxon
    );

  size_q25 =
    bounded_trunc_normal_quantile(
      0.25,
      mu,
      sigma,
      size_lower,
      size_upper_taxon
    );

  size_q50 =
    bounded_trunc_normal_quantile(
      0.50,
      mu,
      sigma,
      size_lower,
      size_upper_taxon
    );

  size_q75 =
    bounded_trunc_normal_quantile(
      0.75,
      mu,
      sigma,
      size_lower,
      size_upper_taxon
    );

  size_q95 =
    bounded_trunc_normal_quantile(
      0.95,
      mu,
      sigma,
      size_lower,
      size_upper_taxon
    );

  size_q975 =
    bounded_trunc_normal_quantile(
      0.975,
      mu,
      sigma,
      size_lower,
      size_upper_taxon
    );

  size_q99 =
    bounded_trunc_normal_quantile(
      0.99,
      mu,
      sigma,
      size_lower,
      size_upper_taxon
    );

  max_ref_q025 =
    nb_max_quantile_conditional_bounded(
      0.025,
      mu,
      sigma,
      log_lambda,
      phi,
      k_ref,
      size_lower,
      size_upper_taxon
    );

  max_ref_q50 =
    nb_max_quantile_conditional_bounded(
      0.50,
      mu,
      sigma,
      log_lambda,
      phi,
      k_ref,
      size_lower,
      size_upper_taxon
    );

  max_ref_q975 =
    nb_max_quantile_conditional_bounded(
      0.975,
      mu,
      sigma,
      log_lambda,
      phi,
      k_ref,
      size_lower,
      size_upper_taxon
    );

  max_ref_rep =
    nb_max_quantile_conditional_bounded(
      uniform_rng(
        1e-12,
        1 - 1e-12
      ),
      mu,
      sigma,
      log_lambda,
      phi,
      k_ref,
      size_lower,
      size_upper_taxon
    );
}
