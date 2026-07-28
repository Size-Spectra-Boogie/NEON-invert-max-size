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

    /*
     * Select the CDF or CCDF representation that is most stable in
     * the relevant tail. The formulas are mathematically equivalent.
     */
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
   *
   * p is a probability on the truncated-normal scale.
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

    /*
     * Interpolate on the CDF scale. Tail-specific calculations avoid
     * subtracting nearly equal probabilities when the whole interval
     * lies far from the latent normal location.
     */
    if (x_upper <= mu) {
      real log_p_normal =
        log_sum_exp(
          log1m(p_safe) +
            normal_lcdf(x_lower | mu, sigma),
          log(p_safe) +
            normal_lcdf(x_upper | mu, sigma)
        );

      p_normal = exp(log_p_normal);
    } else if (x_lower >= mu) {
      real log_survival =
        log_sum_exp(
          log1m(p_safe) +
            normal_lccdf(x_lower | mu, sigma),
          log(p_safe) +
            normal_lccdf(x_upper | mu, sigma)
        );

      p_normal = -expm1(log_survival);
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

    /*
     * The final clamp protects against tiny floating-point excursions
     * outside the modeled support.
     */
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
   *
   * For each event:
   *   N_r ~ NegBinomial2(lambda, phi)
   *
   * The size distribution is a normal truncated to
   * [size_lower, size_upper].
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

    /*
     * Probability of zero individuals across all k_ref events.
     */
    real log_p_zero =
      -total_shape *
      log1p_exp(
        log_lambda -
        log(phi)
      );

    /*
     * Convert a conditional maximum CDF q into the corresponding
     * unconditional CDF value:
     *
     * P(M <= m) =
     *   P(N = 0) +
     *   q * [1 - P(N = 0)].
     */
    real log_target_cdf =
      log_sum_exp(
        log_p_zero,
        log(q) +
        log1m_exp(log_p_zero)
      );

    /*
     * Invert the negative-binomial probability-generating function
     * to obtain the required individual-size CDF value.
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
  int<lower=1> S;       // Number of sites
  int<lower=1> K;       // Number of sampling events

  /*
   * Number of occupied event-by-size-bin cells. Each row represents
   * one positive bin count, rather than one individual.
   */
  int<lower=1> B_obs;

  /*
   * Event-level count data.
   */
  array[K] int<lower=0> n_per_sample;
  array[K] int<lower=1, upper=S> site_id;

  /*
   * Compressed binned-size data.
   *
   * bin_event[b] identifies the sampling event.
   * bin_count[b] is the number of individuals in that bin.
   * bin_lower[b] and bin_upper[b] are the bin edges in millimeters.
   */
  array[B_obs] int<lower=1, upper=K> bin_event;
  array[B_obs] int<lower=1> bin_count;
  vector<lower=0>[B_obs] bin_lower;
  vector<lower=0>[B_obs] bin_upper;

  /*
   * Lower support of the size distribution. Use 0.5 for 1-mm bins
   * whose first biologically valid edge is 0.5 mm.
   */
  real<lower=0> size_lower;

  /*
   * The shared taxon-level endpoint lies between the upper edge of
   * the largest occupied bin across all sites and this multiplier
   * times that edge.
   *
   * Use 1.3 to allow at most a 30% increase above the largest
   * observed taxon-level bin edge.
   */
  real<lower=1> upper_multiplier_max;
  
  /*
   * Standardized number of sampling events used for maximum-size
   * prediction.
   */
  int<lower=1> k_ref;

  /*
   * Numerical/scientific safeguards for the latent normal location
   * and scale. Set mu_upper = 400 and sigma_upper = 200.
   */
  real<lower=1e-6> mu_upper;
  real<lower=1e-6> sigma_upper;

  /*
   * prior_only = 0: posterior sampling
   * prior_only = 1: sample only from the parameter priors
   */
  int<lower=0, upper=1> prior_only;
}

transformed data {
  array[K] int reconstructed_event_count =
    rep_array(0, K);

  array[S] int occupied_bin_rows_by_site =
    rep_array(0, S);

  vector[S] max_bin_edge_site =
    rep_vector(size_lower, S);

  real max_bin_edge_taxon = 
    size_lower;
    
  real size_upper_taxon;

   /*
   * Multinomial coefficient for the vector of bin counts within each
   * event. It is constant with respect to model parameters, but is
   * retained so event-level log_lik is a normalized joint likelihood.
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
      "mu_upper must be greater than size_lower; received ",
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

    int s =
      site_id[j];

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

    occupied_bin_rows_by_site[s] +=
      1;

    bin_multinomial_log_coefficient[j] -=
      lgamma(bin_count[b] + 1.0);

    if (bin_upper[b] > max_bin_edge_site[s]) {
      max_bin_edge_site[s] =
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

  /*
   * A site with no measured individuals has no observed maximum-bin
   * edge from which to define its endpoint interval.
   */
 for (s in 1:S) {
    if (occupied_bin_rows_by_site[s] == 0) {
      reject(
        "Site ",
        s,
        " has no occupied size bins."
      );
    }
  }
  
  max_bin_edge_taxon =
    max(max_bin_edge_site);
    
  /*
   * Fixed taxon-level support cap.
   */
  size_upper_taxon =
    upper_multiplier_max *
    max_bin_edge_taxon;
}

parameters {
  vector<
    lower=log(1e-8),
    upper=log(mu_upper)
  >[S] log_mu_site;

  vector<
    lower=log(1e-8),
    upper=log(sigma_upper)
  >[S] log_sigma_site;

  // Retain count hierarchy
  real alpha_log_lambda;
  real<lower=0, upper=4.5> tau_log_lambda;
  vector[S] z_lambda;

  real<lower=-4, upper=12> log_phi;
}

transformed parameters {
  vector<lower=0>[S] mu =
    exp(log_mu_site);

  vector<lower=0>[S] sigma =
    exp(log_sigma_site);

  vector[S] log_lambda =
    alpha_log_lambda +
    tau_log_lambda * z_lambda;

  real<lower=0> phi =
    exp(log_phi);

  /*
 * Taxon-level upper-support multiplier.
 *
 * upper_fraction = 0:
 *   size_upper_taxon = max_bin_edge_taxon
 *
 * upper_fraction = 1:
 *   size_upper_taxon =
 *     upper_multiplier_max * max_bin_edge_taxon
 */
// real upper_multiplier =
//   1 +
//   (
//     upper_multiplier_max - 1
//   ) *
//   upper_fraction;

/*
 * Shared taxon-level upper support for all sites.
 */
// real size_upper_taxon =
//   max_bin_edge_taxon *
//   upper_multiplier;
}

model {
  log_mu_site ~
    student_t(
      3,
      log(3),
      0.5
    );

  log_sigma_site ~
    normal(
      log(1),
      0.8
    );

  /*
   * Correctly normalized bounded hierarchical priors for site-level
   * log(mu) and log(sigma). The normalizing terms depend on the
   * estimated hyperparameters and therefore cannot be dropped.
   */
  // for (s in 1:S) {
  //   target +=
  //     normal_lpdf(
  //       log_mu_site[s] |
  //       alpha_log_mu,
  //       tau_log_mu
  //     ) -
  //     log_diff_exp(
  //       normal_lcdf(
  //         log(mu_upper) |
  //         alpha_log_mu,
  //         tau_log_mu
  //       ),
  //       normal_lcdf(
  //         log(1e-8) |
  //         alpha_log_mu,
  //         tau_log_mu
  //       )
  //     );
  // 
  //   target +=
  //     normal_lpdf(
  //       log_sigma_site[s] |
  //       alpha_log_sigma,
  //       tau_log_sigma
  //     ) -
  //     log_diff_exp(
  //       normal_lcdf(
  //         log(sigma_upper) |
  //         alpha_log_sigma,
  //         tau_log_sigma
  //       ),
  //       normal_lcdf(
  //         log(1e-8) |
  //         alpha_log_sigma,
  //         tau_log_sigma
  //       )
  //     );
  // }

  /*
   * Endpoint prior. For upper_multiplier_max = 1.3:
   *
   * upper_fraction = 0   -> endpoint at the largest occupied bin edge
   * upper_fraction = 1   -> endpoint at 1.3 times that edge
   *
   * beta(2, 2) has zero density at both boundaries and is centered at
   * the midpoint of the allowable interval.
   */
  // upper_fraction ~
  //   beta(2, 2);

  /*
   * Count hierarchy.
   */
  alpha_log_lambda ~
    normal(
      log(500),
      2.0
    );

  tau_log_lambda ~
    normal(
      0,
      1.5
    );

  z_lambda ~
    std_normal();

  log_phi ~
    normal(
      log(20),
      1.0
    );

  /*
   * Observed-data likelihood.
   */
  if (prior_only == 0) {
  /*
   * Negative-binomial event totals.
   */
  for (j in 1:K) {
    target +=
      neg_binomial_2_log_lpmf(
        n_per_sample[j] |
        log_lambda[site_id[j]],
        phi
      );

    target +=
      bin_multinomial_log_coefficient[j];
  }

  /*
   * Conditional binned-size distribution.
   */
  for (b in 1:B_obs) {
    int j =
      bin_event[b];

    int s =
      site_id[j];

    target +=
      bin_count[b] *
      (
        normal_interval_log_prob(
          bin_lower[b],
          bin_upper[b],
          mu[s],
          sigma[s]
        ) -
        normal_interval_log_prob(
          size_lower,
          size_upper_taxon,
          mu[s],
          sigma[s]
        )
      );
  }
}
}

generated quantities {
  /*
   * Event-level joint log likelihoods.
   */
  vector[K] log_lik;
  vector[K] log_lik_count;
  vector[K] log_lik_size;

  /*
   * Standardized expected counts and zero-count probability across
   * k_ref events.
   */
  vector[S] log_expected_n_ref;
  vector[S] prob_zero_ref;

  /*
   * Summaries of the doubly truncated size distribution.
   */
  vector[S] size_mean;
  vector[S] size_q05;
  vector[S] size_q25;
  vector[S] size_q50;
  vector[S] size_q75;
  vector[S] size_q95;
  vector[S] size_q975;
  vector[S] size_q99;

  /*
   * Conditional quantiles and one predictive realization of the
   * maximum across k_ref events, given at least one individual.
   */
  vector[S] max_ref_q025;
  vector[S] max_ref_q50;
  vector[S] max_ref_q975;
  vector[S] max_ref_rep;

  /*
   * Store the observed maximum-bin edge in the output for convenient
   * endpoint diagnostics.
   */
  vector[S] observed_max_bin_edge =
    max_bin_edge_site;

  log_lik =
    rep_vector(0, K);

  log_lik_count =
    rep_vector(0, K);

  log_lik_size =
    rep_vector(0, K);

  /*
   * Avoid evaluating the observed likelihood in prior-only runs.
   */
  if (prior_only == 0) {
    for (j in 1:K) {
      log_lik_count[j] =
        neg_binomial_2_log_lpmf(
          n_per_sample[j] |
          log_lambda[site_id[j]],
          phi
        );

      log_lik_size[j] =
        bin_multinomial_log_coefficient[j];
    }

    for (b in 1:B_obs) {
      int j =
        bin_event[b];

      int s =
        site_id[j];

      log_lik_size[j] +=
        bin_count[b] *
        (
          normal_interval_log_prob(
            bin_lower[b],
            bin_upper[b],
            mu[s],
            sigma[s]
          ) -
          normal_interval_log_prob(
            size_lower,
            size_upper_taxon,
            mu[s],
            sigma[s]
          )
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
        log_lambda[s] -
        log(phi)
      );

    log_expected_n_ref[s] =
      log(k_ref) +
      log_lambda[s];

    prob_zero_ref[s] =
      exp(log_prob_zero_count);

    size_mean[s] =
      bounded_trunc_normal_mean(
        mu[s],
        sigma[s],
        size_lower,
        size_upper_taxon
      );

    max_ref_q025[s] =
      nb_max_quantile_conditional_bounded(
        0.025,
        mu[s],
        sigma[s],
        log_lambda[s],
        phi,
        k_ref,
        size_lower,
        size_upper_taxon
      );

    max_ref_q50[s] =
      nb_max_quantile_conditional_bounded(
        0.50,
        mu[s],
        sigma[s],
        log_lambda[s],
        phi,
        k_ref,
        size_lower,
        size_upper_taxon
      );

    max_ref_q975[s] =
      nb_max_quantile_conditional_bounded(
        0.975,
        mu[s],
        sigma[s],
        log_lambda[s],
        phi,
        k_ref,
        size_lower,
        size_upper_taxon
      );

    max_ref_rep[s] =
      nb_max_quantile_conditional_bounded(
        uniform_rng(
          1e-12,
          1 - 1e-12
        ),
        mu[s],
        sigma[s],
        log_lambda[s],
        phi,
        k_ref,
        size_lower,
        size_upper_taxon
      );
      
      size_q05[s] =
      bounded_trunc_normal_quantile(
        0.05,
        mu[s],
        sigma[s],
        size_lower,
        size_upper_taxon
      );

    size_q25[s] =
      bounded_trunc_normal_quantile(
        0.25,
        mu[s],
        sigma[s],
        size_lower,
        size_upper_taxon
      );

    size_q50[s] =
      bounded_trunc_normal_quantile(
        0.50,
        mu[s],
        sigma[s],
        size_lower,
        size_upper_taxon
      );

    size_q75[s] =
      bounded_trunc_normal_quantile(
        0.75,
        mu[s],
        sigma[s],
        size_lower,
        size_upper_taxon
      );

    size_q95[s] =
      bounded_trunc_normal_quantile(
        0.95,
        mu[s],
        sigma[s],
        size_lower,
        size_upper_taxon
      );

    size_q975[s] =
      bounded_trunc_normal_quantile(
        0.975,
        mu[s],
        sigma[s],
        size_lower,
        size_upper_taxon
      );

    size_q99[s] =
      bounded_trunc_normal_quantile(
        0.99,
        mu[s],
        sigma[s],
        size_lower,
        size_upper_taxon
      );
  }
}
