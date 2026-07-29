functions {
  /*
   * Log CDF of a GEV distribution.
   *
   * G(y) = exp(-(1 + xi * z)^(-1 / xi))
   * z    = (y - mu) / sigma
   */
  real gev_logcdf(
      real y,
      real mu,
      real sigma,
      real xi
  ) {
    real z = (y - mu) / sigma;

    // Gumbel limit
    if (abs(xi) < 1e-10) {
      return -exp(-z);
    }

    {
      real support = 1 + xi * z;

      if (support <= 0) {
        if (xi > 0) {
          // Below the lower endpoint
          return negative_infinity();
        } else {
          // Above the upper endpoint
          return 0;
        }
      }

      /*
       * log(1 + xi*z) / xi approaches z as xi approaches zero.
       * The expansion improves numerical behavior near xi = 0.
       */
      if (abs(xi * z) < 1e-6) {
        real a =
          z
          - 0.5 * xi * square(z)
          + square(xi) * z * square(z) / 3;

        return -exp(-a);
      }

      return -exp(-log1p(xi * z) / xi);
    }
  }

  /*
   * Interval probability for an event maximum.
   *
   * G is the GEV CDF at reference count n_ref.
   * For count ratio r = n_event / n_ref:
   *
   * H(y) = G(y)^r
   */
  real gev_max_interval_lprob(
      real x_lower,
      real x_upper,
      real mu,
      real sigma,
      real xi,
      real count_ratio
  ) {
    real log_H_upper;
    real log_H_lower;

    if (x_upper <= x_lower) {
      return negative_infinity();
    }

    log_H_upper =
      count_ratio * gev_logcdf(x_upper, mu, sigma, xi);

    log_H_lower =
      count_ratio * gev_logcdf(x_lower, mu, sigma, xi);

    // Interval entirely below the GEV support
    if (log_H_upper == negative_infinity()) {
      return negative_infinity();
    }

    // Lower interval boundary lies below the GEV support
    if (log_H_lower == negative_infinity()) {
      return log_H_upper;
    }

    // Zero-probability interval
    if (log_H_upper <= log_H_lower) {
      return negative_infinity();
    }

    return log_diff_exp(log_H_upper, log_H_lower);
  }

  /*
   * GEV quantile using log(p), avoiding loss of precision when
   * p is close to one.
   */
  real gev_quantile_logp(
      real log_p,
      real mu,
      real sigma,
      real xi
  ) {
    real w = -log_p;
    real a = -log(w);

    // Gumbel limit
    if (abs(xi) < 1e-10) {
      return mu + sigma * a;
    }

    /*
     * Standardized quantile:
     *
     * ((-log(p))^(-xi) - 1) / xi
     *
     * Since a = -log(-log(p)), this equals expm1(xi*a)/xi.
     */
    if (abs(xi * a) < 1e-6) {
      real q_standardized =
        a
        + 0.5 * xi * square(a)
        + square(xi) * a * square(a) / 6;

      return mu + sigma * q_standardized;
    }

    return mu + sigma * expm1(xi * a) / xi;
  }
}

data {
  int<lower=1> E;                         // total sampling events
  int<lower=1> S;                         // number of sites

  array[E] int<lower=1, upper=S> site_id;

  /*
   * Recorded maximum-bin midpoint.
   * For max_bin = 2 and bin_width = 1:
   * latent maximum is in (1.5, 2.5].
   */
  vector<lower=1>[E] max_bin;

  /* data to scale each taxon
  *
  */
  real taxon_center;
  real<lower=1e-8> taxon_spread;
  
  vector[E] bin_lower;
  vector[E] bin_upper;

  // Individuals contributing to each observed event maximum
  array[E] int<lower=1> n_per_sample;

  /*
   * Global reference individuals per event.
   * Site GEV parameters describe an event with n_ref individuals.
   */
  real<lower=1> n_ref;

  // Standardized number of events, for example 20
  int<lower=1> K_ref;

  // Global parameter priors
  real loc_prior_mean;
  real<lower=1e-8> loc_prior_sd;

  real log_scale_prior_mean;
  real<lower=1e-8> log_scale_prior_sd;

  real shape_prior_mean;
  real<lower=1e-8> shape_prior_sd;

  // Priors for among-site standard deviations
  real<lower=1e-8> tau_loc_prior_sd;
  real<lower=1e-8> tau_log_scale_prior_sd;
  real<lower=1e-8> tau_shape_prior_sd;
  
  int<lower=0, upper=1> prior_only;
}

transformed data {
  
  array[S] int events_per_site;
  vector[S] mean_n_site;

  events_per_site = rep_array(0, S);
  mean_n_site = rep_vector(0, S);

  for (e in 1:E) {

    events_per_site[site_id[e]] += 1;
    mean_n_site[site_id[e]] += n_per_sample[e];
  }

  for (s in 1:S) {
    if (events_per_site[s] == 0) {
      reject("Every site must have at least one sampling event.");
    }

    mean_n_site[s] /= events_per_site[s];
  }
}

parameters {
  // Population-level GEV parameters
  real alpha_loc;
  real alpha_log_scale;
  real alpha_shape;

  // Among-site variation
  real<lower=0> tau_loc;
  real<lower=0> tau_log_scale;
  // real<lower=0> tau_shape;

  // Non-centered site effects
  vector[S] z_loc;
  vector[S] z_log_scale;
  // vector[S] z_shape;
}

transformed parameters {
  vector[S] loc_site;
  vector[S] log_scale_site;
  vector<lower=0>[S] scale_site;
  vector[S] shape_site;
  // vector[S] shape_site = rep_vector(0,S);
    // rep_vector(alpha_shape, S);

  for (s in 1:S) {
    loc_site[s] =
      alpha_loc + tau_loc * z_loc[s];

    log_scale_site[s] =
      alpha_log_scale + tau_log_scale * z_log_scale[s];

    scale_site[s] =
      exp(log_scale_site[s]);

    shape_site[s] =
      alpha_shape;// + tau_shape * z_shape[s];
  }
}

model {
  // Population-level priors
  alpha_loc ~ normal(
    loc_prior_mean,
    loc_prior_sd
  );

  alpha_log_scale ~ normal(
    log_scale_prior_mean,
    log_scale_prior_sd
  );

  // alpha_shape ~ normal(
  //   shape_prior_mean,
  //   shape_prior_sd
  // );
  alpha_shape ~ normal(0,0.1);
  
  // Half-normal priors because tau parameters are positive
  tau_loc ~ normal(
    0,
    tau_loc_prior_sd
  );

  tau_log_scale ~ normal(
    0,
    tau_log_scale_prior_sd
  );

  // tau_shape ~ normal(
  //   0,
  //   tau_shape_prior_sd
  // );

  // Non-centered hierarchy
  z_loc ~ std_normal();
  z_log_scale ~ std_normal();
  // z_shape ~ std_normal();

  // Interval-censored GEV likelihood
  if (prior_only == 0) {
  for (e in 1:E) {
    int s = site_id[e];

    real count_ratio =
      n_per_sample[e] / n_ref;

    target += gev_max_interval_lprob(
      bin_lower[e],
      bin_upper[e],
      loc_site[s],
      scale_site[s],
      shape_site[s],
      count_ratio
    );
  }
}
}

generated quantities {
  // Event-level diagnostics
  vector[E] log_lik;
  vector[E] event_max_rep;
  vector[E] event_max_rep_mm;

  /*
   * Standardized site-level outputs:
   * maximum across K_ref site-typical sampling events.
   */
  vector[S] loc_site_mm;
  vector[S] scale_site_mm;
  vector[S] shape_site_mm;

  vector[S] site_ref_q50;
  vector[S] site_ref_q95;
  vector[S] site_ref_q99;
  
  vector[S] site_ref_q50_mm;
  vector[S] site_ref_q95_mm;
  vector[S] site_ref_q99_mm;

  /*
   * One posterior-predictive standardized maximum per site.
   * Summarizing these draws gives the full predictive
   * distribution, including process and parameter uncertainty.
   */
  vector[S] site_ref_rep;
  vector[S] site_ref_rep_mm;

  // Explicit reference-effort information
  vector[S] site_mean_n;
  vector[S] site_ref_total_n;
  vector[S] site_ref_exposure_ratio;

  for (e in 1:E) {
    int s = site_id[e];

    real count_ratio =
      n_per_sample[e] / n_ref;

     if (prior_only == 0) {
    log_lik[e] = gev_max_interval_lprob(
      bin_lower[e],
      bin_upper[e],
      loc_site[s],
      scale_site[s],
      shape_site[s],
      count_ratio
    );
  } else {
    log_lik[e] = 0;
  }

    event_max_rep[e] = gev_quantile_logp(
      log(uniform_rng(1e-12, 1 - 1e-12))
        / count_ratio,
      loc_site[s],
      scale_site[s],
      shape_site[s]
    );
    
    event_max_rep_mm[e] = 
      taxon_center + taxon_spread * event_max_rep[e];
  }

  for (s in 1:S) {
    /*
     * Back transform bin loc scale, & shape:
     */
    loc_site_mm[s] =
      taxon_center + taxon_spread * loc_site[s];

    scale_site_mm[s] =
      taxon_spread * scale_site[s];

    shape_site_mm[s] =
      shape_site[s];
    /*
     * Standardized total sample size:
     *
     * K_ref events * mean individuals per observed event.
     *
     * Relative to the GEV's reference sample size:
     *
     * R_s = K_ref * mean_n_site[s] / n_ref.
     */
    site_mean_n[s] =
      mean_n_site[s];

    site_ref_total_n[s] =
      K_ref * mean_n_site[s];

    site_ref_exposure_ratio[s] =
      site_ref_total_n[s] / n_ref;

    site_ref_q50[s] = gev_quantile_logp(
      log(0.50) / site_ref_exposure_ratio[s],
      loc_site[s],
      scale_site[s],
      shape_site[s]
    );

    site_ref_q95[s] = gev_quantile_logp(
      log(0.95) / site_ref_exposure_ratio[s],
      loc_site[s],
      scale_site[s],
      shape_site[s]
    );

    site_ref_q99[s] = gev_quantile_logp(
      log(0.99) / site_ref_exposure_ratio[s],
      loc_site[s],
      scale_site[s],
      shape_site[s]
    );

    site_ref_rep[s] = gev_quantile_logp(
      log(uniform_rng(1e-12, 1 - 1e-12))
        / site_ref_exposure_ratio[s],
      loc_site[s],
      scale_site[s],
      shape_site[s]
    );
    
    site_ref_q50_mm[s] =
      taxon_center + taxon_spread * site_ref_q50[s];

    site_ref_q95_mm[s] =
      taxon_center + taxon_spread * site_ref_q95[s];

    site_ref_q99_mm[s] =
      taxon_center + taxon_spread * site_ref_q99[s];
    site_ref_rep_mm[s] = 
      taxon_center + taxon_spread * site_ref_rep[s];
  }
}
