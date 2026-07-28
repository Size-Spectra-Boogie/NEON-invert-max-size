functions {
  /*
   * Full Dirichlet-multinomial log probability.
   */
  real dirichlet_multinomial_logpmf(
      array[] int y,
      vector alpha
  ) {
    int B =
      size(y);

    int N =
      0;

    real alpha_total =
      sum(alpha);

    real lp =
      lgamma(alpha_total);

    for (b in 1:B) {
      N +=
        y[b];
    }

    lp +=
      lgamma(N + 1.0) -
      lgamma(N + alpha_total);

    for (b in 1:B) {
      lp +=
        lgamma(
          y[b] +
          alpha[b]
        ) -
        lgamma(
          alpha[b]
        ) -
        lgamma(
          y[b] +
          1.0
        );
    }

    return lp;
  }

  /*
   * Aggregate retained latent-grid mass into observed bins and
   * normalize it to an observed probability vector.
   *
   * log_true_mass includes the latent cell width.
   */
  vector observed_bin_probability(
      vector log_true_mass,
      vector log_retention,
      array[] int grid_bin,
      int B
  ) {
    int G =
      num_elements(
        log_true_mass
      );

    vector[B] log_retained_bin_mass =
      rep_vector(
        negative_infinity(),
        B
      );

    for (g in 1:G) {
      int b =
        grid_bin[g];

      log_retained_bin_mass[b] =
        log_sum_exp(
          log_retained_bin_mass[b],
          log_true_mass[g] +
          log_retention[g]
        );
    }

    return
      softmax(
        log_retained_bin_mass
      );
  }

  /*
   * Aggregate latent fine-grid probabilities into observed-bin
   * intervals without applying retention.
   */
  vector true_bin_probability(
      vector p_true_grid,
      array[] int grid_bin,
      int B
  ) {
    int G =
      num_elements(
        p_true_grid
      );

    vector[B] p_true_bin =
      rep_vector(
        0,
        B
      );

    for (g in 1:G) {
      p_true_bin[
        grid_bin[g]
      ] +=
        p_true_grid[g];
    }

    return p_true_bin;
  }

  /*
   * Quantile from a piecewise-uniform latent fine grid.
   */
  real latent_grid_quantile(
      vector probability,
      vector grid_lower,
      vector grid_upper,
      real q
  ) {
    int G =
      num_elements(
        probability
      );

    real q_safe =
      fmin(
        1 - 1e-12,
        fmax(
          1e-12,
          q
        )
      );

    real cumulative =
      0;

    for (g in 1:G) {
      real next_cumulative =
        cumulative +
        probability[g];

      if (
        q_safe <=
        next_cumulative
      ) {
        real fraction;

        if (
          probability[g] <=
          1e-15
        ) {
          return
            0.5 *
            (
              grid_lower[g] +
              grid_upper[g]
            );
        }

        fraction =
          (
            q_safe -
            cumulative
          ) /
          probability[g];

        fraction =
          fmin(
            1,
            fmax(
              0,
              fraction
            )
          );

        return
          grid_lower[g] +
          fraction *
          (
            grid_upper[g] -
            grid_lower[g]
          );
      }

      cumulative =
        next_cumulative;
    }

    return
      grid_upper[G];
  }
}

data {
  int<lower=1> S;
  int<lower=1> K;
  int<lower=2> B;
  int<lower=B> G;
  int<lower=1> J;

  array[K] int<lower=1, upper=S> site_id;
  array[K, B] int<lower=0> y;

  vector<lower=0>[G] grid_lower;
  vector<lower=0>[G] grid_upper;
  vector<lower=0>[G] grid_mid;
  vector<lower=0>[G] grid_width;

  array[G] int<lower=1, upper=B> grid_bin;

  vector<lower=1e-12, upper=1>[G] retention;

  vector[G] z_size;
  vector[G] q_size;
  matrix[G, J] B_shape;

  /*
   * Set to zero for a one-site taxon and one otherwise.
   */
  int<lower=0, upper=1> use_site_effects;

  /*
   * prior_only = 0 fits the likelihood.
   * prior_only = 1 samples only from the prior.
   */
  int<lower=0, upper=1> prior_only;
}

transformed data {
  vector[G] log_grid_width =
    log(
      grid_width
    );

  vector[G] log_retention =
    log(
      retention
    );

  for (g in 1:G) {
    if (
      grid_upper[g] <=
      grid_lower[g]
    ) {
      reject(
        "grid_upper must exceed grid_lower at grid cell ",
        g
      );
    }
  }

  /*
   * Every observed bin must contain at least one latent grid cell.
   */
  for (b in 1:B) {
    int cells =
      0;

    for (g in 1:G) {
      if (
        grid_bin[g] ==
        b
      ) {
        cells +=
          1;
      }
    }

    if (cells == 0) {
      reject(
        "Observed bin ",
        b,
        " contains no latent grid cells."
      );
    }
  }
}

parameters {
  /*
   * Positive global negative drift in standardized log body size.
   */
  real<lower=0, upper=5> beta_global;

  /*
   * Global quadratic tilt controlling broad dispersion and tail shape.
   */
  real gamma_global;

  /*
   * Taxon-wide smooth shape coefficients.
   */
  vector[J] coef_global_shape;

  /*
   * Partially pooled site location tilts.
   */
  real<lower=0, upper=2> tau_loc;
  vector[S] z_loc;

  /*
   * Partially pooled site dispersion tilts.
   */
  real<lower=0, upper=2> tau_disp;
  vector[S] z_disp;

  /*
   * Strongly regularized site-specific smooth departures.
   */
  real<lower=0, upper=1> tau_site_shape;
  matrix[S, J] z_site_shape;

  /*
   * Shared event-to-event Dirichlet-multinomial concentration.
   */
  real<lower=log(0.5), upper=log(1e5)> log_kappa;
}

model {
  real kappa =
    exp(
      log_kappa
    );

  vector[G] eta_global =
    -beta_global *
    z_size +
    gamma_global *
    q_size +
    B_shape *
    coef_global_shape;

  /*
   * This local array is required by the likelihood but is not written
   * to the posterior output.
   */
  array[S] vector[B] p_observed_site;

  beta_global ~
    normal(
      1,
      0.7
    );

  gamma_global ~
    normal(
      0,
      0.5
    );

  coef_global_shape ~
    normal(
      0,
      0.35
    );

  tau_loc ~
    normal(
      0,
      0.35
    );

  z_loc ~
    std_normal();

  tau_disp ~
    normal(
      0,
      0.25
    );

  z_disp ~
    std_normal();

  tau_site_shape ~
    normal(
      0,
      0.15
    );

  to_vector(
    z_site_shape
  ) ~
    std_normal();

  log_kappa ~
    normal(
      log(50),
      1.0
    );

  /*
   * Construct each site's retained observed-bin distribution locally.
   */
  for (s in 1:S) {
    real location_s =
      0;

    real dispersion_s =
      0;

    vector[J] site_shape_s =
      rep_vector(
        0,
        J
      );

    vector[G] eta_site;

    vector[G] log_true_mass;

    if (
      use_site_effects ==
      1
    ) {
      location_s =
        tau_loc *
        z_loc[s];

      dispersion_s =
        tau_disp *
        z_disp[s];

      site_shape_s =
        tau_site_shape *
        to_vector(
          z_site_shape[s]'
        );
    }

    eta_site =
      eta_global +
      location_s *
      z_size +
      dispersion_s *
      q_size +
      B_shape *
      site_shape_s;

    log_true_mass =
      eta_site +
      log_grid_width;

    p_observed_site[s] =
      observed_bin_probability(
        log_true_mass,
        log_retention,
        grid_bin,
        B
      );
  }

  if (
    prior_only ==
    0
  ) {
    for (event in 1:K) {
      vector[B] alpha_event =
        rep_vector(
          1e-9,
          B
        ) +
        kappa *
        p_observed_site[
          site_id[event]
        ];

      target +=
        dirichlet_multinomial_logpmf(
          y[event],
          alpha_event
        );
    }
  }
}

generated quantities {
  /*
   * Compact scalar and site-level production outputs only.
   */
  real<lower=0> kappa =
    exp(
      log_kappa
    );

  vector[S] size_mean;
  vector[S] size_sd;

  vector[S] size_q50;
  vector[S] size_q75;
  vector[S] size_q95;
  vector[S] size_q99;

  /*
   * Expected proportion of the latent population retained by the mesh.
   */
  vector[S] mean_retention;

  /*
   * Compact shape and upper-grid diagnostics.
   */
  vector[S] decline_fraction;
  vector[S] top_bin_probability;
  array[S] int<lower=0, upper=1> q99_in_top_bin;

  /*
   * Nested scope: all grid-sized variables declared inside this block
   * are temporary and are not written to the output CSV.
   */
  {
    vector[G] eta_global_local =
      -beta_global *
      z_size +
      gamma_global *
      q_size +
      B_shape *
      coef_global_shape;

    for (s in 1:S) {
    real location_s =
      0;

    real dispersion_s =
      0;

    vector[J] site_shape_s =
      rep_vector(
        0,
        J
      );

    vector[G] eta_site;

    vector[G] log_true_mass;

    vector[G] p_true_grid_local;

    vector[B] p_true_bin_local;

    real mean_s;

    real variance_s =
      0;

    real decreasing_count =
      0;

    real cumulative_before_top;

    if (
      use_site_effects ==
      1
    ) {
      location_s =
        tau_loc *
        z_loc[s];

      dispersion_s =
        tau_disp *
        z_disp[s];

      site_shape_s =
        tau_site_shape *
        to_vector(
          z_site_shape[s]'
        );
    }

    eta_site =
      eta_global_local +
      location_s *
      z_size +
      dispersion_s *
      q_size +
      B_shape *
      site_shape_s;

    log_true_mass =
      eta_site +
      log_grid_width;

    p_true_grid_local =
      softmax(
        log_true_mass
      );

    p_true_bin_local =
      true_bin_probability(
        p_true_grid_local,
        grid_bin,
        B
      );

    mean_s =
      dot_product(
        p_true_grid_local,
        grid_mid
      );

    size_mean[s] =
      mean_s;

    for (g in 1:G) {
      variance_s +=
        p_true_grid_local[g] *
        (
          square(
            grid_mid[g] -
            mean_s
          ) +
          square(
            grid_width[g]
          ) /
          12
        );
    }

    size_sd[s] =
      sqrt(
        fmax(
          variance_s,
          0
        )
      );

    size_q50[s] =
      latent_grid_quantile(
        p_true_grid_local,
        grid_lower,
        grid_upper,
        0.50
      );

    size_q75[s] =
      latent_grid_quantile(
        p_true_grid_local,
        grid_lower,
        grid_upper,
        0.75
      );

    size_q95[s] =
      latent_grid_quantile(
        p_true_grid_local,
        grid_lower,
        grid_upper,
        0.95
      );

    size_q99[s] =
      latent_grid_quantile(
        p_true_grid_local,
        grid_lower,
        grid_upper,
        0.99
      );

    mean_retention[s] =
      dot_product(
        p_true_grid_local,
        retention
      );

    for (g in 1:(G - 1)) {
      real density_current =
        p_true_grid_local[g] /
        grid_width[g];

      real density_next =
        p_true_grid_local[g + 1] /
        grid_width[g + 1];

      if (
        density_next <
        density_current
      ) {
        decreasing_count +=
          1;
      }
    }

    decline_fraction[s] =
      decreasing_count /
      (G - 1.0);

    top_bin_probability[s] =
      p_true_bin_local[B];

    cumulative_before_top =
      sum(
        head(
          p_true_bin_local,
          B - 1
        )
      );

    q99_in_top_bin[s] =
      cumulative_before_top <
      0.99;
    }
  }
}
