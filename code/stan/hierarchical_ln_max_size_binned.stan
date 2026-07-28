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
   * The R data builder sets use_site_effects = 0 for a single-site
   * taxon and 1 otherwise. This avoids confounding a one-site effect
   * with the global profile.
   */
  int<lower=0, upper=1> use_site_effects;

  int<lower=0, upper=1> prior_only;
  int<lower=0, upper=1> generate_y_rep;
}

transformed data {
  array[K] int<lower=0> n_event;

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

  for (event in 1:K) {
    int total =
      0;

    for (b in 1:B) {
      total +=
        y[event, b];
    }

    n_event[event] =
      total;
  }

  /*
   * Every observed bin must contain at least one fine-grid cell.
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
   * The term in eta is -beta_global * z_size.
   */
  real<lower=0, upper=5> beta_global;

  /*
   * Global quadratic tilt controlling broad dispersion/tail shape.
   */
  real gamma_global;

  /*
   * Taxon-wide flexible smooth shape.
   */
  vector[J] coef_global_shape;

  /*
   * Partially pooled site location tilt.
   */
  real<lower=0, upper=2> tau_loc;
  vector[S] z_loc;

  /*
   * Partially pooled site dispersion tilt.
   */
  real<lower=0, upper=2> tau_disp;
  vector[S] z_disp;

  /*
   * Strongly regularized site-specific smooth deviations.
   */
  real<lower=0, upper=1> tau_site_shape;
  matrix[S, J] z_site_shape;

  /*
   * Shared event-to-event Dirichlet-multinomial concentration.
   */
  real<lower=log(0.5), upper=log(1e5)> log_kappa;
}

transformed parameters {
  real<lower=0> kappa =
    exp(
      log_kappa
    );

  vector[S] location_tilt;
  vector[S] dispersion_tilt;
  matrix[S, J] site_shape_coefficient;

  vector[G] eta_global =
    -beta_global *
    z_size +
    gamma_global *
    q_size +
    B_shape *
    coef_global_shape;

  array[S] vector[G] p_true_grid;
  array[S] vector[B] p_true_bin;
  array[S] vector[B] p_observed_bin;

  if (
    use_site_effects ==
    1
  ) {
    location_tilt =
      tau_loc *
      z_loc;

    dispersion_tilt =
      tau_disp *
      z_disp;

    site_shape_coefficient =
      tau_site_shape *
      z_site_shape;
  } else {
    location_tilt =
      rep_vector(
        0,
        S
      );

    dispersion_tilt =
      rep_vector(
        0,
        S
      );

    site_shape_coefficient =
      rep_matrix(
        0,
        S,
        J
      );
  }

  for (s in 1:S) {
    vector[J] site_shape =
      to_vector(
        site_shape_coefficient[s]'
      );

    vector[G] eta_site =
      eta_global +
      location_tilt[s] *
      z_size +
      dispersion_tilt[s] *
      q_size +
      B_shape *
      site_shape;

    vector[G] log_true_mass =
      eta_site +
      log_grid_width;

    vector[B] true_bin_mass =
      rep_vector(
        0,
        B
      );

    vector[B] log_retained_bin_mass =
      rep_vector(
        negative_infinity(),
        B
      );

    p_true_grid[s] =
      softmax(
        log_true_mass
      );

    for (g in 1:G) {
      int b =
        grid_bin[g];

      true_bin_mass[b] +=
        p_true_grid[s][g];

      log_retained_bin_mass[b] =
        log_sum_exp(
          log_retained_bin_mass[b],
          log_true_mass[g] +
          log_retention[g]
        );
    }

    p_true_bin[s] =
      true_bin_mass;

    p_observed_bin[s] =
      softmax(
        log_retained_bin_mass
      );
  }
}

model {
  /*
   * A soft biological expectation of decreasing abundance with size.
   * Positive local slopes remain possible through the quadratic and
   * smooth terms.
   */
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
        p_observed_bin[
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
  vector[K] log_lik;

  vector[S] size_mean;
  vector[S] size_sd;

  vector[S] size_q50;
  vector[S] size_q75;
  vector[S] size_q95;
  vector[S] size_q99;

  /*
   * Expected fraction of the latent population retained by the mesh.
   */
  vector[S] mean_retention;

  /*
   * Fraction of adjacent fine-grid density values that decline with
   * increasing size. This is a descriptive posterior diagnostic, not
   * a hard constraint.
   */
  vector[S] decline_fraction;

  /*
   * Probability in the largest modeled observed class and indicator
   * that Q99 lies in that class. These diagnose upper-grid sensitivity.
   */
  vector[S] top_bin_probability;
  array[S] int<lower=0, upper=1> q99_in_top_bin;

  array[K, B] int<lower=0> y_rep;

  for (event in 1:K) {
    if (
      prior_only ==
      0
    ) {
      vector[B] alpha_event =
        rep_vector(
          1e-9,
          B
        ) +
        kappa *
        p_observed_bin[
          site_id[event]
        ];

      log_lik[event] =
        dirichlet_multinomial_logpmf(
          y[event],
          alpha_event
        );
    } else {
      log_lik[event] =
        0;
    }

    if (
      generate_y_rep ==
      1
    ) {
      vector[B] alpha_event =
        rep_vector(
          1e-9,
          B
        ) +
        kappa *
        p_observed_bin[
          site_id[event]
        ];

      vector[B] event_probability =
        dirichlet_rng(
          alpha_event
        );

      y_rep[event] =
        multinomial_rng(
          event_probability,
          n_event[event]
        );
    } else {
      y_rep[event] =
        rep_array(
          0,
          B
        );
    }
  }

  for (s in 1:S) {
    real mean_s =
      dot_product(
        p_true_grid[s],
        grid_mid
      );

    real variance_s =
      0;

    real decreasing_count =
      0;

    real cumulative_before_top =
      sum(
        head(
          p_true_bin[s],
          B - 1
        )
      );

    size_mean[s] =
      mean_s;

    for (g in 1:G) {
      variance_s +=
        p_true_grid[s][g] *
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
        p_true_grid[s],
        grid_lower,
        grid_upper,
        0.50
      );

    size_q75[s] =
      latent_grid_quantile(
        p_true_grid[s],
        grid_lower,
        grid_upper,
        0.75
      );

    size_q95[s] =
      latent_grid_quantile(
        p_true_grid[s],
        grid_lower,
        grid_upper,
        0.95
      );

    size_q99[s] =
      latent_grid_quantile(
        p_true_grid[s],
        grid_lower,
        grid_upper,
        0.99
      );

    mean_retention[s] =
      dot_product(
        p_true_grid[s],
        retention
      );

    for (g in 1:(G - 1)) {
      real density_current =
        p_true_grid[s][g] /
        grid_width[g];

      real density_next =
        p_true_grid[s][g + 1] /
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
      p_true_bin[s][B];

    q99_in_top_bin[s] =
      cumulative_before_top <
      0.99;
  }
}
