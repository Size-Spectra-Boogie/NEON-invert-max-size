

#'
#'
#'

fit_max_model_named = function(df = NULL,
                               mod_name = NULL,
                               rerun = FALSE,
                               overwrite = FALSE,
                               model_types = c('evt','evt_gumbel','efs','efsmm'),
                               ...){
  ## currently, we are saving this locally to avoid pushing to git
  mod_path = paste0(here('ignore/models'),"/",mod_name,'.rds')
  if(any(rerun, !file.exists(mod_path))){
    if(all(file.exists(mod_path),!overwrite)){
      warning('Model file already exists and `overwrite` = FALSE. Set to TRUE to overwrite existing files.')
      break()
    }
    
    fit = fishmax::fit_max_model(
      df,
      model_type = model_types,
      output_dir = here('ignore/models/'),
      refresh = 0
      )
    
    saveRDS(fit, mod_path)
    rm(fit)
    print(paste("Model saved as:",mod_name,".rds"))
  } else{
    print(paste("Model ",mod_name," exists. To overwrite, set `rerun` = TRUE and `overwrite` = TRUE"))
  }
}
fit_max_model_named = purrr::safely(fit_max_model_named)

#'
#'
#'

get_max_wide = function(mod_path = NULL){
  mod = readRDS(mod_path)
  site_taxa = gsub("\\.rds","",lapply(strsplit(mod_path, "/"),"[", 10))
  tab = fishmax::get_max(mod) %>% 
    dplyr::mutate(site_taxa = site_taxa) %>% 
    dplyr::select(site_taxa, everything())
  return(tab)
}

#'
#'
#'

# simulate the sampling from the truncated normal body size distribution
simulation = function(simPars, seed = 1312){
  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Package 'truncnorm' must be installed.")
  }
  
  required_columns <- c("site", "mu", "sigma", "k", "n")
  
  if (!all(required_columns %in% names(simPars))) {
    stop(
      "simPars must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (any(simPars$mu <= 0)) {
    stop("All mu values must be positive.")
  }
  
  if (any(simPars$sigma <= 0)) {
    stop("All sigma values must be positive.")
  }
  
  if (any(simPars$k < 1)) {
    stop("Each site must have at least one sampling event.")
  }
  
  if (!all(lengths(simPars$n) == simPars$k)) {
    stop(
      "For every site, length(n[[site]]) must equal k."
    )
  }
  
  if (any(unlist(simPars$n) < 1)) {
    stop("All event-level sample sizes must be at least 1.")
  }
  
  set.seed(seed)
  
  S <- nrow(simPars)
  K <- sum(simPars$k)
  
  event_table <- vector("list", K)
  size_blocks <- vector("list", K)
  
  event_id <- 0L
  
  for (s in seq_len(S)) {
    for (j in seq_len(simPars$k[s])) {
      event_id <- event_id + 1L
      
      n_j <- as.integer(simPars$n[[s]][j])
      
      x_j <- truncnorm::rtruncnorm(
        n = n_j,
        a = 0,
        b = Inf,
        mean = simPars$mu[s],
        sd = simPars$sigma[s]
      )
      
      # Each block contains all observations from one event.
      size_blocks[[event_id]] <- x_j
      
      event_table[[event_id]] <- data.frame(
        event_id = event_id,
        event_within_site = j,
        site = simPars$site[s],
        site_id = s,
        n_per_sample = n_j,
        stringsAsFactors = FALSE
      )
    }
  }
  
  event_table <- do.call(rbind, event_table)
  
  # Concatenation preserves contiguous observations within events.
  x <- unlist(
    size_blocks,
    use.names = FALSE
  )
  
  event_table$start_idx <-
    cumsum(event_table$n_per_sample) -
    event_table$n_per_sample +
    1L
  
  observation_table <- data.frame(
    event_id = rep(
      event_table$event_id,
      times = event_table$n_per_sample
    ),
    event_within_site = rep(
      event_table$event_within_site,
      times = event_table$n_per_sample
    ),
    site = rep(
      event_table$site,
      times = event_table$n_per_sample
    ),
    site_id = rep(
      event_table$site_id,
      times = event_table$n_per_sample
    ),
    x = x,
    stringsAsFactors = FALSE
  )
  
  stan_data <- list(
    S = as.integer(S),
    K = as.integer(K),
    n_obs = as.integer(length(x)),
    x = as.numeric(x),
    n_per_sample = as.integer(
      event_table$n_per_sample
    ),
    start_idx = as.integer(
      event_table$start_idx
    ),
    site_id = as.integer(
      event_table$site_id
    )
  )
  
  site_truth <- data.frame(
    site_id = seq_len(S),
    site = simPars$site,
    mu_true = simPars$mu,
    sigma_true = simPars$sigma,
    k_observed = simPars$k,
    mean_n_observed = vapply(
      simPars$n,
      mean,
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  
  list(
    stan_data = stan_data,
    observations = observation_table,
    events = event_table,
    site_truth = site_truth
  )
}

# Retention-corrected latent-grid data preparation
#
# This script converts one taxon's site x event x 1-mm size-class data
# into the structure required by:
#
#   03_retention_corrected_size_composition.stan
#
# Observation convention used here:
#
# - The latent biological size range starts at latent_min_size, default
#   0.1 mm.
# - Observed sizeClass == 1 contains all retained individuals from
#   latent_min_size through 1.5 mm.
# - Observed classes b >= 2 represent (b - 0.5, b + 0.5] mm.
# - Retention is evaluated on a fine latent grid and then aggregated
#   into the observed 1-mm classes.
#
# The model is fit to one acceptedTaxonID at a time.
#
# IMPORTANT ABOUT THE MORIN FUNCTION
# The linear predictor supplied by the user is:
#
#   eta = a + b*log10(RL) + c*log10(RL)*log10(M)
#
# A probability must lie in [0, 1]. The default "provided_capped"
# implementation evaluates the supplied empirical expression
#
#   exp(eta) / (1 + exp(1.8))
#
# and caps it at one once the empirical curve reaches complete retention.
# The alternative "logistic" option uses plogis(eta) as a sensitivity
# model. A custom externally calibrated retention function can instead
# be supplied through retention_fun.

morin_retention_probability <- function(
    size_mm,
    a = -2.84,
    b = 5.80,
    c = -3.18,
    mesh_mm = 0.25,
    link = c("provided_capped", "logistic")
) {
  link <- match.arg(link)
  
  if (
    any(!is.finite(size_mm)) ||
    any(size_mm <= 0)
  ) {
    stop("size_mm must contain finite positive values.")
  }
  
  if (
    length(mesh_mm) != 1L ||
    !is.finite(mesh_mm) ||
    mesh_mm <= 0
  ) {
    stop("mesh_mm must be a finite positive scalar.")
  }
  
  relative_length <- size_mm / mesh_mm
  
  eta <-
    a +
    b * log10(relative_length) +
    c *
    log10(relative_length) *
    log10(mesh_mm)
  
  probability <- switch(
    link,
    logistic = stats::plogis(eta),
    provided_capped = pmin(
      1,
      pmax(
        0,
        exp(eta) /
          (
            1 + exp(eta)
          )
      )
    )
  )
  
  as.numeric(probability)
}


make_retention_grid_stan_data <- function(
    bin_data,
    event_data = NULL,
    site_col = "siteID",
    event_cols = "collectYear",
    taxon_col = "acceptedTaxonID",
    size_class_col = "sizeClass",
    count_col = "no_m2",
    latent_min_size = 0.1,
    max_size_class = NULL,
    extra_upper_bins = 2L,
    grid_step = 0.05,
    spline_df = 7L,
    mesh_mm = 0.25,
    morin_a = -2.84,
    morin_b = 5.80,
    morin_c = -3.18,
    morin_link = c("provided_capped", "logistic"),
    retention_fun = NULL,
    prior_only = 0L,
    generate_y_rep = 1L,
    integer_tolerance = 1e-8
) {
  morin_link <- match.arg(morin_link)
  
  if (!is.data.frame(bin_data)) {
    stop("bin_data must be a data.frame.")
  }
  
  required_bin_columns <- unique(
    c(
      site_col,
      event_cols,
      taxon_col,
      size_class_col,
      count_col
    )
  )
  
  missing_columns <- setdiff(
    required_bin_columns,
    names(bin_data)
  )
  
  if (length(missing_columns)) {
    stop(
      "bin_data is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  if (
    length(event_cols) < 1L ||
    any(!nzchar(event_cols))
  ) {
    stop("event_cols must contain at least one valid column name.")
  }
  
  if (
    length(latent_min_size) != 1L ||
    !is.finite(latent_min_size) ||
    latent_min_size <= 0 ||
    latent_min_size >= 1.5
  ) {
    stop(
      "latent_min_size must be a finite positive value below 1.5 mm."
    )
  }
  
  if (
    length(grid_step) != 1L ||
    !is.finite(grid_step) ||
    grid_step <= 0
  ) {
    stop("grid_step must be a finite positive scalar.")
  }
  
  if (
    length(spline_df) != 1L ||
    !is.finite(spline_df) ||
    spline_df < 5L
  ) {
    stop("spline_df must be an integer of at least 5.")
  }
  
  if (
    !prior_only %in% c(0L, 1L) ||
    !generate_y_rep %in% c(0L, 1L)
  ) {
    stop("prior_only and generate_y_rep must each be 0 or 1.")
  }
  
  taxon_values <- unique(
    as.character(
      bin_data[[taxon_col]]
    )
  )
  
  taxon_values <- taxon_values[
    !is.na(taxon_values)
  ]
  
  if (length(taxon_values) != 1L) {
    stop(
      "bin_data must contain exactly one nonmissing taxon. Found: ",
      paste(taxon_values, collapse = ", ")
    )
  }
  
  taxon_id <- taxon_values[1L]
  
  make_key <- function(data, columns) {
    key_parts <- lapply(
      data[columns],
      function(value) {
        value <- as.character(value)
        
        if (anyNA(value)) {
          stop(
            "Missing values are not permitted in key columns: ",
            paste(columns, collapse = ", ")
          )
        }
        
        value
      }
    )
    
    do.call(
      paste,
      c(
        key_parts,
        sep = "\r"
      )
    )
  }
  
  bins <- data.frame(
    site_value = as.character(
      bin_data[[site_col]]
    ),
    size_class_raw = as.numeric(
      bin_data[[size_class_col]]
    ),
    count_raw = as.numeric(
      bin_data[[count_col]]
    ),
    stringsAsFactors = FALSE
  )
  
  for (column in event_cols) {
    bins[[column]] <- bin_data[[column]]
  }
  
  if (
    anyNA(bins$site_value) ||
    anyNA(bins$size_class_raw) ||
    anyNA(bins$count_raw)
  ) {
    stop(
      "site, size class, and count values must not be missing."
    )
  }
  
  if (
    any(!is.finite(bins$size_class_raw)) ||
    any(!is.finite(bins$count_raw)) ||
    any(bins$count_raw < 0)
  ) {
    stop(
      "Size classes must be finite and counts must be finite and nonnegative."
    )
  }
  
  noninteger_class <-
    abs(
      bins$size_class_raw -
        round(bins$size_class_raw)
    ) > integer_tolerance
  
  if (any(noninteger_class)) {
    stop(
      "sizeClass must contain integer 1-mm classes. First invalid value: ",
      bins$size_class_raw[
        which(noninteger_class)[1L]
      ]
    )
  }
  
  noninteger_count <-
    abs(
      bins$count_raw -
        round(bins$count_raw)
    ) > integer_tolerance
  
  if (any(noninteger_count)) {
    stop(
      count_col,
      " must contain integer counts for a Dirichlet-multinomial model. ",
      "First noninteger value: ",
      bins$count_raw[
        which(noninteger_count)[1L]
      ]
    )
  }
  
  bins$size_class <- as.integer(
    round(bins$size_class_raw)
  )
  
  bins$count <- as.integer(
    round(bins$count_raw)
  )
  
  if (any(bins$size_class < 1L)) {
    stop("All size classes must be positive integers beginning at 1.")
  }
  
  aggregate_columns <- c(
    "site_value",
    event_cols,
    "size_class"
  )
  
  bins <- stats::aggregate(
    bins["count"],
    by = bins[aggregate_columns],
    FUN = sum
  )
  
  bins <- bins[
    bins$count > 0L,
    ,
    drop = FALSE
  ]
  
  if (!nrow(bins)) {
    stop("No positive event-by-size-class counts remain.")
  }
  
  observed_max_class <- max(
    bins$size_class
  )
  
  if (is.null(max_size_class)) {
    max_size_class <-
      observed_max_class +
      as.integer(extra_upper_bins)
  }
  
  if (
    length(max_size_class) != 1L ||
    !is.finite(max_size_class) ||
    abs(max_size_class - round(max_size_class)) >
    integer_tolerance
  ) {
    stop("max_size_class must be a finite integer.")
  }
  
  max_size_class <- as.integer(
    round(max_size_class)
  )
  
  if (max_size_class < observed_max_class) {
    stop(
      "max_size_class cannot be below the largest observed class, ",
      observed_max_class,
      "."
    )
  }
  
  B <- max_size_class
  
  # Event roster. Zero-total events can be included when event_data is
  # supplied, although they contain no direct size-composition information.
  if (is.null(event_data)) {
    events <- unique(
      bins[
        c(
          "site_value",
          event_cols
        )
      ]
    )
  } else {
    if (!is.data.frame(event_data)) {
      stop("event_data must be NULL or a data.frame.")
    }
    
    required_event_columns <- unique(
      c(
        site_col,
        event_cols
      )
    )
    
    missing_event_columns <- setdiff(
      required_event_columns,
      names(event_data)
    )
    
    if (length(missing_event_columns)) {
      stop(
        "event_data is missing required columns: ",
        paste(missing_event_columns, collapse = ", ")
      )
    }
    
    event_source <- event_data
    
    if (taxon_col %in% names(event_source)) {
      event_source <- event_source[
        as.character(
          event_source[[taxon_col]]
        ) == taxon_id,
        ,
        drop = FALSE
      ]
    }
    
    events <- data.frame(
      site_value = as.character(
        event_source[[site_col]]
      ),
      stringsAsFactors = FALSE
    )
    
    for (column in event_cols) {
      events[[column]] <- event_source[[column]]
    }
    
    # Only retain sites at which the taxon was observed at least once.
    # Sites with no individuals have no information about size shape.
    positive_sites <- unique(
      bins$site_value
    )
    
    events <- events[
      events$site_value %in%
        positive_sites,
      ,
      drop = FALSE
    ]
    
    events <- unique(events)
  }
  
  if (!nrow(events)) {
    stop("No sampling events remain after constructing the event roster.")
  }
  
  key_columns <- c(
    "site_value",
    event_cols
  )
  
  events$event_key <- make_key(
    events,
    key_columns
  )
  
  bins$event_key <- make_key(
    bins,
    key_columns
  )
  
  missing_positive_events <- setdiff(
    unique(bins$event_key),
    events$event_key
  )
  
  if (length(missing_positive_events)) {
    stop(
      "Some positive event-by-bin observations are absent from event_data."
    )
  }
  
  site_levels <- unique(
    events$site_value
  )
  
  events$site_id <- match(
    events$site_value,
    site_levels
  )
  
  events$event_id <- seq_len(
    nrow(events)
  )
  
  bins$event_id <- events$event_id[
    match(
      bins$event_key,
      events$event_key
    )
  ]
  
  if (anyNA(bins$event_id)) {
    stop("Failed to map one or more positive rows to an event.")
  }
  
  K <- nrow(events)
  S <- length(site_levels)
  
  y <- matrix(
    0L,
    nrow = K,
    ncol = B
  )
  
  for (row in seq_len(nrow(bins))) {
    event <- bins$event_id[row]
    size_class <- bins$size_class[row]
    
    if (size_class > B) {
      stop(
        "Internal error: observed size class exceeds B."
      )
    }
    
    y[event, size_class] <-
      y[event, size_class] +
      bins$count[row]
  }
  
  # Observed size-class intervals. Class 1 deliberately contains the
  # unresolved small sizes from latent_min_size through 1.5 mm.
  observed_bin_lower <- numeric(B)
  observed_bin_upper <- numeric(B)
  
  for (size_class in seq_len(B)) {
    observed_bin_lower[size_class] <-
      if (size_class == 1L) {
        latent_min_size
      } else {
        size_class - 0.5
      }
    
    observed_bin_upper[size_class] <-
      size_class + 0.5
  }
  
  # Construct fine latent cells separately inside each observed class.
  # This guarantees that no fine cell straddles an observed boundary.
  grid_parts <- vector(
    "list",
    B
  )
  
  for (size_class in seq_len(B)) {
    interval_width <-
      observed_bin_upper[size_class] -
      observed_bin_lower[size_class]
    
    number_subcells <- max(
      1L,
      as.integer(
        ceiling(
          interval_width /
            grid_step
        )
      )
    )
    
    edges <- seq(
      observed_bin_lower[size_class],
      observed_bin_upper[size_class],
      length.out = number_subcells + 1L
    )
    
    grid_parts[[size_class]] <- data.frame(
      grid_lower = head(
        edges,
        -1L
      ),
      grid_upper = tail(
        edges,
        -1L
      ),
      grid_bin = size_class
    )
  }
  
  grid <- do.call(
    rbind,
    grid_parts
  )
  
  rownames(grid) <- NULL
  
  grid$grid_mid <-
    0.5 *
    (
      grid$grid_lower +
        grid$grid_upper
    )
  
  grid$grid_width <-
    grid$grid_upper -
    grid$grid_lower
  
  G <- nrow(grid)
  
  # Five-point Gauss-Legendre quadrature for cell-average retention.
  quadrature_nodes <- c(
    -0.906179845938664,
    -0.538469310105683,
    0,
    0.538469310105683,
    0.906179845938664
  )
  
  quadrature_weights <- c(
    0.236926885056189,
    0.478628670499366,
    0.568888888888889,
    0.478628670499366,
    0.236926885056189
  )
  
  custom_retention_supplied <-
    !is.null(retention_fun)
  
  if (!custom_retention_supplied) {
    retention_fun <- function(size_mm) {
      morin_retention_probability(
        size_mm = size_mm,
        a = morin_a,
        b = morin_b,
        c = morin_c,
        mesh_mm = mesh_mm,
        link = morin_link
      )
    }
  } else {
    retention_fun <- match.fun(
      retention_fun
    )
  }
  
  cell_average_retention <- function(
    x_lower,
    x_upper
  ) {
    midpoint <-
      0.5 *
      (
        x_lower +
          x_upper
      )
    
    half_width <-
      0.5 *
      (
        x_upper -
          x_lower
      )
    
    evaluation_sizes <-
      midpoint +
      half_width *
      quadrature_nodes
    
    values <- as.numeric(
      retention_fun(
        evaluation_sizes
      )
    )
    
    if (
      length(values) !=
      length(evaluation_sizes) ||
      any(!is.finite(values)) ||
      any(values < 0) ||
      any(values > 1)
    ) {
      stop(
        "retention_fun must return one finite probability in [0,1] ",
        "for every supplied size."
      )
    }
    
    # The Gauss-Legendre weights sum to two; division by two converts
    # the quadrature integral into an interval average.
    sum(
      quadrature_weights *
        values
    ) / 2
  }
  
  retention <- mapply(
    FUN = cell_average_retention,
    x_lower = grid$grid_lower,
    x_upper = grid$grid_upper
  )
  
  # The empirical curve is positive over the modeled range, but use a
  # tiny numerical floor so log(retention) is always defined in Stan.
  retention <- pmin(
    1,
    pmax(
      1e-12,
      retention
    )
  )
  
  # Standardized log-size terms for location and dispersion tilts.
  log_size <- log(
    grid$grid_mid
  )
  
  weighted_mean <- function(
    value,
    weight
  ) {
    sum(
      value *
        weight
    ) /
      sum(weight)
  }
  
  log_size_mean <- weighted_mean(
    log_size,
    grid$grid_width
  )
  
  log_size_sd <- sqrt(
    weighted_mean(
      (
        log_size -
          log_size_mean
      )^2,
      grid$grid_width
    )
  )
  
  if (
    !is.finite(log_size_sd) ||
    log_size_sd <= 0
  ) {
    stop("Could not standardize the latent size grid.")
  }
  
  z_size <-
    (
      log_size -
        log_size_mean
    ) /
    log_size_sd
  
  q_raw <- z_size^2
  
  q_center <- weighted_mean(
    q_raw,
    grid$grid_width
  )
  
  q_size <- q_raw - q_center
  
  q_sd <- sqrt(
    weighted_mean(
      q_size^2,
      grid$grid_width
    )
  )
  
  if (
    is.finite(q_sd) &&
    q_sd > 0
  ) {
    q_size <- q_size / q_sd
  }
  
  # Smooth basis orthogonal to intercept, z_size, and q_size. This
  # keeps the flexible shape terms distinct from the explicit site
  # location and dispersion tilts.
  raw_basis <- splines::bs(
    log_size,
    df = as.integer(spline_df),
    degree = 3L,
    intercept = TRUE,
    Boundary.knots = range(
      log_size
    )
  )
  
  fixed_design <- cbind(
    intercept = 1,
    z_size = z_size,
    q_size = q_size
  )
  
  projection_coefficients <- qr.solve(
    crossprod(
      fixed_design
    ),
    crossprod(
      fixed_design,
      raw_basis
    )
  )
  
  residual_basis <-
    raw_basis -
    fixed_design %*%
    projection_coefficients
  
  basis_svd <- svd(
    residual_basis,
    nu = min(
      nrow(residual_basis),
      ncol(residual_basis)
    ),
    nv = 0L
  )
  
  tolerance <- max(
    basis_svd$d,
    1
  ) *
    1e-9
  
  keep <- which(
    basis_svd$d >
      tolerance
  )
  
  if (!length(keep)) {
    stop(
      "No smooth basis columns remained after orthogonalization. ",
      "Increase spline_df."
    )
  }
  
  B_shape <-
    basis_svd$u[
      ,
      keep,
      drop = FALSE
    ] *
    sqrt(G)
  
  J <- ncol(
    B_shape
  )
  
  stan_data <- list(
    S = as.integer(S),
    K = as.integer(K),
    B = as.integer(B),
    G = as.integer(G),
    J = as.integer(J),
    site_id = as.integer(
      events$site_id
    ),
    y = unname(
      y
    ),
    grid_lower = as.numeric(
      grid$grid_lower
    ),
    grid_upper = as.numeric(
      grid$grid_upper
    ),
    grid_mid = as.numeric(
      grid$grid_mid
    ),
    grid_width = as.numeric(
      grid$grid_width
    ),
    grid_bin = as.integer(
      grid$grid_bin
    ),
    retention = as.numeric(
      retention
    ),
    z_size = as.numeric(
      z_size
    ),
    q_size = as.numeric(
      q_size
    ),
    B_shape = unname(
      B_shape
    ),
    use_site_effects = as.integer(
      S > 1L
    ),
    prior_only = as.integer(
      prior_only
    ),
    generate_y_rep = as.integer(
      generate_y_rep
    )
  )
  
  total_individuals_by_site <- vapply(
    seq_len(S),
    function(site) {
      sum(
        y[
          events$site_id == site,
          ,
          drop = FALSE
        ]
      )
    },
    numeric(1)
  )
  
  site_mapping <- data.frame(
    site_id = seq_along(
      site_levels
    ),
    siteID = site_levels,
    n_events = as.integer(
      tabulate(
        events$site_id,
        nbins = S
      )
    ),
    total_individuals = as.integer(
      total_individuals_by_site
    ),
    stringsAsFactors = FALSE
  )
  
  event_mapping <- data.frame(
    event_id = events$event_id,
    site_id = events$site_id,
    stringsAsFactors = FALSE
  )
  
  event_mapping[[site_col]] <-
    events$site_value
  
  for (column in event_cols) {
    event_mapping[[column]] <-
      events[[column]]
  }
  
  event_mapping$n_individuals <-
    rowSums(y)
  
  bin_mapping <- data.frame(
    bin_id = seq_len(B),
    sizeClass = seq_len(B),
    lower_mm = observed_bin_lower,
    upper_mm = observed_bin_upper,
    observed_anywhere = colSums(y) > 0,
    total_count = colSums(y),
    stringsAsFactors = FALSE
  )
  
  grid_mapping <- data.frame(
    grid_id = seq_len(G),
    lower_mm = grid$grid_lower,
    upper_mm = grid$grid_upper,
    midpoint_mm = grid$grid_mid,
    width_mm = grid$grid_width,
    observed_bin = grid$grid_bin,
    retention_probability = retention,
    stringsAsFactors = FALSE
  )
  
  metadata <- list(
    acceptedTaxonID = taxon_id,
    latent_min_size = latent_min_size,
    max_size_class = B,
    latent_max_size = max(
      grid$grid_upper
    ),
    observed_max_class = observed_max_class,
    extra_upper_bins = B - observed_max_class,
    mesh_mm = mesh_mm,
    morin_coefficients = c(
      a = morin_a,
      b = morin_b,
      c = morin_c
    ),
    morin_link = if (
      custom_retention_supplied
    ) {
      "custom"
    } else {
      morin_link
    },
    grid_step_target = grid_step,
    total_grid_cells = G,
    smooth_basis_columns = J
  )
  
  list(
    stan_data = stan_data,
    mapping = list(
      site = site_mapping,
      event = event_mapping,
      bins = bin_mapping,
      grid = grid_mapping
    ),
    metadata = metadata
  )
}


# Initialization for retention-corrected latent-grid model
# --------------------------------------------------------
#
# Matches:
#
#   hierarchical_ln_max_size_binned_slim.stan
#
# Initial values are obtained by:
#
# 1. approximately correcting pooled observed bin counts for retention;
# 2. fitting the global negative-drift + quadratic + smooth profile;
# 3. optimizing site-specific location and dispersion tilts;
# 4. estimating a rough shared Dirichlet-multinomial concentration from
#    event-to-event compositional variability.

make_retention_grid_init <- function(
    stan_data,
    chains = 4L,
    seed = 1234L,
    pseudocount = 0.5
) {
  required <- c(
    "S",
    "K",
    "B",
    "G",
    "J",
    "site_id",
    "y",
    "grid_width",
    "grid_bin",
    "retention",
    "z_size",
    "q_size",
    "B_shape",
    "use_site_effects"
  )
  
  missing <- setdiff(
    required,
    names(stan_data)
  )
  
  if (length(missing)) {
    stop(
      "stan_data is missing: ",
      paste(missing, collapse = ", ")
    )
  }
  
  S <- as.integer(
    stan_data$S
  )
  
  K <- as.integer(
    stan_data$K
  )
  
  B <- as.integer(
    stan_data$B
  )
  
  G <- as.integer(
    stan_data$G
  )
  
  J <- as.integer(
    stan_data$J
  )
  
  chains <- as.integer(
    chains
  )
  
  if (
    S < 1L ||
    K < 1L ||
    B < 2L ||
    G < B ||
    J < 1L ||
    chains < 1L
  ) {
    stop("Invalid S, K, B, G, J, or chains in stan_data.")
  }
  
  y <- as.matrix(
    stan_data$y
  )
  
  if (
    nrow(y) != K ||
    ncol(y) != B ||
    any(!is.finite(y)) ||
    any(y < 0)
  ) {
    stop("stan_data$y must be a nonnegative K x B matrix.")
  }
  
  site_id <- as.integer(
    stan_data$site_id
  )
  
  grid_width <- as.numeric(
    stan_data$grid_width
  )
  
  grid_bin <- as.integer(
    stan_data$grid_bin
  )
  
  retention <- as.numeric(
    stan_data$retention
  )
  
  z_size <- as.numeric(
    stan_data$z_size
  )
  
  q_size <- as.numeric(
    stan_data$q_size
  )
  
  B_shape <- as.matrix(
    stan_data$B_shape
  )
  
  if (
    length(site_id) != K ||
    any(site_id < 1L) ||
    any(site_id > S)
  ) {
    stop("site_id must have length K and values between 1 and S.")
  }
  
  if (
    length(grid_width) != G ||
    length(grid_bin) != G ||
    length(retention) != G ||
    length(z_size) != G ||
    length(q_size) != G ||
    nrow(B_shape) != G ||
    ncol(B_shape) != J
  ) {
    stop("Latent-grid dimensions in stan_data are inconsistent.")
  }
  
  if (
    any(grid_width <= 0) ||
    any(retention <= 0) ||
    any(retention > 1) ||
    any(grid_bin < 1L) ||
    any(grid_bin > B)
  ) {
    stop("Invalid grid width, retention probability, or grid-bin index.")
  }
  
  clamp <- function(
    value,
    lower,
    upper
  ) {
    pmin(
      pmax(
        value,
        lower
      ),
      upper
    )
  }
  
  log_sum_exp_r <- function(value) {
    maximum <- max(value)
    
    maximum +
      log(
        sum(
          exp(
            value -
              maximum
          )
        )
      )
  }
  
  aggregate_grid_mass <- function(
    grid_mass
  ) {
    output <- numeric(B)
    
    for (g in seq_len(G)) {
      output[
        grid_bin[g]
      ] <-
        output[
          grid_bin[g]
        ] +
        grid_mass[g]
    }
    
    output
  }
  
  observed_bin_probability <- function(
    eta
  ) {
    log_true_mass <-
      eta +
      log(grid_width)
    
    log_true_mass <-
      log_true_mass -
      log_sum_exp_r(
        log_true_mass
      )
    
    true_grid_probability <-
      exp(log_true_mass)
    
    retained_grid_mass <-
      true_grid_probability *
      retention
    
    observed_bin_mass <-
      aggregate_grid_mass(
        retained_grid_mass
      )
    
    observed_bin_mass /
      sum(observed_bin_mass)
  }
  
  # Cell-width weighted average retention in each observed bin.
  bin_retention <- numeric(B)
  bin_width_total <- numeric(B)
  
  for (g in seq_len(G)) {
    bin_retention[
      grid_bin[g]
    ] <-
      bin_retention[
        grid_bin[g]
      ] +
      grid_width[g] *
      retention[g]
    
    bin_width_total[
      grid_bin[g]
    ] <-
      bin_width_total[
        grid_bin[g]
      ] +
      grid_width[g]
  }
  
  bin_retention <-
    bin_retention /
    bin_width_total
  
  bin_retention <- pmax(
    bin_retention,
    1e-6
  )
  
  site_counts <- matrix(
    0,
    nrow = S,
    ncol = B
  )
  
  for (event in seq_len(K)) {
    site_counts[
      site_id[event],
    ] <-
      site_counts[
        site_id[event],
      ] +
      y[event, ]
  }
  
  global_counts <- colSums(
    site_counts
  )
  
  # Approximate latent true bin probabilities after inverse-retention
  # correction. The pseudocount keeps empty upper bins finite.
  corrected_global_counts <-
    (
      global_counts +
        pseudocount
    ) /
    bin_retention
  
  global_bin_probability <-
    corrected_global_counts /
    sum(corrected_global_counts)
  
  bin_total_width <- numeric(B)
  
  for (g in seq_len(G)) {
    bin_total_width[
      grid_bin[g]
    ] <-
      bin_total_width[
        grid_bin[g]
      ] +
      grid_width[g]
  }
  
  approximate_grid_probability <- numeric(G)
  
  for (g in seq_len(G)) {
    approximate_grid_probability[g] <-
      global_bin_probability[
        grid_bin[g]
      ] *
      grid_width[g] /
      bin_total_width[
        grid_bin[g]
      ]
  }
  
  approximate_log_density <-
    log(
      pmax(
        approximate_grid_probability /
          grid_width,
        1e-12
      )
    )
  
  approximate_log_density <-
    approximate_log_density -
    weighted.mean(
      approximate_log_density,
      w = approximate_grid_probability
    )
  
  # eta = beta*(-z) + gamma*q + B_shape*coef_global_shape.
  design <- cbind(
    negative_z = -z_size,
    q_size = q_size,
    B_shape
  )
  
  weights <-
    sqrt(
      approximate_grid_probability +
        1e-4 / G
    )
  
  weighted_design <-
    design *
    weights
  
  weighted_response <-
    approximate_log_density *
    weights
  
  ridge <- diag(
    c(
      0.05,
      0.10,
      rep(
        0.25,
        J
      )
    )
  )
  
  coefficient_estimate <- tryCatch(
    solve(
      crossprod(
        weighted_design
      ) +
        ridge,
      crossprod(
        weighted_design,
        weighted_response
      )
    ),
    error = function(e) {
      rep(
        0,
        2L + J
      )
    }
  )
  
  beta_global_target <- clamp(
    coefficient_estimate[1L],
    0.05,
    4.0
  )
  
  gamma_global_target <- clamp(
    coefficient_estimate[2L],
    -2.0,
    2.0
  )
  
  coef_global_target <- clamp(
    coefficient_estimate[
      2L +
        seq_len(J)
    ],
    -1.5,
    1.5
  )
  
  eta_global <-
    -beta_global_target *
    z_size +
    gamma_global_target *
    q_size +
    as.vector(
      B_shape %*%
        coef_global_target
    )
  
  fit_site_tilts <- function(site) {
    counts <- site_counts[
      site,
    ]
    
    if (
      sum(counts) <= 0 ||
      stan_data$use_site_effects == 0L
    ) {
      return(
        c(
          location = 0,
          dispersion = 0
        )
      )
    }
    
    objective <- function(parameter) {
      eta <-
        eta_global +
        parameter[1L] *
        z_size +
        parameter[2L] *
        q_size
      
      probability <-
        observed_bin_probability(
          eta
        )
      
      -sum(
        counts *
          log(
            pmax(
              probability,
              1e-12
            )
          )
      )
    }
    
    fit <- tryCatch(
      stats::optim(
        par = c(
          0,
          0
        ),
        fn = objective,
        method = "L-BFGS-B",
        lower = c(
          -2.5,
          -2.0
        ),
        upper = c(
          2.5,
          2.0
        ),
        control = list(
          maxit = 500
        )
      ),
      error = function(e) NULL
    )
    
    if (
      is.null(fit) ||
      any(!is.finite(fit$par))
    ) {
      return(
        c(
          location = 0,
          dispersion = 0
        )
      )
    }
    
    c(
      location = fit$par[1L],
      dispersion = fit$par[2L]
    )
  }
  
  site_tilt_target <- t(
    vapply(
      seq_len(S),
      fit_site_tilts,
      numeric(2)
    )
  )
  
  robust_scale <- function(
    value,
    minimum,
    maximum
  ) {
    if (
      length(value) <= 1L ||
      all(
        abs(value) <
        1e-10
      )
    ) {
      return(minimum)
    }
    
    scale_mad <- stats::mad(
      value,
      center = 0,
      constant = 1.4826
    )
    
    scale_rms <- sqrt(
      mean(
        value^2
      )
    )
    
    clamp(
      max(
        minimum,
        scale_mad,
        scale_rms
      ),
      minimum,
      maximum
    )
  }
  
  tau_loc_target <- robust_scale(
    site_tilt_target[
      ,
      "location"
    ],
    minimum = 0.10,
    maximum = 1.5
  )
  
  tau_disp_target <- robust_scale(
    site_tilt_target[
      ,
      "dispersion"
    ],
    minimum = 0.10,
    maximum = 1.25
  )
  
  z_loc_target <-
    site_tilt_target[
      ,
      "location"
    ] /
    tau_loc_target
  
  z_disp_target <-
    site_tilt_target[
      ,
      "dispersion"
    ] /
    tau_disp_target
  
  if (stan_data$use_site_effects == 0L) {
    z_loc_target[] <- 0
    z_disp_target[] <- 0
  }
  
  # Rough method-of-moments estimate for the common
  # Dirichlet-multinomial concentration.
  inflation_estimates <- numeric(0)
  event_sizes <- numeric(0)
  
  for (site in seq_len(S)) {
    event_index <- which(
      site_id == site &
        rowSums(y) > 0
    )
    
    if (length(event_index) < 2L) {
      next
    }
    
    pooled <- colSums(
      y[
        event_index,
        ,
        drop = FALSE
      ]
    )
    
    pooled_probability <-
      (
        pooled +
          0.5
      ) /
      sum(
        pooled +
          0.5
      )
    
    for (event in event_index) {
      N <- sum(
        y[event, ]
      )
      
      expected <- N *
        pooled_probability
      
      usable <- expected >
        1e-6
      
      degrees_freedom <-
        max(
          1,
          sum(usable) -
            1
        )
      
      pearson <-
        sum(
          (
            y[event, usable] -
              expected[usable]
          )^2 /
            expected[usable]
        )
      
      inflation_estimates <-
        c(
          inflation_estimates,
          max(
            1,
            pearson /
              degrees_freedom
          )
        )
      
      event_sizes <-
        c(
          event_sizes,
          N
        )
    }
  }
  
  if (
    length(inflation_estimates) &&
    length(event_sizes)
  ) {
    F <- stats::median(
      inflation_estimates
    )
    
    N_typical <- stats::median(
      event_sizes
    )
    
    kappa_target <- if (
      F > 1.01 &&
      N_typical > F
    ) {
      (
        N_typical -
          F
      ) /
        (
          F -
            1
        )
    } else {
      100
    }
  } else {
    kappa_target <- 50
  }
  
  kappa_target <- clamp(
    kappa_target,
    2,
    500
  )
  
  set.seed(
    as.integer(seed)
  )
  
  lapply(
    seq_len(
      as.integer(chains)
    ),
    function(chain) {
      beta_global <-
        clamp(
          beta_global_target *
            exp(
              stats::rnorm(
                1,
                0,
                0.03
              )
            ),
          0.02,
          4.5
        )
      
      gamma_global <-
        gamma_global_target +
        stats::rnorm(
          1,
          0,
          0.03
        )
      
      coef_global_shape <-
        coef_global_target +
        stats::rnorm(
          J,
          0,
          0.02
        )
      
      tau_loc <-
        clamp(
          tau_loc_target *
            exp(
              stats::rnorm(
                1,
                0,
                0.03
              )
            ),
          0.03,
          1.8
        )
      
      tau_disp <-
        clamp(
          tau_disp_target *
            exp(
              stats::rnorm(
                1,
                0,
                0.03
              )
            ),
          0.03,
          1.8
        )
      
      location_target_chain <-
        site_tilt_target[
          ,
          "location"
        ] +
        stats::rnorm(
          S,
          0,
          0.01
        )
      
      dispersion_target_chain <-
        site_tilt_target[
          ,
          "dispersion"
        ] +
        stats::rnorm(
          S,
          0,
          0.01
        )
      
      if (stan_data$use_site_effects == 0L) {
        location_target_chain[] <- 0
        dispersion_target_chain[] <- 0
      }
      
      list(
        beta_global = as.numeric(
          beta_global
        ),
        gamma_global = as.numeric(
          gamma_global
        ),
        coef_global_shape = as.numeric(
          coef_global_shape
        ),
        tau_loc = as.numeric(
          tau_loc
        ),
        z_loc = as.numeric(
          location_target_chain /
            tau_loc
        ),
        tau_disp = as.numeric(
          tau_disp
        ),
        z_disp = as.numeric(
          dispersion_target_chain /
            tau_disp
        ),
        tau_site_shape = as.numeric(
          0.05 *
            exp(
              stats::rnorm(
                1,
                0,
                0.03
              )
            )
        ),
        z_site_shape = matrix(
          stats::rnorm(
            S * J,
            0,
            0.01
          ),
          nrow = S,
          ncol = J
        ),
        log_kappa = as.numeric(
          log(kappa_target) +
            stats::rnorm(
              1,
              0,
              0.05
            )
        )
      )
    }
  )
}

#'
#'
#'
fit_ln_named = function(stan_data = NULL, taxaName = NULL, rerun = FALSE, overwrite = FALSE){
  filePath = paste0(here("ignore/models"),"/",taxaName,"_lnbin.rds")
  print(taxaName)
  if(any(rerun, !file.exists(filePath))){
    if(all(file.exists(filePath),!overwrite)){
      warning('Model file already exists and `overwrite` = FALSE. Set to TRUE to overwrite existing files.')
      return(NULL)
    }
    chains = 4L
    # remove unneeded variable for the slim model
    stan_data$generate_y_rep <- NULL
    # this model accounts for a single site internally
    ## set initialization values.
      init_list = make_retention_grid_init(
        stan_data = stan_data,
        chains = chains,
        seed = 1312
      )

      fit = multi_mod_ln_slim$sample(
        data = stan_data,
        seed = 1312,
        chains = chains,
        parallel_chains = chains,
        
        init = init_list,
        
        iter_warmup = 1500,
        iter_sampling = 1000,
        adapt_delta = 0.99,
        max_treedepth = 12,
        refresh = 0
      )
    fit$save_object(file = filePath)
    print(paste0('File saved as: ignore/models/',taxaName,'_lnbin.rds'))
    return(NULL)
  } else{
    print(paste("Model ",taxaName," exists. To overwrite, set `rerun` = TRUE and `overwrite` = TRUE"))
    return(NULL)
  }
}
#'
#'
check_divergences = function(filePath = NULL){
  fit = readRDS(filePath)
  diverge = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent) > 0
  num_diverge = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  return(list(diverge = diverge,
              num_diverge = num_diverge))
}

make_site_index = function(df = NULL){
  taxaName = unique(df$acceptedTaxonID)
  siteYearDf = df %>% 
    named_group_split(siteID, collectYear) %>% 
    map(~.x %>% 
          pmap(~rep(x = ..4, times = ..6)) %>% 
          list %>% 
          unlist %>% 
          as_tibble) %>% 
    bind_rows(.id = 'id')%>% 
    tidyr::separate_wider_delim(id, names = c('siteID','collectYear'), delim = "/", cols_remove = FALSE)
  n_per_sample = siteYearDf %>% 
    summarise(count = n(),.by = 'id') %>% 
    select(count) %>% 
    unlist %>% unname
  start_idx = c(1,(cumsum(n_per_sample)+1)) %>% head(.,-1)
  site_name = as.factor(siteYearDf$siteID)[start_idx]
  site_id = as.character(as.integer(as.factor(siteYearDf$siteID))[start_idx])
  
  return(data.frame(
    acceptedTaxonID = taxaName,
    siteName = site_name,
    siteID = site_id
  ) %>% filter(!duplicated(.))
  )
}

extract_max_nb = function(filePath = NULL, names = NULL, siteIndex = NULL){
  mod = readRDS(filePath)
  vars = paste(c('max_ref_rep', 'size_mean','size_median'), collapse = '|')
  fullSumm = mod$summary()
  filteredSumm = fullSumm[grepl(vars, fullSumm$variable), c('variable','mean', 'median', 'sd', 'q5','q95')]
  filteredSumm$acceptedTaxonID = names
  filteredSumm$siteID = gsub("\\w+\\[(\\d{1,2})\\]", "\\1", filteredSumm$variable)
  summDf = left_join(filteredSumm, data.frame(siteIndex), by = c('acceptedTaxonID','siteID'))
  return(summDf)
}

extract_max_softmax = function(filePath = NULL, name = NULL, siteIndex = NULL){
  mod = readRDS(filePath)
  vars = paste(c('size_mean', 'size_sd','size_q50', 'size_q75','size_q95', 'size_q99'), collapse = '|')
  fullSumm = mod$summary()
  filteredSumm = fullSumm[grepl(vars, fullSumm$variable), c('variable','mean', 'median', 'sd', 'q5','q95')]
  filteredSumm$acceptedTaxonID <- siteIndex$acceptedTaxonID <- name
  filteredSumm$site_id = gsub("\\w+\\[(\\d{1,2})\\]", "\\1", filteredSumm$variable)
  summDf = left_join(filteredSumm, data.frame(siteIndex), by = c('acceptedTaxonID','site_id'))
  return(summDf)
}
###### SPARED(D) CODE ########
# Convert one taxon's event-by-size-bin data into the compressed data
# structure required by hierarchical_nb_binned_lower_trunc.stan.
#
# Expected default input columns:
#   siteID           site identifier
#   collectYear      sampling-event identifier within site
#   acceptedTaxonID  taxon identifier
#   sizeClass        midpoint of a 1-mm size bin
#   no_m2            integer abundance/count in that event and bin
#
# The column `dw` may be present but is not used by this size model.
#
# IMPORTANT:
# - Each function call must contain exactly one taxon.
# - The negative-binomial and multinomial likelihoods require integer
#   values in bin_count_col.
# - When multiple collections can occur within a site-year, pass all
#   columns required to identify an event through event_cols, for example
#   event_cols = c("collectYear", "collectDate", "sampleID").
# - Supply event_data to include sampled events at which the taxon had
#   zero abundance. Without event_data, only positive-occurrence events
#   can be reconstructed.
make_binned_nb_stan_data <- function(
    bin_data,
    event_data = NULL,
    site_col = "siteID",
    event_cols = "collectYear",
    taxon_col = "acceptedTaxonID",
    bin_mid_col = "sizeClass",
    bin_count_col = "no_m2",
    bin_width = 1,
    size_lower = 0.5,
    sigma_floor = NULL,
    k_ref = 20L,
    prior_only = 0L,
    integer_tolerance = 1e-8
) {
  stopifnot(
    is.data.frame(bin_data),
    length(site_col) == 1L,
    length(taxon_col) == 1L,
    length(bin_mid_col) == 1L,
    length(bin_count_col) == 1L,
    length(event_cols) >= 1L,
    length(bin_width) == 1L,
    is.finite(bin_width),
    bin_width > 0,
    length(size_lower) == 1L,
    is.finite(size_lower),
    size_lower >= 0,
    length(k_ref) == 1L,
    is.finite(k_ref),
    k_ref >= 1L,
    prior_only %in% c(0L, 1L),
    integer_tolerance >= 0
  )
  
  if (is.null(sigma_floor)) {
    # Numerical floor equal to 1% of the bin width by default.
    # For 1-mm bins this is 0.01 mm: far below the resolution of
    # the observations, but large enough to prevent exp(log_sigma)
    # from rounding to exactly zero.
    sigma_floor <- bin_width / 100
  }
  
  if (
    length(sigma_floor) != 1L ||
    !is.finite(sigma_floor) ||
    sigma_floor <= 0
  ) {
    stop("sigma_floor must be a finite positive scalar.")
  }
  
  required_bin_cols <- unique(
    c(
      site_col,
      event_cols,
      taxon_col,
      bin_mid_col,
      bin_count_col
    )
  )
  
  missing_bin_cols <- setdiff(
    required_bin_cols,
    names(bin_data)
  )
  
  if (length(missing_bin_cols)) {
    stop(
      "bin_data is missing required columns: ",
      paste(missing_bin_cols, collapse = ", ")
    )
  }
  
  taxon_values <- unique(
    as.character(bin_data[[taxon_col]])
  )
  
  taxon_values <- taxon_values[
    !is.na(taxon_values)
  ]
  
  if (length(taxon_values) != 1L) {
    stop(
      "bin_data must contain exactly one nonmissing taxon. Found: ",
      paste(taxon_values, collapse = ", ")
    )
  }
  
  taxon_id <- taxon_values[1L]
  
  make_key <- function(data, columns) {
    key_parts <- lapply(
      data[columns],
      function(x) {
        x <- as.character(x)
        
        if (anyNA(x)) {
          stop(
            "Missing values are not allowed in key columns: ",
            paste(columns, collapse = ", ")
          )
        }
        
        x
      }
    )
    
    do.call(
      paste,
      c(
        key_parts,
        sep = "\r"
      )
    )
  }
  
  bins <- data.frame(
    site_value = as.character(bin_data[[site_col]]),
    bin_mid = as.numeric(bin_data[[bin_mid_col]]),
    bin_count_raw = as.numeric(bin_data[[bin_count_col]]),
    stringsAsFactors = FALSE
  )
  
  for (column in event_cols) {
    bins[[column]] <- bin_data[[column]]
  }
  
  if (
    anyNA(bins$site_value) ||
    anyNA(bins$bin_mid) ||
    anyNA(bins$bin_count_raw)
  ) {
    stop(
      "bin_data contains missing site, size midpoint, or bin count values."
    )
  }
  
  if (
    any(!is.finite(bins$bin_mid)) ||
    any(!is.finite(bins$bin_count_raw)) ||
    any(bins$bin_count_raw < 0)
  ) {
    stop(
      "Size midpoints must be finite and bin counts must be finite ",
      "and nonnegative."
    )
  }
  
  noninteger <- abs(
    bins$bin_count_raw -
      round(bins$bin_count_raw)
  ) > integer_tolerance
  
  if (any(noninteger)) {
    bad <- which(noninteger)[1L]
    
    stop(
      "`", bin_count_col, "` must contain integer counts for the ",
      "negative-binomial and multinomial likelihoods. First ",
      "noninteger value: ", bins$bin_count_raw[bad]
    )
  }
  
  bins$bin_count <- as.integer(
    round(bins$bin_count_raw)
  )
  
  aggregate_columns <- c(
    "site_value",
    event_cols,
    "bin_mid"
  )
  
  bins <- stats::aggregate(
    bins["bin_count"],
    by = bins[aggregate_columns],
    FUN = sum
  )
  
  bins <- bins[
    bins$bin_count > 0L,
    ,
    drop = FALSE
  ]
  
  if (!nrow(bins)) {
    stop("No positive event-by-bin counts remain after aggregation.")
  }
  
  bins$bin_lower <-
    bins$bin_mid -
    0.5 * bin_width
  
  bins$bin_upper <-
    bins$bin_mid +
    0.5 * bin_width
  
  tolerance <- sqrt(.Machine$double.eps)
  
  below_lower <- bins$bin_lower <
    size_lower - tolerance
  
  if (any(below_lower)) {
    bad <- which(below_lower)[1L]
    
    stop(
      "An occupied bin extends below size_lower. Midpoint = ",
      bins$bin_mid[bad],
      ", lower edge = ",
      bins$bin_lower[bad],
      ", size_lower = ",
      size_lower
    )
  }
  
  bins$bin_lower <- pmax(
    bins$bin_lower,
    size_lower
  )
  
  # Build the event roster. A complete roster is required to represent
  # genuine zero-count events for this taxon.
  if (is.null(event_data)) {
    warning(
      "event_data was not supplied. Events are inferred only from ",
      "positive taxon observations, so true zero-count events are ",
      "not included in the count model."
    )
    
    events <- unique(
      bins[c("site_value", event_cols)]
    )
  } else {
    if (!is.data.frame(event_data)) {
      stop("event_data must be NULL or a data.frame.")
    }
    
    required_event_cols <- unique(
      c(
        site_col,
        event_cols
      )
    )
    
    missing_event_cols <- setdiff(
      required_event_cols,
      names(event_data)
    )
    
    if (length(missing_event_cols)) {
      stop(
        "event_data is missing required columns: ",
        paste(missing_event_cols, collapse = ", ")
      )
    }
    
    event_source <- event_data
    
    # If event_data is taxon-specific and contains the taxon column,
    # retain the current taxon. A general event roster should omit the
    # taxon column and will be used without filtering.
    if (taxon_col %in% names(event_source)) {
      event_source <- event_source[
        as.character(event_source[[taxon_col]]) == taxon_id,
        ,
        drop = FALSE
      ]
    }
    
    events <- data.frame(
      site_value = as.character(event_source[[site_col]]),
      stringsAsFactors = FALSE
    )
    
    for (column in event_cols) {
      events[[column]] <- event_source[[column]]
    }
    
    events <- unique(events)
    
    if (!nrow(events)) {
      stop(
        "No sampling events remain in event_data for taxon ",
        taxon_id,
        "."
      )
    }
    
    if (anyNA(events$site_value)) {
      stop("event_data contains missing site identifiers.")
    }
  }
  
  event_key_columns <- c(
    "site_value",
    event_cols
  )
  
  events$event_key <- make_key(
    events,
    event_key_columns
  )
  
  bins$event_key <- make_key(
    bins,
    event_key_columns
  )
  
  missing_event_keys <- setdiff(
    unique(bins$event_key),
    events$event_key
  )
  
  if (length(missing_event_keys)) {
    stop(
      "Some positive event-by-bin records are absent from event_data."
    )
  }
  
  site_levels <- unique(
    events$site_value
  )
  
  events$site_id <- match(
    events$site_value,
    site_levels
  )
  
  events$event_id <- seq_len(
    nrow(events)
  )
  
  bins$event_id <- events$event_id[
    match(
      bins$event_key,
      events$event_key
    )
  ]
  
  if (anyNA(bins$event_id)) {
    stop("Failed to map one or more occupied bins to an event.")
  }
  
  event_totals <- integer(
    nrow(events)
  )
  
  summed_counts <- tapply(
    bins$bin_count,
    bins$event_id,
    sum
  )
  
  event_totals[
    as.integer(names(summed_counts))
  ] <- as.integer(summed_counts)
  
  bins <- bins[
    order(
      bins$event_id,
      bins$bin_mid
    ),
    ,
    drop = FALSE
  ]
  
  site_total_count <- numeric(
    length(site_levels)
  )
  
  site_event_count <- integer(
    length(site_levels)
  )
  
  for (site in seq_along(site_levels)) {
    use <- events$site_id == site
    
    site_total_count[site] <-
      sum(event_totals[use])
    
    site_event_count[site] <-
      sum(use)
  }
  
  max_bin_edge_by_site <- rep(
    NA_real_,
    length(site_levels)
  )
  
  names(max_bin_edge_by_site) <- site_levels
  
  observed_sites <- unique(
    bins$site_value
  )
  
  for (site_name in observed_sites) {
    max_bin_edge_by_site[site_name] <- max(
      bins$bin_upper[
        bins$site_value == site_name
      ]
    )
  }
  
  stan_data <- list(
    S = as.integer(
      length(site_levels)
    ),
    K = as.integer(
      nrow(events)
    ),
    B_obs = as.integer(
      nrow(bins)
    ),
    n_per_sample = as.integer(
      event_totals
    ),
    site_id = as.integer(
      events$site_id
    ),
    bin_event = as.integer(
      bins$event_id
    ),
    bin_count = as.integer(
      bins$bin_count
    ),
    bin_lower = as.numeric(
      bins$bin_lower
    ),
    bin_upper = as.numeric(
      bins$bin_upper
    ),
    size_lower = as.numeric(
      size_lower
    ),
    sigma_floor = as.numeric(
      sigma_floor
    ),
    k_ref = as.integer(
      k_ref
    ),
    prior_only = as.integer(
      prior_only
    )
  )
  
  site_mapping <- data.frame(
    site_id = seq_along(site_levels),
    siteID = site_levels,
    n_events = site_event_count,
    total_count = site_total_count,
    has_size_data = site_total_count > 0,
    max_observed_bin_edge = unname(
      max_bin_edge_by_site[site_levels]
    ),
    stringsAsFactors = FALSE
  )
  
  event_mapping <- data.frame(
    event_id = events$event_id,
    site_id = events$site_id,
    stringsAsFactors = FALSE
  )
  
  event_mapping[[site_col]] <-
    events$site_value
  
  for (column in event_cols) {
    event_mapping[[column]] <-
      events[[column]]
  }
  
  event_mapping$n_per_sample <-
    event_totals
  
  bin_mapping <- data.frame(
    event_id = bins$event_id,
    site_id = events$site_id[
      bins$event_id
    ],
    stringsAsFactors = FALSE
  )
  
  bin_mapping[[site_col]] <-
    bins$site_value
  
  for (column in event_cols) {
    bin_mapping[[column]] <-
      bins[[column]]
  }
  
  bin_mapping[[taxon_col]] <-
    taxon_id
  
  bin_mapping[[bin_mid_col]] <-
    bins$bin_mid
  
  bin_mapping$bin_lower <-
    bins$bin_lower
  
  bin_mapping$bin_upper <-
    bins$bin_upper
  
  bin_mapping$bin_count <-
    bins$bin_count
  
  metadata <- list(
    acceptedTaxonID = taxon_id,
    total_sites = length(site_levels),
    total_events = nrow(events),
    positive_events = sum(event_totals > 0L),
    total_count = sum(event_totals),
    bin_width = bin_width,
    size_lower = size_lower,
    sigma_floor = sigma_floor,
    k_ref = as.integer(k_ref)
  )
  
  list(
    stan_data = stan_data,
    mapping = list(
      site = site_mapping,
      event = event_mapping,
      bins = bin_mapping
    ),
    metadata = metadata
  )
}

# Generate stable, chain-specific initial values for
# hierarchical_nb_binned_lower_trunc.stan.
#
# The size hierarchy is noncentered:
#
#   log_mu_site[s] =
#     alpha_log_mu + tau_log_mu * z_mu[s]
#
#   log_sigma_site[s] =
#     alpha_log_sigma + tau_log_sigma * z_sigma[s]
#
# Count means are direct site-level parameters and are not pooled.
make_init_binned_nb <- function(
    stan_data,
    chains = 4L,
    seed = 1234L
) {
  required_names <- c(
    "S",
    "K",
    "B_obs",
    "n_per_sample",
    "site_id",
    "bin_event",
    "bin_count",
    "bin_lower",
    "bin_upper",
    "size_lower",
    "sigma_floor"
  )
  
  missing_names <- setdiff(
    required_names,
    names(stan_data)
  )
  
  if (length(missing_names)) {
    stop(
      "stan_data is missing: ",
      paste(missing_names, collapse = ", ")
    )
  }
  
  stopifnot(
    is.list(stan_data),
    length(chains) == 1L,
    is.finite(chains),
    chains >= 1L,
    length(seed) == 1L,
    is.finite(seed)
  )
  
  S <- as.integer(stan_data$S)
  K <- as.integer(stan_data$K)
  B_obs <- as.integer(stan_data$B_obs)
  chains <- as.integer(chains)
  seed <- as.integer(seed)
  
  if (
    S < 1L ||
    K < 1L ||
    B_obs < 1L
  ) {
    stop("S, K, and B_obs must be positive.")
  }
  
  if (S == 1L) {
    warning(
      "S == 1. The initializer will return values for the ",
      "hierarchical model, but a dedicated single-site model is ",
      "usually preferable because the among-site scales cannot be ",
      "learned from one site."
    )
  }
  
  if (
    length(stan_data$n_per_sample) != K ||
    length(stan_data$site_id) != K
  ) {
    stop("n_per_sample and site_id must each have length K.")
  }
  
  if (
    length(stan_data$bin_event) != B_obs ||
    length(stan_data$bin_count) != B_obs ||
    length(stan_data$bin_lower) != B_obs ||
    length(stan_data$bin_upper) != B_obs
  ) {
    stop(
      "bin_event, bin_count, bin_lower, and bin_upper must ",
      "each have length B_obs."
    )
  }
  
  if (
    any(stan_data$site_id < 1L) ||
    any(stan_data$site_id > S)
  ) {
    stop("site_id values must lie between 1 and S.")
  }
  
  if (
    any(stan_data$bin_event < 1L) ||
    any(stan_data$bin_event > K)
  ) {
    stop("bin_event values must lie between 1 and K.")
  }
  
  if (
    any(stan_data$n_per_sample < 0L) ||
    any(stan_data$bin_count <= 0L)
  ) {
    stop(
      "n_per_sample must be nonnegative and bin_count must be positive."
    )
  }
  
  if (
    any(!is.finite(stan_data$bin_lower)) ||
    any(!is.finite(stan_data$bin_upper)) ||
    any(stan_data$bin_upper <= stan_data$bin_lower)
  ) {
    stop(
      "All bin edges must be finite and bin_upper must exceed bin_lower."
    )
  }
  
  if (
    !is.finite(stan_data$size_lower) ||
    stan_data$size_lower < 0 ||
    any(stan_data$bin_lower < stan_data$size_lower)
  ) {
    stop(
      "size_lower must be nonnegative and occupied bins may not ",
      "extend below it."
    )
  }
  
  clamp <- function(x, x_lower, x_upper) {
    pmin(
      pmax(x, x_lower),
      x_upper
    )
  }
  
  # Confirm that the compressed bins reproduce each event total.
  reconstructed <- numeric(K)
  
  for (b in seq_len(B_obs)) {
    event <- stan_data$bin_event[b]
    
    reconstructed[event] <-
      reconstructed[event] +
      stan_data$bin_count[b]
  }
  
  if (
    any(
      reconstructed !=
      as.numeric(stan_data$n_per_sample)
    )
  ) {
    bad <- which(
      reconstructed !=
        as.numeric(stan_data$n_per_sample)
    )[1L]
    
    stop(
      "Binned counts do not reproduce n_per_sample at event ",
      bad,
      ". Reconstructed = ",
      reconstructed[bad],
      "; n_per_sample = ",
      stan_data$n_per_sample[bad]
    )
  }
  
  bin_site <- as.integer(
    stan_data$site_id[
      stan_data$bin_event
    ]
  )
  
  bin_mid <-
    0.5 * (
      stan_data$bin_lower +
        stan_data$bin_upper
    )
  
  bin_width <-
    stan_data$bin_upper -
    stan_data$bin_lower
  
  typical_bin_width <- stats::median(
    bin_width
  )
  
  if (
    !is.finite(typical_bin_width) ||
    typical_bin_width <= 0
  ) {
    typical_bin_width <- 1
  }
  
  sigma_model_floor <-
    as.numeric(stan_data$sigma_floor)
  
  if (
    !is.finite(sigma_model_floor) ||
    sigma_model_floor <= 0
  ) {
    stop("stan_data$sigma_floor must be finite and positive.")
  }
  
  # The optimizer estimates the actual sigma, whose Stan
  # parameterization is:
  #
  # sigma = sigma_model_floor + exp(log_sigma_site).
  sigma_init_floor <-
    sigma_model_floor +
    max(
      1e-4,
      typical_bin_width / 1000
    )
  
  mu_guard <- 400
  sigma_guard <- 200
  
  weighted_bin_moments <- function(
    midpoint,
    width,
    count
  ) {
    total <- sum(count)
    
    if (
      !is.finite(total) ||
      total <= 0
    ) {
      return(
        c(
          mean = NA_real_,
          sd = NA_real_
        )
      )
    }
    
    weighted_mean <-
      sum(count * midpoint) /
      total
    
    weighted_variance <-
      sum(
        count *
          (
            (midpoint - weighted_mean)^2 +
              width^2 / 12
          )
      ) /
      total
    
    c(
      mean = weighted_mean,
      sd = sqrt(
        max(
          weighted_variance,
          0
        )
      )
    )
  }
  
  global_moments <- weighted_bin_moments(
    midpoint = bin_mid,
    width = bin_width,
    count = stan_data$bin_count
  )
  
  global_mean <- global_moments["mean"]
  global_sd <- global_moments["sd"]
  
  if (
    !is.finite(global_mean) ||
    global_mean <= 0
  ) {
    global_mean <- max(
      stan_data$size_lower + 0.5 * typical_bin_width,
      1
    )
  }
  
  if (
    !is.finite(global_sd) ||
    global_sd <= 0
  ) {
    global_sd <- max(
      sigma_init_floor,
      0.25 * global_mean
    )
  }
  
  global_mean <- clamp(
    global_mean,
    1e-4,
    mu_guard
  )
  
  global_sd <- clamp(
    global_sd,
    sigma_init_floor,
    sigma_guard
  )
  
  # Log probability for a normal interval in R. The CDF difference is
  # used first, with a survival-function fallback for the upper tail.
  normal_interval_log_prob_r <- function(
    x_lower,
    x_upper,
    mu,
    sigma
  ) {
    cdf_difference <-
      stats::pnorm(
        x_upper,
        mean = mu,
        sd = sigma
      ) -
      stats::pnorm(
        x_lower,
        mean = mu,
        sd = sigma
      )
    
    bad <- !is.finite(cdf_difference) |
      cdf_difference <= 0
    
    if (any(bad)) {
      cdf_difference[bad] <-
        stats::pnorm(
          x_lower[bad],
          mean = mu,
          sd = sigma,
          lower.tail = FALSE
        ) -
        stats::pnorm(
          x_upper[bad],
          mean = mu,
          sd = sigma,
          lower.tail = FALSE
        )
    }
    
    log(
      pmax(
        cdf_difference,
        .Machine$double.xmin
      )
    )
  }
  
  fit_site_size <- function(site) {
    use <- which(
      bin_site == site
    )
    
    # Sites with no observed individuals are initialized from the pooled
    # taxon moments and will be informed by the hierarchy in Stan.
    if (!length(use)) {
      return(
        c(
          mu = global_mean,
          sigma = global_sd
        )
      )
    }
    
    site_moments <- weighted_bin_moments(
      midpoint = bin_mid[use],
      width = bin_width[use],
      count = stan_data$bin_count[use]
    )
    
    raw_mu <- site_moments["mean"]
    raw_sigma <- site_moments["sd"]
    
    if (
      !is.finite(raw_mu) ||
      raw_mu <= 0
    ) {
      raw_mu <- global_mean
    }
    
    if (
      !is.finite(raw_sigma) ||
      raw_sigma <= 0
    ) {
      raw_sigma <- global_sd
    }
    
    raw_mu <- clamp(
      raw_mu,
      1e-4,
      mu_guard
    )
    
    raw_sigma <- clamp(
      raw_sigma,
      sigma_init_floor,
      sigma_guard
    )
    
    distinct_bins <- nrow(
      unique(
        data.frame(
          x_lower = stan_data$bin_lower[use],
          x_upper = stan_data$bin_upper[use]
        )
      )
    )
    
    if (distinct_bins < 2L) {
      return(
        c(
          mu = raw_mu,
          sigma = max(
            raw_sigma,
            sigma_init_floor,
            0.5 * global_sd
          )
        )
      )
    }
    
    site_max_edge <- max(
      stan_data$bin_upper[use]
    )
    
    local_mu_upper <- min(
      mu_guard,
      max(
        10,
        4 * site_max_edge,
        4 * raw_mu
      )
    )
    
    local_sigma_upper <- min(
      sigma_guard,
      max(
        5,
        2 * site_max_edge,
        4 * raw_sigma
      )
    )
    
    negative_log_likelihood <- function(par) {
      mu <- exp(par[1])
      sigma <- exp(par[2])
      
      log_normalizer <- stats::pnorm(
        stan_data$size_lower,
        mean = mu,
        sd = sigma,
        lower.tail = FALSE,
        log.p = TRUE
      )
      
      log_bin_probability <-
        normal_interval_log_prob_r(
          x_lower = stan_data$bin_lower[use],
          x_upper = stan_data$bin_upper[use],
          mu = mu,
          sigma = sigma
        ) -
        log_normalizer
      
      if (
        !is.finite(log_normalizer) ||
        any(!is.finite(log_bin_probability))
      ) {
        return(
          .Machine$double.xmax^0.25
        )
      }
      
      log_likelihood <- sum(
        stan_data$bin_count[use] *
          log_bin_probability
      )
      
      if (!is.finite(log_likelihood)) {
        return(
          .Machine$double.xmax^0.25
        )
      }
      
      -log_likelihood
    }
    
    fit <- tryCatch(
      stats::optim(
        par = log(
          c(
            raw_mu,
            raw_sigma
          )
        ),
        fn = negative_log_likelihood,
        method = "L-BFGS-B",
        lower = log(
          c(
            1e-4,
            sigma_init_floor
          )
        ),
        upper = log(
          c(
            local_mu_upper,
            local_sigma_upper
          )
        ),
        control = list(
          maxit = 500,
          factr = 1e8
        )
      ),
      error = function(e) NULL
    )
    
    if (
      is.null(fit) ||
      !is.finite(fit$value) ||
      any(!is.finite(fit$par))
    ) {
      return(
        c(
          mu = raw_mu,
          sigma = raw_sigma
        )
      )
    }
    
    estimates <- exp(
      fit$par
    )
    
    c(
      mu = clamp(
        estimates[1],
        1e-4,
        mu_guard
      ),
      sigma = clamp(
        estimates[2],
        sigma_init_floor,
        sigma_guard
      )
    )
  }
  
  site_size_estimates <- t(
    vapply(
      seq_len(S),
      fit_site_size,
      numeric(2)
    )
  )
  
  log_mu_target <- log(
    site_size_estimates[, "mu.mean"]
  )
  
  log_sigma_target <- log(
    pmax(
      site_size_estimates[, "sigma.sd"] -
        sigma_model_floor,
      1e-4
    )
  )
  
  robust_scale <- function(
    x,
    center,
    minimum,
    maximum,
    max_z = 2.5
  ) {
    if (length(x) <= 1L) {
      return(minimum)
    }
    
    mad_scale <- stats::mad(
      x,
      center = center,
      constant = 1.4826
    )
    
    if (!is.finite(mad_scale)) {
      mad_scale <- 0
    }
    
    range_scale <- diff(
      range(x)
    ) / 4
    
    if (!is.finite(range_scale)) {
      range_scale <- 0
    }
    
    deviation_scale <- max(
      abs(x - center)
    ) / max_z
    
    min(
      max(
        minimum,
        mad_scale,
        range_scale,
        deviation_scale
      ),
      maximum
    )
  }
  
  alpha_mu_init <- stats::median(
    log_mu_target
  )
  
  alpha_sigma_init <- stats::median(
    log_sigma_target
  )
  
  tau_mu_init <- robust_scale(
    log_mu_target,
    center = alpha_mu_init,
    minimum = 0.10,
    maximum = 2.0
  )
  
  tau_sigma_init <- robust_scale(
    log_sigma_target,
    center = alpha_sigma_init,
    minimum = 0.10,
    maximum = 2.0
  )
  
  # Direct, unpooled site-level expected event counts.
  lambda_site <- vapply(
    seq_len(S),
    function(site) {
      counts <- stan_data$n_per_sample[
        stan_data$site_id == site
      ]
      
      if (!length(counts)) {
        return(0.1)
      }
      
      max(
        mean(counts),
        0.1
      )
    },
    numeric(1)
  )
  
  log_lambda_target <- log(
    lambda_site
  )
  
  # Site-level method-of-moments estimates for shared NB2 dispersion.
  phi_by_site <- vapply(
    seq_len(S),
    function(site) {
      counts <- stan_data$n_per_sample[
        stan_data$site_id == site
      ]
      
      if (length(counts) < 2L) {
        return(NA_real_)
      }
      
      mean_count <- mean(counts)
      variance_count <- stats::var(counts)
      
      if (
        !is.finite(mean_count) ||
        !is.finite(variance_count) ||
        mean_count <= 0 ||
        variance_count <= mean_count
      ) {
        return(NA_real_)
      }
      
      mean_count^2 /
        (
          variance_count -
            mean_count
        )
    },
    numeric(1)
  )
  
  finite_phi <- phi_by_site[
    is.finite(phi_by_site) &
      phi_by_site > 0
  ]
  
  phi_init <- if (length(finite_phi)) {
    stats::median(
      finite_phi
    )
  } else {
    20
  }
  
  log_phi_init <- clamp(
    log(phi_init),
    -3.95,
    11.95
  )
  
  set.seed(seed)
  
  lapply(
    seq_len(chains),
    function(chain_id) {
      target_log_mu_chain <-
        log_mu_target +
        stats::rnorm(
          S,
          mean = 0,
          sd = 0.01
        )
      
      target_log_sigma_chain <-
        log_sigma_target +
        stats::rnorm(
          S,
          mean = 0,
          sd = 0.01
        )
      
      alpha_log_mu_chain <-
        stats::median(
          target_log_mu_chain
        ) +
        stats::rnorm(
          1,
          mean = 0,
          sd = 0.015
        )
      
      alpha_log_sigma_chain <-
        stats::median(
          target_log_sigma_chain
        ) +
        stats::rnorm(
          1,
          mean = 0,
          sd = 0.015
        )
      
      tau_log_mu_chain <- clamp(
        tau_mu_init *
          exp(
            stats::rnorm(
              1,
              mean = 0,
              sd = 0.02
            )
          ),
        0.02,
        2.40
      )
      
      tau_log_sigma_chain <- clamp(
        tau_sigma_init *
          exp(
            stats::rnorm(
              1,
              mean = 0,
              sd = 0.02
            )
          ),
        0.02,
        2.40
      )
      
      z_mu_chain <-
        (
          target_log_mu_chain -
            alpha_log_mu_chain
        ) /
        tau_log_mu_chain
      
      z_sigma_chain <-
        (
          target_log_sigma_chain -
            alpha_log_sigma_chain
        ) /
        tau_log_sigma_chain
      
      log_lambda_chain <-
        log_lambda_target +
        stats::rnorm(
          S,
          mean = 0,
          sd = 0.02
        )
      
      log_phi_chain <- clamp(
        log_phi_init +
          stats::rnorm(
            1,
            mean = 0,
            sd = 0.03
          ),
        -3.95,
        11.95
      )
      
      list(
        alpha_log_mu = as.numeric(
          alpha_log_mu_chain
        ),
        tau_log_mu = as.numeric(
          tau_log_mu_chain
        ),
        z_mu = as.numeric(
          z_mu_chain
        ),
        alpha_log_sigma = as.numeric(
          alpha_log_sigma_chain
        ),
        tau_log_sigma = as.numeric(
          tau_log_sigma_chain
        ),
        z_sigma = as.numeric(
          z_sigma_chain
        ),
        log_lambda_site = as.numeric(
          log_lambda_chain
        ),
        log_phi = as.numeric(
          log_phi_chain
        )
      )
    }
  )
}

#'
#'
#'

midpoint_resample = function(x){
  round(runif(n = 1, min = x-0.5, max = x+0.5),3)
}
midpoint_resample_vec = Vectorize(midpoint_resample)


#'
#'
#'
make_stanData_taxa = function(df = NULL){
  taxaName = df$acceptedTaxonID
  
  siteYearDf = df %>% 
    named_group_split(siteID, collectYear) %>% 
    map(~.x %>% 
          pmap(~rep(x = ..4, times = ..6)) %>% 
          list %>% 
          unlist %>% 
          midpoint_resample_vec %>% 
          as_tibble) %>% 
    bind_rows(.id = 'id')%>% 
    tidyr::separate_wider_delim(id, names = c('siteID','collectYear'), delim = "/", cols_remove = FALSE)
  
  S = as.integer(count(unique(siteYearDf$siteID)))
  K = as.integer(count(unique(siteYearDf$id)))
  n_obs = nrow(siteYearDf)
  x = unlist(siteYearDf$value)
  n_per_sample = siteYearDf %>% 
    summarise(count = n(),.by = 'id') %>% 
    select(count) %>% 
    unlist %>% unname
  start_idx = c(1,(cumsum(n_per_sample)+1)) %>% head(.,-1)
  site_id = as.integer(as.factor(siteYearDf$siteID))[start_idx]
  k_ref = 20L
  
  stan_data = list(
    S = S,
    K = K,
    n_obs = n_obs,
    x = x,
    n_per_sample = n_per_sample,
    start_idx = start_idx,
    site_id = site_id,
    k_ref = k_ref,
    # set upper mu limit to generous but reasonable
    mu_upper = as.double(400),
    # set upper sigma limit to generous but reasonable
    sigma_upper = as.double(200),
    prior_only = 0
  )
  
  # stan_data$mu_upper <- max(
  #   200,
  #   10 * max(stan_data$x)
  # )
  # 
  # stan_data$sigma_upper <- max(
  #   200,
  #   10 * stats::sd(stan_data$x)
  # )
  
  return(stan_data = stan_data)
}

# set initial values for each model to speed warmup
#'
#'
#'
#'
make_init_list_stable <- function(
    stan_data,
    chains = 4L,
    seed = 1234L,
    boundary_margin = 0.02,
    log_phi_bounds = c(-4, 12)
) {
  stopifnot(
    is.list(stan_data),
    length(stan_data$x) == stan_data$n_obs,
    length(stan_data$n_per_sample) == stan_data$K,
    length(stan_data$site_id) == stan_data$K,
    sum(stan_data$n_per_sample) == stan_data$n_obs,
    all(is.finite(stan_data$x)),
    all(stan_data$x >= 0),
    all(stan_data$n_per_sample >= 0),
    all(stan_data$site_id %in% seq_len(stan_data$S)),
    is.finite(stan_data$mu_upper),
    is.finite(stan_data$sigma_upper),
    stan_data$mu_upper > 0,
    stan_data$sigma_upper > 0,
    boundary_margin > 0,
    boundary_margin < 0.25
  )
  
  S <- as.integer(stan_data$S)
  
  # Stay away from the exact parameter boundaries.
  mu_lower <- 1e-6
  sigma_lower <- 1e-6
  
  mu_upper_init <-
    (1 - boundary_margin) * stan_data$mu_upper
  
  sigma_upper_init <-
    (1 - boundary_margin) * stan_data$sigma_upper
  
  obs_site <- rep(
    stan_data$site_id,
    times = stan_data$n_per_sample
  )
  
  global_mean_x <- mean(stan_data$x)
  global_sd_x <- stats::sd(stan_data$x)
  
  if (!is.finite(global_mean_x) || global_mean_x <= 0) {
    global_mean_x <- min(1, 0.25 * stan_data$mu_upper)
  }
  
  if (!is.finite(global_sd_x) || global_sd_x <= 0) {
    global_sd_x <- max(0.1, 0.25 * global_mean_x)
  }
  
  global_mean_x <- min(
    max(global_mean_x, mu_lower),
    mu_upper_init
  )
  
  global_sd_x <- min(
    max(global_sd_x, sigma_lower),
    sigma_upper_init
  )
  
  # #
  #   * Site-specific MLE for:
  #   *
  #   * X ~ Normal(mu, sigma), conditional on X > 0.
  # *
  #   * Optimization occurs in log(mu), log(sigma), matching
  # * the parameterization used in Stan.
  # #
  fit_site_tnorm <- function(xs) {
    xs <- xs[is.finite(xs) & xs >= 0]
    
    if (length(xs) < 2L) {
      return(
        c(
          mu = global_mean_x,
          sigma = global_sd_x
        )
      )
    }
    
    raw_mu <- mean(xs)
    raw_sigma <- stats::sd(xs)
    
    if (!is.finite(raw_mu) || raw_mu <= 0) {
      raw_mu <- global_mean_x
    }
    
    if (!is.finite(raw_sigma) || raw_sigma <= 0) {
      raw_sigma <- global_sd_x
    }
    
    raw_mu <- min(
      max(raw_mu, mu_lower),
      mu_upper_init
    )
    
    raw_sigma <- min(
      max(raw_sigma, sigma_lower),
      sigma_upper_init
    )
    
    negative_log_likelihood <- function(par) {
      mu <- exp(par[1])
      sigma <- exp(par[2])
      
      if (
        !is.finite(mu) ||
        !is.finite(sigma) ||
        sigma <= 0
      ) {
        return(.Machine$double.xmax^0.25)
      }
      
      log_normalizer <- stats::pnorm(
        0,
        mean = mu,
        sd = sigma,
        lower.tail = FALSE,
        log.p = TRUE
      )
      
      log_likelihood <-
        sum(
          stats::dnorm(
            xs,
            mean = mu,
            sd = sigma,
            log = TRUE
          )
        ) -
        length(xs) * log_normalizer
      
      if (!is.finite(log_likelihood)) {
        return(.Machine$double.xmax^0.25)
      }
      
      -log_likelihood
    }
    
    fit <- tryCatch(
      stats::optim(
        par = log(c(raw_mu, raw_sigma)),
        fn = negative_log_likelihood,
        method = "L-BFGS-B",
        lower = log(c(mu_lower, sigma_lower)),
        upper = log(
          c(mu_upper_init, sigma_upper_init)
        ),
        control = list(
          maxit = 500,
          factr = 1e8
        )
      ),
      error = function(e) NULL
    )
    
    if (
      is.null(fit) ||
      !is.finite(fit$value) ||
      any(!is.finite(fit$par))
    ) {
      return(
        c(
          mu = raw_mu,
          sigma = raw_sigma
        )
      )
    }
    
    estimates <- exp(fit$par)
    
    c(
      mu = min(
        max(estimates[1], mu_lower),
        mu_upper_init
      ),
      sigma = min(
        max(estimates[2], sigma_lower),
        sigma_upper_init
      )
    )
  }
  
  site_size_estimates <- t(
    vapply(
      seq_len(S),
      function(s) {
        fit_site_tnorm(
          stan_data$x[obs_site == s]
        )
      },
      numeric(2)
    )
  )
  
  site_mu <- site_size_estimates[, "mu"]
  site_sigma <- site_size_estimates[, "sigma"]
  
  # Mean count per standardized event at each site.
  site_lambda <- vapply(
    seq_len(S),
    function(s) {
      ns <- stan_data$n_per_sample[
        stan_data$site_id == s
      ]
      
      if (!length(ns)) {
        return(
          max(mean(stan_data$n_per_sample), 0.1)
        )
      }
      
      max(mean(ns), 0.1)
    },
    numeric(1)
  )
  
  log_mu_init <- log(site_mu)
  log_sigma_init <- log(site_sigma)
  log_lambda_init <- log(site_lambda)
  
  # #
  #   * Robust center. This prevents one unusually large site
  # * from determining the initial population-level center.
  # #
  alpha_mu_init <- stats::median(log_mu_init)
  alpha_sigma_init <- stats::median(log_sigma_init)
  alpha_lambda_init <- stats::median(log_lambda_init)
  
  robust_scale <- function(
    x,
    center,
    minimum,
    maximum,
    max_z = 2.5
  ) {
    mad_scale <- stats::mad(
      x,
      center = center,
      constant = 1.4826
    )
    
    if (!is.finite(mad_scale)) {
      mad_scale <- 0
    }
    
    range_scale <- diff(range(x)) / 4
    
    if (!is.finite(range_scale)) {
      range_scale <- 0
    }
    
    # #
    #   * Ensure no initial group effect is excessively far from
    # * the hierarchical center.
    # #
    deviation_scale <-
      max(abs(x - center)) / max_z
    
    value <- max(
      minimum,
      mad_scale,
      range_scale,
      deviation_scale
    )
    
    min(value, maximum)
  }
  
  tau_mu_init <- robust_scale(
    log_mu_init,
    center = alpha_mu_init,
    minimum = 0.10,
    maximum = 2.0
  )
  
  tau_sigma_init <- robust_scale(
    log_sigma_init,
    center = alpha_sigma_init,
    minimum = 0.10,
    maximum = 2.0
  )
  
  tau_lambda_init <- robust_scale(
    log_lambda_init,
    center = alpha_lambda_init,
    minimum = 0.30,
    maximum = 4.0
  )
  
  z_lambda_init <-
    (log_lambda_init - alpha_lambda_init) /
    tau_lambda_init
  
  # #
  #   * Method-of-moments initial phi values:
  #   *
  #   * Var(N) = mean(N) + mean(N)^2 / phi
  # *
  #   * so phi = mean(N)^2 / [Var(N) - mean(N)].
  # #
  phi_by_site <- vapply(
    seq_len(S),
    function(s) {
      ns <- stan_data$n_per_sample[
        stan_data$site_id == s
      ]
      
      if (length(ns) < 2L) {
        return(NA_real_)
      }
      
      mean_n <- mean(ns)
      var_n <- stats::var(ns)
      
      if (
        !is.finite(mean_n) ||
        !is.finite(var_n) ||
        mean_n <= 0 ||
        var_n <= mean_n
      ) {
        return(NA_real_)
      }
      
      mean_n^2 / (var_n - mean_n)
    },
    numeric(1)
  )
  
  finite_phi <- phi_by_site[
    is.finite(phi_by_site) &
      phi_by_site > 0
  ]
  
  phi_init <- if (length(finite_phi)) {
    stats::median(finite_phi)
  } else {
    # Moderately weak overdispersion if data cannot identify it.
    50
  }
  
  phi_init <- min(
    max(phi_init, exp(log_phi_bounds[1] + 0.1)),
    exp(log_phi_bounds[2] - 0.1)
  )
  
  clamp <- function(x, lower, upper) {
    pmin(pmax(x, lower), upper)
  }
  
  log_mu_lower <- log(mu_lower)
  log_mu_upper <- log(mu_upper_init)
  
  log_sigma_lower <- log(sigma_lower)
  log_sigma_upper <- log(sigma_upper_init)
  
  log_phi_lower <- log_phi_bounds[1] + 0.05
  log_phi_upper <- log_phi_bounds[2] - 0.05
  
  set.seed(seed)
  
  lapply(
    seq_len(chains),
    function(chain_id) {
      # Small chain-specific perturbations
      log_mu_chain <- clamp(
        log_mu_init +
          stats::rnorm(S, 0, 0.01),
        log_mu_lower,
        log_mu_upper
      )
      
      log_sigma_chain <- clamp(
        log_sigma_init +
          stats::rnorm(S, 0, 0.01),
        log_sigma_lower,
        log_sigma_upper
      )
      
      list(
        alpha_log_mu =
          alpha_mu_init +
          stats::rnorm(1, 0, 0.02),
        
        alpha_log_sigma =
          alpha_sigma_init +
          stats::rnorm(1, 0, 0.02),
        
        log_mu_site =
          log_mu_chain,
        
        log_sigma_site =
          log_sigma_chain,
        
        tau_log_mu =
          tau_mu_init *
          exp(stats::rnorm(1, 0, 0.02)),
        
        tau_log_sigma =
          tau_sigma_init *
          exp(stats::rnorm(1, 0, 0.02)),
        
        alpha_log_lambda =
          alpha_lambda_init +
          stats::rnorm(1, 0, 0.03),
        
        tau_log_lambda =
          tau_lambda_init *
          exp(stats::rnorm(1, 0, 0.02)),
        
        z_lambda =
          z_lambda_init +
          stats::rnorm(S, 0, 0.02),
        
        log_phi =
          clamp(
            log(phi_init) +
              stats::rnorm(1, 0, 0.03),
            log_phi_lower,
            log_phi_upper
          )
      )
    }
  )
}

make_init_list_single_site <- function(
    stan_data,
    chains = 4L,
    seed = 1234L,
    boundary_margin = 0.02,
    log_phi_bounds = c(-4, 12),
    sigma_relative_floor = 0.01
) {
  stopifnot(
    is.list(stan_data),
    stan_data$S == 1L,
    length(stan_data$x) == stan_data$n_obs,
    length(stan_data$n_per_sample) == stan_data$K,
    length(stan_data$site_id) == stan_data$K,
    all(stan_data$site_id == 1L),
    sum(stan_data$n_per_sample) == stan_data$n_obs,
    all(is.finite(stan_data$x)),
    all(stan_data$x >= 0),
    all(stan_data$n_per_sample >= 0),
    is.finite(stan_data$mu_upper),
    is.finite(stan_data$sigma_upper),
    stan_data$mu_upper > 0,
    stan_data$sigma_upper > 0,
    boundary_margin > 0,
    boundary_margin < 0.25,
    length(log_phi_bounds) == 2L,
    log_phi_bounds[1] < log_phi_bounds[2]
  )
  
  x <- stan_data$x
  counts <- stan_data$n_per_sample
  
  clamp <- function(x, lower, upper) {
    pmin(
      pmax(x, lower),
      upper
    )
  }
  
  
  # Stay comfortably inside the bounds used by Stan.
  
  mu_lower_init <- 1e-6
  
  mu_upper_init <-
    (1 - boundary_margin) *
    stan_data$mu_upper
  
  sigma_upper_init <-
    (1 - boundary_margin) *
    stan_data$sigma_upper
  
  observed_mean <- mean(x)
  observed_sd <- stats::sd(x)
  
  if (
    !is.finite(observed_mean) ||
    observed_mean <= 0
  ) {
    observed_mean <- min(
      1,
      0.25 * stan_data$mu_upper
    )
  }
  
  #  A small data-scale-dependent floor prevents a very narrow
  #  sample from initializing sigma almost exactly at zero.
  #  This affects initialization only, not the posterior support.
  
  sigma_lower_init <- max(
    1e-6,
    sigma_relative_floor * observed_mean
  )
  
  sigma_lower_init <- min(
    sigma_lower_init,
    0.25 * sigma_upper_init
  )
  
  if (
    !is.finite(observed_sd) ||
    observed_sd <= 0
  ) {
    observed_sd <- max(
      sigma_lower_init,
      0.10 * observed_mean
    )
  }
  
  raw_mu <- clamp(
    observed_mean,
    mu_lower_init,
    mu_upper_init
  )
  
  raw_sigma <- clamp(
    observed_sd,
    sigma_lower_init,
    sigma_upper_init
  )
  
  
  # Fit the same positive-truncated normal used in Stan.
  
  # This provides initialization for the latent normal mu and
  # sigma, which are not generally identical to the observed
  # mean and standard deviation after truncation.
  
  negative_log_likelihood <- function(par) {
    mu <- exp(par[1])
    sigma <- exp(par[2])
    
    if (
      !is.finite(mu) ||
      !is.finite(sigma) ||
      sigma <= 0
    ) {
      return(.Machine$double.xmax^0.25)
    }
    
    log_normalizer <- stats::pnorm(
      0,
      mean = mu,
      sd = sigma,
      lower.tail = FALSE,
      log.p = TRUE
    )
    
    log_likelihood <-
      sum(
        stats::dnorm(
          x,
          mean = mu,
          sd = sigma,
          log = TRUE
        )
      ) -
      length(x) * log_normalizer
    
    if (!is.finite(log_likelihood)) {
      return(.Machine$double.xmax^0.25)
    }
    
    -log_likelihood
  }
  
  
  # With fewer than two observations, the truncated-normal
  # scale cannot be estimated empirically.
  
  size_fit <- if (length(x) >= 2L) {
    tryCatch(
      stats::optim(
        par = log(
          c(
            raw_mu,
            raw_sigma
          )
        ),
        fn = negative_log_likelihood,
        method = "L-BFGS-B",
        lower = log(
          c(
            mu_lower_init,
            sigma_lower_init
          )
        ),
        upper = log(
          c(
            mu_upper_init,
            sigma_upper_init
          )
        ),
        control = list(
          maxit = 500,
          factr = 1e8
        )
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  
  if (
    is.null(size_fit) ||
    !is.finite(size_fit$value) ||
    any(!is.finite(size_fit$par))
  ) {
    mu_init <- raw_mu
    sigma_init <- raw_sigma
  } else {
    estimates <- exp(size_fit$par)
    
    mu_init <- clamp(
      estimates[1],
      mu_lower_init,
      mu_upper_init
    )
    
    sigma_init <- clamp(
      estimates[2],
      sigma_lower_init,
      sigma_upper_init
    )
  }
  
  
  # Expected count per equal-effort sampling event.
  
  lambda_init <- mean(counts)
  
  if (
    !is.finite(lambda_init) ||
    lambda_init <= 0
  ) {
    lambda_init <- 0.1
  }
  
  log_lambda_init <-
    log(max(lambda_init, 0.1))
  
  
  # Method-of-moments initialization for NB2 dispersion:
  # Var(N) = lambda + lambda^2 / phi
  # phi = lambda^2 / [Var(N) - lambda].
  
  # A single event cannot estimate phi, and underdispersion
  # relative to Poisson does not yield a finite NB2 estimate.
  
  if (length(counts) >= 2L) {
    mean_count <- mean(counts)
    variance_count <- stats::var(counts)
    
    if (
      is.finite(mean_count) &&
      is.finite(variance_count) &&
      mean_count > 0 &&
      variance_count > mean_count
    ) {
      phi_init <-
        mean_count^2 /
        (variance_count - mean_count)
    } else {
      phi_init <- 20
    }
  } else {
    phi_init <- 20
  }
  
  log_phi_lower <-
    log_phi_bounds[1] +
    0.05
  
  log_phi_upper <-
    log_phi_bounds[2] -
    0.05
  
  log_phi_init <- clamp(
    log(phi_init),
    log_phi_lower,
    log_phi_upper
  )
  
  log_mu_lower <-
    log(mu_lower_init)
  
  log_mu_upper <-
    log(mu_upper_init)
  
  log_sigma_lower <-
    log(sigma_lower_init)
  
  log_sigma_upper <-
    log(sigma_upper_init)
  
  set.seed(seed)
  
  lapply(
    seq_len(chains),
    function(chain_id) {
      
      # The Stan parameters are vectors of length one, so
      # log_mu_site, log_sigma_site, and log_lambda must each
      # be supplied as length-one numeric vectors.
      
      list(
        log_mu_site = c(
          clamp(
            log(mu_init) +
              stats::rnorm(1, 0, 0.01),
            log_mu_lower,
            log_mu_upper
          )
        ),
        
        log_sigma_site = c(
          clamp(
            log(sigma_init) +
              stats::rnorm(1, 0, 0.01),
            log_sigma_lower,
            log_sigma_upper
          )
        ),
        
        log_lambda = c(
          log_lambda_init +
            stats::rnorm(1, 0, 0.02)
        ),
        
        log_phi = clamp(
          log_phi_init +
            stats::rnorm(1, 0, 0.03),
          log_phi_lower,
          log_phi_upper
        )
      )
    }
  )
}
#'
#'
#'
make_init_list <- function(
    stan_data,
    chains = 4L,
    seed = 1234L
) {
  stopifnot(
    is.list(stan_data),
    length(stan_data$x) == stan_data$n_obs,
    length(stan_data$n_per_sample) == stan_data$K,
    length(stan_data$site_id) == stan_data$K,
    sum(stan_data$n_per_sample) == stan_data$n_obs
  )
  
  S <- stan_data$S
  
  # Site identity for each individual size observation
  obs_site <- rep(
    stan_data$site_id,
    times = stan_data$n_per_sample
  )
  
  # Global fallbacks
  global_mean_x <- mean(stan_data$x)
  global_sd_x <- stats::sd(stan_data$x)
  
  if (!is.finite(global_sd_x) || global_sd_x <= 0) {
    global_sd_x <- max(global_mean_x * 0.25, 0.1)
  }
  
  global_mean_n <- mean(stan_data$n_per_sample)
  
  # Empirical size summaries by site
  site_mu <- vapply(
    seq_len(S),
    function(s) {
      xs <- stan_data$x[obs_site == s]
      
      if (length(xs) > 0L) {
        max(mean(xs), 1e-3)
      } else {
        max(global_mean_x, 1e-3)
      }
    },
    numeric(1)
  )
  
  site_sigma <- vapply(
    seq_len(S),
    function(s) {
      xs <- stan_data$x[obs_site == s]
      
      if (length(xs) >= 2L) {
        sx <- stats::sd(xs)
        
        if (is.finite(sx) && sx > 0) {
          return(max(sx, 0.05))
        }
      }
      
      max(global_sd_x, 0.05)
    },
    numeric(1)
  )
  
  # Mean equal-effort count per event at each site
  # Zero-count events are retained here.
  site_lambda <- vapply(
    seq_len(S),
    function(s) {
      ns <- stan_data$n_per_sample[
        stan_data$site_id == s
      ]
      
      if (length(ns) > 0L) {
        max(mean(ns), 0.1)
      } else {
        max(global_mean_n, 0.1)
      }
    },
    numeric(1)
  )
  
  log_site_mu <- log(site_mu)
  log_site_sigma <- log(site_sigma)
  log_site_lambda <- log(site_lambda)
  
  # Across-site initial centers
  alpha_mu_init <- mean(log_site_mu)
  alpha_sigma_init <- mean(log_site_sigma)
  alpha_lambda_init <- mean(log_site_lambda)
  
  safe_sd <- function(x, fallback) {
    sx <- stats::sd(x)
    
    if (!is.finite(sx) || sx <= 0) {
      fallback
    } else {
      sx
    }
  }
  
  # Initial among-site heterogeneity
  tau_mu_init <- min(
    max(safe_sd(log_site_mu, 0.20), 0.10),
    1.50
  )
  
  tau_sigma_init <- min(
    max(safe_sd(log_site_sigma, 0.20), 0.10),
    1.50
  )
  
  # Broader because counts may differ by hundreds-fold among sites
  tau_lambda_init <- min(
    max(safe_sd(log_site_lambda, 0.75), 0.30),
    2.50
  )
  
  # Standardized non-centered site effects
  z_mu_init <-
    (log_site_mu - alpha_mu_init) /
    tau_mu_init
  
  z_sigma_init <-
    (log_site_sigma - alpha_sigma_init) /
    tau_sigma_init
  
  z_lambda_init <-
    (log_site_lambda - alpha_lambda_init) /
    tau_lambda_init
  
  set.seed(seed)
  
  lapply(
    seq_len(chains),
    function(chain_id) {
      list(
        alpha_log_mu =
          alpha_mu_init +
          stats::rnorm(1, 0, 0.01),
        
        alpha_log_sigma =
          alpha_sigma_init +
          stats::rnorm(1, 0, 0.01),
        
        alpha_log_lambda =
          alpha_lambda_init +
          stats::rnorm(1, 0, 0.02),
        
        tau_log_mu =
          tau_mu_init,
        
        tau_log_sigma =
          tau_sigma_init,
        
        tau_log_lambda =
          tau_lambda_init,
        
        z_mu =
          z_mu_init +
          stats::rnorm(S, 0, 0.01),
        
        z_sigma =
          z_sigma_init +
          stats::rnorm(S, 0, 0.01),
        
        z_lambda =
          z_lambda_init +
          stats::rnorm(S, 0, 0.01),
        
        log_phi =
          log(20) +
          stats::rnorm(1, 0, 0.02)
      )
    }
  )
}

#'
#'
#'
fit_negbin_named = function(stan_data = NULL, taxaName = NULL, rerun = FALSE, overwrite = FALSE){
  filePath = paste0(here("ignore/models"),"/",taxaName,"_negbin.rds")
  print(taxaName)
  if(any(rerun, !file.exists(filePath))){
    if(all(file.exists(filePath),!overwrite)){
      warning('Model file already exists and `overwrite` = FALSE. Set to TRUE to overwrite existing files.')
      return(NULL)
    }
    chains = 4L
    
    # add a modifier for only single site data
    if(stan_data$S == 1L){
      
      init_list = make_init_list_single_site(
        stan_data = stan_data,
        chains = chains,
        seed = 1312
      )
      fit = single_mod$sample(
        data = stan_data,
        seed = 1312,
        chains = chains,
        parallel_chains = chains,
        
        init = init_list,
        
        iter_warmup = 1500,
        iter_sampling = 1000,
        adapt_delta = 0.99,
        max_treedepth = 12,
        refresh = 0
      )
    } else{
      
      init_list = make_init_list_stable(
        stan_data = stan_data,
        chains = chains,
        seed = 1312
      )
      
      fit = mod$sample(
        data = stan_data,
        seed = 1312,
        chains = chains,
        parallel_chains = chains,
        
        init = init_list,
        
        iter_warmup = 1500,
        iter_sampling = 1000,
        adapt_delta = 0.99,
        max_treedepth = 12,
        refresh = 0
      )
    }
    fit$save_object(file = filePath)
    print(paste0('File saved as: ignore/models/',taxaName,'_negbin.rds'))
    return(NULL)
  } else{
    print(paste("Model ",taxaName," exists. To overwrite, set `rerun` = TRUE and `overwrite` = TRUE"))
    return(NULL)
  }
}

# fit_negbin_named = purrr::safely(fit_negbin_named)