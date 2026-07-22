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
  
  # /*
  #   * Site-specific MLE for:
  #   *
  #   * X ~ Normal(mu, sigma), conditional on X > 0.
  # *
  #   * Optimization occurs in log(mu), log(sigma), matching
  # * the parameterization used in Stan.
  # */
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
  
  # /*
  #   * Robust center. This prevents one unusually large site
  # * from determining the initial population-level center.
  # */
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
    
    # /*
    #   * Ensure no initial group effect is excessively far from
    # * the hierarchical center.
    # */
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
  
  # /*
  #   * Method-of-moments initial phi values:
  #   *
  #   * Var(N) = mean(N) + mean(N)^2 / phi
  # *
  #   * so phi = mean(N)^2 / [Var(N) - mean(N)].
  # */
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
check_divergences = function(filePath = NULL){
  fit = readRDS(filePath)
  diverge = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent) > 0
  num_diverge = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  return(list(diverge = diverge,
              num_diverge = num_diverge))
}

###### SPARED(D) CODE ########
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