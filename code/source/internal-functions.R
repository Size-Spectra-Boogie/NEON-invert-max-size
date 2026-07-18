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
    mu_upper = 200,
    # set upper sigma limit to generous but reasonable
    sigma_upper = 100,
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
    adapt_delta = 0.95,
    max_treedepth = 12,
    refresh = 0
  )
  fit$save_object(file = filePath)
  print(paste0('File saved as: ignore/models/',taxaName,'_negbin.rds'))
  return(NULL)
  } else{
    print(paste("Model ",taxaName," exists. To overwrite, set `rerun` = TRUE and `overwrite` = TRUE"))
  return(NULL)
    }
}

# fit_negbin_named = purrr::safely(fit_negbin_named)

# set initial values across all models

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

make_init_list_stable <- function(
    stan_data,
    chains = 4L,
    seed = 1234L
) {
  stopifnot(
    is.list(stan_data),
    length(stan_data$x) == stan_data$n_obs,
    length(stan_data$n_per_sample) == stan_data$K,
    length(stan_data$site_id) == stan_data$K,
    sum(stan_data$n_per_sample) == stan_data$n_obs,
    is.finite(stan_data$mu_upper),
    is.finite(stan_data$sigma_upper)
  )
  
  S <- stan_data$S
  
  obs_site <- rep(
    stan_data$site_id,
    times = stan_data$n_per_sample
  )
  
  global_mean_x <- mean(stan_data$x)
  global_sd_x <- stats::sd(stan_data$x)
  
  if (!is.finite(global_sd_x) || global_sd_x <= 0) {
    global_sd_x <- max(0.1, 0.25 * global_mean_x)
  }
  
  site_mu <- vapply(
    seq_len(S),
    function(s) {
      xs <- stan_data$x[obs_site == s]
      
      value <- if (length(xs)) {
        mean(xs)
      } else {
        global_mean_x
      }
      
      min(
        max(value, 1e-6),
        0.95 * stan_data$mu_upper
      )
    },
    numeric(1)
  )
  
  site_sigma <- vapply(
    seq_len(S),
    function(s) {
      xs <- stan_data$x[obs_site == s]
      
      value <- if (length(xs) >= 2L) {
        stats::sd(xs)
      } else {
        global_sd_x
      }
      
      if (!is.finite(value) || value <= 0) {
        value <- global_sd_x
      }
      
      min(
        max(value, 1e-4),
        0.95 * stan_data$sigma_upper
      )
    },
    numeric(1)
  )
  
  site_lambda <- vapply(
    seq_len(S),
    function(s) {
      ns <- stan_data$n_per_sample[
        stan_data$site_id == s
      ]
      
      if (length(ns)) {
        max(mean(ns), 0.1)
      } else {
        max(mean(stan_data$n_per_sample), 0.1)
      }
    },
    numeric(1)
  )
  
  log_mu_init <- log(site_mu)
  log_sigma_init <- log(site_sigma)
  log_lambda_init <- log(site_lambda)
  
  safe_sd <- function(x, fallback) {
    value <- stats::sd(x)
    
    if (!is.finite(value) || value <= 0) {
      fallback
    } else {
      value
    }
  }
  
  alpha_mu_init <- mean(log_mu_init)
  alpha_sigma_init <- mean(log_sigma_init)
  alpha_lambda_init <- mean(log_lambda_init)
  
  tau_mu_init <- min(
    max(safe_sd(log_mu_init, 0.2), 0.05),
    2
  )
  
  tau_sigma_init <- min(
    max(safe_sd(log_sigma_init, 0.2), 0.05),
    2
  )
  
  tau_lambda_init <- min(
    max(safe_sd(log_lambda_init, 0.75), 0.2),
    3
  )
  
  z_lambda_init <-
    (log_lambda_init - alpha_lambda_init) /
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
        
        log_mu_site =
          log_mu_init +
          stats::rnorm(S, 0, 0.005),
        
        log_sigma_site =
          log_sigma_init +
          stats::rnorm(S, 0, 0.005),
        
        tau_log_mu = tau_mu_init,
        tau_log_sigma = tau_sigma_init,
        
        alpha_log_lambda =
          alpha_lambda_init +
          stats::rnorm(1, 0, 0.02),
        
        tau_log_lambda = tau_lambda_init,
        
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
