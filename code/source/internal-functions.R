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