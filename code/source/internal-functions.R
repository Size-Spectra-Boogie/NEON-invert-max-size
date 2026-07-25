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

make_binned_nb_stan_data <- function(
    bin_data,
    event_data = NULL,
    site_col = "siteID",
    event_col = "collectYear",
    taxon_col = "acceptedTaxonID",
    bin_mid_col = "sizeClass",
    bin_count_col = "no_m2",
    bin_width = 1,
    size_lower = 0.5,
    upper_multiplier_max = 1.3,
    k_ref = 20L,
    mu_upper = 400,
    sigma_upper = 200,
    prior_only = 0L,
    integer_tolerance = 1e-8,
    allow_zero_size_sites = FALSE
) {
  stopifnot(
    is.data.frame(bin_data),
    is.numeric(bin_width),
    length(bin_width) == 1L,
    is.finite(bin_width),
    bin_width > 0,
    is.numeric(size_lower),
    length(size_lower) == 1L,
    is.finite(size_lower),
    size_lower >= 0,
    is.numeric(upper_multiplier_max),
    length(upper_multiplier_max) == 1L,
    is.finite(upper_multiplier_max),
    upper_multiplier_max > 1,
    is.numeric(mu_upper),
    is.finite(mu_upper),
    mu_upper > size_lower,
    is.numeric(sigma_upper),
    is.finite(sigma_upper),
    sigma_upper > 0,
    k_ref >= 1L,
    prior_only %in% c(0L, 1L)
  )
  
  required_bin_columns <- c(
    site_col,
    event_col,
    taxon_col,
    bin_mid_col,
    bin_count_col
  )
  
  missing_bin_columns <- setdiff(
    required_bin_columns,
    names(bin_data)
  )
  
  if (length(missing_bin_columns)) {
    stop(
      "bin_data is missing required columns: ",
      paste(missing_bin_columns, collapse = ", ")
    )
  }
  
  /*
    * Each call should contain one taxon.
  */
    taxon_values <- unique(
      as.character(bin_data[[taxon_col]])
    )
  
  taxon_values <- taxon_values[
    !is.na(taxon_values)
  ]
  
  if (length(taxon_values) != 1L) {
    stop(
      "bin_data must contain exactly one taxon. Found: ",
      paste(taxon_values, collapse = ", ")
    )
  }
  
  taxon_id <- taxon_values[1L]
  
  bins <- data.frame(
    site = as.character(bin_data[[site_col]]),
    event = as.character(bin_data[[event_col]]),
    taxon = as.character(bin_data[[taxon_col]]),
    bin_mid = as.numeric(bin_data[[bin_mid_col]]),
    bin_count_raw = as.numeric(bin_data[[bin_count_col]]),
    stringsAsFactors = FALSE
  )
  
  if (
    anyNA(bins$site) ||
    anyNA(bins$event) ||
    anyNA(bins$taxon) ||
    anyNA(bins$bin_mid) ||
    anyNA(bins$bin_count_raw)
  ) {
    stop(
      "bin_data contains missing site, event, taxon, ",
      "sizeClass, or no_m2 values."
    )
  }
  
  if (any(!is.finite(bins$bin_mid))) {
    stop("All sizeClass values must be finite.")
  }
  
  if (
    any(!is.finite(bins$bin_count_raw)) ||
    any(bins$bin_count_raw < 0)
  ) {
    stop("All no_m2 values must be finite and nonnegative.")
  }
  
  /*
    * Negative-binomial and multinomial likelihoods require
  * integer-valued counts.
  */
    noninteger <- abs(
      bins$bin_count_raw -
        round(bins$bin_count_raw)
    ) > integer_tolerance
  
  if (any(noninteger)) {
    bad <- which(noninteger)[1L]
    
    stop(
      "`", bin_count_col, "` contains non-integer values. ",
      "The current Stan model uses negative-binomial and ",
      "multinomial count likelihoods and therefore requires ",
      "integer counts. First non-integer value: ",
      bins$bin_count_raw[bad],
      ". Use original integer counts or revise the likelihood."
    )
  }
  
  bins$bin_count <- as.integer(
    round(bins$bin_count_raw)
  )
  
  /*
    * Aggregate duplicate site-event-bin rows.
  */
    bins <- stats::aggregate(
      bin_count ~ site + event + taxon + bin_mid,
      data = bins,
      FUN = sum
    )
  
  /*
    * Stan only needs occupied event-bin cells.
  */
    bins <- bins[
      bins$bin_count > 0L,
      ,
      drop = FALSE
    ]
  
  if (!nrow(bins)) {
    stop("No positive bin counts remain after aggregation.")
  }
  
  bins$bin_lower <-
    bins$bin_mid -
    0.5 * bin_width
  
  bins$bin_upper <-
    bins$bin_mid +
    0.5 * bin_width
  
  tolerance <- sqrt(.Machine$double.eps)
  
  if (
    any(
      bins$bin_lower <
      size_lower - tolerance
    )
  ) {
    bad <- which(
      bins$bin_lower <
        size_lower - tolerance
    )[1L]
    
    stop(
      "An occupied bin extends below size_lower. ",
      "Taxon = ", taxon_id,
      "; site = ", bins$site[bad],
      "; event = ", bins$event[bad],
      "; midpoint = ", bins$bin_mid[bad],
      "; lower edge = ", bins$bin_lower[bad],
      "; size_lower = ", size_lower
    )
  }
  
  /*
    * Remove negligible floating-point differences at the
  * lower boundary.
  */
    bins$bin_lower <-
    pmax(
      bins$bin_lower,
      size_lower
    )
  
  /*
    * Construct the complete event table.
  *
    * event_data should contain every sampled site-event
  * combination, including events in which this taxon was absent.
  */
    if (is.null(event_data)) {
      warning(
        "event_data was not supplied. Sampling events are inferred ",
        "only from rows where the taxon was observed. True zero-count ",
        "events will therefore be omitted from the negative-binomial ",
        "count model."
      )
      
      events <- unique(
        bins[c("site", "event")]
      )
    } else {
      required_event_columns <- c(
        site_col,
        event_col
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
      
      /*
        * If event_data is taxon-specific, retain only the current taxon.
      * If it is a general sampling-event roster without a taxon column,
      * all listed site-event combinations are retained.
      */
        if (taxon_col %in% names(event_source)) {
          event_source <- event_source[
            as.character(event_source[[taxon_col]]) ==
              taxon_id,
            ,
            drop = FALSE
          ]
        }
      
      events <- unique(
        data.frame(
          site = as.character(
            event_source[[site_col]]
          ),
          event = as.character(
            event_source[[event_col]]
          ),
          stringsAsFactors = FALSE
        )
      )
      
      if (!nrow(events)) {
        stop(
          "No sampling events remain in event_data for taxon ",
          taxon_id,
          "."
        )
      }
      
      if (
        anyNA(events$site) ||
        anyNA(events$event)
      ) {
        stop(
          "event_data contains missing site or event identifiers."
        )
      }
    }
  
  /*
    * Site and event identifiers may contain repeated numeric or
  * character values, so use a composite key.
  */
    event_key <- function(site, event) {
      paste(
        site,
        event,
        sep = "\r"
      )
    }
  
  events$key <- event_key(
    events$site,
    events$event
  )
  
  bins$key <- event_key(
    bins$site,
    bins$event
  )
  
  missing_events <- setdiff(
    unique(bins$key),
    events$key
  )
  
  if (length(missing_events)) {
    stop(
      "Some positive bin records do not have matching rows ",
      "in event_data."
    )
  }
  
  /*
    * Retain the order in which sites first appear in the event table.
  */
    site_levels <- unique(
      events$site
    )
  
  events$site_id <- match(
    events$site,
    site_levels
  )
  
  events$event_id <-
    seq_len(nrow(events))
  
  bins$event_id <- events$event_id[
    match(
      bins$key,
      events$key
    )
  ]
  
  if (anyNA(bins$event_id)) {
    stop(
      "Failed to match one or more bin records to sampling events."
    )
  }
  
  /*
    * Event totals are the sums of the size-bin counts.
  *
    * Events absent from bins retain a total of zero.
  */
    event_totals <-
    integer(nrow(events))
  
  summed_counts <- tapply(
    bins$bin_count,
    bins$event_id,
    sum
  )
  
  event_totals[
    as.integer(names(summed_counts))
  ] <- as.integer(summed_counts)
  
  /*
    * Sort compressed rows for easier inspection.
  * Stan does not require this ordering.
  */
    bins <- bins[
      order(
        bins$event_id,
        bins$bin_mid
      ),
      ,
      drop = FALSE
    ]
  
  /*
    * Number of observed individuals at each site.
  */
    site_total_counts <- tapply(
      event_totals,
      events$site_id,
      sum
    )
  
  site_total_counts_complete <-
    integer(length(site_levels))
  
  site_total_counts_complete[
    as.integer(names(site_total_counts))
  ] <- as.integer(site_total_counts)
  
  zero_size_sites <- which(
    site_total_counts_complete == 0L
  )
  
  if (
    length(zero_size_sites) &&
    !allow_zero_size_sites
  ) {
    stop(
      "The following sites have no measured individuals for taxon ",
      taxon_id,
      ": ",
      paste(
        site_levels[zero_size_sites],
        collapse = ", "
      ),
      ". Remove those sites, or set allow_zero_size_sites = TRUE ",
      "after removing the corresponding zero-size-site rejection ",
      "from the Stan transformed-data block."
    )
  }
  
  /*
    * Taxon-level endpoint diagnostics.
  */
    max_bin_edge_taxon <-
    max(bins$bin_upper)
  
  size_upper_taxon_cap <-
    upper_multiplier_max *
    max_bin_edge_taxon
  
  /*
    * Site-level observed maximum-bin edges are diagnostic only.
  * The Stan likelihood uses one shared taxon-level endpoint.
  */
    max_bin_edge_by_site <- tapply(
      bins$bin_upper,
      bins$site,
      max
    )
  
  site_mapping <- data.frame(
    site_id = seq_along(site_levels),
    siteID = site_levels,
    total_count = site_total_counts_complete,
    has_size_data =
      site_total_counts_complete > 0L,
    max_observed_bin_edge =
      as.numeric(
        max_bin_edge_by_site[
          match(
            site_levels,
            names(max_bin_edge_by_site)
          )
        ]
      ),
    stringsAsFactors = FALSE
  )
  
  event_mapping <- data.frame(
    event_id = events$event_id,
    site_id = events$site_id,
    siteID = events$site,
    collectYear = events$event,
    n_per_sample = event_totals,
    stringsAsFactors = FALSE
  )
  
  bin_mapping <- data.frame(
    event_id = bins$event_id,
    siteID = bins$site,
    collectYear = bins$event,
    acceptedTaxonID = bins$taxon,
    sizeClass = bins$bin_mid,
    bin_lower = bins$bin_lower,
    bin_upper = bins$bin_upper,
    bin_count = bins$bin_count,
    stringsAsFactors = FALSE
  )
  
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
    
    upper_multiplier_max = as.numeric(
      upper_multiplier_max
    ),
    
    k_ref = as.integer(
      k_ref
    ),
    
    mu_upper = as.numeric(
      mu_upper
    ),
    
    sigma_upper = as.numeric(
      sigma_upper
    ),
    
    prior_only = as.integer(
      prior_only
    )
  )
  
  metadata <- list(
    acceptedTaxonID =
      taxon_id,
    
    max_bin_edge_taxon =
      max_bin_edge_taxon,
    
    upper_multiplier_max =
      upper_multiplier_max,
    
    size_upper_taxon_cap =
      size_upper_taxon_cap,
    
    observed_total =
      sum(event_totals),
    
    observed_events =
      sum(event_totals > 0L),
    
    total_events =
      length(event_totals),
    
    total_sites =
      length(site_levels)
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

#'
#'
#'
# Create stable, chain-specific initial values for the hierarchical
# negative-binomial + binned doubly truncated-normal Stan model.
#
# The function expects the stan_data object returned by
# make_binned_nb_stan_data().
make_init_list_binned_nb <- function(
    stan_data,
    chains = 4L,
    seed = 1234L,
    boundary_margin = 0.02,
    log_phi_bounds = c(-4, 12),
    upper_fraction_center = 0.5
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
    "upper_multiplier_max",
    "mu_upper",
    "sigma_upper"
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
    is.finite(seed),
    boundary_margin > 0,
    boundary_margin < 0.25,
    length(log_phi_bounds) == 2L,
    all(is.finite(log_phi_bounds)),
    log_phi_bounds[1] < log_phi_bounds[2],
    upper_fraction_center > 0,
    upper_fraction_center < 1
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
    stop("S, K, and B_obs must all be positive.")
  }
  
  if (
    length(stan_data$n_per_sample) != K ||
    length(stan_data$site_id) != K
  ) {
    stop(
      "n_per_sample and site_id must each have length K."
    )
  }
  
  if (
    length(stan_data$bin_event) != B_obs ||
    length(stan_data$bin_count) != B_obs ||
    length(stan_data$bin_lower) != B_obs ||
    length(stan_data$bin_upper) != B_obs
  ) {
    stop(
      "bin_event, bin_count, bin_lower, and bin_upper ",
      "must each have length B_obs."
    )
  }
  
  if (
    any(!is.finite(stan_data$n_per_sample)) ||
    any(stan_data$n_per_sample < 0) ||
    any(stan_data$n_per_sample != round(stan_data$n_per_sample))
  ) {
    stop("n_per_sample must contain nonnegative integer counts.")
  }
  
  if (
    any(!is.finite(stan_data$site_id)) ||
    any(stan_data$site_id < 1L) ||
    any(stan_data$site_id > S) ||
    any(stan_data$site_id != round(stan_data$site_id))
  ) {
    stop("site_id must contain integers between 1 and S.")
  }
  
  if (
    any(!is.finite(stan_data$bin_event)) ||
    any(stan_data$bin_event < 1L) ||
    any(stan_data$bin_event > K) ||
    any(stan_data$bin_event != round(stan_data$bin_event))
  ) {
    stop("bin_event must contain integers between 1 and K.")
  }
  
  if (
    any(!is.finite(stan_data$bin_count)) ||
    any(stan_data$bin_count <= 0) ||
    any(stan_data$bin_count != round(stan_data$bin_count))
  ) {
    stop("bin_count must contain positive integer counts.")
  }
  
  if (
    any(!is.finite(stan_data$bin_lower)) ||
    any(!is.finite(stan_data$bin_upper)) ||
    any(stan_data$bin_upper <= stan_data$bin_lower)
  ) {
    stop(
      "All bin limits must be finite and bin_upper must exceed ",
      "bin_lower."
    )
  }
  
  if (
    !is.finite(stan_data$size_lower) ||
    stan_data$size_lower < 0 ||
    any(stan_data$bin_lower < stan_data$size_lower)
  ) {
    stop(
      "size_lower must be nonnegative and no occupied bin may ",
      "extend below it."
    )
  }
  
  if (
    !is.finite(stan_data$upper_multiplier_max) ||
    stan_data$upper_multiplier_max <= 1
  ) {
    stop("upper_multiplier_max must be greater than one.")
  }
  
  if (
    !is.finite(stan_data$mu_upper) ||
    !is.finite(stan_data$sigma_upper) ||
    stan_data$mu_upper <= 0 ||
    stan_data$sigma_upper <= 0
  ) {
    stop("mu_upper and sigma_upper must be finite and positive.")
  }
  
  clamp <- function(x, lower, upper) {
    pmin(
      pmax(x, lower),
      upper
    )
  }
  
  # Confirm that compressed bin counts reproduce event totals.
  reconstructed_counts <- numeric(K)
  
  for (b in seq_len(B_obs)) {
    j <- stan_data$bin_event[b]
    
    reconstructed_counts[j] <-
      reconstructed_counts[j] +
      stan_data$bin_count[b]
  }
  
  if (
    any(
      reconstructed_counts !=
      as.numeric(stan_data$n_per_sample)
    )
  ) {
    bad <- which(
      reconstructed_counts !=
        as.numeric(stan_data$n_per_sample)
    )[1L]
    
    stop(
      "Compressed bin counts do not reproduce n_per_sample at event ",
      bad,
      ". Reconstructed = ",
      reconstructed_counts[bad],
      "; n_per_sample = ",
      stan_data$n_per_sample[bad]
    )
  }
  
  # Map every occupied bin row to its site.
  bin_site <- as.integer(
    stan_data$site_id[
      stan_data$bin_event
    ]
  )
  
  occupied_rows_by_site <- tabulate(
    bin_site,
    nbins = S
  )
  
  if (any(occupied_rows_by_site == 0L)) {
    empty_sites <- which(
      occupied_rows_by_site == 0L
    )
    
    stop(
      "The current Stan model requires at least one occupied size ",
      "bin at every site. Sites without size data: ",
      paste(empty_sites, collapse = ", ")
    )
  }
  
  bin_mid <-
    0.5 * (
      stan_data$bin_lower +
        stan_data$bin_upper
    )
  
  bin_width <-
    stan_data$bin_upper -
    stan_data$bin_lower
  
  max_bin_edge_taxon <-
    max(stan_data$bin_upper)
  
  upper_multiplier_init <-
    1 +
    (
      stan_data$upper_multiplier_max - 1
    ) *
    upper_fraction_center
  
  size_upper_taxon_init <-
    max_bin_edge_taxon *
    upper_multiplier_init
  
  if (size_upper_taxon_init <= stan_data$size_lower) {
    stop(
      "The provisional taxon upper endpoint must exceed size_lower."
    )
  }
  
  # Stay comfortably inside Stan parameter boundaries.
  mu_lower_init <- 1e-6
  sigma_lower_absolute <- 1e-6
  
  mu_upper_init <-
    (1 - boundary_margin) *
    stan_data$mu_upper
  
  sigma_upper_init <-
    (1 - boundary_margin) *
    stan_data$sigma_upper
  
  if (
    mu_upper_init <= mu_lower_init ||
    sigma_upper_init <= sigma_lower_absolute
  ) {
    stop(
      "mu_upper or sigma_upper is too close to the lower ",
      "initialization safeguard."
    )
  }
  
  # Stable scalar subtraction on the log scale:
  # log(exp(log_x) - exp(log_y)), assuming log_x > log_y.
  logspace_subtract <- function(log_x, log_y) {
    if (
      !is.finite(log_x) &&
      is.infinite(log_x) &&
      log_x < 0
    ) {
      return(-Inf)
    }
    
    if (
      !is.finite(log_y) &&
      is.infinite(log_y) &&
      log_y < 0
    ) {
      return(log_x)
    }
    
    if (
      !is.finite(log_x) ||
      !is.finite(log_y) ||
      log_y >= log_x
    ) {
      return(-Inf)
    }
    
    log_x +
      log1p(
        -exp(log_y - log_x)
      )
  }
  
  # R equivalent of normal_interval_log_prob() in the Stan model.
  normal_interval_log_prob_r <- function(
    lower,
    upper,
    mu,
    sigma
  ) {
    if (
      !is.finite(mu) ||
      !is.finite(sigma) ||
      sigma <= 0
    ) {
      return(
        rep(-Inf, length(lower))
      )
    }
    
    vapply(
      seq_along(lower),
      function(i) {
        lo <- lower[i]
        up <- upper[i]
        
        if (
          !is.finite(lo) ||
          !is.finite(up) ||
          up <= lo
        ) {
          return(-Inf)
        }
        
        if (up <= mu) {
          log_cdf_upper <- stats::pnorm(
            up,
            mean = mu,
            sd = sigma,
            log.p = TRUE
          )
          
          log_cdf_lower <- stats::pnorm(
            lo,
            mean = mu,
            sd = sigma,
            log.p = TRUE
          )
          
          value <- logspace_subtract(
            log_cdf_upper,
            log_cdf_lower
          )
        } else if (lo >= mu) {
          log_ccdf_lower <- stats::pnorm(
            lo,
            mean = mu,
            sd = sigma,
            lower.tail = FALSE,
            log.p = TRUE
          )
          
          log_ccdf_upper <- stats::pnorm(
            up,
            mean = mu,
            sd = sigma,
            lower.tail = FALSE,
            log.p = TRUE
          )
          
          value <- logspace_subtract(
            log_ccdf_lower,
            log_ccdf_upper
          )
        } else {
          outside_probability <-
            stats::pnorm(
              lo,
              mean = mu,
              sd = sigma
            ) +
            stats::pnorm(
              up,
              mean = mu,
              sd = sigma,
              lower.tail = FALSE
            )
          
          outside_probability <- clamp(
            outside_probability,
            0,
            1 - .Machine$double.eps
          )
          
          value <-
            log1p(
              -outside_probability
            )
        }
        
        # Direct-probability fallback for rare floating-point ties.
        if (!is.finite(value)) {
          interval_probability <-
            stats::pnorm(
              up,
              mean = mu,
              sd = sigma
            ) -
            stats::pnorm(
              lo,
              mean = mu,
              sd = sigma
            )
          
          value <- log(
            max(
              interval_probability,
              .Machine$double.xmin
            )
          )
        }
        
        value
      },
      numeric(1)
    )
  }
  
  # Approximate moments from grouped bins, including the within-bin
  # variance of a uniform distribution. These are starting values only.
  weighted_bin_moments <- function(
    mid,
    width,
    count
  ) {
    total_count <- sum(count)
    
    if (
      !is.finite(total_count) ||
      total_count <= 0
    ) {
      return(
        c(
          mean = NA_real_,
          sd = NA_real_
        )
      )
    }
    
    mean_value <-
      sum(
        count * mid
      ) /
      total_count
    
    variance_value <-
      sum(
        count *
          (
            (mid - mean_value)^2 +
              width^2 / 12
          )
      ) /
      total_count
    
    c(
      mean = mean_value,
      sd = sqrt(
        max(
          variance_value,
          0
        )
      )
    )
  }
  
  global_moments <- weighted_bin_moments(
    mid = bin_mid,
    width = bin_width,
    count = stan_data$bin_count
  )
  
  global_mean <- global_moments["mean"]
  global_sd <- global_moments["sd"]
  
  typical_bin_width <-
    stats::median(bin_width)
  
  if (
    !is.finite(typical_bin_width) ||
    typical_bin_width <= 0
  ) {
    typical_bin_width <- 1
  }
  
  sigma_data_floor <-
    max(
      sigma_lower_absolute,
      typical_bin_width / sqrt(12)
    )
  
  if (
    !is.finite(global_mean) ||
    global_mean <= 0
  ) {
    global_mean <-
      0.5 *
      (
        stan_data$size_lower +
          size_upper_taxon_init
      )
  }
  
  if (
    !is.finite(global_sd) ||
    global_sd <= 0
  ) {
    global_sd <-
      max(
        sigma_data_floor,
        0.20 * global_mean
      )
  }
  
  global_mean <- clamp(
    global_mean,
    mu_lower_init,
    mu_upper_init
  )
  
  global_sd <- clamp(
    global_sd,
    sigma_data_floor,
    sigma_upper_init
  )
  
  # Estimate provisional latent mu and sigma for one site by maximizing
  # the exact grouped-bin likelihood conditional on the provisional
  # taxon-level endpoint.
  fit_site_binned_tnorm <- function(site) {
    use <- which(
      bin_site == site
    )
    
    site_moments <- weighted_bin_moments(
      mid = bin_mid[use],
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
      mu_lower_init,
      mu_upper_init
    )
    
    raw_sigma <- clamp(
      raw_sigma,
      sigma_data_floor,
      sigma_upper_init
    )
    
    # A site occupying only one distinct bin provides weak information
    # about sigma. Avoid initializing it at an artificial near-zero MLE.
    distinct_bins <- nrow(
      unique(
        data.frame(
          lower = stan_data$bin_lower[use],
          upper = stan_data$bin_upper[use]
        )
      )
    )
    
    if (distinct_bins < 2L) {
      return(
        c(
          mu = raw_mu,
          sigma = clamp(
            max(
              raw_sigma,
              0.5 * global_sd,
              sigma_data_floor
            ),
            sigma_data_floor,
            sigma_upper_init
          )
        )
      )
    }
    
    # For initialization, prevent the optimizer from wandering far
    # beyond the modeled size support even though Stan permits it.
    mu_optimizer_upper <- min(
      mu_upper_init,
      max(
        2 * size_upper_taxon_init,
        2 * raw_mu,
        stan_data$size_lower + typical_bin_width
      )
    )
    
    sigma_optimizer_upper <- min(
      sigma_upper_init,
      max(
        size_upper_taxon_init - stan_data$size_lower,
        2 * raw_sigma,
        typical_bin_width
      )
    )
    
    mu_optimizer_upper <- max(
      mu_optimizer_upper,
      raw_mu
    )
    
    sigma_optimizer_upper <- max(
      sigma_optimizer_upper,
      raw_sigma
    )
    
    negative_log_likelihood <- function(par) {
      mu_value <- exp(par[1])
      sigma_value <- exp(par[2])
      
      log_normalizer <-
        normal_interval_log_prob_r(
          lower = stan_data$size_lower,
          upper = size_upper_taxon_init,
          mu = mu_value,
          sigma = sigma_value
        )[1L]
      
      log_bin_probability <-
        normal_interval_log_prob_r(
          lower = stan_data$bin_lower[use],
          upper = stan_data$bin_upper[use],
          mu = mu_value,
          sigma = sigma_value
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
      
      log_likelihood <-
        sum(
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
            mu_lower_init,
            sigma_data_floor
          )
        ),
        upper = log(
          c(
            mu_optimizer_upper,
            sigma_optimizer_upper
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
    
    estimates <- exp(fit$par)
    
    c(
      mu = clamp(
        estimates[1],
        mu_lower_init,
        mu_upper_init
      ),
      sigma = clamp(
        estimates[2],
        sigma_data_floor,
        sigma_upper_init
      )
    )
  }
  
  site_size_estimates <- t(
    vapply(
      seq_len(S),
      fit_site_binned_tnorm,
      numeric(2)
    )
  )
  
  site_mu <-
    site_size_estimates[, "mu"]
  
  site_sigma <-
    site_size_estimates[, "sigma"]
  
  log_mu_init <-
    log(site_mu)
  
  log_sigma_init <-
    log(site_sigma)
  
  # Site expected event counts.
  site_lambda <- vapply(
    seq_len(S),
    function(site) {
      event_counts <-
        stan_data$n_per_sample[
          stan_data$site_id == site
        ]
      
      if (!length(event_counts)) {
        return(0.1)
      }
      
      value <- mean(event_counts)
      
      if (
        !is.finite(value) ||
        value <= 0
      ) {
        value <- 0.1
      }
      
      max(value, 0.1)
    },
    numeric(1)
  )
  
  log_lambda_init <-
    log(site_lambda)
  
  # Robust hierarchical centers.
  alpha_mu_init <-
    stats::median(log_mu_init)
  
  alpha_sigma_init <-
    stats::median(log_sigma_init)
  
  alpha_lambda_init <-
    stats::median(log_lambda_init)
  
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
    
    range_scale <-
      diff(range(x)) / 4
    
    if (!is.finite(range_scale)) {
      range_scale <- 0
    }
    
    deviation_scale <-
      max(abs(x - center)) /
      max_z
    
    value <- max(
      minimum,
      mad_scale,
      range_scale,
      deviation_scale
    )
    
    min(
      value,
      maximum
    )
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
    (
      log_lambda_init -
        alpha_lambda_init
    ) /
    tau_lambda_init
  
  # Method-of-moments NB2 dispersion estimates by site:
  #
  # Var(N) = lambda + lambda^2 / phi
  # phi = lambda^2 / [Var(N) - lambda].
  phi_by_site <- vapply(
    seq_len(S),
    function(site) {
      event_counts <-
        stan_data$n_per_sample[
          stan_data$site_id == site
        ]
      
      if (length(event_counts) < 2L) {
        return(NA_real_)
      }
      
      mean_count <- mean(event_counts)
      variance_count <- stats::var(event_counts)
      
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
    stats::median(finite_phi)
  } else {
    # Prior median when event replication cannot identify dispersion.
    20
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
    log(sigma_data_floor)
  
  log_sigma_upper <-
    log(sigma_upper_init)
  
  # Keep hierarchical-scale parameters away from Stan bounds.
  tau_size_lower <- 0.02
  tau_size_upper <- 2.40
  tau_lambda_lower <- 0.05
  tau_lambda_upper <- 4.40
  
  set.seed(seed)
  
  lapply(
    seq_len(chains),
    function(chain_id) {
      # Perturb the site size parameters slightly on the log scale.
      log_mu_chain <- clamp(
        log_mu_init +
          stats::rnorm(
            S,
            mean = 0,
            sd = 0.01
          ),
        log_mu_lower,
        log_mu_upper
      )
      
      log_sigma_chain <- clamp(
        log_sigma_init +
          stats::rnorm(
            S,
            mean = 0,
            sd = 0.01
          ),
        log_sigma_lower,
        log_sigma_upper
      )
      
      alpha_mu_chain <-
        stats::median(log_mu_chain) +
        stats::rnorm(
          1,
          mean = 0,
          sd = 0.015
        )
      
      alpha_sigma_chain <-
        stats::median(log_sigma_chain) +
        stats::rnorm(
          1,
          mean = 0,
          sd = 0.015
        )
      
      tau_mu_chain <- clamp(
        tau_mu_init *
          exp(
            stats::rnorm(
              1,
              mean = 0,
              sd = 0.02
            )
          ),
        tau_size_lower,
        tau_size_upper
      )
      
      tau_sigma_chain <- clamp(
        tau_sigma_init *
          exp(
            stats::rnorm(
              1,
              mean = 0,
              sd = 0.02
            )
          ),
        tau_size_lower,
        tau_size_upper
      )
      
      # Preserve the provisional site log-lambdas after perturbing the
      # non-centered hierarchical parameters.
      target_log_lambda <-
        log_lambda_init +
        stats::rnorm(
          S,
          mean = 0,
          sd = 0.02
        )
      
      alpha_lambda_chain <-
        stats::median(target_log_lambda) +
        stats::rnorm(
          1,
          mean = 0,
          sd = 0.02
        )
      
      tau_lambda_chain <- clamp(
        tau_lambda_init *
          exp(
            stats::rnorm(
              1,
              mean = 0,
              sd = 0.02
            )
          ),
        tau_lambda_lower,
        tau_lambda_upper
      )
      
      z_lambda_chain <-
        (
          target_log_lambda -
            alpha_lambda_chain
        ) /
        tau_lambda_chain
      
      z_lambda_chain <- clamp(
        z_lambda_chain,
        -3,
        3
      )
      
      # Scalar endpoint fraction. Perturb on the logit scale and keep
      # it away from the hard boundaries zero and one.
      upper_fraction_chain <- clamp(
        stats::plogis(
          stats::qlogis(
            upper_fraction_center
          ) +
            stats::rnorm(
              1,
              mean = 0,
              sd = 0.15
            )
        ),
        0.05,
        0.95
      )
      
      list(
        alpha_log_mu =
          as.numeric(alpha_mu_chain),
        
        alpha_log_sigma =
          as.numeric(alpha_sigma_chain),
        
        log_mu_site =
          as.numeric(log_mu_chain),
        
        log_sigma_site =
          as.numeric(log_sigma_chain),
        
        tau_log_mu =
          as.numeric(tau_mu_chain),
        
        tau_log_sigma =
          as.numeric(tau_sigma_chain),
        
        upper_fraction =
          as.numeric(upper_fraction_chain),
        
        alpha_log_lambda =
          as.numeric(alpha_lambda_chain),
        
        tau_log_lambda =
          as.numeric(tau_lambda_chain),
        
        z_lambda =
          as.numeric(z_lambda_chain),
        
        log_phi =
          as.numeric(
            clamp(
              log_phi_init +
                stats::rnorm(
                  1,
                  mean = 0,
                  sd = 0.03
                ),
              log_phi_lower,
              log_phi_upper
            )
          )
      )
    }
  )
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

extract_max = function(filePath = NULL, names = NULL, siteIndex = NULL){
  mod = readRDS(filePath)
  vars = paste(c('max_ref_rep', 'size_mean','size_median'), collapse = '|')
  fullSumm = mod$summary()
  filteredSumm = fullSumm[grepl(vars, fullSumm$variable), c('variable','mean', 'median', 'sd', 'q5','q95')]
  filteredSumm$acceptedTaxonID = names
  filteredSumm$siteID = gsub("\\w+\\[(\\d{1,2})\\]", "\\1", filteredSumm$variable)
  summDf = left_join(filteredSumm, data.frame(siteIndex), by = c('acceptedTaxonID','siteID'))
  return(summDf)
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