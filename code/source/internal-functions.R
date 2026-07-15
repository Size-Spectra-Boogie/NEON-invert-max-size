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

fit_max_model_named = function(df = NULL, mod_name = NULL, rerun = FALSE, overwrite = FALSE){
  ## currently, we are saving this locally to avoid pushing to git
  mod_path = paste0(here('ignore/models'),"/",mod_name,'.rds')
  if(any(rerun, !file.exists(mod_path))){
    if(all(file.exists(mod_path),!overwrite)){
      warning('Model file already exists and `overwrite` = FALSE. Set to TRUE to overwrite existing files.')
      break()
    }
    fit = fishmax::fit_max_model(df, model_type = c('efs','efsmm'), output_dir = here('ignore/models/'))
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
