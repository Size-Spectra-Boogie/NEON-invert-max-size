source(here::here('code/00_init.R'))
i_am('code/03_estimate-max-size.R')

macroDW <- readRDS(file = here("data/macro_dw_raw.rds"))

# clean macros data.frame and aggregate abundances within sites by taxonID
macroDWList = macroDW %>% 
  dplyr::select(siteID, collectDate, sampleID, acceptedTaxonID,sizeClass, dw, no_m2) %>% 
  dplyr::mutate(collectDate = as.Date(collectDate)) %>% 
  dplyr::mutate(collectMonth = as.Date(
    paste0(
      lubridate::year(collectDate),
      "-",
      lubridate::month(collectDate),
      "-01"),
    format = "%Y-%m-%d"),
    collectYear = lubridate::year(collectDate)) %>% 
  dplyr::summarise(no_m2 = sum(no_m2), .by = c(siteID, collectYear, acceptedTaxonID, sizeClass, dw)) %>% 
  named_group_split(siteID, acceptedTaxonID)

## How many groups only have 1 size class measured?

macroDWList %>% keep(~nrow(.) == 1) %>% length()
# 2k taxa-site combinations

## Which has the most observations?
max_r = map_int(macroDWList, ~nrow(.x)) %>% max();max_r


## practice workflow on a single, well-sampled taxa
x = Filter(function(x) nrow(x) == max_r, macroDWList) %>% purrr::flatten_df()

x %>% 
  ggplot()+
  geom_histogram(aes(x = sizeClass, y = ..density.., weight = no_m2))
# In each sampling year, we sample the length classes
# to account for the binning structure, we resample from a uniform
# distribution +/- 0.5mm from the length class


x_samp = x %>% 
  slice_sample(n = 100, by = collectYear, weight_by = no_m2, replace = TRUE)

x_list = x_samp %>% 
  named_group_split(collectYear) %>% 
  map(~.x %>%
        select(sizeClass) %>%
        mutate(sizeClass_est = midpoint_resample_vec(sizeClass)) %>% 
        select(sizeClass_est) %>% 
        unlist)

# fit the model with the resampled data sets
# x_mod = fishmax::fit_max_model(x_list)
# saveRDS(x_mod, here('data/models/CHISP8_PRLA.rds'))
x_mod = readRDS(here('data/models/CHISP8_PRLA.rds'))
# get the estimates for each model
get_max(x_mod)

# plot all estimates
plot_max(x_mod)

## try with another species with less sampling

x2 = macroDWList[[1]]

x2_samp = x2 %>% 
  slice_sample(n = 1000, by = collectYear, weight_by = no_m2, replace = TRUE)

x2_list = x2_samp %>% 
  named_group_split(collectYear) %>% 
  map(~.x %>%
        select(sizeClass) %>%
        mutate(sizeClass_est = midpoint_resample_vec(sizeClass)) %>% 
        select(sizeClass_est) %>% 
        unlist)

# x2_mod = fishmax::fit_max_model(x2_list)
# saveRDS(x2_mod, here('data/models/ABLSP_ARIK.rds'))
x2_mod = readRDS(here('data/models/ABLSP_ARIK.rds'))

get_max(x2_mod)

plot_max(x2_mod)

## With resampling of 100 body sizes for each year
# model           max_fit.50% max_lwr.10% max_upr.90%
# 1    EVT (GEV)   10.355262    7.759288   15.011650
# 2 EVT (Gumbel)    9.729430    7.438284   13.445614
# 3          EFS    9.310764    7.585413   12.833442
# 4        EFSmm    8.722388    8.402264    9.091755


## with resampling 1000 body sizes for each year
# model           max_fit.50% max_lwr.10% max_upr.90%
# 1    EVT (GEV)   10.504437    7.898656   15.446824
# 2 EVT (Gumbel)   10.015032    7.729730   13.894558
# 3          EFS    9.331321    7.636065   12.748386
# 4        EFSmm    9.502151    9.387443    9.618906

## try with another species with only only 3 size classes observed (the minimum required)

x3 = macroDWList[[7]]

x3_samp = x3 %>% 
  slice_sample(n = 100, by = collectYear, weight_by = no_m2, replace = TRUE)

x3_list = x3_samp %>% 
  named_group_split(collectYear) %>% 
  map(~.x %>%
        select(sizeClass) %>%
        mutate(sizeClass_est = midpoint_resample_vec(sizeClass)) %>% 
        select(sizeClass_est) %>% 
        unlist)

# x3_mod = fishmax::fit_max_model(x3_list)
# saveRDS(x3_mod, here('data/models/AEOSP2_ARIK.rds'))
x3_mod = readRDS(here('data/models/AEOSP2_ARIK.rds'))

get_max(x3_mod)

plot_max(x3_mod)
