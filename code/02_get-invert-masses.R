here::i_am('code/02_get-invert-masses.R')
source(here::here('code/00_init.R'))

# 1) Load in macro data
# Do we need to update the data file with new releases? set rerun_data = TRUE in 00_init.R
# else load in previously pulled data
source(here('code/01_get-clean-inverts.R'))


# 2) Add LW coefficients, estimate dry weights  ------------------------------------
coeff <- read.csv(here("data/macro_lw_coeffs.csv"))

# add length weight coefficients by taxon
MN.lw <- LW_coef(
  x = macros,
  lw_coef = coeff,
  percent = TRUE
)

# questionable measurements ####
# filter out individuals that were "damaged" and measurement was affected
# this is a flag which is added by NEON
MN.no.damage <- MN.lw %>%
  filter(!str_detect(
    sampleCondition,
    "measurement"
  )) %>%
  est_dw(fieldData = fieldData)

# 3) filter out NA values in dw
macro_dw <- MN.no.damage %>%
  filter(!is.na(dw), !is.na(no_m2))

# what % of columns are maintained?
nrow(macro_dw) / nrow(MN.no.damage)

saveRDS(macro_dw, file = here("data/macro_dw_raw.rds"))
