here::i_am('code/00_init.R')

## load in the libraries and set modeling and plotting specifications
library(here)
library(magrittr)
library(neonstore)
library(magrittr) 
library(arrow)
library(aws.s3)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(junkR)
library(cmdstanr)
library(fishmax)
library(viridis)

theme_set(theme_minimal())

rerun_data = FALSE

## source the internal function files
purrr::walk(list.files(here('code/source/'), pattern = '*.R', full.names = TRUE), \(x) source(x))
