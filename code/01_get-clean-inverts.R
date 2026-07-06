here::i_am('code/01_get-clean-inverts.R')
source(here::here('code/00_init.R'))

if(rerun_data){
# code to import and clean invertebrate data from s3 bucket
# only needs to be run upon new data releases

## set the s3 bucket specifications
bucket = "neonsizedata-392202703749-us-east-2-an"
prefix = "neonstore/parquet"

## connect to the s3 bucket
s3 = arrow::s3_bucket(bucket = paste0(bucket,"/",prefix,"/"))
remote = neon_remote_db(bucket = s3)

## read the invertebrate files in and save to data

neon_remote(
  table = "inv_taxonomyProcessed-basic",
  product = "DP1.20120.001",
  db = remote
)  %>%
  dplyr::collect() %>% 
  dplyr::select(siteID, collectDate, sampleID, 
                acceptedTaxonID,scientificName,
                morphospeciesID,class, subclass,
                order, family, genus, specificEpithet,
                sizeClass, estimatedTotalCount,
                sampleCondition
                ) %>% 
  saveRDS(here('data/NEON-macros.rds'))
  macros <<- readRDS(here('data/NEON-macros.rds'))
  
  macros %>% 
    dplyr::select(acceptedTaxonID, class, subclass, order, family, genus, specificEpithet, scientificName) %>% 
    dplyr::distinct() %>% 
    saveRDS(here('data/NEON-macros-tax.rds'))
  
neon_remote(
    table = "inv_fieldData",
    product = "DP1.20120.001",
    db = remote) %>% 
    collect() %>% 
   select(siteID, sampleID, benthicArea) %>% 
  saveRDS(here('data/NEON-macros-field.rds'))

  fieldData <<- readRDS(here('data/NEON-macros-field.rds'))

  print("Invertebrate data were updated.")
} else{
  macros <<- readRDS(here('data/NEON-macros.rds'))
  print("Invertebrate data were not updated. Loaded previously cleaned data.")
}