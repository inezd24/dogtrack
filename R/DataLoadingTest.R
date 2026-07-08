# ==============================================================================
# Dogtrack R package update and examination
# Package : dogtrack
# Purpose : Separately apply pipeline to examine each function (including non-
#           exported functions) and make any necessary changes. 
# Created : 17 June 2026
# Updated : 18 June 2026 [see below]
# ==============================================================================


# ==============================================================================

# Last updated:
if (exists("last_updated")) {
  cat(paste0("File last updated on: ", last_updated))
} else {
  cat("File not yet updated. Update with today's date.\n")
  last_updated <- Sys.Date()
  cat(paste0("File last updated on: ", last_updated))
}

# Emtpy environment [optional]
rm(list = ls())

# Set working directory
setwd("U:/Data/GPS")

# Import packages
library(dogtrack)
library(tidyverse)
library(ggplot2)


# =================== STEP 1: LOAD METADATA ====================================

# Set col aliases
column_aliases = list(household_id = 'Household ID',
                      dog_id = 'Dog ID',
                      field_site = 'Field site')

# Load metadata
metadata_uganda <- load_metadata(metadata_path = "MetaData.xlsx",
                                  col_aliases = column_aliases,
                                  sheet = 'Uganda')
metadata_chad <- load_metadata(metadata_path = "MetaData.xlsx",
                                 col_aliases = column_aliases,
                                 sheet = 'Chad')


# =================== LOADING PREFIX MAPS ======================================

# Load prefix map Uganda
pfx_map_uganda <- build_prefix_map(metadata_uganda)
cat(sprintf("  %d dogs across %d sites: %s\n\n", 
            nrow(metadata_uganda), nrow(pfx_map_uganda), 
            paste(pfx_map_uganda$field_site, collapse = ", ")))
# 499 dogs across 3 sites: Masaka, Arua, Soroti

# Load prefix map Chad
pfx_map_chad <- build_prefix_map(metadata_chad)
cat(sprintf("  %d dogs across %d sites: %s\n\n", 
            nrow(metadata_chad), nrow(pfx_map_chad), 
            paste(pfx_map_chad$field_site, collapse = ", ")))
# 484 dogs across 4 sites: SARH, Bongor, N'Djamena (TN), N'Djamena (TP)



# =================== LOADING STATIC DATA ======================================

# Load static data Uganda
static_data_uganda <- load_static_tests(data_dir = "U:/Data/GPS/Uganda", 
                                        prefix_map = pfx_map_uganda,
                                        file_pattern = "^static",
                                        device_type = 'columbus',
                                        col_aliases = NULL) 
# load_static_tests(): loaded 46288 GPS fixes from 9 static test session(s) 
# across 9 site/number combination(s).


# Load static data Chad
static_data_chad <- load_static_tests(data_dir = "U:/Data/GPS/Chad", 
                                      prefix_map = pfx_map_chad,
                                      file_pattern = "^static",
                                      device_type = 'mix',
                                      col_aliases = NULL) 


# ====================== LOAD DOG DATA =========================================

# Load dog data Uganda
dogdata_uganda <- load_dog_data(data_dir = "U:/Data/GPS/Uganda",
                                metadata = metadata_uganda,
                                device_type = 'columbus',
                                col_aliases = NULL) 
# Info: 2355966 GPS fixes from 263 dogs.
# after metadata join -- 237 dogs, 2085914 GPS fixes retained.
# CSV files with no metadata match (excluded): 1 
# 6 file(s) are empty and will be SKIPPED:
# U:/Data/GPS/Uganda/UAR046-01F.CSV
# U:/Data/GPS/Uganda/UAR046-01H.CSV
# U:/Data/GPS/Uganda/UAU017-01F.CSV
# U:/Data/GPS/Uganda/UAU017-01L.CSV
# U:/Data/GPS/Uganda/UAU017-01M.CSV
# U:/Data/GPS/Uganda/UAU017-01Y.CSV

# Load dog data Chad
dogdata_chad <- load_dog_data(data_dir = "U:/Data/GPS/Chad",
                              metadata = metadata_chad,
                              device_type = 'mix',
                              col_aliases = NULL) # 296 files
# found 296 dog CSV file(s) after session selection.                                 
# device types detected -- catlogger: 206, columbus: 90 
# 59 file(s) disambiguated via `version` (multiple files shared a base_key).
# loaded 39528451 GPS fixes from 288 dogs.
# after metadata join -- 286 dogs, 39150750 GPS fixes retained.
# CSV files with no metadata match (excluded): 2  TNR029-01 and TNR056-02
#  9 file(s) are empty and will be SKIPPED:
#U:/Data/GPS/Chad/Bongor/Bongor rural chien/TBR007-01-102-not fully sure about the ID/25191901.CSV
#U:/Data/GPS/Chad/Bongor/Bongor rural humains observation conjointe/TBR007-01-P-89/24125101.CSV
#U:/Data/GPS/Chad/Bongor/Bongor rural humains observation conjointe/TBR007-01-P-89/24132226.CSV
#U:/Data/GPS/Chad/Bongor/Bongor rural humains observation conjointe/TBR007-01-P-89/25195718.CSV
#U:/Data/GPS/Chad/Bongor/Bongor rural humains observation conjointe/TBR007-01-P-89/25200037.CSV
#U:/Data/GPS/Chad/Bongor/Bongor rural humains observation conjointe/TBR007-01-P-89/25204001.CSV
#U:/Data/GPS/Chad/Bongor/Bongor rural humains observation conjointe/TBR007-01-P-89/26051919.CSV
#U:/Data/GPS/Chad/Bongor/Bongor rural humains observation conjointe/TBR007-01-P-89/26071731.CSV
#U:/Data/GPS/Chad/Bongor/Bongor rural humains observation conjointe/TBR007-01-P-89/26072706.CSV
# In load_dog_data(data_dir = "U:/Data/GPS/Chad", metadata = metadata_chad,  :
# load_dog_data(): 58 file(s) do not match the expected naming convention 
# ([Country][Location][Setting]-[house]-[individual]) and are EXCLUDED:


# =================== DERIVE SITE THRESHOLDS ===================================

# Uganda
thresholds_uganda <- derive_site_thresholds(static_data = static_data_uganda, 
                                     prefix_map = pfx_map_uganda, 
                                     dog_data = dogdata_uganda,
                                     height_overrides = list('Arua' = 1600,
                                                             'Masaka' = 1600,
                                                             'Soroti' = 1500),
                                     baseline_dist_noise_m = 50)
# derive_site_thresholds(): site 'Soroti' has no static tests; 
# height_threshold from override, dist_noise_m from baseline.

# Chad
thresholds_chad <- derive_site_thresholds(static_data = static_data_chad, 
                                          prefix_map = pfx_map_chad, 
                                          dog_data = dogdata_chad,
                                          height_overrides = list('SARH' = 450,
                                                                  'Bongor' = 450,
                                                                  "N'Djamena" = 450),
                                          baseline_dist_noise_m = 50)
# derive_site_thresholds(): site 'SARH' has no static tests; 
# height_threshold from override, dist_noise_m from baseline.


# ==================== FLAGGING & SUMMARY: UGANDA ==============================


# Remove individual who we are uncertain of. 

# HDOP
dogdata_uganda <- filter_hdop(dogdata_uganda) 
# 2346397 fixes | flagged: 1100 (0.05%) | NA HDOP: 1 | threshold: HDOP > 4.9

# Speed
dogdata_uganda <- filter_speed(dogdata_uganda)
# 2346397 fixes | flagged: 3811 (0.162%) | NA speed: 94750 | 
# delta_t below noise floor (< 8s): 94427
# of 3811 speed-flagged fixes -- 6 also HDOP-flagged, 3805 HDOP-clean (investigate)

# Angle
dogdata_uganda <- filter_angle(dogdata_uganda, thresholds_uganda)
# 2346397 fixes | evaluable triplets: 46818 (2.0%) | flagged: 10425 
# (0.444% of all, 22.267% of evaluable) | angle threshold: > 150 deg
#Warning message:
 # In filter_angle(dogdata_uganda, thresholds_uganda) :
 # filter_angle(): 60 fix(es) with duplicate datetime within a dog track -- 
 # these fixes receive flag_angle = FALSE.

# Height
dogdata_uganda <- filter_height(dogdata_uganda, thresholds_uganda)
# filter_height(): 2346397 fixes | flagged: 444744 (18.954%) | NA HEIGHT: 1
# of 444744 height-flagged -- 187 also HDOP-flagged, 444557 HDOP-clean

# Summarise 
flag_summary(dogdata_uganda)


# ==================== FLAGGING & SUMMARY: CHAD ================================

# HDOP
dogdata_chad <- filter_hdop(dogdata_chad) 
# 39150750 fixes | flagged: 4263 (0.01%) | NA HDOP: 0 | threshold: HDOP > 4.9

# Speed
dogdata_chad <- filter_speed(dogdata_chad)
# 39150750 fixes | flagged: 947 (0.002%) | NA speed: 38563212 | 
# delta_t below noise floor (< 8s): 129283
# of 947 speed-flagged fixes -- 0 also HDOP-flagged, 947 HDOP-clean (investigate)

# Angle
dogdata_chad <- filter_angle(dogdata_chad, thresholds_chad)
# filter_angle(): 39150750 fixes | evaluable triplets: 4172 (0.0%) | 
# flagged: 234 (0.001% of all, 5.609% of evaluable) | angle threshold: > 150 deg
# Warning message:
  # In filter_angle(dogdata_chad, thresholds_chad) :
  # filter_angle(): 38433643 fix(es) with duplicate datetime within a dog track 
  # -- these fixes receive flag_angle = FALSE.

# Height
dogdata_chad <- filter_height(dogdata_chad, thresholds_chad)
# filter_height(): 39150750 fixes | flagged: 73189 (0.187%) | NA HEIGHT: 0

# Summarise
flag_summary(dogdata_chad)

# =================== DATA EXAMINATION:UGANDA ==================================


# Effect of household
households <- dogdata_uganda %>% 
  group_by(field_site, setting, household_id, dog_id) %>% 
  summarise(obs = n())
ao <- aov(obs ~ household_id, data = households)
summary(ao)

# Effect of setting
overview <- dogdata_uganda %>% 
  group_by(field_site, setting, dog_id) %>% 
  summarise(obs=n())
ao <- aov(obs ~ setting, data = overview)
summary(ao)
ggplot(overview, aes(x=setting, y=obs, color=setting)) +
  geom_boxplot() +
  theme_minimal()
ggsave("setting.png", height = 16, width = 16)

# Effect of field sites
overview <- dogdata_uganda %>% 
  group_by(field_site, setting, dog_id) %>% 
  summarise(obs=n())
ao <- aov(obs ~ field_site, data = overview)
summary(ao)
ggplot(overview, aes(x=field_site, y=obs, color=setting)) +
  geom_boxplot() +
  theme_minimal() +
  theme(text = element_text(size = 25))
ggsave("fieldsite_setting.png", height = 16, width = 16)


# =================== DATA EXAMINATION:CHAD ====================================


# Effect of household
households <- dogdata_chad %>% 
  group_by(field_site, setting, household_id, dog_id) %>% 
  summarise(obs = n())
ao <- aov(obs ~ household_id, data = households)
summary(ao) # very slight effect of household

# Effect of setting
overview <- dogdata_chad %>% 
  group_by(field_site, setting, dog_id) %>% 
  summarise(obs=n())
ao <- aov(obs ~ setting, data = overview)
summary(ao)
ggplot(overview, aes(x=setting, y=obs, color=setting)) +
  geom_boxplot() +
  theme_minimal()
ggsave("setting_chad.png", height = 16, width = 16)

# Effect of field sites
overview <- dogdata_chad %>% 
  group_by(field_site, setting, dog_id) %>% 
  summarise(obs=n())
ao <- aov(obs ~ field_site, data = overview)
summary(ao)
ggplot(overview, aes(x=field_site, y=obs, color=setting)) +
  geom_boxplot() +
  theme_minimal() +
  theme(text = element_text(size = 25))
ggsave("fieldsite_setting_chad.png", height = 16, width = 16)

# Effect of device type on satellite number
ao <- aov(SAT ~ device_type, data = dogdata_chad)
summary(ao) # very significant differences
ggplot(dogdata_chad, aes(x=as.numeric(SAT), fill=device_type)) +
  geom_density() +
  theme_minimal() +
  theme(text = element_text(size = 25))
