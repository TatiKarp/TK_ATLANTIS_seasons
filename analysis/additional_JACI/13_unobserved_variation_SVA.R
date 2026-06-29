## Run SVA analysis

library(sva)
library(dplyr)
library(edgeR)


# location 
conf <- config::get()

source(file.path(conf$script_path, "R/utils_DE.R"))

# load data
master.Table <- read.csv(file.path(conf$data_path, "Season/Season_new_date/ATLANTIS_master_table_seasons_new.csv"), header = TRUE) %>%
  dplyr::select(-c("VISDAT_c", "month", "season",  "season_2", "Date", "X", "Customer.ID",
                   "PAT_NO.y", "Requestion.y", "RQN", "DayYear.Maarten", "difference", "season.Maaike", "season_2_Maaike"))

# 1. create 2 groups 
master.Table <- master.Table %>%
  mutate(new_seasons = as.factor(if_else(((DayYear.Maaike >= 135) & (DayYear.Maaike < 319)), "summer_autumn", "winter_spring")))
table(master.Table$new_seasons)
table(master.Table$new_seasons, master.Table$asthma.status)

# change the order of the factors
master.Table$asthma.status <- factor(master.Table$asthma.status, levels = c("H","A"))

# 2. expression data
expression.data <- read.csv(file.path(conf$data_path, "Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"), header = TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table$GenomeScan_ID)) %>%
  as.matrix()

DGEL <- edgeR::DGEList(expression.data)
# Get your normalized log-CPM matrix
logCPM <- cpm(DGEL, log = TRUE, prior.count = 3)

# Define full model (with season) and null model
mod  <- model.matrix(~new_seasons + age + gender + smoking.status + asthma.status, data = master.Table)
mod0 <- model.matrix(~age + gender + smoking.status + asthma.status, data = master.Table)  # null: no season


# 3. Estimate surrogate variables
svobj <- svaseq(logCPM, mod, mod0)

# 4. Add SVs to your design matrix and rerun edgeR
mod_sv <- cbind(mod, svobj$sv)




