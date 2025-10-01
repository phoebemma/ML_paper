library(trainomeMetaData)
library(dplyr)
library(tidyverse)
library(reliefdata)
library(stringi)

source("R/Trainome_functions.R")
# Load the relief gene counts
Relief_genes <- extract_rsem_gene_counts("data/gene_counts/Relief_RSEM_outputs/")



# extraction sequence is that number following "R" on the sample name
seq_df <- as.data.frame(colnames(Relief_genes[, -1])) %>%
        mutate(seq_id = as.double(stri_extract_first_regex(colnames(Relief_genes[, -1]), "\\d+")))




# rename the first column to seq_sample_id to match the other datasets
colnames(seq_df)[colnames(seq_df) == "colnames(Relief_genes[, -1])"] <- "seq_sample_id"



# Extract the 

# Get sequence information from the excel file obtained from Kristian
Relief_metadata <- readxl::read_excel("Relief_sampleIDs.xlsx")  %>%
        dplyr::select(-c( weight))

# The subject exists as numbers, change to characters
Relief_metadata$subject <- as.character(Relief_metadata$subject)
# Change the "sample" column name to seq_id
colnames(Relief_metadata)[colnames(Relief_metadata) == "sample"] <- "seq_id"
# Rename subject to participant
colnames(Relief_metadata)[colnames(Relief_metadata) == "subject"] <- "participant"


# change the relief volume participant to character, as it would be needed later to get condition
# relief_volume$participant <- as.character(relief_volume$participant)

# Add study name
 Relief_metadata$study <- "ReLiEf"
 relief_volume <- relief_volume %>%
         mutate(participant = as.character(participant))

# extract the muscle thickness
thickness <- reliefdata::relief_thickness %>%
        mutate(time = case_when(time == "pre" ~ "PreExc",
                                time == "post" ~ "PostExc")) %>%
        filter(time == "PreExc" | time == "PostExc") %>%
        # some participants had multiple measurements
        group_by(participant, leg, time) %>%
        summarise(avg_VLthick = mean(VLthick)) %>%
        ungroup() %>%
        dplyr::select(participant, leg, time,avg_VLthick) %>%
        inner_join(relief_volume, by = c("participant", "leg")) %>%
        mutate(volume = if_else(condition == "low", 1, 3))


thickness_pre <- thickness %>%
        filter(time == "PreExc")


thickness_post <- thickness %>%
        filter(time == "PostExc")

thickness_merged <- thickness_pre %>%
        inner_join(thickness_post,by = c("participant", "leg", "volume") ) %>%
        mutate(pct_muscle_change = ((avg_VLthick.y - avg_VLthick.x) / avg_VLthick.x) * 100) %>%
        dplyr::select(participant, leg, pct_muscle_change, volume) 




Relief_metadf <- Relief_metadata %>%
        inner_join(seq_df, by = "seq_id") %>%
        right_join(relief_participants, by =  "participant") %>%
        mutate(BMI = round(weight/((height/100)^2), digits = 1) ) %>%
        dplyr::select(study, participant, age, sex, time_rep, seq_sample_id, BMI) %>%
        #Extract leg column
        mutate( leg = case_when(time_rep == "t1rnaL" ~ "L",
                                       time_rep == "t2rnaL" ~ "L",
                                       time_rep == "t3rnaL" ~ "L",
                                       time_rep == "t1rnaR"  ~ "R",
                                       time_rep == "t2rnaR"  ~ "R",
                                       time_rep == "t3rnaR"  ~ "R"),
                # Add time column
                time = case_when(time_rep == "t1rnaL" ~ "PreExc",
                                 time_rep == "t2rnaL" ~ "MidExc",
                                 time_rep == "t3rnaL" ~ "PostExc",
                                 time_rep == "t1rnaR"  ~ "PreExc",
                                 time_rep == "t2rnaR"  ~ "MidExc",
                                 time_rep == "t3rnaR"  ~ "PostExc") ) %>%
        dplyr::select(study, participant, age, sex, time,leg, seq_sample_id, BMI) %>%
        dplyr::filter(time == "PreExc") %>%
        full_join(thickness_merged, by = c("participant", "leg")) %>%
        dplyr::select(-time) %>%
        #some participants had samples collected only on one leg at baseline
        #fill in their missing values
        group_by(participant) %>%
        # carry the last observed value
        mutate(across(c(study,sex, seq_sample_id, age, BMI), ~ ifelse(is.na(.), zoo::na.locf(., na.rm= F,
                                                                                       fromLast = F), .))) %>%
        # carry next observed value
        mutate(across(c(study, sex, seq_sample_id, age, BMI), ~ ifelse(is.na(.), zoo::na.locf(., na.rm = FALSE, fromLast = TRUE), .))) %>%
        drop_na() %>%
        ungroup()%>%
        mutate_if(is.numeric,round, digits = 2) 

Relief_metadf$condition <- "RM10"
        
# remove duplicate samples
Relief_metadf <- Relief_metadf[!duplicated(Relief_metadf),]
   
length(unique(Relief_metadf$participant))
 

saveRDS(Relief_metadf, "data/Relief/Relief_metadata.RDS")
