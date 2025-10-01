library(AOData)
library(dplyr)
library(tidyverse)
library(trainomeMetaData)
library(stringi)
source("R/Trainome_functions.R")

# download_ome(download = "alphaomega_gene_rsem")
# 
 # AO_gene <- readRDS("ome-data/alphaomega_gene_rsem.rds")
 # AO_gene <- AO_gene$expected_count
#AO_thickness <- AOData::thickness

#extract gene expression data

alpha_omega_counts <- extract_rsem_gene_counts("data/gene_counts/Alpha_Omega_RSEM_outputs/")


# Extract the sequence ID to match the extraction sequence given in the metadata
# extraction sequence is that number following "s" on the sample name
seq_df <- as.data.frame(colnames(alpha_omega_counts [, -1])) %>%
        mutate(seq_id = as.double(stri_extract_first_regex(colnames(alpha_omega_counts [, -1]), "\\d+")))



# rename the first column to seq_sample_id to match the other datasets
colnames(seq_df)[colnames(seq_df) == "colnames(alpha_omega_counts[, -1])"] <- "seq_sample_id"


# read the document containing sequence ids
A_O_seq_list <- readxl::read_excel("Alpha_Omega_sample_list_transcriptomics.xlsx") %>%
        filter(Tissue == "muscle") %>%
        mutate(time = case_when(time_rep == "T1rna1" ~ "PreExc",
                                time_rep == "T4rna1" ~ "PostExc")) %>%
        dplyr::select(subject, Tissue, sample, time)
AOData::thickness


AO_thickness <- AOData::thickness %>%
        mutate(time = case_when(time == "T1" ~ "PreExc",
                                time == "T2" ~ "PreTrain",
                                time == "T4" ~ "PostExc")) %>%
                       filter(time == "PreExc" | time == "PostExc") %>%
        dplyr::select(-VIthick)


thickness_pre <- AO_thickness %>%
        filter(time == "PreExc")

thickness_post <- AO_thickness %>%
        filter(time == "PostExc")

thickness_merged <- thickness_pre %>%
        inner_join(thickness_post,by = c("participant", "leg", "condition") ) %>%
        mutate(pct_muscle_change = ((VLthick.y - VLthick.x) / VLthick.x) * 100) %>%
        dplyr::select(participant, leg, pct_muscle_change, condition) %>%
        drop_na(pct_muscle_change)

length(unique(thickness_merged$participant))


sequenced_samples <- seq_df %>%
        inner_join(A_O_seq_list, by = c("seq_id" = "sample"))%>%
        inner_join(AOData::idkeys, by = c("subject" = "participant"))%>%
        # rename column to match other datasets
        rename(BMI = bmi, participant = subject) %>%
        inner_join(AOData::seq_samples, by = c("participant", "seq_id" = "extraction_seq" ) ) %>%
        rename(time = time.x) %>%
        filter(time == "PreExc") %>%
        dplyr::select(participant, seq_sample_id,  leg,  age, sex,BMI, condition) %>%
        right_join(thickness_merged, by = c("participant", "leg", "condition")) %>%
        #some participants had samples collected only on one leg at baseline
        #fill in their missing values
        group_by(participant) %>%
        # carry the last observed value
        mutate(across(c(sex, seq_sample_id,  age, BMI), ~ ifelse(is.na(.), zoo::na.locf(., na.rm= F,
                                                                                             fromLast = F), .))) %>%
        # carry next observed value
        mutate(across(c( sex, seq_sample_id, age, BMI), ~ ifelse(is.na(.), zoo::na.locf(., na.rm = FALSE, fromLast = TRUE), .)))%>%
        drop_na() %>%
        ungroup()%>%
        mutate_if(is.numeric,round, digits = 2) %>%
 
        
               mutate(participant = as.character(participant))
     
length(unique(sequenced_samples$participant))     


# check those with duplicates and see if the match to different conditions
# duplicates <- sequenced_samples %>%
#         group_by(seq_sample_id) %>%
#         filter(n() > 1) %>%
#         arrange(seq_sample_id)

# Add the study name and volume of exercise

sequenced_samples$study <- "Alpha/Omega"
sequenced_samples$volume <- 3

sequenced_samples["sex"][sequenced_samples["sex"] == "m" ] <- "male"

sequenced_samples["sex"][sequenced_samples["sex"] == "f" ] <- "female"


hist(sequenced_samples$age)
hist(sequenced_samples$BMI)
hist(sequenced_samples$pct_muscle_change)

# remove duplicate rows
sequenced_samples <- sequenced_samples[!duplicated(sequenced_samples),]


saveRDS(sequenced_samples, "data/Alpha_Omega/Alpha_Omega_metadata.RDS")
