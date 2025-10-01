library(trainomeMetaData)
library(dplyr)
library(tidyverse)



copd_metadata <- trainomeMetaData::copd_participants %>%
        # calculate BMI from weight and height of each participant
        mutate(BMI = round(weight/((height/100)^2), digits = 1)) %>%
        # Join with COPD samples data
        inner_join(trainomeMetaData::copd_samples, by = c("study", "participant", "sex", "treatment")) %>%
        # remove the unneeded columns
        dplyr::select(-c(sample.weight, total.rna, rqi)) %>%
        
        # select only pre-exercise and postexercise data
        filter(time == "PreExc" | time == "PostExc") %>%
        # select only those with sequenced samples
      #  drop_na(seq_sample_id) %>%
        # extract the Vastuls lateralis measurements
        inner_join(trainomeMetaData::copd_thickness %>%
                           dplyr::select(-c(RFthick, Qthick)), by = c("study", "participant", "sex", "leg", "treatment", "condition",
                                                                      "time"))%>%
        dplyr::select(-c(diagnosis, treatment, height, weight))

 
                                                                      

# extract the pre and post exercise
copd_pre <- copd_metadata %>%
        filter(time == "PreExc")

copd_post <- copd_metadata %>%
        filter(time == "PostExc")


# merge them
copd_merged <- copd_pre %>%
        inner_join(copd_post, by = c("study", "participant", "leg", "age" , "sex", 
                                     "BMI", "condition" ))%>%
        rename(seq_sample_id = seq_sample_id.x,
               seq_sample_id_post = seq_sample_id.y) %>%
        # extract the percentage change of the values of interest
        mutate(pct_muscle_change = ((VLthick.y - VLthick.x) / VLthick.x) * 100) %>%
        mutate_if(is.numeric,round, digits = 2) %>%
        dplyr::select("study", "participant", "sex", "seq_sample_id", 
                      "condition",  "age", "BMI", "leg" ,  "pct_muscle_change") %>%
        group_by(participant) %>%
        # carry the last observed value
        mutate(across(c( seq_sample_id), ~ ifelse(is.na(.), zoo::na.locf(., na.rm= F,
                                                                                   fromLast = F), .))) %>%
        # carry next observed value
        mutate(across(c(  seq_sample_id), ~ ifelse(is.na(.), zoo::na.locf(., na.rm = FALSE, fromLast = TRUE), .)))%>%
        drop_na() %>%
        ungroup()%>%
        mutate_if(is.numeric,round, digits = 2) 


# specify the volume
copd_merged["volume"] <- 3

hist(copd_merged$age)
hist(copd_merged$BMI)
hist(copd_merged$pct_muscle_change)


# remove duplicate samples
copd_merged <- copd_merged[!duplicated(copd_merged),]

saveRDS(copd_merged, "data/COPD/copd_metadata.RDS")

length(unique(copd_merged$participant))
