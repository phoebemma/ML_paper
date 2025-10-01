library(trainomeMetaData)
library(dplyr)
library(tidyverse)
library(contratraindata)

contra_thickness <- contratraindata::ct_thickness %>%
        filter(time == "t1" | time == "t4") %>%
        pivot_wider(names_from = c(leg, time), values_from = thickness)%>%
        mutate(pct_muscle_change_R = ((R_t4 - R_t1) / R_t1) * 100,
               pct_muscle_change_L = ((L_t4 - L_t1) / L_t1) * 100) %>%
        pivot_longer(cols = starts_with("pct"), names_to = "leg", values_to = "pct_muscle_change")%>%
        # extract the right and left legs into a column
        mutate(leg = case_when(str_detect(leg, "_L") ~ "L",
                               str_detect(leg, "_R") ~ "R")) %>%
        dplyr::select(participant, leg, pct_muscle_change)%>%
        # rename participant to match the names from trainometadata
        mutate(participant = gsub("ID", "FP", participant)) %>%
        inner_join(trainomeMetaData::ct_samples, by = c("participant", "leg")) %>%
        filter(time == "t1") %>%
        dplyr::select(study, participant,  sex, leg, pct_muscle_change, condition, seq_sample_id ) 

#hist(contra_thickness$pct_muscle_change)


# Extract participant informtaion and merge to data we have
contra_participants <- contratraindata::ct_participants %>%
        mutate(BMI = round(weight/((height/100)^2), digits = 1)) %>%
        dplyr::select(participant, group, age, sex, BMI) %>%
        # inner_join(contra_thickness, by = "participant") %>%
        # # select the pre and post exercise data

        inner_join(contratraindata::ct_seqsamples, by = c("participant")) %>%
        dplyr::filter(time == "t1" ) %>%
        # rename participant to match the names from trainometadata
        mutate(participant = gsub("ID", "FP", participant)) %>%
        full_join(contra_thickness, by = c("participant",  "sex", "leg")) %>%
        rename(seq_sample_id = seq_sample_id.x) %>%
          dplyr::select(study, participant, age, sex, BMI, leg,pct_muscle_change,condition, seq_sample_id ) %>%
        # create new column called volume from the condition column
        # this is done to match the volume and condition from the other datasets
        mutate(volume = case_when(condition == "set0" ~ 0,
                                  condition == "set6" ~ 6,
                                  condition == "set3" ~ 3)) %>%
        group_by(participant) %>%
        # carry the last observed value
        mutate(across(c( seq_sample_id, age, BMI), ~ ifelse(is.na(.), zoo::na.locf(., na.rm= F,
                                                                                             fromLast = F), .))) %>%
        # carry next observed value
        mutate(across(c(  seq_sample_id, age, BMI), ~ ifelse(is.na(.), zoo::na.locf(., na.rm = FALSE, fromLast = TRUE), .)))%>%
        drop_na() %>%
        ungroup()%>%
        mutate_if(is.numeric,round, digits = 2) 
       
# remove duplicate samples
contra_participants <- contra_participants[!duplicated(contra_participants),] 

contra_participants$condition <- "RM10"
hist(contra_participants$BMI)

hist(contra_participants$pct_muscle_change)

length(unique(contra_participants$participant))

saveRDS(contra_participants, "data/ContraTrain/Contratrain_metadata.RDS")
