library(trainomeMetaData)
library(dplyr)
library(tidyverse)


# volume metadata


Vol_thick <- trainomeMetaData::vol_mr 
# Change to match full dataset
Vol_thick["time"][Vol_thick["time"] == "pre"] <- "w0"
Vol_thick["time"][Vol_thick["time"] == "post"] <- "w12" 
Vol_thick_pre <- Vol_thick %>%
        filter(time == "w0")
Vol_thick_post <- Vol_thick %>%
        filter(time == "w12")

Vol_thick_merged <- Vol_thick_pre %>%
        inner_join(Vol_thick_post, by = c("study", "participant", "leg", "condition")) %>%
        mutate(pct_muscle_change = ((csa.y - csa.x) / csa.x) * 100)  %>%
        dplyr::select(study, participant, leg, condition, pct_muscle_change)
        

Vol_metadata <- trainomeMetaData::vol_samples %>%
       inner_join(trainomeMetaData::vol_participants, by = c("study", "participant", "sex")) %>%
        # calculate BMI for each participant
        mutate(BMI = round(weight/((height/100)^2), digits = 1)) %>%
       dplyr::select(study, participant, sex, seq_sample_id,  condition, time, age, BMI, leg) %>%
        # select only pre and post exercise data
        subset(time == "w0") %>%
        inner_join(Vol_thick_merged, by = c("study", "participant", "leg", "condition")) %>%
        mutate(volume = ifelse(condition == "single", 1, 3)) %>%
        group_by(participant) %>%
        # carry the last observed value
        mutate(across(c(sex, seq_sample_id, time, age, BMI), ~ ifelse(is.na(.), zoo::na.locf(., na.rm= F,
                                                                                        fromLast = F), .))) %>%
        # carry next observed value
        mutate(across(c( sex, seq_sample_id,time, age, BMI), ~ ifelse(is.na(.), zoo::na.locf(., na.rm = FALSE, fromLast = TRUE), .)))%>%
        drop_na() %>%
        ungroup()%>%
        mutate_if(is.numeric,round, digits = 2) %>%
        dplyr::select("study", "participant", "sex", "seq_sample_id",
                      "condition",  "age", "BMI", "leg" ,   "volume" , "pct_muscle_change" ) 
        
#change the time variables to PreExc and PostExc
# Vol_metadata["time"][Vol_metadata["time"] == "w0"] <- "PreExc"
#Vol_metadata["time"][Vol_metadata["time"] == "w12"] <- "PostExc"
# Make the conditions all RM10 as that is what they are
         
Vol_metadata$condition <- "RM10"    

# remove duplicate samples
Vol_metadata <- Vol_metadata[!duplicated(Vol_metadata),]


hist(Vol_metadata$age)

length(unique(Vol_metadata$participant))

saveRDS(Vol_metadata, "data/Volume/Volume_metadata.RDS")
