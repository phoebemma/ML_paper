source("R/Trainome_functions.R")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(trainomeMetaData)
library(edgeR)

alpha_omega_counts <- extract_rsem_gene_counts("data/gene_counts/Alpha_Omega_RSEM_outputs/")
# the metadata sequence sample ids contain only numbers, match the gene expression to it
# 
# names(alpha_omega_counts)[-1] <- gsub("[^0-9]", "", names(alpha_omega_counts)[-1])


# load the alpha omega metadata 
alpha_omega_Pre_exc_metadata <- readRDS("data/Alpha_Omega/Alpha_Omega_metadata.RDS") %>%
        dplyr::select(study, participant, age, sex, leg, seq_sample_id, BMI, pct_muscle_change )

alpha_omega <- alpha_omega_counts %>% 
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, all_of(alpha_omega_Pre_exc_metadata$seq_sample_id))



## Keep nonzero rows
nonzero <- alpha_omega %>%
  dplyr::filter(rowSums(alpha_omega[, -1]) != 0)


nonzero <- nonzero %>%
        dplyr::filter(filterByExpr(nonzero[,-1], min.total.count = 126)) # at least, one and half counts per participant




# Create dge lists and calculate norm factors
dge_ao   <- DGEList(nonzero[,-1])


dge_ao  <- calcNormFactors(dge_ao, method = "TMM")


# get the normalised counts

normalised_ao_counts <- as.data.frame(cpm(dge_ao, normalized.lib.sizes = T))

normalised_ao_counts$gene_name <- nonzero$gene_name

# Make the transcript_ID the first column

normalised_ao_counts <- normalised_ao_counts %>%
  dplyr::select(gene_name, everything())


saveRDS(normalised_ao_counts, "data/gene_counts/cpm_normalised_PreExc_AO_counts.RDS")






# Volume
Vol_counts <- extract_rsem_gene_counts("data/gene_counts/Volume_RSEM_outputs_new/")
#remove everything before the . in sample_id. 
colnames(Vol_counts) <- gsub(".*?\\.", "", colnames(Vol_counts) )


# load the volume metadata 
volume_metadata <- readRDS("data/Volume/Volume_metadata.RDS") %>%
        dplyr::select(study, participant, age, sex, leg, seq_sample_id, BMI, pct_muscle_change )



volume <- Vol_counts %>% 
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, all_of(volume_metadata$seq_sample_id))


## Keep nonzero rows
nonzero_vol <- volume %>%
  dplyr::filter(rowSums(volume[,-1]) != 0)


nonzero_vol <- nonzero_vol %>%
        dplyr::filter(filterByExpr(nonzero_vol[,-1], min.total.count = 75)) # at least, one and half counts per participant

# Create dge lists and calculate norm factors
dge_vol   <- DGEList(nonzero_vol[,-1])


dge_vol  <- calcNormFactors(dge_vol, method = "TMM")


# get the normalised counts

normalised_vol_counts <- as.data.frame(cpm(dge_vol, normalized.lib.sizes = T))

normalised_vol_counts$gene_name <- nonzero_vol$gene_name

# Make the transcript_ID the first column

normalised_vol_counts <- normalised_vol_counts %>%
  dplyr::select(gene_name, everything())


saveRDS(normalised_vol_counts, "data/gene_counts/cpm_normalised_PreExc_vol_counts.RDS")



# Contratrain data

Ct_counts <- extract_rsem_gene_counts("data/gene_counts/Contratrain_RSEM_outputs_new/")


# load the volume metadata 
ct_metadata <- readRDS("data/ContraTrain/Contratrain_metadata.RDS") %>%
        dplyr::select(study, participant, age, sex, leg, seq_sample_id, BMI, pct_muscle_change )




Contratrain <- Ct_counts %>% 
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, all_of(ct_metadata$seq_sample_id))


## Keep nonzero rows
nonzero_ct <- Contratrain %>%
  dplyr::filter(rowSums(Contratrain[,-1]) != 0)


nonzero_ct <- nonzero_ct %>%
        dplyr::filter(filterByExpr(nonzero_ct[,-1], min.total.count = 96))  # at least, one and half counts per participant


# Create dge lists and calculate norm factors
dge_ct   <- DGEList(nonzero_ct[,-1])


dge_ct  <- calcNormFactors(dge_ct, method = "TMM")


# get the normalised counts

normalised_ct_counts <- as.data.frame(cpm(dge_ct, normalized.lib.sizes = T))

normalised_ct_counts$gene_name <- nonzero_ct$gene_name

# Make the transcript_ID the first column

normalised_ct_counts <- normalised_ct_counts %>%
  dplyr::select(gene_name, everything())


saveRDS(normalised_ct_counts, "data/gene_counts/cpm_normalised_PreExc_contratrain_counts.RDS")




# COPD


Copd_counts <- extract_rsem_gene_counts("data/gene_counts/COPD_RSEM_outputs_new/")

# load COPD metadata
copd_metadata <- readRDS("data/COPD/copd_metadata.RDS") %>%
        dplyr::select(study, participant, age, sex, leg, seq_sample_id, BMI, pct_muscle_change )


Copd <- Copd_counts %>% 
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, all_of(copd_metadata$seq_sample_id))


## Keep nonzero rows
nonzero_copd <- Copd %>%
  dplyr::filter(rowSums(Copd[,-1]) != 0)



nonzero_copd <- nonzero_copd %>%
        dplyr::filter(filterByExpr(nonzero_copd[,-1], min.total.count = 155)) 


# Create dge lists and calculate norm factors
dge_copd   <- DGEList(nonzero_copd[,-1])


dge_copd  <- calcNormFactors(dge_copd, method = "TMM")


# get the normalised counts

normalised_copd_counts <- as.data.frame(cpm(dge_copd, normalized.lib.sizes = T))

normalised_copd_counts$gene_name <- nonzero_copd$gene_name

# Make the transcript_ID the first column

normalised_copd_counts <- normalised_copd_counts %>%
  dplyr::select(gene_name, everything())


saveRDS(normalised_copd_counts, "data/gene_counts/cpm_normalised_PreExc_COPD_counts.RDS")







# Relief

Relief_counts <- extract_rsem_gene_counts("data/gene_counts/Relief_RSEM_outputs/")


# Load the Relief metadata
Relief_metadata <- readRDS("data/Relief/Relief_metadata.RDS") %>%
        dplyr::select(study, participant, age, sex, leg, seq_sample_id, BMI, pct_muscle_change )


Relief <- Relief_counts %>% 
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, all_of(Relief_metadata$seq_sample_id))


## Keep nonzero rows
nonzero_Relief <- Relief %>%
  dplyr::filter(rowSums(Relief[,-1]) != 0)



nonzero_Relief <- nonzero_Relief %>%
        dplyr::filter(filterByExpr(nonzero_Relief[,-1], min.total.count = 117)) 

# Create dge lists and calculate norm factors
dge_Relief   <- DGEList(nonzero_Relief[,-1])


dge_Relief  <- calcNormFactors(dge_Relief, method = "TMM")


# get the normalised counts

normalised_Relief_counts <- as.data.frame(cpm(dge_Relief, normalized.lib.sizes = T))

normalised_Relief_counts$gene_name <- nonzero_Relief$gene_name

# Make the transcript_ID the first column

normalised_Relief_counts <- normalised_Relief_counts %>%
  dplyr::select(gene_name, everything())


saveRDS(normalised_Relief_counts, "data/gene_counts/cpm_normalised_PreExc_Relief_counts.RDS")






###### Batch effects correction
# Merge the normalised gene counts into one

all_normalised_transcripts <- normalised_copd_counts %>%
  full_join(normalised_vol_counts, by = "gene_name") %>%
  full_join(normalised_ct_counts, by = "gene_name") %>%
  full_join(normalised_ao_counts, by = "gene_name") %>%
  full_join(normalised_Relief_counts, by = "gene_name") %>%
  drop_na() 




# merge the metadata
all_metadata <- rbind(copd_metadata, volume_metadata ) %>%
        rbind(ct_metadata) %>%
        rbind(alpha_omega_Pre_exc_metadata) %>%
        rbind(Relief_metadata) %>%
        mutate(across(c("age"), round, 0),
               participant = paste0(study,"_", participant ))

length(unique(all_metadata$participant))
hist(all_metadata$age)
hist(all_metadata$BMI)
hist(all_metadata$pct_muscle_change)



# Some of the sequence samples are in duplicates
# These is from studies that took only one sample at baseline, but two at postexc


meta_df <- all_metadata %>%
        distinct(seq_sample_id, .keep_all = T)
# Check if everything matches except the transcript_id
match(colnames(all_normalised_transcripts), meta_df$seq_sample_id)

all_transcripts_reordered <- all_normalised_transcripts[ , c("gene_name",  meta_df$seq_sample_id)] 

all_transcripts_reordered[,-1] <- round(all_transcripts_reordered[,-1], 2)

# recheck# recheckround()
match(colnames(all_transcripts_reordered),  meta_df$seq_sample_id)


# Define the batch, which in this case is the styd
batch <- factor(meta_df$study)

# correct batch effect
corrected_counts <- removeBatchEffect(all_transcripts_reordered[,-1], batch = batch)
corrected_counts$gene_name <- unlist(all_transcripts_reordered$gene_name)

#boxplot(corrected_counts, main = "batch corrected data")

saveRDS(corrected_counts, "data/gene_counts/batch_corrected_RSEM_counts")

