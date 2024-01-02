library(ontologyIndex)
library(plotly)
library(dplyr)
library(stringr)
hpo <- get_OBO("../data/hp.obo")
data_tf <- read.csv2("../output/save/data_tf_mat_tnn.csv", na = "")
data_non_hpo <- data_tf[,!substr(names(data_tf), 0, 2) == "HP"] %>% mutate_at("ID_record_id", as.integer)
data_hpo <- data_tf[,substr(names(data_tf), 0, 2) == "HP" | names(data_tf) == "ID_record_id"] %>% mutate_at("ID_record_id", as.integer)
names(data_hpo) <- str_replace(names(data_hpo), "\\.", "\\:")
names <- names(data_hpo)
for(i in 2:length(names))
{
  if(length(tm <- get_term_property(hpo,  "name", names[i])) == 1)
  {
    names[i] <- tm[1]
  }
}

names(data_hpo) <- names

# Variable response
data_non_hpo$RS_result[data_non_hpo$ID_record_id %in% c(21,31,41,70,72,156,175,209,218,517,519,520,553,578,601,693,997,1218,1360,1462)] <- T
data_non_hpo$RS_result[!data_non_hpo$ID_record_id %in% c(21,31,41,70,72,156,175,209,218,517,519,520,553,578,601,693,997,1218,1360,1462)] <- F

genotypes <- data_tnn[,c("record_id", "gene_name")] %>%
  rename(ID_record_id = record_id) %>% 
  mutate_at("ID_record_id", as.integer) %>% 
  filter(!duplicated(ID_record_id)) %>%
  inner_join(data_non_hpo, by = "ID_record_id")

data_non_hpo <- genotypes
rm(genotypes)

data_non_hpo$gene_name[data_non_hpo$gene_name == ""] <- NA
data_non_hpo$gene_name[!data_non_hpo$RS_result] <- NA
data_non_hpo <- data_non_hpo %>% rename(Result = gene_name) %>% select(-RS_result, -ID_center)

data_non_hpo <- data_non_hpo %>%
  inner_join(data_hpo, by = "ID_record_id")

data <- data_non_hpo[,-1]
data$CS_age <- data$CS_age < 35
data$CS_stage <- data$`Stage 4 chronic kidney disease` | data$`Stage 5 chronic kidney disease`
data <- data %>% 
  select(-Retinopathy) %>% 
  select(-CS_origins_asiatique, -CS_origins_caucasian, -CS_origins_subsahafricaorcaribea, -CS_origins_northafrican, -CS_origins_otherspecified) %>%
  select(-`Stage 1 chronic kidney disease`, -`Stage 2 chronic kidney disease`, -`Stage 3 chronic kidney disease`, -`Stage 4 chronic kidney disease`, -`Stage 5 chronic kidney disease`, -`Status post organ transplantation`) %>%
  select(-CS_prev_test)
data <- data %>% 
  rename('ES before 35 yo' = CS_age, 'Male gender' = CS_gender_male, 'Consanguinity' = CS_consanguinite, 'First degree history' = CS_fam_first, 'Second degree history' = CS_fam_second,
         'Biopsy performed' = CS_prev_biopsy, 'Biological TMA' = CS_tma_bio, 'Histological TMA' = CS_tma_histo,
         'CKD stage IV-V' = CS_stage)
cols <- colSums(data[,-1]) > 1
cols <- c(T, cols)
data <- data[,cols]

write.csv2(data, "../output/data_shiny.csv", row.names = F)
write.csv2(data, "app/data.csv", row.names = F)