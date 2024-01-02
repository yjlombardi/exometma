library(ontologyIndex)
library(dplyr)
hpo <- get_OBO("../data/hp.obo")

#### LOAD AND TRANSFORM DATA ####
load_and_transform <- function(data_final)
{  
  data_nonhpo <- data_final[,substr(names(data_final), 0, 2) != "HP"]
  data_hpo <- data_final[,substr(names(data_final), 0, 2) == "HP"]
  
  # Rename or remove some ambiguous HPO terms
  # Cakut 
  data_hpo <- data_hpo %>% rename(HP_0000110 = HP_0010935)
  
  # Abnormality of the liver
  data_hpo$HP_0001395 <- data_hpo$HP_0001395 | data_hpo$HP_0001392
  data_hpo <- data_hpo %>% select(-HP_0001392)
  
  # Retinopathy & hypertensive retinopathy
  data_hpo$HP_0001095 <- data_hpo$HP_0001095 | (data_hpo$HP_0000488 & data_hpo$HP_0000822)
  data_hpo$HP_0000488[data_hpo$HP_0001095] <- F
  
  # Hypertension
  data_hpo <- data_hpo %>% select(-HP_0000822)
  
  # Nephrolithiasis
  data_hpo$HP_0000787 <- data_hpo$HP_0000787 | data_hpo$HP_0034368
  data_hpo <- data_hpo %>% select(-HP_0034368)
  
  # Polycystic kidney dysplasia
  data_hpo$HP_0000113 <- data_hpo$HP_0000113 | data_hpo$HP_0005562
  data_hpo <- data_hpo %>% select(-HP_0005562)
  
  # Simple renal cyst (by approximation: if renal cyst or simple renal cyst, and not multiple and not polycystic)
  data_hpo$HP_0012581 <- data_hpo$HP_0000107 | data_hpo$HP_0012581
  data_hpo$HP_0012581[data_hpo$HP_0000113] <- F
  data_hpo <- data_hpo %>% select(-HP_0000107)
  
  # Left ventricular hypertrophy
  data_hpo$HP_0001712 <- data_hpo$HP_0001712 | data_hpo$HP_0001639
  data_hpo <- data_hpo %>% select(-HP_0001639)
  
  # Young adult onset
  data_hpo <- data_hpo %>% select(-HP_0011462)
  
  # Family history
  data_hpo <- data_hpo %>% select(-HP_0032316)
  
  # Microangiopathic hemolytic anemia
  data_hpo <- data_hpo %>% select(-HP_0001937)
  
  # Abnormal renal arteriole morphology
  data_hpo <- data_hpo %>% select(-HP_0033889)
  
  # Nephrosclerosis
  data_hpo <- data_hpo %>% select(-HP_0009741)
  
  # Bilateral renal atrophy
  data_hpo$HP_0012586[data_hpo$HP_0008717] <- F
  
  names(data_hpo) <- str_replace(names(data_hpo), "_", "\\:")
  names <- names(data_hpo)
  matrix <- get_term_descendancy_matrix(hpo, names)
  dec <- colSums(matrix) == 0
  which <- dec
  unique <- names(which)[which == T]
  unique_nm <- unique
  
  for(i in 1:length(unique))
  {
    if(length(tm <- get_term_property(hpo,  "name", unique[i])) == 1)
    {
      unique_nm[i] <- tm[1]
    }
  }
  
  for(i in 1:length(unique))
  {
    row <- matrix[unique[i],]
    row <- row[row == T]
    if(length(row) == 0)
    {
      next
    }
    else
    {
      data_hpo[,unique[i]] <- rowSums(data_hpo[,c(unique[i], names(row))]) > 0
    }
  }
  
  data_hpo <- data_hpo[,unique]
  
  data_nonhpo$CS_tma_bio <- data_nonhpo$tma %in% c("1:Biological+HistoNA", "2:Biological+Histo-", "4:Biological+Histo+")
  data_nonhpo$CS_tma_histo <- data_nonhpo$tma %in% c("3:Biological-Histo+", "4:Biological+Histo+")
  data_nonhpo <- data_nonhpo %>% select(-tma) %>% select(-CS_early_gout)
  data_nonhpo <- data_nonhpo %>% select(-CS_fam_brother, -CS_fam_sister, -CS_fam_father, -CS_fam_mother, -CS_fam_child)
  data_prepared <- cbind(data_nonhpo, data_hpo)
  
  return(data_prepared)
}

data_tf_tnn <- load_and_transform(data_final_tnn)
gc()
write_csv2(data_tf_tnn, "../output/save/data_tf_mat_tnn.csv", na = "")

#### SAVE AND CLEAN ####
rm(data_tf_tnn)
gc()