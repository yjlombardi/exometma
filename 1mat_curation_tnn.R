library(dplyr)

library(readr)
library(dplyr)
library(stringr)
library(ontologyIndex)
hpo <- get_OBO("../data/hp.obo")

#### SELECTION ####
# Nombre de lignes
data <- data_tnn[,!colSums(is.na(data_tnn) | data_tnn == "") == nrow(data_tnn)]
print(nrow(data)[1])

# Nombre d'entrées
data <- data[is.na(data$redcap_repeat_instance),]
print(nrow(data)[1])

# On sélectionne les entrées avec biomnis_tag (= prélèvement reçu au laboratoire)
data <- data[data$biomnis_tag != "",]
print(nrow(data)[1])

# On sélectionne les entrées avec date de prélèvement disponible
data <- data[data$sample_date_collec != "",]
data <- data[order(data$sample_date_collec),] # Pour être sûr d'enlever les entrées dupliquées chronologiquement parlant
print(nrow(data)[1])

# On selectionne les entrées qui ne concernent pas un donneur
data <- data[data$donor_screened == 0 | is.na(data$donor_screened),]
print(nrow(data)[1])

# On sélectionne les entrées non dupliquées
data <- data[!duplicated(data$ipp_exome),]
data$name <- trimws(toupper(data$name))

data <- data[!duplicated(data[,c("dob_exome", "name")]),]
print(nrow(data)[1])

# On sélectionne uniquement les entrées d'intérêt
inclusion <- read.csv2("../data/tma_inclusion_liaison.csv") %>% mutate_at("v0_rec", as.integer) %>% rename(record_id = v0_rec) 
data <- inner_join(data %>% mutate_at("record_id", as.integer), inclusion, by = "record_id", keep = F)
print(nrow(data)[1])

#### CORRECTION ####
## Des resultats manquants
result_tnn_filled <- read.csv2("../data/result_tnn_filled.csv")
# Version temporaire : on considère T ceux qui sont à 1, et F ceux qui sont à 0 ou NA
data$findings_ok[data$record_id %in% result_tnn_filled[result_tnn_filled$result == 1 & !is.na(result_tnn_filled$result), "record_id"]] <- 1
data$findings_ok[data$record_id %in% result_tnn_filled[result_tnn_filled$result == 0 | is.na(result_tnn_filled$result), "record_id"]] <- 0

#### RENOMMAGE ####
data$ID_record_id <- as.integer(data$record_id) # Record identifyer
data$ID_center <- "TNN" # Center code
data$CS_age <- as.numeric(difftime(data$sample_date_collec, data$dob_exome, units = "days"))/365.25 # Age at ES (Custom)
data$CS_gender_male <- data$sex == 0 # Male gender (Custom)
data$CS_origins_caucasian <- !is.na(data$origins.factor) & data$origins.factor == "Caucasien" # Origine caucasienne (Custom)
data$CS_origins_northafrican <- !is.na(data$origins.factor) & data$origins.factor %in% c("Algérien", "Maroc", "Tunisie") # Origine nord africaine (Custom)
data$CS_origins_subsahafricaorcaribea <- !is.na(data$origins.factor) & data$origins.factor %in% c("Africain", "Antilles") # Origine subsaharienne ou caribéenne (Custom)
data$CS_origins_asiatique <- !is.na(data$origins.factor) & data$origins.factor == "Asiatique" # Origine asiatique (Custom)
data$CS_origins_otherspecified <- !is.na(data$origins.factor) & data$origins.factor %in% c("Turquie", "Autre") # Origine autre (Custom)
data$CS_consanguinite <- !is.na(data$consanguinity) & data$consanguinity == 1
data$HP_0011462 <- !is.na(data$inclusion_criteria_exomic___0) & data$inclusion_criteria_exomic___0 == 1 | 
  (!is.na(data$inclusion_criteria_exomic___3) & data$inclusion_criteria_exomic___3 == 1) | (!is.na(data$facteur_prescp_2022___0) & data$facteur_prescp_2022___0 == 1) | 
  (!is.na(data$fmg2025) & data$fmg2025 == 1) |   (!is.na(data$factor_ex_test___0) & data$factor_ex_test___0 == 1) | data$CS_age < 35 # Young adult onset
data$HP_0000790 <- data$avoid_biopsy___0 == 1 | data$avoid_biopsy___2 == 1 # Hematuria
data$HP_0000093 <- data$avoid_biopsy___1 == 1 | data$avoid_biopsy___2 == 1 # Proteinuria
data$HP_0012623 <- data$ckd_stages.factor == "Stade I, pas dinsuffisance rénale" # Stage I CKD
data$HP_0012624 <- data$ckd_stages.factor == "Stade II, DFG > 60 ml/min [HP:0012624]" # Stage II CKD
data$HP_0012625 <- data$ckd_stages.factor == "Stade III, DFG 60-30 ml/min [HP:0012625]" # Stage III CKD
data$HP_0012626 <- data$ckd_stages.factor == "Stade IV, DFG  30- 15 ml/min [HP:0012626]" # Stage IV CKD
data$HP_0003774 <- data$ckd_stages.factor %in% c("Stade V non dialysé/greffé [HP:0003774]", "Stade V dialysé [HP:0003774]", "Patient transplanté [HP:0032444]") # Stage V CKD
data$HP_0032444 <- data$ckd_stages.factor == "Patient transplanté [HP:0032444]" # Status post organ transplantation
data$HP_0100820 <- !is.na(data$nephropathy_class.factor) & data$nephropathy_class.factor == "Glomérulaire (HP:0031263)" # Glomerulopathy
data$HP_0005576 <- !is.na(data$nephropathy_class.factor) & data$nephropathy_class.factor == "Tubulo-interstitielle fibrosante (HP:0005576)" # Tubulointerstitial fibrosis
data$HP_0009741 <- !is.na(data$nephropathy_class.factor) & data$nephropathy_class.factor == "Vasculaire (HP:0009741)" # Nephrosclerosis
data$HP_0000107 <- (!is.na(data$nephropathy_class.factor) & data$nephropathy_class.factor == "Maladie kystique ou trouble du développement rénal (HP:0000107)") |
  (!is.na(data$renal_cysts) & data$renal_cysts == 1) # Renal cyst
data$HP_0000124 <- !is.na(data$nephropathy_class.factor) & data$nephropathy_class.factor == "Tubulopathie sans I.Rénale (bartter/giltelman/Acidoses...) (HP:0000124)" # Renal tubular dysfunction
data$HP_0000085 <- !is.na(data$horseshoe) & data$horseshoe == 1 # Horseshoe kidney
data$HP_0012581 <- !is.na(data$isolated_cysts) & data$isolated_cysts == 1 # Simple renal cyst
data$HP_0005562 <- !is.na(data$cyst_bilat) & data$cyst_bilat == 1 # Multiple renal cyst
data$HP_0012586 <- !is.na(data$small_kidney) & data$small_kidney == 1 | 
  (!is.na(data$left_kidney_lenght) & data$left_kidney_lenght == 0 & !is.na(data$right_kidney_lenght) & data$right_kidney_lenght == 0) # Bilateral renal atrophy
data$HP_0000113 <- !is.na(data$left_kidney_mx) & data$left_kidney_mx %in% c(2,3) & !is.na(data$right_kidney_mx) & data$right_kidney_mx %in% c(2,3) # Polycystic kidney dysplasia
data$HP_0008717 <- (!is.na(data$left_kidney_lenght) & data$left_kidney_lenght == 0 & !is.na(data$right_kidney_lenght) & data$right_kidney_lenght != 0) |
  (!is.na(data$left_kidney_lenght) & data$left_kidney_lenght != 0 & !is.na(data$right_kidney_lenght) & data$right_kidney_lenght == 0) # Unilateral renal atrophy
data$HP_0001407 <- !is.na(data$hepatic_cysts) & data$hepatic_cysts == 1 # Hepatic cysts
data$HP_0000822 <- !is.na(data$hypertension) & data$hypertension == 1 # Hypertension
data$HP_0001997 <- !is.na(data$gout) & data$gout == 1 # Gout
data$CS_early_gout <- !is.na(data$gout) & data$gout == 1 # Gout with early onset (Custom)
data$HP_0002907 <- !is.na(data$hematuria_micro) & data$hematuria_micro == 1 # Microscopic hematuria
data$CS_prev_test <- !is.na(data$genetic_test_before) & data$genetic_test_before == 1 # Previous genetic test (Custom)
data$CS_prev_biopsy <- !is.na(data$renal_biopsy) & data$renal_biopsy == 1 # Previous renal biopsy (Custom)
data$HP_0000097 <- !is.na(data$pbr_first_findings___0) & data$pbr_first_findings___0 == 1 # Focal segmental glomerulosclerosis
data$HP_0033889 <- !is.na(data$pbr_first_findings___1) & data$pbr_first_findings___1 == 1 # Abnormal renal arteriole morphology
data$HP_0032599 <- !is.na(data$pbr_first_findings___2) & data$pbr_first_findings___2 == 1 # Abnormal renal tubular epithelial morphology
data$HP_0001937 <- !is.na(data$pbr_first_findings___3) & data$pbr_first_findings___3 == 1 |
  (!is.na(data$tma.x) & data$tma.x == 1) | (!is.na(data$tma_current) & data$tma_current == 1) # Microangiopathic hemolytic anemia
data$HP_0030949 <- !is.na(data$pbr_first_findings___4) & data$pbr_first_findings___4 == 1 # Glomerular deposits
data$HP_0000794 <- !is.na(data$pbr_first_findings___5) & data$pbr_first_findings___5 == 1 # IgA deposition in the glomerulus
data$HP_0012577 <- !is.na(data$pbr_first_findings___6) & data$pbr_first_findings___6 == 1 # Thin glomerular basement membrane
data$HP_0032583 <- !is.na(data$pbr_first_findings___8) & data$pbr_first_findings___8 == 1 # Renal glomerular foam cells
data$HP_0010935 <- !is.na(data$cakut_suspect) & data$cakut_suspect == 1 # Abnormality of the upper urinary tract
data$HP_0000076 <- !is.na(data$vur) & data$vur == 1 # Vesicoureteral reflux
data$HP_0000787 <- !is.na(data$lithiasis_suspected) & data$lithiasis_suspected == 1 # Nephrolithiasis
data$HP_0000707 <- !is.na(data$atcd_neurological) & data$atcd_neurological == 1 # Abnormality of the nervous system
data$HP_0002960 <- !is.na(data$past_hist_aid) & data$past_hist_aid == 1 # Autoimmunity
data$CS_fam_first <- data$facteur_prescp_2022___1 == 1 | (!is.na(data$first_degree_np) & data$first_degree_np == 1) |
  (!is.na(data$factor_ex_test___1) & data$factor_ex_test___1 == 1) # First degree history (Custom)
data$CS_fam_mother <- !is.na(data$first_degree_ascent___0) & data$first_degree_ascent___0 == 1 # First degree history: mother (Custom)
data$CS_fam_father <- !is.na(data$first_degree_ascent___1) & data$first_degree_ascent___1 == 1 # First degree history: father (Custom)
data$CS_fam_brother <- !is.na(data$first_degree_ascent___2) & data$first_degree_ascent___2 == 1 # First degree history: brother (Custom)
data$CS_fam_sister <- !is.na(data$first_degree_ascent___3) & data$first_degree_ascent___3 == 1 # First degree history: sister (Custom)
data$CS_fam_child <- !is.na(data$first_degree_ascent___4) & data$first_degree_ascent___4 == 1 # First degree history: child (Custom)
data$CS_fam_second <- !is.na(data$nd_pedigree) & data$nd_pedigree == 1 # Second degree history (Custom)
data$HP_0032316 <- data$inclusion_criteria_exomic___1 == 1 | data$inclusion_criteria_exomic___3 == 1 | data$CS_fam_first | data$CS_fam_second # Family history
data$RS_result <- data$findings_ok == 1 # Le résultat
data_struct <- data[,substr(names(data), 0, 2) %in% c("ID", "CS", "HP", "RS")]

#### CHARGEMENT DES CODES ISSUS DES DONNES NON STRUCTUREES ####
data_unstruct_source <- read.csv2("../data/data_tnn_hpos.csv", na.strings = "")
data_unstruct <- data.frame(ID_record_id = as.integer(data_unstruct_source$record_id))
# HPOS
for(i in 1:nrow(data_unstruct_source))
{
  currow <- na.omit(as.character(data_unstruct_source[i,paste("hpo", seq(1:20), sep = "")]))
  if((l = length(currow)) == 0)
  {
    next
  }
  for(j in 1:l)
  {
    curname <- paste(str_replace(currow[j], "\\:", "_"), sep = "")
    data_unstruct[i,curname] <- T
  }
}
data_unstruct[is.na(data_unstruct)] <- F

#### MISE EN COMMUN ####
join <- left_join(data_struct, data_unstruct, by = "ID_record_id", suffix = c("_ST", "_US"))

## 0 : phénotype absent | 1 : structuré seulement | 2 : non strcuturé seulement | 3 : dans les deux
# Colonnes disponibles dans les 2 sources
to_adj <- names(join)[grepl("ST", names(join), fixed = TRUE)]
if(length(to_adj) > 0)
{
  for(i in 1:length(to_adj))
  {
    curname <- substring(to_adj[i], 0, 10)
    curst <- to_adj[i]
    curus <- paste(curname, "_US", sep = "")
    join[,curname] <- as.integer(ifelse(join[,curst] == 0 & join[,curus] == 0, 0,
                                        ifelse(join[,curst] == 1 & join[,curus] == 0, 1,
                                               ifelse(join[,curst] == 0 & join[,curus] == 1, 2, 3))))
  }
}

# Colonnes disponibles uniquement dans 1 source
cols_struct <- names(data_struct)[!names(data_struct) %in% names(data_unstruct)]
cols_unstruct <- names(data_unstruct)[!names(data_unstruct) %in% names(data_struct)]
cols_struct <- cols_struct[substr(cols_struct, 0, 2) == 'HP']
cols_unstruct <- cols_unstruct[substr(cols_unstruct, 0, 2) == 'HP']
join[,cols_struct] <- join[,cols_struct] %>% mutate_all(function(x){as.integer(ifelse(x, 1, 0))})
join[,cols_unstruct] <- join[,cols_unstruct] %>% mutate_all(function(x){as.integer(ifelse(x, 2, 0))})
to_rm <- which(grepl( "ST", names(join), fixed = TRUE) | grepl( "US", names(join), fixed = TRUE))
if(length(to_rm) > 0)
{
  data_final_tnn <- join[,-to_rm]
} else {
  data_final_tnn <- join
}

# Restriction de la population
inclusion <- read.csv2("../data/tma_inclusion_liaison.csv") %>% 
  mutate_at("v0_rec", as.integer) %>% 
  rename(ID_record_id = v0_rec) %>%
  filter(tma != '0:No')

data_final_tnn <- inner_join(data_final_tnn, inclusion, by = "ID_record_id")
sums <- colSums(data_final_tnn[,substr(names(data_final_tnn), 0, 2) == 'HP'])
names <- names(sums)[sums == 0]
data_final_tnn <- data_final_tnn[,!names(data_final_tnn) %in% names]
data_final_tnn[,-1] <- data_final_tnn[,-1] %>%
  mutate_if(is.integer, function(x){x > 0})
data_final_tnn$tma <- as.factor(data_final_tnn$tma)

#### NETTOYAGE ####
rm(data)
rm(join)
rm(data_struct)
rm(data_unstruct)
rm(data_unstruct_source)
rm(result_known_tnn_2_empty)
rm(result_known_tnn_2_filled)
rm(result_known_tnn_empty)
rm(result_known_tnn_filled)
rm(result_tnn_empty)
rm(result_tnn_filled)
rm(c_good)
rm(curname)
rm(currow)
rm(i)
rm(j)
rm(l)
rm(ipp_app)
rm(to_adj)
rm(to_rm)
rm(curst)
rm(curus)
gc()


