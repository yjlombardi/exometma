library(dplyr, warn.conflicts = F)
library(Hmisc, warn.conflicts = F)

hospit <- read.csv2("../data/exome_hospit.csv")
hospit2 <- unique(data_tnn$record_id[data_tnn$dna_during_hospit.factor == 'Yes' & !is.na(data_tnn$dna_during_hospit.factor)])
hospit <- unique(c(hospit$record_id, hospit2))

out <- read.csv("../data/out.csv")
names(out)[1] <- "v0_rec"
out$v0_rec <- as.integer(out$v0_rec)
out <- out[,!names(out) == "patient_num"]
data_tnn_clean$v0_rec <- as.integer(data_tnn_clean$v0_rec)

data_bio <- left_join(data_tnn_clean[data_tnn_clean$v2_dt <= "2022-12-31",], out, by = "v0_rec", keep = FALSE)
label(data_bio$platelets) <- "Platelets, G/L"
label(data_bio$hemoglobin) <- "Hemoglobin, g/dL"
label(data_bio$haptoglobin) <- "Haptoglobin, g/L"
label(data_bio$schizo_quant) <- "Schizocytes, %"
label(data_bio$schizo) <- "Schizocytes, presence"
data_bio$ldh <- as.numeric(unname(data_bio$ldh))
label(data_bio$ldh) <- "Lactate dehydrogenase, U/L"

data_bio$rosner <- data_bio$rosner >= 15
data_bio$anticardiolipine <- data_bio$anticardiolipine >= 15
data_bio$antib2gp1 <- data_bio$antib2gp1 >= 20
data_bio$antib2gp1_qual <- data_bio$antib2gp1_qual == 1
w <- is.na(data_bio$antib2gp1) & is.na(data_bio$antib2gp1_qual)
data_bio$antib2gp1 <- !is.na(data_bio$antib2gp1) & data_bio$antib2gp1 | !is.na(data_bio$antib2gp1_qual) & data_bio$antib2gp1_qual
data_bio$antib2gp1[w] <- NA
data_bio$adamts13 <- data_bio$adamts13 <= 20

label(data_bio$c3) <- "C3, g/l"
label(data_bio$c4) <- "C4, g/l"
label(data_bio$antib2gp1) <- "Anti-beta 2 glycoprotein 1, presence"
label(data_bio$anticardiolipine) <- "Anti-cardiolipin, presence"
label(data_bio$rosner) <- "Pathological rosner index"
label(data_bio$adamts13) <- "Pathological ADAMTS 13 activity"

verif_tma <- read.csv2("../output/verif_mat_maj_short_230606.csv")[,c("v0_rec", "certain_tma_manual", "certain_tmabio_manual")]
data_bio <- left_join(data_bio, verif_tma, by = "v0_rec", keep = F)
rm(verif_tma)

## MAT bio
data_bio$tmabio <- NA

# Ceux qui ont schizo+
w1 <- !is.na(data_bio$schizo) & data_bio$schizo

# Ceux qui ont plq < 100 et hapto < 0.2
w2 <- !is.na(data_bio$platelets) & !is.na(data_bio$haptoglobin) & data_bio$platelets < 100 & data_bio$haptoglobin < 0.2

# Ceux qui ont une réponse positive pour MAT bio
w3 <- data_bio$v28_tma

# Ceux qui ont été vérifiés à la main
w4 <- !is.na(data_bio$certain_tmabio_manual)

# Ceux qui ont été vérifié à la main et sont positifs
w5 <- !is.na(data_bio$certain_tmabio_manual) & data_bio$certain_tmabio_manual

# Adjudication
data_bio$tmabio <- (w3 & (w1 | w2) & !w4) | w5

## MAT histo
data_bio$tmahisto <- NA

# Ceux qui ont histo+
ww1 <- !is.na(data_bio$v17_pbrres4) & data_bio$v17_pbrres4 == 'Presence'

# Ceux qui ont été vérifié à la main et sont négatifs
ww2 <- !is.na(data_bio$certain_tma_manual) & !data_bio$certain_tma_manual

# Ceux qui ont histo-
ww3 <- (!is.na(data_bio$v17_pbrres4) & data_bio$v17_pbrres4 != 'Presence') | is.na(data_bio$v17_pbrres4)

# Adjudication
data_bio$tmahisto <- ww1 & !(ww2 | ww3)
data_bio$tmahisto[data_bio$v17_pbrres4 == 'No biopsy performed'] <- NA

## MAT
data_bio$tma <- NA
data_bio$tma[!data_bio$tmabio & (is.na(data_bio$tmahisto) | !data_bio$tmahisto)] <- "0:No"
data_bio$tma[data_bio$tmabio & is.na(data_bio$tmahisto)] <- "1:Biological+HistoNA"
data_bio$tma[data_bio$tmabio & !is.na(data_bio$tmahisto) & !data_bio$tmahisto] <- "2:Biological+Histo-"
data_bio$tma[!data_bio$tmabio & !is.na(data_bio$tmahisto) & data_bio$tmahisto] <- "3:Biological-Histo+"
data_bio$tma[data_bio$tmabio & !is.na(data_bio$tmahisto) & data_bio$tmahisto] <- "4:Biological+Histo+"
data_bio$tma <- as.factor(data_bio$tma)

data_bio$tmagrp <- data_bio$tma != "0:No"

## Hospit
data_bio$hospit <- data_bio$v0_rec %in% hospit

## Adjudication diag LM
data_bio$v43_ex_pos <- F
data_bio$v43_ex_pos[data_bio$v0_rec %in% c(
  21,
  31,
  41,
  70,
  72,
  156,
  175,
  209,
  218,
  517,
  519,
  520,
  553,
  578,
  601,
  693,
  997,
  1218,
  1360,
  1462
  
)] <- T

## Tables
data_bio[,c("tma", "platelets", "hemoglobin", "ldh", "haptoglobin", "schizo", "schizo_quant")] %>%
  tbl_strata(
    strata = c(tma),
    .tbl_fun =
      ~ .x %>%
      tbl_summary(type = list(where(is.numeric) ~ "continuous2"),

        missing = "no") %>%
      add_n( )
  )
#
data_bio[,c("tma", "platelets", "hemoglobin", "ldh", "haptoglobin", "schizo", "schizo_quant")] %>%
  tbl_summary(by=tma) %>%
  add_p()

data_bio[data_bio$tma != '0:No',c("tma", "platelets", "hemoglobin", "ldh", "haptoglobin", "schizo", "schizo_quant", "hospit")] %>%
  tbl_summary(by=hospit) %>%
  add_p()

data_bio[data_bio$hospit,c("tma", "platelets", "hemoglobin", "ldh", "haptoglobin", "schizo", "schizo_quant")] %>%
  tbl_summary(by=tma) %>%
  add_p()

data_bio[!data_bio$hospit,c("tma", "platelets", "hemoglobin", "ldh", "haptoglobin", "schizo", "schizo_quant")] %>%
  tbl_summary(by=tma) %>%
  add_p()

data_bio$pretest_probability <- as.numeric(data_bio$pretest_probability)
data_bio$v4_age <- as.numeric(data_bio$v4_age)

droplevels(data_bio[data_bio$tma != '0:No',]) %>%
  tbl_summary(by=v43_ex_pos) %>%
  add_p()