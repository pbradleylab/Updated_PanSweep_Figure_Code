#Determine Genes to BLAST:
library(readr)
library(tidyverse)


#Paths:
Path_tables <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/PanSweep_Results/Most_Abd/PanSweep_Analysis_Output_2025-01-29/PanSweep_Analysis_TablesOnly.rds"
Path_base <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/BLAST/"
################################################################################
#Filter out genes for blasting form B theta '100196'
Tables <- read_rds(Path_tables)
eggNOG <- Tables$uhgp_90_eggNOG
Bt_genes <- eggNOG %>% filter(Species_id == "100196")

Not_lineage <- c("Erysipelotrichia", "Negativicutes", "unclassified Lachnospiraceae", "Bacilli", "Escherichia", "Ruminococcaceae",
                 "Lachnoclostridium", "Firmicutes", "Clostridia", "Blautia", "Oscillospiraceae", "Sutterellaceae", "Planctomycetes",
                 "Eubacteriaceae", "Paenibacillaceae", "Streptococcus oralis", "Dorea", "Clostridiaceae")

#Select 5 genes from each condition:
set.seed(42)
F_NA <- Bt_genes %>% filter(is.na(Predicted_taxonomic_group), Lineage_Shared == FALSE) %>% sample_n(5)
T_NA <- Bt_genes %>% filter(is.na(Predicted_taxonomic_group), Lineage_Shared == TRUE) %>% sample_n(5)
F_B <- Bt_genes %>% filter(Predicted_taxonomic_group == "Bacteroidaceae", Lineage_Shared == FALSE) %>% sample_n(5)
T_B <- Bt_genes %>% filter(Predicted_taxonomic_group == "Bacteroidaceae", Lineage_Shared == TRUE) %>% sample_n(5)
T_Not_Lin <- Bt_genes %>% filter(Predicted_taxonomic_group %in% Not_lineage, Lineage_Shared == TRUE) %>% sample_n(5)
F_Not_Lin <- Bt_genes %>% filter(Predicted_taxonomic_group %in% Not_lineage, Lineage_Shared == FALSE) %>% sample_n(5)

Blast_Genes <- rbind(F_NA, T_NA, F_B, T_B, T_Not_Lin, F_Not_Lin) 

write_tsv(Blast_Genes, paste0(Path_base, "eggNOG_Genes_to_BLAST.tsv"))
write_csv(Blast_Genes, paste0(Path_base, "eggNOG_Genes_to_BLAST.csv"))
################################################################################
#
Blast_Genes <- read_tsv(paste0(Path_base, "eggNOG_Genes_to_BLAST.tsv"))








