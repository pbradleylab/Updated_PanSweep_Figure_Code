#eggNOG Taxa Analysis - Core vs acc
library(tidyverse)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(pillar)
library(readr)

###############################################################################
#Paths:
Base_path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/"
Genome_data_Path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/Le_Genome_Pangenome/MGYG-HGUT-00060/pan-genome/"
path_uhgp_90_cluster <- "/home/majernik14la/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Parquet/UHGP_Parquet/uhgp_90_cluster"

Path_core <- paste0(Genome_data_Path, "core_genes.faa")
Path_acc <- paste0(Genome_data_Path, "accessory_genes.faa")
Path_all <- paste0(Genome_data_Path, "pan-genome.faa")

#Save files:
path_UHGP_90_genes <- paste0(Base_path, "UHGP_90_Le_gene.tsv")
path_core_lst <- paste0(Base_path, "Le_gene_Core_acc_all.tsv")
###############################################################################
#New Paths:
Base_path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/"

Path_UHGP_90_eggNog <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Arrow_DB_Build/Making_UHGP_Parquet_DB/uhgp_db/uhgp_90_eggNOG_Fix_v2.tsv"
path_UHGP_90_genes <- paste0(Base_path, "UHGP_90_Le_gene.tsv")
path_core_lst <- paste0(Base_path, "Le_gene_Core_acc_all.tsv")

Path_core_eggNOG_summ <- paste0(Base_path, "core_eggNOG_summ.tsv")
Path_acc_eggNOG_summ <- paste0(Base_path, "acc_eggNOG_summ.tsv")
###############################################################################
Base_path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/"
Path_core_eggNOG_summ <- paste0(Base_path, "core_eggNOG_summ.tsv")
Path_acc_eggNOG_summ <- paste0(Base_path, "acc_eggNOG_summ.tsv")
###############################################################################
Base_path2 <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Genome_Analysis/"
G_meta_path <- paste0(Base_path2, "genomes-all_metadata.tsv")
Path_3 <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/Core_acc_Analysis/"
G_Pangenome_Path <- paste0(Path_3, "locus/g_100060/genes_presence-absence_locus.csv")
###############################################################################
Base_path3 <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/Le_Genome_Pangenome/"
###############################################################################
#Read in Core genes:
Core_genes <- read_lines(Path_core) %>%
  str_subset("^>") %>%
  str_extract("^>([^ ]+)") %>%
  str_remove("^>")
#Read in acc genes:
Acc_genes <- read_lines(Path_acc) %>%
  str_subset("^>") %>%
  str_extract("^>([^ ]+)") %>%
  str_remove("^>")
#Read in all genes:
all_genes <- read_lines(Path_all) %>%
  str_subset("^>") %>%
  str_extract("^>([^ ]+)") %>%
  str_remove("^>")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create Le only gene_id to uhgp-90_id converter:
Gut_Trans <- gsub("GUT_GENOME","", all_genes) %>%
  gsub("_", "", .)
  
prqt_uhgp_90 <- arrow::open_dataset(path_uhgp_90_cluster, format = "parquet")

all_genes_uhgp_90 <- prqt_uhgp_90 %>%
  filter(gene_id %in% Gut_Trans) %>%
  dplyr::collect() %>%
  dplyr::select(-c(group)) %>%
  rename(cluster_id_n = cluster_id)

final_all_genes_uhgp_90 <- all_genes_uhgp_90 %>%
  map_dfc(function(x) as.character(x)) %>%
  map_dfc(function(x) str_pad(x, 11, "left", pad = "0"))

Gut_Trans2 <- all_genes %>% as.tibble() %>% rename(gene_id = value) %>%
  mutate("gene_id_n" = gsub("GUT_GENOME","", gene_id) %>%
           gsub("_", "", .))

UHGP_90_Le_genes <- final_all_genes_uhgp_90 %>% mutate("cluster_id" =
                                 gsub('^(.{6})(.*)$','GUT_GENOME\\1_\\2',
                                      final_all_genes_uhgp_90$cluster_id_n)) %>%
  rename(gene_id_n = gene_id) %>%
  left_join(Gut_Trans2, by = "gene_id_n")
#Create a list of all UHGP genes that belong to Le:
write_tsv(UHGP_90_Le_genes, file = path_UHGP_90_genes)
length(Core_genes) <-length(all_genes)
length(Acc_genes) <- length(all_genes)
Le_genes_core_acc_all <- data.frame(all = all_genes, core = Core_genes, acc = Acc_genes)
write_tsv(Le_genes_core_acc_all, file = path_core_lst)
################################################################################
#Call in UHGP-90.tsv by read in lines:
UHGP_90_eggNOG <-read_tsv(Path_UHGP_90_eggNog, col_names = FALSE)
UHGP_90_Le_genes <- read_tsv(path_UHGP_90_genes)
Le_Core_Acc <- read_tsv(path_core_lst)

#Find the core and accessory cluster ID:
UHGP_90_core_id <- UHGP_90_Le_genes %>%
  filter(gene_id %in% Le_Core_Acc$core)

UHGP_90_acc_id <- UHGP_90_Le_genes %>%
  filter(gene_id %in% Le_Core_Acc$acc)

UHGP_90_Le_eggNOG_Core <- UHGP_90_eggNOG %>%
  filter(X1 %in% UHGP_90_core_id$cluster_id) %>%
  select(X1, X5)

UHGP_90_Le_eggNOG_acc <- UHGP_90_eggNOG %>%
  filter(X1 %in% UHGP_90_acc_id$cluster_id) %>%
  select(X1, X5)

#Summarise table to get counts for taxa:
core_eggNOG_summ <- UHGP_90_Le_eggNOG_Core %>%
  group_by(X5) %>%
  summarise(count = n(),
            X1 = paste(X1, collapse = ",")) %>%
  rename(cluster_id = X1, Putative_Taxa = X5)

acc_eggNOG_summ <- UHGP_90_Le_eggNOG_acc %>%
  group_by(X5) %>%
  summarise(count = n(),
            X1 = paste(X1, collapse = ",")) %>%
  rename(cluster_id = X1, Putative_Taxa = X5)

write_tsv(core_eggNOG_summ, Path_core_eggNOG_summ)
write_tsv(acc_eggNOG_summ, Path_acc_eggNOG_summ)
################################################################################
#Identify suspect taxa
acc_eggNOG_summ <- read_tsv(Path_acc_eggNOG_summ)
core_eggNOG_summ <- read_tsv(Path_core_eggNOG_summ)

Too_General <- c("Bacteria", "unclassified Bacteria", "Viruses", "dsDNA viruses, no RNA stage" )
  
Exspected_Taxa <-c("unclassified Lachnospiraceae", "Bacillota", "Clostridia", "Lachnospirales", "Lachnospiraceae", "Lachnospira",
                   "Blautia", "Coprococcus", "Dorea", "Lachnospira", "Oribacterium", "Roseburia", "Ruminococcus", "Ruminococcaceae",
                   "unclassified Clostridiales", "Butyrivibrio", "Pseudobutyrivibrio", "Lachnoanaerobaculum", "Lachnoclostridium",
                   "Oribacterium", "Eubacteriaceae")
Acc_contamination <- acc_eggNOG_summ %>%
  filter(! Putative_Taxa %in% Exspected_Taxa) %>%
  filter(!Putative_Taxa %in% Too_General)

Core_contamination <- core_eggNOG_summ %>%
  filter(! Putative_Taxa %in% Exspected_Taxa) %>%
  filter(!Putative_Taxa %in% Too_General)

Core_non <- core_eggNOG_summ %>%
  filter(Putative_Taxa %in% Exspected_Taxa) %>%
  filter(!Putative_Taxa %in% Too_General)

Acc_non <- acc_eggNOG_summ %>%
  filter(Putative_Taxa %in% Exspected_Taxa) %>%
  filter(!Putative_Taxa %in% Too_General)

Core_contamination_t <- sum(Core_contamination$count)
Core_non_t <- sum(Core_non$count)
Acc_contamination_t <-sum(Acc_contamination$count)
Acc_non_t <- sum(Acc_non$count)
Total_Core <- sum(Core_contamination_t, Core_non_t)
Total_acc <- sum(Acc_contamination_t, Acc_non_t)

Perent_con_Core <- Core_contamination_t/Total_Core * 100
Percent_con_Acc <- Acc_contamination_t/Total_acc *100

Core_too_Gen <- core_eggNOG_summ %>%
  filter(Putative_Taxa %in% Too_General)

Acc_too_Gen <- acc_eggNOG_summ %>%
  filter(Putative_Taxa %in% Too_General)

Acc_core <- bind_rows(
  Acc_contamination %>% mutate(Core_acc_gene = "accessory", Contaminate = "True"),
  Acc_non %>% mutate(Core_acc_gene = "accessory", Contaminate = "False"),
  Core_contamination %>% mutate(Core_acc_gene = "core", Contaminate = "True"),
  Core_non %>% mutate(Core_acc_gene = "core", Contaminate = "False"),
  Acc_too_Gen %>% mutate(Core_acc_gene = "accessory", Contaminate = "Too_General"),
  Core_too_Gen %>% mutate(Core_acc_gene = "core", Contaminate = "Too_General")
)

write_tsv(Acc_core, paste0(Base_path, "accessory_core_contaminate_genes.tsv"))
################################################################################
#ID if genes are in MAGs or Isolates:
Acc_core <- read_tsv(paste0(Base_path, "accessory_core_contaminate_genes.tsv"))
Genome_metadata <- read_tsv(G_meta_path)
Gene_metadata <- read_csv(G_Pangenome_Path)
Le_only_Gut_Trans <- read_tsv( path_UHGP_90_genes) 

Isolate_Genomes <- Genome_metadata %>%
  filter(Genome_type == "Isolate") %>%
  filter(MGnify_accession == "MGYG-HGUT-00060")

Acc_core_ex <- Acc_core %>%
  separate_longer_delim(cluster_id, delim = ",")

Isolate_Genomes_Genes <- Gene_metadata %>%
  select(any_of(Isolate_Genomes$Genome)) %>%
  pivot_longer(col = any_of(Isolate_Genomes$Genome), names_to = "Genome", values_to = "Genes") %>% 
  separate_longer_delim(Genes, delim = "\t") %>%
  inner_join(Le_only_Gut_Trans, by = c("Genes" = "gene_id")) %>%
  select(-c(gene_id_n, cluster_id_n)) %>%
  inner_join(Acc_core_ex, by = "cluster_id")

acc_genes_isolate <- Isolate_Genomes_Genes %>%
  filter(Core_acc_gene == "accessory") %>%
  distinct(cluster_id, .keep_all = TRUE) %>%
  group_by(Contaminate) %>%
  summarise(n = n(), .groups = 'drop') %>%
  bind_rows(summarise(., Contaminate = "Total", n = sum(n)))

                 
core_genes_isolate <- Isolate_Genomes_Genes %>%
  filter(Core_acc_gene == "core") %>%
  distinct(cluster_id, .keep_all = TRUE) %>%
  group_by(Contaminate) %>%
  summarise(n = n(), .groups = 'drop') %>%
  bind_rows(summarise(., Contaminate = "Total", n = sum(n)))


write_tsv(acc_genes_isolate, paste0(Base_path, "Acc_genes_Isolate_count.tsv"))
write_tsv(core_genes_isolate, paste0(Base_path, "Core_genes_Isolate_count.tsv"))

################################################################################


 