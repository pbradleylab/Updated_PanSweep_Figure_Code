library(readr)
library(tidyverse)
library(ggplot2)
library(viridis)

#Paths:
path_for_Species_Abundance <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/Sp_Abundance_Files/species_marker_read_counts.tsv"
path_for_midas_genes <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/cirrhosis-analysis/data/processed/midas_snv_merge/genes/"
Base_path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Genome_Analysis/"
G_meta_path <- paste0(Base_path, "genomes-all_metadata.tsv")
path_for_sample_metadata <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/cirrhosis-analysis/data/manual/phylogenize/cirrhosis-metadata.tsv"
Graph_Directory <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Genome_Metadata <- read_tsv(G_meta_path)
Species_Abundance <- read_tsv(path_for_Species_Abundance)
sample_meta <- read_tsv(path_for_sample_metadata) %>% select(c("sample","env"))
################################################################################
#Functions:
separate_taxonomy_with_s <- function(inpt, taxa_col){
  inpt <- inpt %>%
    separate(taxa_col, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
    mutate(across(c("d", "p", "c", "o", "f", "g", "s"), ~ gsub("[dpcofgs]__","", .))) %>%
    rename_with(~ case_when(
      . == "d" ~ "Domain",
      . == "p" ~ "Phylum",
      . == "c" ~ "Class",
      . == "o" ~ "Order",
      . == "f" ~ "Family",
      . == "g" ~ "Genus",
      . == "s" ~ "Species",
      TRUE ~ .
    ))
  return(inpt)
}
###############################################################################
#Graph the number of genomes types vs Species abundance:
Gene_Meta_Sp <- Genome_Metadata %>% select(c("Genome_type", "MGnify_accession"))

Gene_Meta_Sp_Iso_MAG <- Gene_Meta_Sp %>%
  group_by(MGnify_accession, Genome_type) %>%
  summarise(count = n())

Gene_Meta_Sp_count <- Gene_Meta_Sp %>%
  group_by(MGnify_accession) %>%
  summarise(count = n())

Gene_Meta_Sp_Iso_MAG <- Gene_Meta_Sp_Iso_MAG %>% 
  group_by(MGnify_accession) %>%
  mutate(total = sum(count)) %>%
  mutate(percent = count / total *100)

Gene_Meta_MAG_per <- Gene_Meta_Sp_Iso_MAG %>% filter(Genome_type == "MAG")

Gene_Meta_all <- Gene_Meta_Sp_count %>% left_join(Gene_Meta_MAG_per, by = "MGnify_accession") %>% select("MGnify_accession", "count.x", "percent") %>%
  mutate(percent = replace_na(percent, 0)) %>% rename(percent_MAG = percent, Num_Genomes = count.x) %>%
  mutate(MGnify_accession = gsub("MGYG-HGUT-", 1, MGnify_accession)) %>% rename(species_id = MGnify_accession)

Lg_Sp_abd <- Species_Abundance %>% pivot_longer(cols = 2:ncol(Species_Abundance),
                                                names_to = "Sample",
                                                values_to = "Relative_Abd") %>%
  mutate(species_id = as.character(species_id))
Sp_abd_Geneome_counts <- Lg_Sp_abd %>% left_join(Gene_Meta_all, by = "species_id")

Sp_abd_Geneome_counts_env <- Sp_abd_Geneome_counts %>% left_join(sample_meta, by = c("Sample" = "sample"))

Mn_Me_Sp_abd <- Sp_abd_Geneome_counts_env %>%
  group_by(species_id, env) %>%
  mutate(Avg = mean(Relative_Abd), Med = median(Relative_Abd))

Log_Sp_abd <- Sp_abd_Geneome_counts_env %>%
  group_by(species_id, env) %>%
  mutate(Log = log(Relative_Abd + 1)) %>%
  mutate(L_Avg = mean(Log))

Plot5 <- ggplot(data = Log_Sp_abd, aes(x = Num_Genomes, y = L_Avg, shape = factor(env), color = percent_MAG)) +
  geom_point() +
  scale_color_viridis_c(option = "viridis") +
  labs(title = "Number of genomes vs Mean(Log + 1) Species Abundance in Cirrhosis Dataset", x = "Number of Genomes in UHGG", y = "Log(Relative Abundance + 1)", color = "Precent MAG", shape = "Case vs Control")
Plot5

ggsave("Genome_Num_vs_Species_Abd_Log_Avg.tiff", plot = Plot5, device = "tiff", width = 7.5,height=6, dpi=300, units = "in")
################################################################################

