#Genome Analysis
library(readr)
library(tidyverse)
library(ggplot2)

################################################################################
#Paths#
Base_path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Genome_Analysis/"
Base_sv_path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Genome_Analysis/"
Base_path2 <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/"

G_G100060_path <- paste0(Base_path,"Sig_Genes_Pangenomes_100060.csv")
G_G100271_path <- paste0(Base_path,"Sig_Genes_Pangenomes_100271.csv")
G_G100078_path <- paste0(Base_path,"Sig_Genes_Pangenomes_100078.csv")
G_G102528_path <- paste0(Base_path,"Sig_Genes_Pangenomes_102528.csv")
G_G100217_path <- paste0(Base_path,"Sig_Genes_Pangenomes_100217.csv")
G_G101380_path <- paste0(Base_path,"Sig_Genes_Pangenomes_101380.csv")

Egg_load <- paste0(Base_path2, "uhgp_90_eggNOG.tsv")

G_meta_path <- paste0(Base_path, "genomes-all_metadata.tsv")
#Made in code
Gene_Genome_Type_Counts_path <- paste0(Base_path,"Gene_Genome_Type_Counts.tsv")
Gene_Genome_Type_Counts_with_Metadata_path <- paste0(Base_path,"Gene_Genome_Type_Counts_with_Metadata.csv")
Gene_Genome_Type_Counts_with_Metadata_eggNOG_path <- paste0(Base_path,"Gene_Genome_Type_Counts_with_Metadata_eggNOG.csv")
genes_with_genome_meta_Cleaned_path <- paste0(Base_path, "genes_with_genome_meta_Cleaned.rds")

################################################################################
#Load in eggNOG data:
uhgp_90_eggNOG <- read_tsv(file = Egg_load)
#Load in Significant genes by genome
G_G100060 <- read_csv(file = G_G100060_path)
G_G100271 <- read_csv(file = G_G100271_path)
G_G100078 <- read_csv(file = G_G100078_path)
G_G102528 <- read_csv(file = G_G102528_path)
G_G100217 <- read_csv(file = G_G100217_path)
G_G101380 <- read_csv(file = G_G101380_path)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Now we will clean, select, and concatenate:
DFg_n <- c("G_G102528", "G_G100078", "G_G100271", "G_G100060", "G_G100217", "G_G101380")
DFs_g <- list(G_G102528, G_G100078, G_G100271, G_G100060, G_G100217, G_G101380) %>% set_names(DFg_n)
#Create gene_id column:
Dfs_g <- lapply(DFs_g, function(x){
  x <- x %>%
    mutate(Gene_id = paste0("UHGG", UHGG_Gene)) %>% select("Gene_id", "No. isolates")
  return(x)
}) %>% setNames(DFg_n)
Genes_isolates <- Dfs_g %>% bind_rows()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Load eggnog data from analysis
Tmp <- read_tsv(file = Egg_load)
Gene_Genome_Analysis <- Genes_isolates %>% 
  left_join(Tmp %>% select(c("Species_id", "Gene_id", "Species", 
                             "Lineage_Shared", "cor_max_species")), by = "Gene_id")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Isolate or MAG for each Gene:
DFg_l <- list("G_G102528", "G_G100078", "G_G100271", "G_G100060", "G_G100217", "G_G101380")
DFs_g <- list(G_G102528, G_G100078, G_G100271, G_G100060, G_G100217, G_G101380) %>% set_names(DFg_l)

DF_g_l <- DFs_g %>%
  lapply(function(x){
    x %>% 
      select(-c(2:15)) %>%
      #2:15 needed for metadata removal
      pivot_longer(!UHGG_Gene, names_to = "all_Genomes", values_to = "pres_Genes")
  }) %>% setNames(DFg_l)

DF_g_l_p <- DF_g_l %>%
  lapply(function(x){
    x %>%
      mutate(present_Genomes = if_else(!is.na(x$pres_Genes), x$all_Genomes, NA))
  }) %>% setNames(DFg_l)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Now lets add the metadata:
all_genome_meta <- read_tsv(G_meta_path)
Needed_genome_Meta <- all_genome_meta %>% select(c("Genome", "Genome_type", "Contamination"))

genes_with_genome_meta <- DF_g_l_p %>%
  lapply(function(x){
    left_join(x, Needed_genome_Meta, by = c("present_Genomes" = "Genome"))
  })
Gene_Genome_I_M_Counts <- genes_with_genome_meta %>%
  lapply( function(x){
    x%>%
      group_by(UHGG_Gene, Genome_type) %>%
      summarize(count = n())
  })

Gene_Genome_I_M_Counts_unlist <- Gene_Genome_I_M_Counts %>% enframe(name =  "Species", value = "Data") %>% unnest(Data) %>% select(-c("Species"))
Gene_Genome_Type_Counts <- Gene_Genome_I_M_Counts_unlist %>%
  mutate(UHGG_Gene = paste0("Gut_Genome", UHGG_Gene)) %>%
  na.omit(cols = "Genome_type") %>% pivot_wider(names_from = Genome_type, values_from = count) %>% rename(Gene = UHGG_Gene)

write_tsv(Gene_Genome_Type_Counts, Gene_Genome_Type_Counts_path)
###########################################################
#Make Gene_Genome_Type_Counts With Additional Metadata:
Gene_Genome_Type_Counts_with_Meta <- Gene_Genome_Type_Counts %>% 
  mutate("Num" = gsub('Gut_Genome', "UHGG", Gene)) %>%
  left_join(uhgp_90_eggNOG, by=c("Num" = "Gene_id")) %>%
  group_by(Species_id) %>%
  select(Gene, Isolate, MAG, Lineage_Shared, Species_id, Species, cor_max_species)

write_csv(Gene_Genome_Type_Counts_with_Meta, Gene_Genome_Type_Counts_with_Metadata_path)
###########################################################
Total_M_I_Counts <- all_genome_meta %>%
  group_by(MGnify_accession, Genome_type) %>%
  summarize(count = n())
###########################################################
genes_with_genome_meta_Cleaned <- genes_with_genome_meta %>%
  lapply( function(x){
    x%>%
      select(-c("all_Genomes", "pres_Genes"))
  })
write_rds(genes_with_genome_meta_Cleaned, genes_with_genome_meta_Cleaned_path)
############################################################
Gene_Genome_Type_Counts_with_Meta_eggNOG <- Gene_Genome_Type_Counts %>% 
  mutate("Num" = gsub('Gut_Genome', "UHGG", Gene)) %>%
  left_join(uhgp_90_eggNOG, by=c("Num" = "Gene_id")) %>%
  group_by(Species_id)
write_csv(Gene_Genome_Type_Counts_with_Meta_eggNOG, Gene_Genome_Type_Counts_with_Metadata_eggNOG_path )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Create list on contaminate and non-contaminate genes:
Con_Genes <- Gene_Genome_Analysis$Gene_id[!Gene_Genome_Analysis$Lineage_Shared] %>% gsub("UHGG", "", .)
Non_Genes <- Gene_Genome_Analysis$Gene_id[Gene_Genome_Analysis$Lineage_Shared] %>% gsub("UHGG", "", .)
#Pre-clean duplicates out:
#Note: Genomes should never be shared between species

DeRep_genes_with_genome_meta <- genes_with_genome_meta %>%
  lapply(function(x){
    x %>%
      distinct(present_Genomes, .keep_all = TRUE)
  })  
#Create lists of Contamination
Cont_val_total <- c()
ContValues_con_l <- DeRep_genes_with_genome_meta %>%
  lapply(function(x){
    Cont_val <- x %>%
      filter(UHGG_Gene %in% Con_Genes) %>%
      pull(Contamination)
    Cont_val_total <- rbind(Cont_val_total, Cont_val)
    return(Cont_val_total)
  })
ContValues_non_l <- DeRep_genes_with_genome_meta %>%
  lapply(function(x){
    Cont_val <- x %>%
      filter(x$UHGG_Gene %in% Non_Genes) %>%
      pull(Contamination)
    Cont_val_total <- rbind(Cont_val_total, Cont_val)
    return(Cont_val_total)
  })
#Extract all of the values into a single vector:
ContValues_con <- ContValues_con_l %>% unlist() %>% as.vector()
ContValues_con <- ContValues_con[which(!is.na(ContValues_con))]
con_Df <- cbind("Contamination" = ContValues_con,"Gene" = rep("Contaminate", times = length(ContValues_con)))
ContValues_non <- ContValues_non_l %>% unlist() %>% as.vector()
ContValues_non <- ContValues_non[which(!is.na(ContValues_non))]
non_Df <- cbind("Contamination" = ContValues_non,"Gene" = rep("Non-contaminate", times = length(ContValues_non)))
ContValues <- rbind.data.frame(con_Df, non_Df)
ContValues$Contamination <- as.numeric(ContValues$Contamination)

 # ggplot(ContValues, aes(x = Contamination, fill = Gene)) + 
 #   geom_histogram(alpha = 0.5, position = "identity") +
 #   labs(title = "Percent Contamination in Geneomes",
 #        x = "Percent Contamination",
 #        y = "Number of Genomes")
################################################################################ 
 # ggplot(ContValues, aes(Contamination, fill = Gene)) +
 #   geom_histogram(alpha = 0.5, aes(y =  ..density..), position = 'identity')+
 #   labs(title = "Percent Contamination in Geneomes",
 #        x = "Percent Contamination",
 #        y = "Number of Genomes")
 
 ContValues$Gene <- factor(ContValues$Gene, levels = c("Non-contaminate", "Contaminate")) 
 
p3 <- ggplot(ContValues, aes(Contamination, fill = Gene, color = Gene)) +
  geom_density(alpha = 0.7, aes(y =  ..density..), position = 'identity')+
   labs(title = "Percent Contamination in Genomes that Contain Significant Genes",
        x = "Percent Contamination",
        y = "Number of Genomes") +
   scale_fill_manual(name = "Lineage Test",
                     values = c("Contaminate" = "orange", "Non-contaminate" = "deepskyblue4"),
                     labels = c("Pass", "Fail")
                     ) +
   scale_color_manual(name = "Lineage Test",
     values = c("Contaminate" = "orange3", "Non-contaminate" = "deepskyblue4"),
     labels = c("Pass", "Fail")
   )
 
 ggsave(paste0(Base_sv_path,"Percent_Contamination_in_Genomes.tiff"), 
        plot = p3, device = "tiff", width = 7.5,height=5, dpi=300, units = "in")
 ###############################################################################
 Genome_summart <- ContValues %>% 
   group_by(Gene) %>%
   summarise(
     min = min(Contamination),
     max = max(Contamination),
     mean = mean(Contamination),
     median = median(Contamination)
   )
 ##############################################################################
 #Make Gene_Genome_Type_Counts With Additional Metadata & BLAST Results:

  BLST <- read_csv(file = paste0(Base_path, "BLAST_RSLT.csv"))
 G_G_Egg_BLST <-  Gene_Genome_Type_Counts_with_Meta_eggNOG %>%
   left_join(BLST, by = c("Num" = "Gene_id"))

 write_csv(G_G_Egg_BLST, file = paste0(Base_path, "Gene_Genome_Analysis_with_EggNOG_and_BLAST.csv"))