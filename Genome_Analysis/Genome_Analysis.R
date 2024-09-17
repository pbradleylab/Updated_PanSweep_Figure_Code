#Genome Analysis
library(readr)
library(tidyverse)

################################################################################
#Paths#
Base_path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/MIDAS_Genome_analysis/"
Base_sv_path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/MIDAS_Genome_analysis"


G_G100060_path <- paste0(Base_path,"Sig_Genes_Pangenomes_100060.csv")
G_G100271_path <- paste0(Base_path,"Sig_Genes_Pangenomes_100271.csv")
G_G102580_path <- paste0(Base_path,"Sig_Genes_Pangenomes_102580.csv")
G_G103694_path <- paste0(Base_path,"Sig_Genes_Pangenomes_103694.csv")

Egg_load <- paste0(Base_path, "uhgp_90_eggNOG.tsv")

G_meta_path <- paste0(Base_path, "genomes-all_metadata.tsv")

Gene_Genome_Type_Counts_path <- paste0(Base_path,"Gene_Genome_Type_Counts.tsv")
genes_with_genome_meta_Cleaned_path <- paste0(Base_path, "genes_with_genome_meta_Cleaned.rds")

################################################################################
#Load in Significant genes by genome
G_G100060 <- read_csv(file = G_G100060_path)
G_G100271 <- read_csv(file = G_G100271_path)
G_G102580 <- read_csv(file = G_G102580_path)
G_G103694 <- read_csv(file = G_G103694_path)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Now we will clean, select, and concatenate:
DFg_n <- c("G_G103694", "G_G102580", "G_G100271", "G_G100060")
DFs_g <- list(G_G103694, G_G102580, G_G100271, G_G100060) %>% set_names(DFg_n)
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
DFg_l <- list("G_G103694", "G_G102580", "G_G100271", "G_G100060")
DFs_g <- list(G_G103694, G_G102580, G_G100271, G_G100060) %>% set_names(DFg_l)

DF_g_l <- DFs_g %>%
  lapply(function(x){
    x %>% 
      select(-c(2:15)) %>%
      
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Remove genes that are in nearly every genome:
remove <- c("180526_01821", "055467_01854")
Con_Genes <- Gene_Genome_Analysis$Gene_id[!Gene_Genome_Analysis$Lineage_Shared] %>% gsub("UHGG", "", .) %>% setdiff(remove)
Non_Genes <- Gene_Genome_Analysis$Gene_id[Gene_Genome_Analysis$Lineage_Shared] %>% gsub("UHGG", "", .)

DeRep_genes_with_genome_meta <- genes_with_genome_meta %>%
  lapply(function(x){
    x %>%
      distinct(present_Genomes, .keep_all = TRUE)
  })  

Cont_val_total <- c()
ContValues_con_l <- DeRep_genes_with_genome_meta %>%
  lapply(function(x){
    Cont_val <- x %>%
      filter(x$UHGG_Gene %in% Con_Genes) %>%
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
#Make a second one with only non-contaminate genes:
Non_Cont_Values_DF <-as.data.frame(non_Df)
Non_Cont_Values_DF$Contamination <- as.numeric(Non_Cont_Values_DF$Contamination)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Make a plot that has 
p3 <- Non_Cont_Values_DF %>%
  ggplot(aes(x=Contamination, fill=Gene)) +
  geom_histogram() +
  geom_vline(aes(xintercept = ContValues_con[1], color = "Contaminate"), linetype = "dashed") +
  geom_vline(aes(xintercept = ContValues_con[2], color = "Contaminate"), linetype = "dashed") +
  scale_color_manual(name = "Gene", values = c(Contaminate = "blue")) +
  labs(title = "Percent Contamination in Geneomes",
       x = "Percent Contamination",
       y = "Number of Genomes",
       fill = NULL)

p3
ggsave(paste0(Base_path,"Percent_Contamination_in_Genomes.tiff"), 
       plot = p3, device = "tiff", width = 7.5,height=5, dpi=300, units = "in")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
