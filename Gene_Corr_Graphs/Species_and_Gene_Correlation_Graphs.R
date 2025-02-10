                        #Graphs for Figure 1 and Supplemental#
library(tidyverse)
library(readr)
library(ggplot2)
library(patchwork)
                                    #Paths#

Base_sv_path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Genome_Analysis/"
Base_path2 <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/"
Base_path3 <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/"

path_for_genomes_all_metadata <- paste0(Base_path3, "Metadata/Genome_metadata/metadata.tsv")
path_for_species_abundance <- paste0(Base_path3, "PanSweep_Analysis_Updated/data/species/species_marker_read_counts.tsv")
path_for_sample_metadata <- paste0(Base_path3, "Sample_meta.tsv")
path_for_gene_reads <- paste0(Base_path2, "Sig_Gene_reads.rds")
path_to_cor_org_species <- paste0(Base_path2, "Species_id_cor_and_orig.rds")
Egg_load <- paste0(Base_path2, "uhgp_90_eggNOG.tsv")
path_phylo_md2 <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/cirr_md.tsv"

#Saving#
Graph_Save_Dir <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/Review_Response/Figure1/"
Supp_Graph_Save_Dir <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Gene_Corr_Graphs/Sup"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        #Load in species id and names of species#
genome_metadata <- read_tsv(file= path_for_genomes_all_metadata)

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

meta_genome_sep_taxa <- genome_metadata %>% separate_taxonomy_with_s("Lineage") %>%
  mutate(species_id = as.character(species_id))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                    #Load in metadata transformation Function# 
library(PanSweep)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                          #Load in Species abundance#
phylo_md2 <- read_tsv(path_phylo_md2, show_col_types = FALSE)
merge_fn_counts = base::sum
Species_Abd <-read_tsv(path_for_species_abundance) %>%
  merge_columns_tbl(md=phylo_md2, fn=merge_fn_counts) %>%
  pivot_longer(!species_id, names_to = "subject", values_to = "species_count") %>% 
  pivot_wider(names_from = "species_id", values_from = "species_count") %>% column_to_rownames("subject")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                              #Load in Gene Reads#
SGR <- read_rds(file = path_for_gene_reads)
SIG_Gene_Reads <- lapply(names(SGR), function(x){
  SGR[[x]] %>% as.data.frame() %>% rownames_to_column("gene_id") %>%
    pivot_longer(!gene_id, names_to = "subject", values_to = "gene_r") %>% 
    pivot_wider(names_from = "gene_id", values_from = "gene_r") %>% column_to_rownames("subject")
}) %>% set_names(names(SGR))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                            #Load in  Sample Metadata#
sample_meta <- read_tsv(path_for_sample_metadata)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                          #Create Plots for figure 1#
sp_reads_mtx <- t(Species_Abd)
colnames(sp_reads_mtx) <- gsub("-", "", colnames(sp_reads_mtx))
lel_reads_mtx <- t(SIG_Gene_Reads$`100060`)
cols_both <- intersect(colnames(sp_reads_mtx), colnames(lel_reads_mtx))

Va_con_gf <- tibble(cols = cols_both, rank.species = dplyr::percent_rank(sp_reads_mtx["101444", cols_both]), rank.gene=dplyr::percent_rank(lel_reads_mtx["UHGG047117_02379", cols_both]), 
                    gene.zero = lel_reads_mtx["UHGG047117_02379", cols_both] == 0) %>% left_join(sample_meta, by = c("cols" = "subject")) %>% select(-c("sample")) %>%
  ggplot(aes(x=rank.species, y=rank.gene, shape = interaction(env, gene.zero))) + 
  geom_point(size=5, color = "orange") +
  geom_smooth(method='lm', color = "gray28", aes(group = 1)) +
  scale_shape_manual(values = c("control.TRUE" = 1, "control.FALSE" = 16, "case.TRUE" = 2, "case.FALSE" = 17)) +
  theme(  legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 20))
                                                                  
Le_con_gf <- tibble(cols = cols_both, rank.species = dplyr::percent_rank(sp_reads_mtx["100060", cols_both]), rank.gene=dplyr::percent_rank(lel_reads_mtx["UHGG047117_02379", cols_both]), 
                    gene.zero = lel_reads_mtx["UHGG047117_02379", cols_both] == 0) %>% left_join(sample_meta, by = c("cols" = "subject")) %>% select(-c("sample")) %>%
  ggplot(aes(x=rank.species, y=rank.gene, shape = interaction(env, gene.zero))) + 
  geom_point(size=5, color = "deepskyblue4") +
  geom_smooth(method='lm', color = "gray28", aes(group = 1)) +
  scale_shape_manual(values = c("control.TRUE" = 1, "control.FALSE" = 16, "case.TRUE" = 2, "case.FALSE" = 17)) +
  theme(  legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 20))

Va_ncon_gf <- tibble(cols = cols_both, rank.species = dplyr::percent_rank(sp_reads_mtx["101444", cols_both]), rank.gene=dplyr::percent_rank(lel_reads_mtx["UHGG152466_01649", cols_both]), 
                     gene.zero = lel_reads_mtx["UHGG152466_01649", cols_both] == 0) %>% left_join(sample_meta, by = c("cols" = "subject")) %>% select(-c("sample")) %>%
  ggplot(aes(x=rank.species, y=rank.gene, shape = interaction(env, gene.zero))) + 
  geom_point(size=5, color = "orange") +
  geom_smooth(method='lm', color = "gray28", aes(group = 1)) +
  scale_shape_manual(values = c("control.TRUE" = 1, "control.FALSE" = 16, "case.TRUE" = 2, "case.FALSE" = 17)) +
  theme(  legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 20))

Le_ncon_gf <- tibble(cols = cols_both, rank.species = dplyr::percent_rank(sp_reads_mtx["100060", cols_both]), rank.gene=dplyr::percent_rank(lel_reads_mtx["UHGG152466_01649", cols_both]), 
                     gene.zero = lel_reads_mtx["UHGG152466_01649", cols_both] == 0) %>% left_join(sample_meta, by = c("cols" = "subject")) %>% select(-c("sample")) %>%
  ggplot(aes(x=rank.species, y=rank.gene, shape = interaction(env, gene.zero))) + 
  geom_point(size=5, color = "deepskyblue4") +
  geom_smooth(method='lm', color = "gray28", aes(group = 1)) +
  scale_shape_manual(values = c("control.TRUE" = 1, "control.FALSE" = 16, "case.TRUE" = 2, "case.FALSE" = 17)) +
  theme(  legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 20))

setwd(Graph_Save_Dir)
ggsave("Final_update_V.atypica_and_Contaminate_Gene.tiff", plot = Va_con_gf, device = "tiff", width = 7.5,height=6, dpi=300, units = "in")
ggsave("Final_update_L.eligens_and_Contaminate_Gene.tiff", plot = Le_con_gf, device = "tiff", width = 7.5,height=6, dpi=300, units = "in")
ggsave("Final_update_V.atypica_and_Non-Contaminate_Gene.tiff", plot = Va_ncon_gf, device = "tiff", width = 7.5,height=6, dpi=300, units = "in")
ggsave("Final_update_L.eligens_and_Non-Contaminate_Gene.tiff", plot = Le_ncon_gf, device = "tiff", width = 7.5,height=6, dpi=300, units = "in")
#################################################################################
                              #Make Supplemental Graphs#
setwd(Supp_Graph_Save_Dir)
        #Load in all gene max correlation species and originating species#
Max_org_species <- read_rds(path_to_cor_org_species)
Species_for_graph <- Max_org_species %>% filter(mark == "max")


                 #Separate out subjects and combine with species#
Species_Abd_join <- Species_Abd %>% rownames_to_column("subject")
Slct_Gns100060 <- SIG_Gene_Reads$'100060' %>% rownames_to_column("subject")
Slct_Gns100078 <- SIG_Gene_Reads$'100078' %>% rownames_to_column("subject")
Slct_Gns100271 <- SIG_Gene_Reads$'100271' %>% rownames_to_column("subject")
Slct_Gns102528 <- SIG_Gene_Reads$'102528' %>% rownames_to_column("subject")
Slct_Gns100217 <- SIG_Gene_Reads$'100217' %>% rownames_to_column("subject")
Slct_Gns101380 <- SIG_Gene_Reads$'101380' %>% rownames_to_column("subject")

Graph_tbl100060 <- Slct_Gns100060 %>% left_join(Species_Abd_join, by = 'subject') %>% 
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>% 
  left_join(sample_meta, by = c("subject" = "subject"))
Graph_tbl100078 <- Slct_Gns100078 %>% left_join(Species_Abd_join, by = 'subject') %>%
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("subject" = "subject"))  
Graph_tbl100271 <- Slct_Gns100271 %>% left_join(Species_Abd_join, by = 'subject') %>%
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("subject" = "subject"))  
Graph_tbl102528 <- Slct_Gns102528 %>% left_join(Species_Abd_join, by = 'subject') %>%
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("subject" = "subject"))
Graph_tbl100217 <- Slct_Gns100217 %>% left_join(Species_Abd_join, by = 'subject') %>% 
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>% 
  left_join(sample_meta, by = c("subject" = "subject"))
Graph_tbl101380 <- Slct_Gns101380 %>% left_join(Species_Abd_join, by = 'subject') %>%
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("subject" = "subject"))  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Make a vector of gene names by species#
Gene_100060 <- Species_for_graph %>% filter(Species == 100060)
Gene_V_100060 <- Gene_100060$Gene
                    #Make graphs for L. eligens (100060)#
po60 <-map(Gene_V_100060, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
   p<-ggplot(Graph_tbl100060, aes(x = !!sym("100060"), y = !!sym(g), shape = factor(env))) +
     geom_point(size = 2, color = "deepskyblue4") +
     geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
     guides(
       color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
       shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
     ) +
     #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
     #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
     labs(#title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("L. eligens")),
          x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
          y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("L. eligens"))
            ) +
     theme(legend.position = "none",
           axis.text.x = element_text(size = 5), 
           axis.text.y = element_text(size = 5),
           plot.title = element_text(size = 5),
           axis.title.x = element_text(size = 5),
           axis.title.y = element_text(size = 5)
           )
   return(p)
   #ggsave(filename = paste0(g, "_Org_Sp_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                #Make a vector of gene names by species#
Gene_100078 <- Species_for_graph %>% filter(Species == 100078)
Gene_V_100078 <- Gene_100078$Gene
                  #Make graphs for Lachnospira rogosae (100078)#
po78 <-map(Gene_V_100078, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
  p <- ggplot(Graph_tbl100078, aes(x = !!sym("100078"), y = !!sym(g), shape = factor(env))) +
    geom_point(size = 2, color = "deepskyblue4") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(#title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("Lachnospira rogosae")),
         x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
         y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("Lachnospira rogosae"))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
          )
  return(p)
  #ggsave(filename = paste0(g, "_Org_Sp_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Make a vector of gene names by species#
Gene_100271 <- Species_for_graph %>% filter(Species == 100271)
Gene_V_100271 <- Gene_100271$Gene
#Make graphs for Roseburia inulinivorans (100271)#
po71 <- map(Gene_V_100271, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
  p <- ggplot(Graph_tbl100271, aes(x = !!sym("100271"), y = !!sym(g), shape = factor(env))) +
    geom_point(size = 2, color = "deepskyblue4") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(#title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("Roseburia inulinivorans")),
         x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
         y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("Roseburia inulinivorans"))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
          )
  return(p)
  #ggsave(filename = paste0(g, "_Org_Sp_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Make a vector of gene names by species#
Gene_102528 <- Species_for_graph %>% filter(Species == 102528)
Gene_V_102528 <- Gene_102528$Gene
#Make graphs for "Anaerostipes hadrus" (102528)#
po28 <-map(Gene_V_102528, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
  p <- ggplot(Graph_tbl102528, aes(x = !!sym("102528"), y = !!sym(g), shape = factor(env))) +
    geom_point(size = 2, color = "deepskyblue4") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(#title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("Anaerostipes hadrus")),
         x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
         y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("Anaerostipes hadrus"))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
          )
  return(p)
  #ggsave(filename = paste0(g, "_Org_Sp_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Gene_100217 <- Species_for_graph %>% filter(Species == 100217)
Gene_V_100217 <- Gene_100217$Gene
#Make graphs for "Acetatifactor sp900066565" (100217)#
po17 <-map(Gene_V_100217, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
  p<-ggplot(Graph_tbl100217, aes(x = !!sym("100217"), y = !!sym(g), shape = factor(env))) +
    geom_point(size = 2, color = "deepskyblue4") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(#title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("L. eligens")),
      x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
      y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("Acetatifactor sp900066565"))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
    )
  return(p)
  #ggsave(filename = paste0(g, "_Org_Sp_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Gene_101380 <- Species_for_graph %>% filter(Species == 101380)
Gene_V_101380 <- Gene_101380$Gene
#Make graphs for "Faecalicatena gnavus" (101380)#
po80 <-map(Gene_V_101380, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
  p <- ggplot(Graph_tbl101380, aes(x = !!sym("101380"), y = !!sym(g), shape = factor(env))) +
    geom_point(size = 2, color = "deepskyblue4") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(#title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("Lachnospira rogosae")),
      x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
      y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("Faecalicatena gnavus"))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
    )
  return(p)
  #ggsave(filename = paste0(g, "_Org_Sp_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})
################################################################################
#Make correlating species to gene graphs
setwd(Supp_Graph_Save_Dir)

#Separate out subjects and combine with species#
Species_Abd_join <- Species_Abd %>% rownames_to_column("subject")
Slct_Gns100060 <- SIG_Gene_Reads$'100060' %>% rownames_to_column("subject")
Slct_Gns100078 <- SIG_Gene_Reads$'100078' %>% rownames_to_column("subject")
Slct_Gns100271 <- SIG_Gene_Reads$'100271' %>% rownames_to_column("subject")
Slct_Gns102528 <- SIG_Gene_Reads$'102528' %>% rownames_to_column("subject")
Slct_Gns100217 <- SIG_Gene_Reads$'100217' %>% rownames_to_column("subject")
Slct_Gns101380 <- SIG_Gene_Reads$'101380' %>% rownames_to_column("subject")

Graph_tbl100060 <- Slct_Gns100060 %>% left_join(Species_Abd_join, by = 'subject') %>% 
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>% 
  left_join(sample_meta, by = c("subject" = "subject"))
Graph_tbl100078 <- Slct_Gns100078 %>% left_join(Species_Abd_join, by = 'subject') %>%
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("subject" = "subject"))  
Graph_tbl100271 <- Slct_Gns100271 %>% left_join(Species_Abd_join, by = 'subject') %>%
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("subject" = "subject"))  
Graph_tbl102528 <- Slct_Gns102528 %>% left_join(Species_Abd_join, by = 'subject') %>%
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("subject" = "subject"))
Graph_tbl100217 <- Slct_Gns100217 %>% left_join(Species_Abd_join, by = 'subject') %>% 
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>% 
  left_join(sample_meta, by = c("subject" = "subject"))
Graph_tbl101380 <- Slct_Gns101380 %>% left_join(Species_Abd_join, by = 'subject') %>%
  mutate(across(-c("subject"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("subject" = "subject"))  

#Make the plots for the correlated species vs gene reads
pc60 <- map2(Gene_100060$Species_Cor, Gene_100060$Gene, function(.x, .y){
  #Plot_name <- paste(g, "_Plot_orig")
  Sp_Name <- meta_genome_sep_taxa %>% filter(species_id == .x) %>% select("Species") %>% as.character()
  p <- ggplot(Graph_tbl100060, aes(x = !!sym(.x), y = !!sym(.y), shape = factor(env))) +
    geom_point(size = 2, color = "orange") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(#title = bquote(.(.y) ~"Gene Counts versus Relative Abundance" ~ italic(.(Sp_Name))),
      x = bquote(log[10] ~ "Gene count + 1 of" ~ .(.y)),
      y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic(.(Sp_Name)))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
          )
  return(p)
  #ggsave(filename = paste0(.y, "_max_correlation_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})

pc71 <- map2(Gene_100271$Species_Cor, Gene_100271$Gene, function(.x, .y){
  #Plot_name <- paste(g, "_Plot_orig")
  Sp_Name <- meta_genome_sep_taxa %>% filter(species_id == .x) %>% select("Species") %>% as.character()
  p <- ggplot(Graph_tbl100271, aes(x = !!sym(.x), y = !!sym(.y), shape = factor(env))) +
    geom_point(size = 2, color = "orange") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(
      x = bquote(log[10] ~ "Gene count + 1 of" ~ .(.y)),
      y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic(.(Sp_Name)))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
          )
  return(p)
  #ggsave(filename = paste0(.y, "_max_correlation_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})

pc78 <- Cor_plt_100078 <- map2(Gene_100078$Species_Cor, Gene_100078$Gene, function(.x, .y){
  #Plot_name <- paste(g, "_Plot_orig")
  Sp_Name <- meta_genome_sep_taxa %>% filter(species_id == .x) %>% select("Species") %>% as.character()
 p <- ggplot(Graph_tbl100078, aes(x = !!sym(.x), y = !!sym(.y), shape = factor(env))) +
    geom_point(size = 2, color = "orange") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(
      x = bquote(log[10] ~ "Gene count + 1 of" ~ .(.y)),
      y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic(.(Sp_Name)))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
          )
 return(p)
 #ggsave(filename = paste0(.y, "_max_correlation_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})

pc28 <- map2(Gene_102528$Species_Cor, Gene_102528$Gene, function(.x, .y){
  Sp_Name <- meta_genome_sep_taxa %>% filter(species_id == .x) %>% select("Species") %>% as.character()
  p <- ggplot(Graph_tbl102528, aes(x = !!sym(.x), y = !!sym(.y), shape = factor(env))) +
    geom_point(size = 2, color = "orange") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(
         x = bquote(log[10] ~ "Gene count + 1 of" ~ .(.y)),
         y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic(.(Sp_Name)))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
          )
  return(p)
  #ggsave(filename = paste0(.y, "_max_correlation_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})
pc17 <- map2(Gene_100217$Species_Cor, Gene_100217$Gene, function(.x, .y){
  #Plot_name <- paste(g, "_Plot_orig")
  Sp_Name <- meta_genome_sep_taxa %>% filter(species_id == .x) %>% select("Species") %>% as.character()
  p <- ggplot(Graph_tbl100217, aes(x = !!sym(.x), y = !!sym(.y), shape = factor(env))) +
    geom_point(size = 2, color = "orange") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(#title = bquote(.(.y) ~"Gene Counts versus Relative Abundance" ~ italic(.(Sp_Name))),
      x = bquote(log[10] ~ "Gene count + 1 of" ~ .(.y)),
      y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic(.(Sp_Name)))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
    )
  return(p)
  #ggsave(filename = paste0(.y, "_max_correlation_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})
pc80 <- Cor_plt_101380 <- map2(Gene_101380$Species_Cor, Gene_101380$Gene, function(.x, .y){
  #Plot_name <- paste(g, "_Plot_orig")
  Sp_Name <- meta_genome_sep_taxa %>% filter(species_id == .x) %>% select("Species") %>% as.character()
  p <- ggplot(Graph_tbl101380, aes(x = !!sym(.x), y = !!sym(.y), shape = factor(env))) +
    geom_point(size = 2, color = "orange") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(
      x = bquote(log[10] ~ "Gene count + 1 of" ~ .(.y)),
      y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic(.(Sp_Name)))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 5)
    )
  return(p)
  #ggsave(filename = paste0(.y, "_max_correlation_Graph.tiff"), plot = p, device = "tiff", width = 2 ,height= 2.0, dpi=300, units = "in")
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(patchwork)
#Make all corrplot figures
pg1 <- list(
  po60[[1]],
  pc60[[1]],
  po60[[2]],
  pc60[[2]],
  po60[[3]],
  pc60[[3]],
  po60[[4]],
  pc60[[4]],
  po60[[5]],
  pc60[[5]],
  po60[[6]],
  pc60[[6]],
  po60[[7]],
  pc60[[7]],
  po60[[8]],
  pc60[[8]],
  po60[[9]],
  pc60[[9]],
  po60[[10]],
  pc60[[10]]
)

Pg1 <- patchwork::wrap_plots(pg1, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg1.pdf", plot = Pg1, width = 7.5, height = 10, units = "in")

pg2 <- list(
  po60[[11]],
  pc60[[11]],
  po60[[12]],
  pc60[[12]],
  po60[[13]],
  pc60[[13]],
  po60[[14]],
  pc60[[14]],
  po60[[15]],
  pc60[[15]],
  po60[[16]],
  pc60[[16]],
  po60[[17]],
  pc60[[17]],
  po60[[18]],
  pc60[[18]],
  po60[[19]],
  pc60[[19]],
  po60[[20]],
  pc60[[20]]
)

Pg2 <- patchwork::wrap_plots(pg2, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg2.pdf", plot = Pg2, width = 7.5, height = 10, units = "in")

pg3 <- list(
  po60[[21]],
  pc60[[21]],
  po60[[22]],
  pc60[[22]],
  po60[[23]],
  pc60[[23]],
  po60[[24]],
  pc60[[24]],
  po60[[25]],
  pc60[[25]],
  po60[[26]],
  pc60[[26]],
  po60[[27]],
  pc60[[27]],
  po60[[28]],
  pc60[[28]],
  po60[[29]],
  pc60[[29]],
  po60[[30]],
  pc60[[30]]
)

Pg3 <- patchwork::wrap_plots(pg3, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg3.pdf", plot = Pg3, width = 7.5, height = 10, units = "in")

pg4 <- list(
  po60[[31]],
  pc60[[31]],
  po60[[32]],
  pc60[[32]],
  po60[[33]],
  pc60[[33]],
  po60[[34]],
  pc60[[34]],
  po60[[35]],
  pc60[[35]],
  po60[[36]],
  pc60[[36]],
  po60[[37]],
  pc60[[37]],
  po60[[38]],
  pc60[[38]],
  po60[[39]],
  pc60[[39]],
  po60[[40]],
  pc60[[40]]
)

Pg4 <- patchwork::wrap_plots(pg4, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg4.pdf", plot = Pg4, width = 7.5, height = 10, units = "in")

pg5 <- list(
  po60[[41]],
  pc60[[41]],
  po60[[42]],
  pc60[[42]],
  po60[[43]],
  pc60[[43]],
  po60[[44]],
  pc60[[44]],
  po60[[45]],
  pc60[[45]],
  po60[[46]],
  pc60[[46]],
  po28[[1]],
  po28[[1]]
)

Pg5 <- patchwork::wrap_plots(pg5, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg5.pdf", plot = Pg5, width = 7.5, height = 10, units = "in")

pg6 <- list(
  po78[[1]],
  pc78[[1]],
  po78[[2]],
  pc78[[2]],
  po78[[3]],
  pc78[[3]],
  po78[[4]],
  pc78[[4]],
  po78[[5]],
  pc78[[5]],
  po78[[6]],
  pc78[[6]],
  po78[[7]],
  pc78[[7]],
  po78[[8]],
  pc78[[8]],
  po78[[9]],
  pc78[[9]]
)

Pg6 <- patchwork::wrap_plots(pg6, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg6.pdf", plot = Pg6, width = 7.5, height = 10, units = "in")

pg7 <- list(
  po71[[1]],
  pc71[[1]],
  po71[[2]],
  pc71[[2]],
  po71[[3]],
  pc71[[3]],
  po71[[4]],
  pc71[[4]],
  po71[[5]],
  pc71[[5]],
  po71[[6]],
  pc71[[6]],
  po71[[7]],
  pc71[[7]],
  po71[[8]],
  pc71[[8]],
  po71[[9]],
  pc71[[9]],
  po71[[10]],
  pc71[[10]]
)

Pg7 <- patchwork::wrap_plots(pg7, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg7.pdf", plot = Pg7, width = 7.5, height = 10, units = "in")

pg8 <- list(
  po71[[11]],
  pc71[[11]],
  po71[[12]],
  pc71[[12]],
  po71[[13]],
  pc71[[13]],
  po71[[14]],
  pc71[[14]],
  po71[[15]],
  pc71[[15]],
  po71[[16]],
  pc71[[16]],
  po71[[17]],
  pc71[[17]],
  po71[[18]],
  pc71[[18]],
  po71[[19]],
  pc71[[19]],
  po71[[20]],
  pc71[[20]]
)

Pg8 <- patchwork::wrap_plots(pg8, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg8.pdf", plot = Pg8, width = 7.5, height = 10, units = "in")

pg9 <- list(
  po71[[21]],
  pc71[[21]],
  po71[[22]],
  pc71[[22]],
  po80[[1]],
  pc80[[1]],
  po80[[2]],
  pc80[[2]],
  po80[[3]],
  pc80[[3]],
  po80[[4]],
  pc80[[4]],
  po80[[5]],
  pc80[[5]],
  po17[[1]],
  pc17[[1]],
  po17[[2]],
  pc17[[2]],
  po17[[3]],
  pc17[[3]]
)

Pg9 <- patchwork::wrap_plots(pg9, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg9.pdf", plot = Pg9, width = 7.5, height = 10, units = "in")