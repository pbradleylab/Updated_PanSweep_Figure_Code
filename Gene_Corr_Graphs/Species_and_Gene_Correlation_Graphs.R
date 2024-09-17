                        #Graphs for Figure 1 and Supplemental#
library(tidyverse)
library(readr)
library(ggplot2)
                                    #Paths#
path_for_genomes_all_metadata <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/UHGP_files/metadata.tsv"

path_for_Species_Abundance <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/Sp_Abundance_Files/species_marker_read_counts.tsv"

path_for_sample_metadata <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/cirrhosis-analysis/data/manual/phylogenize/cirrhosis-metadata.tsv"

path_for_gene_reads <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/Species_GeneC_Graph_Files/Sig_Gene_reads.rds"
Graph_Save_Dir <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/Species_GeneC_Graph_Files"
path_to_cor_org_species <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/MIDAS_Cirrhisis_Analysis/Species_GeneC_Graph_Files/Species_id_cor_and_orig.rds"

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
meta_genome_sep_taxa %>% filter(species_id == 101444) %>% select("Species") %>% print()
meta_genome_sep_taxa %>% filter(species_id == 100060) %>% select("Species") %>% print()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                          #Load in Species abundance#
Species_Abd <-read_tsv(path_for_Species_Abundance) %>% 
  pivot_longer(!species_id, names_to = "run", values_to = "species_count") %>% 
  pivot_wider(names_from = "species_id", values_from = "species_count") %>% column_to_rownames("run")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                              #Load in Gene Reads#
Sig_Gene_Reads <- read_rds(file = path_for_gene_reads)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                            #Load in  Sample Metadata#
sample_meta <- read_tsv(path_for_sample_metadata) %>% select(c("sample","env"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                            #Create Figure 1 Graphs#
SpAb_VA <- Species_Abd %>% select('101444')%>% rownames_to_column("Run")
colnames(SpAb_VA) <- c("Run","V.atypica")
SpAb_Le <- Species_Abd %>% select('100060') %>% rownames_to_column("Run")
colnames(SpAb_Le) <- c("Run","Le")

Slct_Gns <- Sig_Gene_Reads$'100060' %>% rownames_to_column("Run")

Graph_tbl <- Slct_Gns %>% left_join(SpAb_VA, by = 'Run') %>% left_join(SpAb_Le, by = 'Run')

log_Plot <- Graph_tbl %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10)))

log_grph_df_meta <- log_Plot %>% left_join(sample_meta, by = c("Run" = "sample"))



Cont_Gene <- log_grph_df_meta %>% select(c("Run", "Le", "V.atypica", "UHGG047117_02376", "env"))
Org_Gene <- log_grph_df_meta %>% select(c("Run", "Le", "V.atypica", "UHGG152466_01649", "env"))

Plot1 <- ggplot(Cont_Gene, aes(x = V.atypica, y = UHGG047117_02376, shape = factor(env))) +
  geom_point(size = 5, color = "orange") +
  geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
  guides(
    color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
    shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
  ) +
  #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
  scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  labs(
    x = NULL,
    y = NULL) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20))





Plot2 <- ggplot(Cont_Gene, aes(x = Le, y = UHGG047117_02376, shape = factor(env))) +
  geom_point(size = 5, color = "deepskyblue4") +
  geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
  guides(
    color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
    shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
  ) +
  #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
  scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  labs(
    x = NULL,
    y = NULL) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20))

Plot3 <- ggplot(Org_Gene, aes(x = V.atypica, y = UHGG152466_01649, shape = factor(env))) +
  geom_point(size = 5, color = "orange") +
  geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
  guides(
    color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
    shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
  ) +
  #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
  scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  labs(
    x = NULL,
    y = NULL) +
  theme( legend.position = "none",
         axis.text.x = element_text(size = 20), 
         axis.text.y = element_text(size = 20))

Plot4 <- ggplot(Org_Gene, aes(x = Le, y = UHGG152466_01649, shape = factor(env))) +
  geom_point(size = 5, color = "deepskyblue4") +
  geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
  guides(
    color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
    shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
  ) +
  #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
  scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  labs(
    x = NULL,
    y = NULL) +
  theme(  legend.position = "none",
          axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 20))

Plot1
Plot2
Plot3
Plot4
setwd(Graph_Save_Dir)
ggsave("2No_Title_V.atypica_and_Contaminate_Gene.tiff", plot = Plot1, device = "tiff", width = 7.5,height=6, dpi=300, units = "in")
ggsave("2No_Title_L.eligens_and_Contaminate_Gene.tiff", plot = Plot2, device = "tiff", width = 7.5,height=6, dpi=300, units = "in")
ggsave("2No_Title_V.atypica_and_Non-Contaminate_Gene.tiff", plot = Plot3, device = "tiff", width = 7.5,height=6, dpi=300, units = "in")
ggsave("2No_Title_L.eligens_and_Non-Contaminate_Gene.tiff", plot = Plot4, device = "tiff", width = 7.5,height=6, dpi=300, units = "in")
#################################################################################
                              #Make Supplemental Graphs#
        #Load in all gene max correlation species and originating species#
Max_org_species <- read_rds(path_to_cor_org_species)
Species_for_graph <- Max_org_species %>% filter(mark == "max")


                 #Separate out runs and combine with species#
Species_Abd_join <- Species_Abd %>% rownames_to_column("Run")
Slct_Gns100060 <- Sig_Gene_Reads$'100060' %>% rownames_to_column("Run")
Slct_Gns100271 <- Sig_Gene_Reads$'100271' %>% rownames_to_column("Run")
Slct_Gns103694 <- Sig_Gene_Reads$'103694' %>% rownames_to_column("Run")
Slct_Gns102580 <- Sig_Gene_Reads$'102580' %>% rownames_to_column("Run")

Graph_tbl100060 <- Slct_Gns100060 %>% left_join(Species_Abd_join, by = 'Run') %>% 
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>% 
  left_join(sample_meta, by = c("Run" = "sample"))
Graph_tbl100271 <- Slct_Gns100271 %>% left_join(Species_Abd_join, by = 'Run') %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("Run" = "sample"))  
Graph_tbl103694 <- Slct_Gns103694 %>% left_join(Species_Abd_join, by = 'Run') %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("Run" = "sample"))  
Graph_tbl102580 <- Slct_Gns102580 %>% left_join(Species_Abd_join, by = 'Run') %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("Run" = "sample"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Make a vector of gene names by species#
Gene_100060 <- Species_for_graph %>% filter(Species == 100060)
Gene_V_100060 <- Gene_100060$Gene
                    #Make graphs for L. eligens (100060)#
plot_100060 <- map(Gene_V_100060, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
    ggplot(Graph_tbl100060, aes(x = !!sym("100060"), y = !!sym(g), shape = factor(env))) +
     geom_point(size = 5, color = "deepskyblue4") +
     geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
     guides(
       color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
       shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
     ) +
     #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
     scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
     labs(title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("L. eligens")),
          x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
          y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("L. eligens"))
            ) +
     theme(legend.position = "none",
           axis.text.x = element_text(size = 20), 
           axis.text.y = element_text(size = 20))
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                #Make a vector of gene names by species#
Gene_100271 <- Species_for_graph %>% filter(Species == 100271)
Gene_V_100271 <- Gene_100271$Gene
                  #Make graphs for Roseburia inulinivorans (100271)#
plot_100271 <- map(Gene_V_100271, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
  ggplot(Graph_tbl100271, aes(x = !!sym("100271"), y = !!sym(g), shape = factor(env))) +
    geom_point(size = 5, color = "deepskyblue4") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("Roseburia inulinivorans")),
         x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
         y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("Roseburia inulinivorans"))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 20))
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Make a vector of gene names by species#
Gene_103694 <- Species_for_graph %>% filter(Species == 103694)
Gene_V_103694 <- Gene_103694$Gene
#Make graphs for "Agathobacter faecis" (103694)#
plot_103694 <- map(Gene_V_103694, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
  ggplot(Graph_tbl103694, aes(x = !!sym("103694"), y = !!sym(g), shape = factor(env))) +
    geom_point(size = 5, color = "deepskyblue4") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("Agathobacter faecis")),
         x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
         y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("Agathobacter faecis"))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 20))
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Make a vector of gene names by species#
Gene_102580 <- Species_for_graph %>% filter(Species == 102580)
Gene_V_102580 <- Gene_102580$Gene
#Make graphs for "Lachnospira sp000436475" (102580)#
plot_102580 <- map(Gene_V_102580, function(g){
  #Plot_name <- paste(g, "_Plot_orig")
  ggplot(Graph_tbl102580, aes(x = !!sym("102580"), y = !!sym(g), shape = factor(env))) +
    geom_point(size = 5, color = "deepskyblue4") +
    geom_smooth(method = "loess", color = "gray28", aes(group = 1)) +
    guides(
      color = guide_legend(override.aes = list(linetype = 0)),  # Combined legend title
      shape = guide_legend(override.aes = list(linetype = 0))  # Combined legend title
    ) +
    #scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")) +
    #scale_y_continuous(limits = c(-0.1, 0.5), breaks = seq(-0.1, 0.5, 0.1), labels = c("-0.1", "0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
    labs(title = bquote(.(g) ~"Gene Counts versus Relative Abundance" ~ italic("Lachnospira sp000436475")),
         x = bquote(log[10] ~ "Gene count + 1 of" ~ .(g)),
         y = bquote(log[10] ~ "Relative Abundance + 1 of" ~ italic("Lachnospira sp000436475"))
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 20))
})
################################################################################
#Make contaminate gene graphs
Test <- map2(Gene_100060$Gene, Gene_100060$Species_Cor, ~.x + .y {
  
})