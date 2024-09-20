                        #Graphs for Figure 1 and Supplemental#
library(tidyverse)
library(readr)
library(ggplot2)
                                    #Paths#
path_for_genomes_all_metadata <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Metadata/Genome_metadata/metadata.tsv"

path_for_Species_Abundance <-"~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/PanSweep_Analysis_Updated/data/species/species_marker_read_counts.tsv"

path_for_sample_metadata <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Metadata/Sample_Metadata/cirrhosis-metadata_mod.tsv"

path_for_gene_reads <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Sig_Gene_reads.rds"

path_to_cor_org_species <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Species_id_cor_and_orig.rds"

Base_path2 <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/"
Egg_load <- paste0(Base_path2, "uhgp_90_eggNOG.tsv")

Graph_Save_Dir <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Gene_Corr_Graphs/Figure1"

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
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5), labels = c("0","0.5", "1.0","1.5", "2.0","2.5")) +
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
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5), labels = c("0","0.5", "1.0","1.5", "2.0","2.5")) +
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
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5), labels = c("0","0.5", "1.0","1.5", "2.0","2.5")) +
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
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5), labels = c("0","0.5", "1.0","1.5", "2.0","2.5")) +
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
setwd(Supp_Graph_Save_Dir)
        #Load in all gene max correlation species and originating species#
Max_org_species <- read_rds(path_to_cor_org_species)
Species_for_graph <- Max_org_species %>% filter(mark == "max")


                 #Separate out runs and combine with species#
Species_Abd_join <- Species_Abd %>% rownames_to_column("Run")
Slct_Gns100060 <- Sig_Gene_Reads$'100060' %>% rownames_to_column("Run")
Slct_Gns100078 <- Sig_Gene_Reads$'100078' %>% rownames_to_column("Run")
Slct_Gns100271 <- Sig_Gene_Reads$'100271' %>% rownames_to_column("Run")
Slct_Gns102528 <- Sig_Gene_Reads$'102528' %>% rownames_to_column("Run")

Graph_tbl100060 <- Slct_Gns100060 %>% left_join(Species_Abd_join, by = 'Run') %>% 
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>% 
  left_join(sample_meta, by = c("Run" = "sample"))
Graph_tbl100078 <- Slct_Gns100078 %>% left_join(Species_Abd_join, by = 'Run') %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("Run" = "sample"))  
Graph_tbl100271 <- Slct_Gns100271 %>% left_join(Species_Abd_join, by = 'Run') %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("Run" = "sample"))  
Graph_tbl102528 <- Slct_Gns102528 %>% left_join(Species_Abd_join, by = 'Run') %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("Run" = "sample"))
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
po271 <- map(Gene_V_100271, function(g){
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
po2528 <-map(Gene_V_102528, function(g){
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
################################################################################
#Make correlating species to gene graphs
setwd(Supp_Graph_Save_Dir)

#Separate out runs and combine with species#
Species_Abd_join <- Species_Abd %>% rownames_to_column("Run")
Slct_Gns100060 <- Sig_Gene_Reads$'100060' %>% rownames_to_column("Run")
Slct_Gns100078 <- Sig_Gene_Reads$'100078' %>% rownames_to_column("Run")
Slct_Gns100271 <- Sig_Gene_Reads$'100271' %>% rownames_to_column("Run")
Slct_Gns102528 <- Sig_Gene_Reads$'102528' %>% rownames_to_column("Run")

Graph_tbl100060 <- Slct_Gns100060 %>% left_join(Species_Abd_join, by = 'Run') %>% 
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>% 
  left_join(sample_meta, by = c("Run" = "sample"))
Graph_tbl100078 <- Slct_Gns100078 %>% left_join(Species_Abd_join, by = 'Run') %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("Run" = "sample"))  
Graph_tbl100271 <- Slct_Gns100271 %>% left_join(Species_Abd_join, by = 'Run') %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("Run" = "sample"))  
Graph_tbl102528 <- Slct_Gns102528 %>% left_join(Species_Abd_join, by = 'Run') %>%
  mutate(across(-c("Run"), ~ log(.x +1, base = 10))) %>%
  left_join(sample_meta, by = c("Run" = "sample"))

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

pc271 <- map2(Gene_100271$Species_Cor, Gene_100271$Gene, function(.x, .y){
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

pc2528 <- map2(Gene_102528$Species_Cor, Gene_102528$Gene, function(.x, .y){
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
  pc60[[25]]
)

Pg3 <- patchwork::wrap_plots(pg3, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg3.pdf", plot = Pg3, width = 7.5, height = 10, units = "in")

pg4 <- list(
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
  pc78[[8]]
)

Pg4 <- patchwork::wrap_plots(pg4, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg4.pdf", plot = Pg4, width = 7.5, height = 10, units = "in")

pg5 <- list(
  po271[[1]],
  pc271[[1]],
  po271[[2]],
  pc271[[2]],
  po271[[3]],
  pc271[[3]],
  po271[[4]],
  pc271[[4]],
  po271[[5]],
  pc271[[5]],
  po271[[6]],
  pc271[[6]],
  po271[[7]],
  pc271[[7]],
  po2528[[1]],
  pc2528[[1]]
)

Pg5 <- patchwork::wrap_plots(pg5, nrow = 5, ncol = 4, byrow = TRUE)
ggsave(filename ="Pg5.pdf", plot = Pg5, width = 7.5, height = 10, units = "in")