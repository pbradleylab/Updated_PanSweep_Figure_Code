      #Creating Save files from analysis for Figures and Supplemental#
#Run PanSweep Analysis code
#Paths#
Gene_Reads_Paths <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/PanSweep_Analysis_Updated/data/midas_snv_merge/genes/"
Wrt_Gene_Reads_Path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Sig_Gene_reads.rds"
Wrt_Species_id_cor_and_orig <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Species_id_cor_and_orig.rds"
Egg_Save <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/uhgp_90_eggNOG.tsv"
#Create file for gene reads:
species_set <- unique(Genes_intrest_extr$Species_id)

Gene_reads <- map(species_set, ~ read_tsv(paste0(Gene_Reads_Paths,
                                                 .x,
                                                 "/",
                                                 .x,
                                                 ".genes_reads.tsv"),
                                          col_types=cols())) %>% setNames(species_set)
Sig_Gene_Reads <- map(species_set,  ~dplyr::filter(Gene_reads[[.x]], gene_id %in% Genes_intrest_extr$Gene_id) %>% pivot_longer(!gene_id, names_to = "run", values_to = "Gene_count") %>% pivot_wider(names_from = "gene_id", values_from = "Gene_count") %>% column_to_rownames("run")) %>% setNames(species_set)

#Save as RDS for later loading:
write_rds(Sig_Gene_Reads, file = Wrt_Gene_Reads_Path)

#Save Max correlation and originating species information for each gene
write_rds(Species_Cor_DF, file = Wrt_Species_id_cor_and_orig)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#For geneome Analysis:
Tmp <- uhgp_90_eggNOG
write_tsv(Tmp, file = Egg_Save)