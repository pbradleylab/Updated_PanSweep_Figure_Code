      #Creating Save files from analysis for Figures and Supplemental#
#Run PanSweep Analysis code
#Paths#
Gene_Reads_Paths <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/PanSweep_Analysis_Updated/data/midas_snv_merge/genes/"
Wrt_Gene_Reads_Path <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Sig_Gene_reads.rds"
Wrt_Species_id_cor_and_orig <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Species_id_cor_and_orig.rds"
Egg_Save <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/uhgp_90_eggNOG.tsv"
#Create file for gene reads:
species_set_corr <- unique(Genes_intr_extr$Species_id)
#Read in and merge gene counts:
Gene_reads <- purrr::map(species_set_corr, ~ {
  gpath <- file.path(path_to_read_counts, .x, paste0(.x, ".genes_reads.tsv"))
  if (is_compressed) {
    gpath <- paste0(gpath, ".lz4")
    tsv <- read_tsv_arrow(gpath)
  } else {
    tsv <- read_tsv(gpath, col_types=cols())
  }
  merge_columns_tbl(tsv, md=phylo_md2, fn=merge_fn_counts)
}) %>% setNames(species_set_corr)

Sig_Gene_reads <- purrr::map(species_set_corr, ~ {
  dplyr::filter(Gene_reads[[.x]], gene_id %in% Genes_intr_extr$Gene_id) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
}) %>% setNames(species_set_corr)
#Save as RDS for later loading:
write_rds(Sig_Gene_reads, file = Wrt_Gene_Reads_Path)

#Save Max correlation and originating species information for each gene
write_rds(Species_Cor_DF, file = Wrt_Species_id_cor_and_orig)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#For geneome Analysis:
Tmp <- uhgp_90_eggNOG
write_tsv(Tmp, file = Egg_Save)