loadData_Path = "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/PanSweep_Analysis_Updated/PanSweep_Analysis_Output_2024-09-24/PanSweep_Analysis_Output.rds"
Save_Pth <- "~/Documents/Bradley_Lab/MIDAS_Analysis_Main_Folder/PanSweep_Updated_Data_Analysis/Updated_PanSweep_Figure_Code/Gene_Corr_Graphs/"
loadData <- read_rds(loadData_Path)

Test <- loadData$N.Sp_corr %>%
  .$"100060" %>%
  scores() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename("Gene_id" = "rowname") %>%
  left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id")

gp <- Test %>% ggplot(aes(x = NMDS1, y = NMDS2, color = Predicted_taxonomic_group)) +
  geom_point() +
  geom_text_repel(label = Test$Gene_id) +
  annotate("text", x = 0.3, y = 0.4, label = paste("stress", loadData$N.Stress$"100060"), size = 5, color = "black") +
  xlim(-0.55, 0.7) +
  ylim(-0.35, 0.4) +
  theme(panel.background = element_rect(fill = "#e5ecf6"))
gp

ggsave(paste0(Save_Pth,"fig_2.tiff"), plot = gp, device = "tiff", width = 16,height=7, dpi=300, units = "in")