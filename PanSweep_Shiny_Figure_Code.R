#PanSweep Shiny Application#
library(shiny)
library(umap)
library(vegan)
library(tidyverse)
library(plotly)
library(DT)
library(knitr)
library(kableExtra)
library(ggrepel)

PanSweep_Shiny <- function(loadData_Path){
  loadData <- read_rds(loadData_Path)
if(interactive()){
  
  ui <- navbarPage( "PanSweep",
                    
                    tabPanel("Analysis Report",
                             sidebarPanel(selectInput(inputId = "report",
                                                      label = "Choose Report:",
                                                      choices = c("Overall Report", "Species Report", "UHGP-90 Repeat ids", "UHGP-50 Repeat ids"))
                             ),
                             mainPanel(DTOutput("AnaR")),
                    ),
                    
                    tabPanel("eggNOG & Correlation Report",
                             mainPanel(DTOutput("eNR"))
                    ),
                    tabPanel("Ordination & Heatmap", 
                             fluidRow(column(2,
                                             selectInput(inputId = "ord_plt",
                                                         label = "Choose analysis:",
                                                         choices = c("UMAP","NMDS","PCoA", "NMDS_Zoom", "NMDS_NoM")),
                                             sliderInput("n_n",
                                                         label = "Number of n_neighbors:",
                                                         min = 2,
                                                         max = 10,                                            
                                                         value = 2),
                                             sliderInput("min_dist",
                                                         label = "min_dist:",
                                                         min = 0.1,
                                                         max = 0.9,                                            
                                                         value = 0.1,
                                                         step = 0.1),
                                             selectInput(inputId = "species_c",
                                                         label = "Choose species:",
                                                         choices = names(loadData$N.Sp_corr)),
                                             actionButton("reset", "Reset")
                                             
                             ),
                             column(10, plotOutput("ordination_plot"), 
                                    plotlyOutput("corrMax2"), 
                                    tableOutput("mtData")
                             ),
                             
                             
                             )
                    ),
                    tabPanel("NMDS", 
                             sidebarPanel(
                               selectInput(inputId = "species_c2",
                                           label = "Choose species:",
                                           choices = names(loadData$N.Sp_corr))
                             ),
                             mainPanel(plotlyOutput("NMDS"),
                                       plotOutput("stressPlot"))
                    )
  )
  
  server <- function(input, output, session) {

    clickData <- reactiveValues(x = vector(), y = vector(), cd = vector())
    #Builds up click data:
    observeEvent(event_data("plotly_click"), {
      new_clickData <- event_data("plotly_click")
      if (!is.null(new_clickData)){
        clickData$x <- c(clickData$x, new_clickData$x)
        clickData$y <- c(clickData$y, new_clickData$y)
        clickData$cd <- c(clickData$cd, new_clickData$customdata)
      }
    })
    #Resets click data:
    observeEvent(input$reset,{
      clickData$x <- vector()
      clickData$y <- vector()
      clickData$cd <- vector()
    })

    
    observeEvent(input$species_c, {
      req(loadData)
        updateSliderInput(session = session, "n_n", value = 2, min = 2, 
                          max = loadData$M.Sp_corr%>%
                            .[[paste0(input$species_c, sep='')]] %>%
                            {nrow(.)/3} %>%
                            ceiling()
        )
    }
      )
    
    output$ordination_plot <-  renderPlot({
      if(input$ord_plt == "UMAP"){
          n_n = input$n_n
          p <- loadData$U.Sp_corr %>%
            .[[paste0(input$species_c, sep='')]] %>%
            .[[paste0(input$n_n, sep='')]] %>%
            .[[paste0(input$min_dist, sep='')]] %>%
            .$"layout" %>%
            as.data.frame()%>%
            rownames_to_column() %>%
            rename("Gene_id" = "rowname") %>%
            left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id") %>%
            plot_ly(., x = ~V1, y = ~V2, type = 'scatter', mode = 'markers', 
                    color= ~replace(.$Predicted_taxonomic_group, is.na(.$Predicted_taxonomic_group), "NA"),
                    text = paste(.$Gene_id), 
                    customdata = ~paste(.$Gene_id)) %>%  
            layout(showlegend = TRUE,
                   plot_bgcolor = "#e5ecf6", 
                   xaxis = list( 
                     title = "0"),  
                   yaxis = list( 
                     title = "1"),
                   annotations = list(text = ~paste("n_neighbors:", input$n_n, "min_dist", input$min_dist), showarrow=FALSE ), ##MAKE PRERDY##
                   theme(plot.title.position = element_text(vjust = 0.5))
            )
          event_register(p, 'plotly_click')
        
        
      }
      else if(input$ord_plt == "NMDS"){
        # top_lft <- c("UHGG063307_00097", "UHGG228006_01743", "UHGG148794_01632", "UHGG258864_02598", "UHGG000216_02074")
        # NM <- c("UHGG210928_01151", "UHGG047117_02382", "UHGG047117_02376", "UHGG047117_02379", "UHGG041746_02390", "UHGG047117_02378",
        #         "UHGG047117_02377", "UHGG055467_01854", "UHGG047117_02383", "UHGG047117_02381", "UHGG047117_02380", "UHGG192308_01194", 
        #         "UHGG029873_02068", "UHGG032185_00468", "UHGG047117_02375")
        dt_p <-  loadData$N.Sp_corr %>%
            .[[paste0(input$species_c, sep='')]] %>%
            scores() %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            rename("Gene_id" = "rowname") %>%
            left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id")
        
          gp <- TEst %>% ggplot(aes(x = "NMDS1", y = "NMDS2", color = "Predicted_taxonomic_group")) +
                  geom_point() +
                  geom_text_repel(label = "Gene_id") +
          annotate("text", x = 0, y = 0, label = paste("stress", loadData$N.Stress[[paste0(input$species_c2, sep='')]]), size = 5, color = "black")
          gp
      }      
      else if(input$ord_plt == "NMDS_Zoom"){
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
          annotate("text", x = 0.3, y = 0.4, label = paste("stress", loadData$N.Stress$"100060"), size = 5, color = "black") 
         # xlim(-0.55, 0.7) +
        #  ylim(-0.35, 0.4)
        gp
      }
      else if(input$ord_plt == "NMDS_NoM"){
        loadData$N.Sp_corr %>%
          .[[paste0(input$species_c, sep='')]] %>%
          scores() %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          rename("Gene_id" = "rowname") %>%
          left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id") %>%
          plot_ly(x = ~NMDS1, y = ~NMDS2, type = "scatter", mode = "markers",  color= ~replace(.$Predicted_taxonomic_group, is.na(.$Predicted_taxonomic_group), "NA"),
                  text = paste(.$Gene_id), customdata = ~paste(.$Gene_id), showlegend = TRUE, legendgroup = "markers")%>%
          #add_text(textfont = list(color = "black"), textposition = ~Text_Position, showlegend = FALSE) %>%
          layout(plot_bgcolor = "#e5ecf6",
                 #xaxis = list(range = c(-0.25, 0.8)),
                 #yaxis = list(range = c(-0.11, 0.17)),
                 annotations = list(text = ~paste("stress", loadData$N.Stress[[paste0(input$species_c2, sep='')]]), showarrow=FALSE)
          )
      }
      
      else if(input$ord_plt == "PCoA"){
          loadData$P.Sp_corr %>%
            .[[paste0(input$species_c, sep='')]] %>%
            .$points %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            rename("Gene_id" = "rowname") %>%
            left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id") %>%
            plot_ly(x = ~V1, y = ~V2, type = 'scatter', mode = 'markers', 
                    color= ~replace(.$Predicted_taxonomic_group, is.na(.$Predicted_taxonomic_group), "NA"),
                    text = paste(.$Gene_id), 
                    customdata = ~paste(.$Gene_id)) %>%
            layout(showlegend = TRUE, plot_bgcolor = "#e5ecf6")
      }
    })
    
    output$eNR <- renderDT({
      loadData$Analysis_output$uhgp_90_eggNOG %>%
        select("Gene_id", "Species_id", "Species","Lineage_Shared", "cor_max_species", "Sp_rank", "Family_max_rank", "Fdrs", "cluster_id", "Predicted_taxonomic_group", 
               "Predicted_protein_name", "eggNOG_free_text_description") %>%
        mutate(Fdrs = format(signif(Fdrs, 3), scientific = TRUE))%>%
        filter(Gene_id %in% c("UHGG210793_02030", "UHGG047117_02380", "UHGG258864_02598", "UHGG152466_01649", "UHGG047117_02377")) %>%
        rename("Gene ID" = "Gene_id", "Species ID" = "Species_id", "Lineage Shared" = "Lineage_Shared", 
               "Species with Max Correlation Value" = "cor_max_species", "Species Correlation Rank" = "Sp_rank", "Family Correlation Rank" = "Family_max_rank",
               "Cluster ID" = "cluster_id","EggNOG Predicted taxonomic group" = "Predicted_taxonomic_group",
               "Predicted protein name" = "Predicted_protein_name", 
               "eggNOG free text description" = "eggNOG_free_text_description") %>%
        arrange("Species ID", desc(Fdrs), "Lineage Shared")
    })
    
    output$AnaR <- renderDT({
      if (input$report == "Overall Report"){
        loadData$Analysis_output$Analysis_report %>% datatable(colnames = NULL)
      } else if (input$report == "UHGP-90 Repeat ids"){
        loadData$Analysis_output$UHGP_90_cluster_id_summ
      } else if (input$report == "UHGP-50 Repeat ids"){
        loadData$Analysis_output$UHGP_50_cluster_id_summ
      } else if (input$report == "Species Report"){
        loadData$Analysis_output$Num_Sig_Genes_per_sp %>%
          rename("Species Id" = "Species_id") %>%
          rename("Number of Significant Genes in Species" = "n")
      }
    })
      
    
    output$NMDS <- renderPlotly({
        loadData$N.Sp_corr %>%
          .[[paste0(input$species_c2, sep='')]] %>%
          scores() %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          rename("Gene_id" = "rowname") %>%
          left_join(., loadData$Analysis_output$uhgp_90_eggNOG, "Gene_id") %>%
          plot_ly(x = ~NMDS1, y = ~NMDS2, type = "scatter", mode = "markers",  color= ~replace(.$Predicted_taxonomic_group, is.na(.$Predicted_taxonomic_group), "NA"),
                  text = paste(.$Gene_id))%>%
          layout(showlegend = TRUE,
                 annotations = list(text = ~paste("stress", loadData$N.Stress[[paste0(input$species_c2, sep='')]]), showarrow=FALSE))
    })
    
    output$stressPlot <- renderPlot({
        loadData$N.Sp_corr %>%
          .[[paste0(input$species_c2, sep='')]] %>%
          stressplot()
    })
    
    output$corrMax2 <- renderPlotly({
        corMax <- loadData$M.Sp_corr%>%
          .[[paste0(input$species_c, sep='')]] %>%
          as.data.frame() %>%
          mutate(across(everything(), ~  1 - .)) %>%
          as.matrix()
        
        if (length(clickData$cd > 0)){
          
          highlight_row <- isolate(clickData$cd)
          
          heatmap_plot <- plot_ly(x = rownames(corMax), y = colnames(corMax),
                                  z = corMax, zmin = 0, zmax = 1,
                                  type = "heatmap",
                                  colorscale = "Greys",
                                  colorbar = list(title = "Greys", y = 0.45, len = 0.45)
          )
          
          HiLite_mtx <-matrix(NA, nrow = nrow(corMax), ncol = ncol(corMax))
          rownames(HiLite_mtx) <- rownames(corMax)
          colnames(HiLite_mtx) <- colnames(corMax)
          HiLite_mtx[,paste0(highlight_row, sep='')] <- corMax[,paste0(highlight_row, sep='')]
          HiLite_mtx[paste0(highlight_row, sep=''),] <- corMax[paste0(highlight_row, sep=''),]
          
          heatmap_plot <- heatmap_plot %>% 
            add_trace(
              z = HiLite_mtx,
              type = "heatmap",
              colorscale = "Viridis",
              colorbar = list(title = "Viridis", y = 1, len = 0.45)
            ) %>%
            layout(
              title = "Jaccard Similarity"
            )
        } 
        else{
          heatmap_plot <-
            plot_ly(x = rownames(corMax), y = colnames(corMax),
                    z = corMax, zmin = 0, zmax = 1,
                    type = "heatmap") %>%
            layout(
              title = "Jaccard Similarity"
            )
        }
    })
    output$mtData <- function()({
        if(length(clickData$cd) > 0){ 
          
          Test_c <- isolate(clickData$cd)
          
          loadData$Analysis_output$uhgp_90_eggNOG %>%
            select("Gene_id", "Species_id", "Species","Lineage_Shared", "cor_max_species", "Fdrs", "cluster_id", "Predicted_taxonomic_group", 
                   "Predicted_protein_name", "eggNOG_free_text_description") %>%
            mutate(Fdrs = format(signif(Fdrs, 3), scientific = TRUE))%>%
            filter(Gene_id %in% Test_c) %>%
            rename("Gene ID" = "Gene_id", "Species ID" = "Species_id", "Lineage Shared" = "Lineage_Shared", 
                   "Species with Max Correlation Value" = "cor_max_species", "Cluster ID" = "cluster_id","EggNOG Predicted taxonomic group" = "Predicted_taxonomic_group",
                   "Predicted protein name" = "Predicted_protein_name", 
                   "eggNOG free text description" = "eggNOG_free_text_description") %>%
            knitr::kable("html") %>%
            kable_styling("striped", full_width = F)
          
        }
        else if (length(clickData$x) >0){ #x is the print out from the histogram and then print nothing if nothing else
          
          Test_c <- c(isolate(clickData$x), isolate(clickData$y))
          
          loadData$Analysis_output$uhgp_90_eggNOG %>%
            select("Gene_id", "Species_id", "Species","Lineage_Shared", "cor_max_species", "Fdrs", "cluster_id", "Predicted_taxonomic_group", 
                   "Predicted_protein_name", "eggNOG_free_text_description") %>%
            mutate(Fdrs = format(signif(Fdrs, 3), scientific = TRUE))%>%
            filter(Gene_id %in% Test_c) %>%
            rename("Gene ID" = "Gene_id", "Species ID" = "Species_id", "Lineage Shared" = "Lineage_Shared", 
                   "Species with Max Correlation Value" = "cor_max_species", "Cluster ID" = "cluster_id","EggNOG Predicted taxonomic group" = "Predicted_taxonomic_group",
                   "Predicted protein name" = "Predicted_protein_name", 
                   "eggNOG free text description" = "eggNOG_free_text_description") %>%
            knitr::kable("html") %>%
            kable_styling("striped", full_width = F)
    }
      })
  }
  
  shinyApp(ui = ui, server = server)}
}