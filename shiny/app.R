# set up environment ------------------------------------------------------

library(shiny)
library(DESeq2)
library(tidyverse)
library(plotly)
library(DT)
# read in input data ------------------------------------------------------

dds <- readRDS("input_data/dds_sex_tissue.RDS")
ds <- readRDS("input_data/ds_sex_tissue.RDS")
metadata <- read_tsv("input_data/metadata.tsv")
annotation <- read_tsv("input_data/Amblyomma-americanum_combined.tsv")

# process input data ------------------------------------------------------

# process metadata for diffex vars and other user-relevant information
# TODO: consider just writing this metadata out and reading this in already formatted
metadata <- metadata %>%
  filter(tissue != "cell_line") %>% # filter out samples that will make it so model is not full rank
  filter(!study_title %in% "Amblyomma americanum RNA-seq following E. coli treatment") %>% # rm e. coli study, too batched
  filter(!library_name %in% c("Um", "Uf", "Im", "If")) %>% # filter e chaf exposed samples
  select(library_name, experiment_title, study_title, sex, host_meal, tissue, blood_meal_hour, blood_meal_hour_range, total_spots) %>%
  group_by(library_name) %>%
  mutate(total_spots = sum(total_spots)) %>%
  mutate(blood_meal_hour_range = factor(blood_meal_hour_range, levels = c("0", "12_48", "72_144", "168_264", "none"))) %>%
  distinct() %>%
  mutate(sex_tissue = paste0(sex, "_x_", tissue))

# process annotation information so it matches the default EVM annot names 
annotations <- annotation %>%
  mutate(gene = gsub("Amblyomma-americanum_", "", gene_name), .before = gene_name,
         gene = gsub("-", "_", gene),
         gene = gsub(".model", ".TU", gene)) 

# perform variance stabilize transformation for PCA and some other plots
vsd <- vst(ds, blind = FALSE)

# ui logic ----------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Differential Expression Explorer"),
  tabsetPanel(
    tabPanel("PCA Plot", 
             sidebarLayout(
               sidebarPanel(
                 # Add a drop-down menu to select the coloring variable
                 selectInput("color_variable", "Select Coloring Variable", 
                             choices = c("sex", "tissue", "blood_meal_hour_range", "study_title")),
               ),
               mainPanel(
                 plotlyOutput("pca_plot")
               )
             )),
    tabPanel("DE Analysis", 
             sidebarLayout(
               sidebarPanel(
                 # Filtering options (log2FC, p-value, basemean)
                 numericInput("log2FoldChange_filter", "Filter by log2FC", min = -5, max = 5, value = 2),
                 numericInput("padj_filter", "Filter by p-value", min = 0, max = 1, value = 0.05),
                 numericInput("basemean_filter", "Filter by basemean", min = 0, max = 10000000, value = 100),
                 selectInput("contrast_variable", "Select Contrast", 
                             choices = c("sex_tissue"),
                             multiple = FALSE),
                 selectInput("condition1", "Condition 1", 
                             choices =  c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"), 
                             multiple = FALSE),
                 selectInput("condition2", "Condition 2", 
                             choices =  c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"), 
                             multiple = FALSE),
                 actionButton("get_diff_res", "Get Differential Results")
               ),
               mainPanel(
                 plotlyOutput("ma_plot"),
                 DTOutput("gene_table")
               )
             ))
  )
)

# server logic ------------------------------------------------------------


server <- function(input, output, session) {

  # PCA plot
  output$pca_plot <- renderPlotly({
    selected_variable <- input$color_variable
    pcad <- plotPCA(vsd, intgroup = c("sex_tissue"), returnData = TRUE) %>%
      mutate(library_name = vsd$library_name) %>%
      left_join(metadata)
    percentVar <- round(100 * attr(pcad, "percentVar"))
    
    pca_plot <- ggplot(pcad, aes(PC1, PC2, color = !!sym(selected_variable), text = library_name)) +
      geom_point(size = 3) +
      labs(x = paste0("PC1: ", percentVar[1], "% variance"),
           y = paste0("PC2: ", percentVar[2], "% variance"),
           color = input$color_variable) +
      coord_fixed() +
      theme_classic()
    
    ggplotly(pca_plot)
  })
  
  # differential expression results (MA plot, table)
  observe({
    updateSelectInput(session, "condition1", choices = c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"))
    updateSelectInput(session, "condition2", choices = c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"))
  })
  
  # Function to convert DESeqResults object to a data frame
  diff_results <- eventReactive(input$get_diff_res, {
    selected_contrasts <- input$contrast_variable
    selected_conditions1 <- input$condition1
    selected_conditions2 <- input$condition2
    
    # selected_contrasts <- "sex_tissue"
    # selected_conditions1 <- "female_x_salivary_gland"
    # selected_conditions2 <- "female_x_whole"
    # Extract differential expression results based on selected contrasts and conditions
    ds_results <- results(ds, contrast = c(selected_contrasts, selected_conditions1, selected_conditions2))
    
    # Convert DESeqResults object to a data frame
    ds_results_df <- ds_results %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      left_join(annotations)

    # demarcate significant
    ds_results_significant <- ds_results_df %>%
      mutate(significant = ifelse(abs(log2FoldChange) >= input$log2FoldChange_filter & 
                                  padj <= input$padj_filter & 
                                  baseMean >= input$basemean_filter, "significant", "not significant"))
    
    # Return the filtered data frame
    ds_results_significant
  })
  
  # Render the MA plot
  output$ma_plot <- renderPlotly({
    ma_plot <- ggplot(diff_results(), aes(x = log2FoldChange, y = -log10(padj), 
                                          color = significant, label = gene)) +
      geom_point() +
      labs(x = "log2FC", y = "-log10(padj)") +
      theme_classic()
    
    ggplotly(ma_plot)
  })
  
  # Sample gene table
  output$gene_table <- renderDT({
    ds_results <- diff_results() %>%
      filter(significant == "significant")
    
    ds_results %>%
      datatable(options = list(pageLength = 10))
  })
}

# run server --------------------------------------------------------------

shinyApp(ui = ui, server = server)
