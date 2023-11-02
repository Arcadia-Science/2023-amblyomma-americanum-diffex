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
vsda <- assay(vsd)
colnames(vsda) <- vsd$library_name

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
                 plotlyOutput("pca_plot"),
                 DTOutput("metadata_table")
               )
             )),
    tabPanel("DE Analysis", 
             sidebarLayout(
               sidebarPanel(
                 # Filtering options (log2FC, p-value, basemean)
                 numericInput("log2FoldChange_filter", "Filter by log2FC", min = -5, max = 5, value = 2),
                 numericInput("padj_filter", "Filter by p-value", min = 0, max = 1, value = 0.05),
                 numericInput("basemean_filter", "Filter by basemean", min = 0, max = 10000000, value = 20),
                 selectInput("contrast_variable", "Select Contrast", 
                             choices = c("sex_tissue"),
                             multiple = FALSE),
                 selectInput("condition1", "Condition 1", 
                             choices =  c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"), 
                             multiple = FALSE),
                 selectInput("condition2", "Condition 2", 
                             choices =  c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"), 
                             multiple = FALSE),
                 actionButton("get_diff_res", "Get Differential Results"),
                 downloadButton('download_data', 'Download Filtered Results'),
                 helpText(
                   "In this tab, you can upload a CSV file with your own information.",
                   "A join will be executed on the 'gene_name' column, so please ensure your file has a column named 'gene_name' with gene identifiers formatted as follows: ",
                   "'Amblyomma-americanum_evm.model.contig-XXXXX-X.X' (for example, Amblyomma-americanum_evm.model.contig-91628-1.2)",
                   "Once the file is uploaded and processed, you can select a column from your file to color the points in the output plot.",
                   "The output plot is a volcano plot with log2(Fold Change) on the x-axis and -log10(Adjusted p Value (BH)) on the y-axis,",
                   "where each point represents a gene. The coloring of points is based on the column you select."
                 ),
                 fileInput("file", "Choose CSV File", 
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                 radioButtons("color_by", "Color Points By", 
                              choices = c("Significance" = "significance", "Uploaded Variable" = "uploaded"),
                              selected = "significance"),
                 uiOutput("select_column_ui"), # This will show up only if a file is uploaded
                 actionButton("redraw_plots", "Re-draw Plots with New Color"),
               ),
               mainPanel(
                 plotlyOutput("volcano_plot"),
                 plotlyOutput("ma_plot"),
                 DTOutput("gene_table")
               )
             )),
    tabPanel("Gene",
             sidebarLayout(
               sidebarPanel(
                 selectizeInput("selected_gene", "Enter Gene Name:", 
                             #choices = rownames(vsda), multiple = FALSE),
                             choices = NULL, multiple = FALSE, options = NULL),
                 actionButton("plot_gene", "Plot Gene")
               ),
               mainPanel(
                 plotlyOutput("gene_boxplot")
               )
             )),
  )
)

# server logic ------------------------------------------------------------


server <- function(input, output, session) {

  # overview (PCA plot, metadata table) -------------------------------------

  # PCA plot
  output$pca_plot <- renderPlotly({
    selected_variable <- input$color_variable
    pcad <- plotPCA(vsd, intgroup = c("sex_tissue"), returnData = TRUE) %>%
      mutate(library_name = vsd$library_name) %>%
      left_join(metadata)
    percentVar <- round(100 * attr(pcad, "percentVar"))
    
    pca_plot <- ggplot(pcad, aes(PC1, PC2, color = !!sym(selected_variable), 
                                 label = library_name)) +
      geom_point(size = 3) +
      labs(x = paste0("PC1: ", percentVar[1], "% variance"),
           y = paste0("PC2: ", percentVar[2], "% variance"),
           color = input$color_variable) +
      coord_fixed() +
      theme_classic()
    
    ggplotly(pca_plot)
  })
  
  output$metadata_table <- renderDT({
    metadata %>%
      datatable(options = list(pageLength = 10))
  })
  
  # differential expression results (volcano plot, table) -------------
  redraw_trigger <- reactiveVal(0)
  
  observeEvent(input$redraw_plots, {
    redraw_trigger(redraw_trigger() + 1)
  })
  
  observe({
    updateSelectInput(session, "condition1", choices = c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"))
    updateSelectInput(session, "condition2", choices = c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"))
  })
  
  # extract deseq2 results and return df
  diff_results <- eventReactive(input$get_diff_res, {
    selected_contrasts <- input$contrast_variable
    selected_conditions1 <- input$condition1
    selected_conditions2 <- input$condition2
    
    # extract differential expression results based on selected contrasts and conditions
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

  # process the user-provided uploaded file
  uploaded_data <- reactive({
    req(input$file)
    read_csv(input$file$datapath)
  })

  # dynamic UI for selecting a column from the uploaded data
  output$select_column_ui <- renderUI({
    req(input$file)  # Only show if a file is uploaded
    selectInput("color_column", "Select Column to Color By", choices = names(uploaded_data()))
  })
  
  # Combine the diff_results() data with the uploaded_data()
  combined_data <- reactive({
    req(input$color_by)  # Ensure a color_by selection is made
    diff_results_df <- diff_results()
    if (input$color_by == "uploaded") {
      req(input$color_column)  # Ensure a color column is selected
      req(uploaded_data())  # Ensure a file is uploaded
      # Join the uploaded data with the differential expression results data
      joined_data <- left_join(diff_results_df, uploaded_data(), by = "gene_name")
      joined_data
    } else {
      # If not coloring by uploaded data, use the original data
      diff_results_df
    }
  })
  
  # volcano plot
  output$volcano_plot <- renderPlotly({
    redraw_trigger()  # This line will cause this function to re-run each time the button is clicked    
    combined_df <- combined_data()  # Get the combined data
    color_var <- if(input$color_by == "significance") {
      combined_df$significant
    } else {
      combined_df[[input$color_column]]
    }
    
    volcano_plot <- ggplot(combined_df, aes(x = log2FoldChange, y = -log10(padj), 
                                                color = color_var, label = gene)) +
      geom_point(alpha = 0.5) +
      labs(x = "log2(Fold Change)", y = "-log10(Adjusted p Value (BH))") +
      theme_classic()
    
    ggplotly(volcano_plot)
  })
  
  # MA plot
  output$ma_plot <- renderPlotly({
    redraw_trigger()  # This line will cause this function to re-run each time the button is clicked
    combined_df <- combined_data()  # Get the combined data
    color_var <- if(input$color_by == "significance") {
      combined_df$significant
    } else {
      combined_df[[input$color_column]]
    }
    
    ma_plot <- ggplot(combined_df, aes(x = log2(baseMean), y = log2FoldChange, 
                                          color = color_var, label = gene)) +
      geom_point(alpha = 0.5) +
      labs(x = "log2(Mean Count)", y = "log2(Fold Change)") +
      theme_classic()
    
    ggplotly(ma_plot)
  })
  
  # significant gene table
  output$gene_table <- renderDT({
    # diff_results_df <- diff_results()
    # # determine what to color by (default significance, otherwise user-uploaded column)
    # color_var <- if(input$color_by == "significance") {
    #   data$significant
    # } else {
    #   req(uploaded_data())
    #   joined_data <- left_join(data, uploaded_data(), by = "gene_name")
    #   joined_data[[input$color_column]]
    # }
    ds_results <- diff_results() %>%
      filter(significant == "significant")
    
    ds_results %>%
      datatable(options = list(pageLength = 10))
  })
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste0('filtered_results',
             '_log2FC', input$log2FoldChange_filter,
             '_basemean', input$basemean_filter,
             '_padj', input$padj_filter,
             '_', input$condition1,
             '_', input$condition2,
             '.tsv')
    },
    content = function(file) {
      ds_results <- diff_results() %>%
        filter(significant == "significant")
      write_tsv(ds_results, file)
    }
  )
  
  # gene-specific visualization ------------------------------------------
  # update the choices for selected_gene -- using server-side selectize speeds up gene selection 
  # bc we have 30k genes to choose from
  observe({
    gene_choices <- rownames(vsda)
    updateSelectizeInput(session, "selected_gene", choices = gene_choices, server = TRUE)
  })

  observe({
    if (!is.null(input$selected_gene)) {
      gene_conditions <- colnames(vsd)
      updateSelectInput(session, "gene_condition", choices = gene_conditions, selected = gene_conditions)
    } else {
      updateSelectInput(session, "gene_condition", choices = NULL, selected = NULL)
    }
  })

  # plot gene count boxplot
  output$gene_boxplot <- renderPlotly({
    selected_gene <- input$selected_gene
    
    if (!is.null(selected_gene)) {
      gene_data <- data.frame(gene_count = vsda[selected_gene, ]) %>%
        rownames_to_column("library_name") %>%
        left_join(metadata, by = "library_name")

      gene_boxplot <- ggplot(gene_data, aes(x = sex_tissue, y = gene_count)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(alpha = 0.5, aes(label = library_name)) +
        labs(x = "Condition", y = "Normalized Count") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))

      ggplotly(gene_boxplot)
    }
  })
}


# run server --------------------------------------------------------------

shinyApp(ui = ui, server = server)
