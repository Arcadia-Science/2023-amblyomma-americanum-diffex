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
annotations <- read_tsv("input_data/Amblyomma-americanum_combined.tsv")

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
annotations <- annotations %>%
  mutate(gene = gsub("Amblyomma-americanum_", "", gene_name), .before = gene_name,
         gene = gsub("-", "_", gene),
         gene = gsub(".model", ".TU", gene)) %>%
  # define a consensus annotation. Use eggnog first, then kegg if eggnog is blank, then deepsig, or nothing
  mutate(consensus_annotation = ifelse(!is.na(egg_Description), egg_Description,
                                       ifelse(!is.na(KO_definition), KO_definition, 
                                              ifelse(!is.na(deepsig_feature), deepsig_feature, "unknown function")))) %>%
  mutate(combined_gene = paste0(gene_name, "; ", consensus_annotation))

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
                 helpText(
                   "Optional: upload a CSV file with your own information.",
                   "A join will be executed on the 'gene_name' column, so please ensure your file has a column named 'gene_name' with gene identifiers formatted as follows: ",
                   "'Amblyomma-americanum_evm.model.contig-XXXXX-X.X' (for example, Amblyomma-americanum_evm.model.contig-91628-1.2)",
                   "Once the file is uploaded, you can select a column from your file to color the points in the output plot.",
                   "If you do not upload a file, points will be colored by significance of differential expression results, inferrered from the filtering criteria specified above."
                 ),
                 fileInput("file", "Choose CSV File", 
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                 radioButtons("color_by", "Color Points By", 
                              choices = c("Significance" = "significance", "Uploaded Variable" = "uploaded"),
                              selected = "significance"),
                 uiOutput("select_column_ui"), # This will show up only if a file is uploaded
                 actionButton("get_diff_res", "Get Differential Results"),
                 downloadButton('download_data', 'Download Filtered Results')
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
                                choices = NULL, multiple = FALSE, options = NULL),
                 actionButton("plot_gene", "Plot Gene")
               ),
               mainPanel(
                 plotlyOutput("gene_boxplot")
               )
             )),
    tabPanel("Expression by Condition",
             sidebarLayout(
               sidebarPanel(
                 helpText(
                   "This tab reports on expression per sample.",
                   "It uses counts that have been normalized for sequencing depth and undergone variance stabilizing transformation (VST).",
                   "VST outputs transformed data on the log2 scale.",
                   "With this method, negative values indicate a count of less than 1, e.g. that the gene is not detected in the given sample."
                 ),
                 selectInput("condition_input", "Select Condition:",
                             choices = unique(metadata$sex_tissue)),
                 numericInput("expression_threshold", "Expression Threshold:", value = 0, min = -5),
                 actionButton("view_expression", "View Expression")
               ),
               mainPanel(
                 tableOutput("expression_table"),
                 plotlyOutput("expression_plot")
               )
             ))
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
  observe({
    updateSelectInput(session, "condition1", choices = c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"))
    updateSelectInput(session, "condition2", choices = c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"))
  })
  
  # Process the user-provided uploaded file
  uploaded_data <- reactive({
    req(input$file)
    read_csv(input$file$datapath)
  })
  
  # Dynamic UI for selecting a column from the uploaded data
  output$select_column_ui <- renderUI({
    req(input$file)  # Only show if a file is uploaded
    selectInput("color_column", "Select Column to Color By", choices = names(uploaded_data()))
  })
  
  # Extract deseq2 results and return df
  diff_results <- eventReactive(input$get_diff_res, {
    selected_contrasts <- input$contrast_variable
    selected_conditions1 <- input$condition1
    selected_conditions2 <- input$condition2
    
    # Extract differential expression results based on selected contrasts and conditions
    ds_results <- results(ds, contrast = c(selected_contrasts, selected_conditions1, selected_conditions2))
    ds_results_df <- ds_results %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      left_join(annotations)
    
    # If a file is uploaded and "Uploaded Variable" is selected for coloring
    if (!is.null(input$file) && input$color_by == "uploaded") {
      req(input$color_column)
      uploaded_data <- read_csv(input$file$datapath)
      ds_results_df <- left_join(ds_results_df, uploaded_data(), by = "gene_name")
    }
    
    # Demarcate significant
    ds_results_significant <- ds_results_df %>%
      mutate(significant = ifelse(abs(log2FoldChange) >= input$log2FoldChange_filter & 
                                  padj <= input$padj_filter & 
                                  baseMean >= input$basemean_filter, "significant", "not significant"))
    
    # Return the filtered data frame
    ds_results_significant
  })
  
  # volcano plot
  output$volcano_plot <- renderPlotly({
    diff_results_df <- diff_results()  # Get the combined data
    color_var <- if(input$color_by == "significance") {
      diff_results_df$significant
    } else {
      diff_results_df[[input$color_column]]
    }
    
    volcano_plot <- ggplot(diff_results_df, aes(x = log2FoldChange, y = -log10(padj), 
                                                color = color_var, 
                                                label = combined_gene)) +
      geom_point(alpha = 0.5) +
      labs(x = "log2(Fold Change)", y = "-log10(Adjusted p Value (BH))") +
      theme_classic()
    
    ggplotly(volcano_plot)
  })
  
  # MA plot
  output$ma_plot <- renderPlotly({
    diff_results_df <- diff_results()  # Get the combined data
    color_var <- if(input$color_by == "significance") {
      diff_results_df$significant
    } else {
      diff_results_df[[input$color_column]]
    }
    
    ma_plot <- ggplot(diff_results_df, aes(x = log2(baseMean), y = log2FoldChange, 
                                          color = color_var, label = combined_gene)) +
      geom_point(alpha = 0.5) +
      labs(x = "log2(Mean Count)", y = "log2(Fold Change)") +
      theme_classic()
    
    ggplotly(ma_plot)
  })
  
  # significant gene table
  output$gene_table <- renderDT({
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
  
  # expression by condition -------------------------------------------------
  expression_data <- reactive({
    req(input$view_expression)  # Ensure the button is clicked
    req(input$condition_input)  # Ensure a contrast is selected
    # Filter the expression data for the selected contrast
    # Assuming vsda holds the expression data where rows are genes and columns are samples
    selected_samples <- metadata %>% filter(sex_tissue == input$condition_input) %>% pull(library_name)
    condition_expression <- vsda[, selected_samples]
    expression_df <- condition_expression %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(cols = -gene, names_to = "library_name", values_to = "vst")
    expression_df
  })
  
  output$expression_table <- renderTable({
    expression_data() %>%
      filter(vst >= input$expression_threshold)
  })
  
  output$expression_plot <- renderPlotly({
    expression_df <- expression_data()
    expression_plt <- ggplot(expression_df, aes(x = library_name, y = vst)) +
      geom_jitter() +
      labs(x = "Sample", y = "Normalized Counts (VST)") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    expression_plt
    ggplotly(expression_plt)
  })
  
}


# run server --------------------------------------------------------------

shinyApp(ui = ui, server = server)
