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
                   "This tab reports on expression per condition.",
                   "It uses counts that have been normalized for sequencing depth and undergone variance stabilizing transformation (VST).",
                   "VST outputs transformed data on the log2 scale.",
                   "With this method, negative values indicate a count of less than 1, e.g. that the gene is not detected in the given sample.",
                   "Using this information, we calculate genes that are never expressed in samples in that condition,",
                   "that are sometimes expressed in samples in that conditions, or that are always expressed.",
                   "For genes that are always expressed (e.g. expressed in all samples), you can further by the percentile of the genes expression.",
                   "You can use the minimum count/percentile (the lowest count observed in all samples in the group) or the mean count/percentile."
                 ),
                 selectInput("condition_input", "Select Condition:", choices = unique(metadata$sex_tissue)),
                 radioButtons("expression_metric", "Expression Metric:", choices = c("Minimum" = "min", "Mean" = "mean")),
                 sliderInput("expression_percentile", "Expression Percentile:", min = 0, max = 100, value = 50),
                 numericInput("expression_threshold", "Expression Threshold:", value = 0, min = -6),
                 actionButton("view_expression", "View Expression"),
                 downloadButton('download_never', 'Download Never Expressed Table'),
                 downloadButton('download_sometimes', 'Download Sometimes Expressed Table'),
                 downloadButton('download_always', 'Download Always Expressed Table')
               ),
               mainPanel(
                 plotOutput("expression_plot"),
                 HTML("<h3>Never Expressed Genes</h3>"),
                 DTOutput("never_expression_table"),
                 HTML("<h3>Sometimes Expressed Genes</h3>"),
                 DTOutput("sometimes_expression_table"),
                 HTML("<h3>Always Expressed Genes</h3>"),
                 DTOutput("always_expression_table")
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
    # format per-sample expression data as a long data frame
    expression_df <- condition_expression %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(cols = -gene, names_to = "library_name", values_to = "vst")
    
    # create a df of genes that are never expressed in the input condition
    not_expressed_df <- expression_df %>%
      group_by(gene) %>%
      summarize(max_vst = max(vst)) %>%
      filter(max_vst < 0)
    
    # create a df of genes that are sometimes expressed in the input condition
    sometimes_expressed_df <- expression_df %>%
      group_by(gene) %>%
      summarize(max_vst = max(vst),
                min_vst = min(vst)) %>%
      mutate(sometimes_expression = ifelse(max_vst > 0 & min_vst < 0, "sometimes", "other")) %>%
      filter(sometimes_expression == "sometimes") %>%
      select(-sometimes_expression)
    
    # create a df of genes that are always expressed and define magnitude of expression via percentiles
    expressed_df <- expression_df %>%
      filter(!gene %in% not_expressed_df$gene) %>%
      filter(!gene %in% sometimes_expressed_df$gene)
    
    # calculate percentile of expression for genes that are always expressed
    expressed_df_percentiles <- quantile(x = expressed_df$vst, probs = seq(0, 1, by = 0.01))
    
    # define a function for finding the percentile in a group
    find_percentile <- function(value, percentiles) {
      sum(value >= percentiles)
    }
    
    # calculate the minimum and mean expression for each gene across samples
    expressed_df_summary <- expressed_df %>%
      group_by(gene) %>%
      summarize(mean_vst = mean(vst),
                min_vst = min(vst)) %>%
      rowwise() %>%
      mutate(min_percentile = find_percentile(min_vst, expressed_df_percentiles) - 1,
             mean_percentile = find_percentile(mean_vst, expressed_df_percentiles) - 1)
    
    # join the summarized count info with the long-format df
    always_expressed_df <- left_join(expressed_df, expressed_df_summary)
    
    # determine the expression metric and percentile to use for filtering
    expressed_df_summary$expression_metric <- if(input$expression_metric == "min"){
      expressed_df_summary$min_percentile
    } else {
      expressed_df_summary$mean_percentile
    }
    
    expression_percentile <- input$expression_percentile
    
    # filter the always expressed genes by the chosen percentile
    always_expressed_filtered <- expressed_df_summary %>%
      filter(expression_metric >= expression_percentile)
    
    return(list(always_expressed_df = always_expressed_df,
                always_expressed_filtered = always_expressed_filtered,
                not_expressed_df = not_expressed_df,
                sometimes_expressed_df = sometimes_expressed_df))
  })

  # plot the distribution of expression values for genes that are always expressed
  output$expression_plot <- renderPlot({
    expression_list <- expression_data()
    if(input$expression_metric == "min"){
      expression_plt <- ggplot(expression_list$always_expressed_df, 
                               aes(x = vst, 
                                   fill = ifelse(min_percentile > 50, TRUE, FALSE)))
    } else {
      expression_plt <- ggplot(expression_list$always_expressed_df, 
                               aes(x = vst, 
                                   fill = ifelse(mean_percentile > 50, TRUE, FALSE))) }
    expression_plt + 
      geom_histogram(binwidth = 0.1) +
      labs(title = "Distribution of Normalized Count Values",
           x = "Normalized Counts (VST)",
           y = "Frequency",
           fill = "Greater than input percentile") +
      theme_classic()
  })
  
  output$never_expression_table <- renderDT({
    expression_data()$not_expressed_df%>%
      datatable(options = list(pageLength = 10))
  })
  
  output$sometimes_expression_table <- renderDT({
    expression_data()$sometimes_expressed_df %>%
      datatable(options = list(pageLength = 10))
  })

  output$always_expression_table <- renderDT({
    expression_data()$always_expressed_filtered %>%
      select(-expression_metric) %>% # rm internal col used for filtering
      datatable(options = list(pageLength = 10))
  })
  
  # include download handlers for tables
  output$download_never <- downloadHandler(
    filename = function() {
      paste0('never_expressed_', input$condition_input, '.csv')
    },
    content = function(file) {
      write_csv(expression_data()$not_expressed_df, file)
    }
  )
  
  output$download_sometimes <- downloadHandler(
    filename = function() {
      paste0('sometimes_expressed_', input$condition_input, '.csv')
    },
    content = function(file) {
      write_csv(expression_data()$sometimes_expressed_df, file)
    }
  )
  
  output$download_always <- downloadHandler(
    filename = function() {
      paste0('always_expressed_', input$condition_input, '.csv')
    },
    content = function(file) {
      write_csv(expression_data()$expressed_df, file)
    }
  )
}


# run server --------------------------------------------------------------

shinyApp(ui = ui, server = server)
