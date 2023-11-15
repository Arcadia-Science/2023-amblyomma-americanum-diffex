# set up environment ------------------------------------------------------

library(shiny)
library(shinyBS)
library(DESeq2)
library(tidyverse)
library(plotly)
library(DT)

# set run options ---------------------------------------------------------

options(shiny.host = "0.0.0.0")
options(shiny.port = 6626)

# process input data ------------------------------------------------------

# read in and process annotation data
annotations <- read_tsv("input_data/Amblyomma-americanum_combined.tsv.gz")

# process annotation information so it matches the default EVM annot names
annotations <- annotations %>%
  mutate(gene = gsub("Amblyomma-americanum_", "", gene_name), .before = gene_name,
         gene = gsub("-", "_", gene),
         gene = gsub(".model", ".TU", gene)) %>%
  # define a consensus annotation. Use eggnog first, then kegg if eggnog is blank, then deepsig, or nothing.
  # also append deepsig feature to eggnog and kegg definitions if deepsig feature is not an NA
  mutate(consensus_annotation = ifelse(!is.na(egg_Description), paste0(egg_Description, "; ", ifelse(!is.na(deepsig_feature), deepsig_feature, "no deepsig annot")),
                                       ifelse(!is.na(KO_definition), paste0(KO_definition, "; ", ifelse(!is.na(deepsig_feature), deepsig_feature, "no deepsig annot")),
                                              ifelse(!is.na(deepsig_feature), deepsig_feature, "unknown function")))) %>%
  mutate(combined_gene = paste0(gene_name, "; ", consensus_annotation))

# ui logic ----------------------------------------------------------------

ui <- fluidPage(
  # custom header so model selection is reactive & at the top of the app
  tags$header(
    tags$div(
      class = "header",
      tags$h1("Differential Expression Explorer for", tags$i("Amblyomma americanum")),
      tags$p("This interactive tool allows you to explore differential expression results for two differential expression models.",
             "Use the dropdown menu to select a model, and navigate through the tabs to view different visualizations and analyses.",
             "Because this analysis is built from publicly available data, we had to get creative with our experimental design.",
             "We built two models for differential expression that compare different subsets of samples.",
             "The first model, 'sex_tissue,' combines the RNA-seq sample's sex and tissue of origin and compares these values.",
             "This model takes advantage of all 20 publicly available samples which did not have strong batch effects.",
             "The second model, 'sex_tissue_blood_meal_hour,' further incorporates the time in the blood meal in the comparison.",
             "We only had enough replicates to include 12 samples in this model."),
      selectInput("selected_model", "Select Model",
                  choices = c("sex_tissue", "sex_tissue_blood_meal_hour"))
    )
  ),
  
  # CSS to style the header
  tags$style(
    HTML(
      ".header {
        background-color: #333;
        color: white;
        padding: 20px;
        text-align: left;
      }"
    )
  ),
  
  mainPanel(
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
                   helpText(
                     "Hover over each box for guidelines of selection filtering values and interpretting results."
                   ),
                   # Filtering options (log2FC, p-value, basemean)
                   numericInput("log2FoldChange_filter", "Filter by log2FC", min = -5, max = 5, value = 1),
                   bsTooltip(id = "log2FoldChange_filter",
                             title = "A threshold of 1 or -1 (corresponding to 2-fold up or down changes) is often used as a starting point. This threshold can help filter out genes with small changes in expression that may not be biologically meaningful."),
                   numericInput("padj_filter", "Filter by adjusted p-value", min = 0, max = 1, value = 0.05),
                   bsTooltip(id = "padj_filter",
                             title = "A threshold of 0.05 is most common, although values ranging from 0.01 to 0.1 are also used. Lower values are more stringent and result in fewer false positivies but possibly more false negatives."),
                   numericInput("basemean_filter", "Filter by base mean count", min = 5, max = 10000000, value = 10),
                   bsTooltip(id = "basemean_filter", 
                             title = "A commonly-used threshold is a base mean count of 10 or higher and this value should not drop below 5 without strong justification. This threshold helps to filter out genes with very low expression, which may be more prone to statistical noise."),
                   selectInput("condition1", "Condition 1", choices =  NULL, multiple = FALSE),
                   selectInput("condition2", "Condition 2", choices =  NULL, multiple = FALSE),
                   helpText(
                     "Note: Differentially expressed genes with a positive log2FC value are more highly expressed in condition 1."
                   ),
                   helpText(
                     "Optional: upload a CSV file with your own information.",
                     "A join will be executed on the 'gene_name' column, so please ensure your file has a column named 'gene_name' with gene identifiers formatted as",
                     "'Amblyomma-americanum_evm.model.contig-XXXXX-X.X'",
                     "Once the file is uploaded, you can select a column from your file to color the points in the output plot.",
                     "If you do not upload a file, points will be colored by significance of differential expression results, inferrered from the filtering criteria specified above."
                   ),
                   fileInput("file", "Choose CSV File", 
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                   radioButtons("color_by", "Color Points By", 
                                choices = c("Significance" = "significant", "Uploaded Variable" = "uploaded"),
                                selected = "significant"),
                   uiOutput("select_column_ui"), # this will show up only if a file is uploaded
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
                   helpText(
                     "This tab plots the normalized gene counts for each sample, separated by condition.",
                     "It uses counts that have been normalized for sequencing depth and undergone variance stabilizing transformation (VST).",
                     "VST outputs transformed data on the log2 scale.",
                     "With this method, negative values indicate a count of less than 1, e.g. that the gene is not detected in the given sample."
                   ),
                   helpText(
                     "NOTE to paste genes instead of scrolling, click on the box, backspace once, and paste your gene."
                   ),
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
                   # make a tag to reverse slider highlighting so it goes from user number up to 100
                   tags$style( type = "text/css", "
#reverseSlider .irs-bar {
    border-top: 1px solid #ddd;
    border-bottom: 1px solid #ddd;
    background: linear-gradient(to bottom, #DDD -50%, #FFF 150%);
}
#reverseSlider .irs-bar-edge {
    border: 1px solid #ddd;
    background: linear-gradient(to bottom, #DDD -50%, #FFF 150%);
    border-right: 0;
}
#reverseSlider .irs-line {
    background: #428bca;
    border: 1px solid #428bca;
}
"),
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
                   selectInput("condition_input", "Select Condition:", choices = NULL, multiple = FALSE),
                   radioButtons("expression_metric", "Expression Metric:", choices = c("Minimum" = "min", "Mean" = "mean")),
                   div(id = "reverseSlider", 
                       sliderInput(inputId = "expression_percentile", 
                                  label = "Expression Percentile",
                                  min = 0, max = 100,
                                  value = 50)
                   ),
                   actionButton("view_expression", "View Expression"),
                   downloadButton('download_never', 'Download Never Expressed Table'),
                   downloadButton('download_sometimes', 'Download Sometimes Expressed Table'),
                   downloadButton('download_always', 'Download Always Expressed Table')
                 ),
                 mainPanel(
                   plotOutput("expression_plot"),
                   HTML("<h3>Always Expressed Genes</h3>"),
                   DTOutput("always_expression_table"),
                   HTML("<h3>Search for Expression of Genes</h3>"),
                   DTOutput("search_expression_table")
                 )
               ))
    )
  )
)


# server logic ------------------------------------------------------------


server <- function(input, output, session) {

  # Reactive expression to load and process the selected model
  selected_model_data <- reactive({
    # load in selected model data ----------------
    ds <- readRDS(paste0("input_data/ds_", input$selected_model, ".RDS"))
    dds <- readRDS(paste0("input_data/dds_", input$selected_model, ".RDS"))
    
    # perform variance stabilize transformation for normalization
    vsd <- vst(ds, blind = FALSE)
    # extract the VST counts as a matrix
    vsda <- assay(vsd)
    # reset colnams to sample names
    colnames(vsda) <- vsd$library_name
    # switch the gene names to the ones used by noveltree and PC
    tmp <- annotations %>% select(gene, gene_name)
    vsda <- vsda %>% 
      as.data.frame() %>% # convert from matrix to data frame
      rownames_to_column("gene") %>% # add rownames as column
      left_join(tmp, by = "gene") %>% # join to df with both annotation types
      mutate(gene_name = ifelse(is.na(gene_name), gene, gene_name)) %>% # replace NA gene_names with gene 
      select(-gene) %>% # remove gene column
      column_to_rownames("gene_name") %>% # put gene_name as new rowname
      as.matrix() # convert back to matrix
    
    
    # read in metadata and process ----------------
    metadata <- read_tsv("input_data/metadata.tsv") %>%
      filter(tissue != "cell_line") %>% # filter out samples that will make it so model is not full rank
      filter(!study_title %in% "Amblyomma americanum RNA-seq following E. coli treatment") %>% # rm e. coli study, too batched
      filter(!library_name %in% c("Um", "Uf", "Im", "If")) %>% # filter e chaf exposed samples
      select(library_name, run_accession, bioproject, experiment_title,
             study_title, sex, host_meal, tissue, blood_meal_hour, 
             blood_meal_hour_range, total_spots, publication_url, 
             publication_doi, publication_title) %>%
      group_by(library_name) %>%
      mutate(total_spots = sum(total_spots)) %>%
      mutate(blood_meal_hour_range = factor(blood_meal_hour_range, levels = c("0", "12_48", "72_144", "168_264", "none"))) %>%
      distinct()
    # format the contrast column conditionally
    if(input$selected_model == "sex_tissue"){
      # process metadata for diffex vars
      metadata <- metadata %>%
        mutate(sex_tissue = paste0(sex, "_x_", tissue))
    } else {
      metadata <- metadata %>%
        mutate(sex_tissue_blood_meal_hour = paste0(sex, "_x_", tissue, "_x_", blood_meal_hour_range))
      metadata_tally <- metadata %>%
        group_by(sex_tissue_blood_meal_hour) %>%
        tally() %>%
        filter(n > 1) # keep groups with more than one sample
      metadata <- metadata %>%
        filter(sex_tissue_blood_meal_hour %in% metadata_tally$sex_tissue_blood_meal_hour)
    }
    
    list(ds = ds, dds = dds, vsd = vsd, vsda = vsda, metadata = metadata)
  })
  
  # Update the condition choices based on the selected model
  observe({
    model_data <- selected_model_data()
    if (input$selected_model == "sex_tissue") {
      updateSelectInput(session, "condition1",
                        choices = c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"))
      updateSelectInput(session, "condition2",
                        choices = c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"))
      updateSelectInput(session, "condition_input",
                        choices = c("female_x_salivary_gland", "male_x_whole", "female_x_whole", "female_x_midgut"))
    } else {
      updateSelectInput(session, "condition1",
                        choices = c("female_x_salivary_gland_x_72_144", "male_x_whole_x_72_144", "female_x_midgut_x_72_144", "female_x_salivary_gland_x_12_48", "female_x_whole_x_72_144"))
      updateSelectInput(session, "condition2",
                        choices = c("female_x_salivary_gland_x_72_144", "male_x_whole_x_72_144", "female_x_midgut_x_72_144", "female_x_salivary_gland_x_12_48", "female_x_whole_x_72_144"))
      updateSelectInput(session, "condition_input",
                        choices = c("female_x_salivary_gland_x_72_144", "male_x_whole_x_72_144", "female_x_midgut_x_72_144", "female_x_salivary_gland_x_12_48", "female_x_whole_x_72_144"))
    }
  })
  
  # create a function to get model data
  
  # overview (PCA plot, metadata table) -------------------------------------

  # PCA plot
  output$pca_plot <- renderPlotly({
    model_data <- selected_model_data()
    metadata <- model_data$metadata
    vsd <- model_data$vsd

    pcad <- plotPCA(vsd, intgroup = c(input$selected_model), returnData = TRUE) %>%
      mutate(library_name = vsd$library_name) %>%
      left_join(metadata)
    percentVar <- round(100 * attr(pcad, "percentVar"))
    
    pca_plot <- ggplot(pcad, aes(PC1, PC2, color = !!sym(input$color_variable), 
                                 label = library_name)) +
      geom_point(size = 3) +
      labs(x = paste0("PC1: ", percentVar[1], "% variance"),
           y = paste0("PC2: ", percentVar[2], "% variance"),
           color = input$color_variable) +
      coord_fixed() +
      theme_classic()
    
    ggplotly(pca_plot) %>%
      # control position of legend to not crowd plot
      layout(legend = list(orientation = "h", y = -0.25))
  })
  
  output$metadata_table <- renderDT({
    model_data <- selected_model_data()
    model_data$metadata %>%
      datatable(options = list(pageLength = 10))
  })
  
  # differential expression results (volcano plot, table) -------------

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
    model_data <- selected_model_data()
    ds <- model_data$ds

    selected_contrast <- input$selected_model
    selected_condition1 <- input$condition1
    selected_condition2 <- input$condition2
    
    # Extract differential expression results based on selected contrasts and conditions
    ds_results <- results(ds, contrast = c(selected_contrast, selected_condition1, selected_condition2))
    ds_results_df <- ds_results %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      left_join(annotations)
    
    # if a file is uploaded and "Uploaded Variable" is selected for coloring, join it to the ds_results_df data.frame
    if (!is.null(input$file) && input$color_by == "uploaded") {
      req(input$color_column)
      uploaded_data <- read_csv(input$file$datapath)
      ds_results_df <- left_join(ds_results_df, uploaded_data(), by = "gene_name")
    }
    
    # demarcate significant
    ds_results_significant <- ds_results_df %>%
      mutate(significant = ifelse(abs(log2FoldChange) >= input$log2FoldChange_filter & 
                                  padj <= input$padj_filter & 
                                  baseMean >= input$basemean_filter, "significant", "not significant"))
    
    # Return the filtered data frame
    ds_results_significant
  })
  
  # volcano plot
  output$volcano_plot <- renderPlotly({
    diff_results_df <- diff_results()

    if(input$color_by == "significant") {
       color_var <- diff_results_df$significant
       color_name <- "Meets thresholds"
    } else {
       color_var <- diff_results_df[[input$color_column]]
       color_name <- input$color_column
    }
  
    volcano_plot <- ggplot(diff_results_df, aes(x = log2FoldChange, 
                                                y = -log10(padj), 
                                                color = color_var,
                                                text = str_wrap(combined_gene, 50))) +
      geom_point(alpha = 0.5) +
      labs(x = "log2(Fold Change)", 
           y = "-log10(Adjusted p Value (BH))", 
           color = color_name,
           label = "gene") +
      geom_hline(yintercept = -log10(input$padj_filter), linetype = 2) +
      geom_vline(xintercept = input$log2FoldChange_filter, linetype = 2) +
      geom_vline(xintercept = -input$log2FoldChange_filter, linetype = 2) +
      theme_classic()
    
    ggplotly(volcano_plot, tooltip = c("text", "x", "y"))
  })
  
  # MA plot
  output$ma_plot <- renderPlotly({
    diff_results_df <- diff_results()
    
    if(input$color_by == "significant") {
      color_var <- diff_results_df$significant
      color_name <- "Meets thresholds"
    } else {
      color_var <- diff_results_df[[input$color_column]]
      color_name <- input$color_column
    }
    
    ma_plot <- ggplot(diff_results_df, aes(x = log2(baseMean), 
                                           y = log2FoldChange, 
                                           color = color_var, 
                                           text = str_wrap(combined_gene, 50))) +
      geom_point(alpha = 0.5) +
      labs(x = "log2(Mean Count)", y = "log2(Fold Change)", color = color_name) +
      geom_hline(yintercept = input$log2FoldChange_filter, linetype = 2) +
      geom_hline(yintercept = -input$log2FoldChange_filter, linetype = 2) +
      geom_vline(xintercept = log2(input$basemean_filter), linetype = 2) +
      theme_classic()
    
    ggplotly(ma_plot,  tooltip = c("text", "x", "y"))
  })
  
  # significant gene table
  output$gene_table <- renderDT({
    diff_results() %>%
      filter(significant == "significant") %>%
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
    model_data <- selected_model_data()
    vsda <- model_data$vsda
    gene_choices <- rownames(vsda)
    updateSelectizeInput(session, "selected_gene", choices = gene_choices, server = TRUE)
  })

  # plot gene count boxplot
  output$gene_boxplot <- renderPlotly({
    model_data <- selected_model_data()
    metadata <- model_data$metadata
    vsda <- model_data$vsda
    selected_gene <- input$selected_gene
    
    if (!is.null(selected_gene)) {
      gene_data <- data.frame(gene_count = vsda[selected_gene, ]) %>%
        rownames_to_column("library_name") %>%
        left_join(metadata, by = "library_name")

      gene_boxplot <- ggplot(gene_data, aes(x = !!sym(input$selected_model), y = gene_count)) +
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
    model_data <- selected_model_data()
    metadata <- model_data$metadata
    vsda <- model_data$vsda
    
    req(input$view_expression)  # ensure the button is clicked
    req(input$condition_input)  # ensure a contrast is selected
    # Filter the expression data for the selected contrast
    # Assuming vsda holds the expression data where rows are genes and columns are samples
    selected_samples <- metadata %>% filter(!!sym(input$selected_model) == input$condition_input) %>% pull(library_name)
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
      mutate(sometimes_expression = ifelse(max_vst >= 0 & min_vst < 0, "sometimes", "other")) %>%
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
                                   fill = ifelse(min_percentile >= input$expression_percentile, TRUE, FALSE)))
    } else {
      expression_plt <- ggplot(expression_list$always_expressed_df, 
                               aes(x = vst, 
                                   fill = ifelse(mean_percentile >= input$expression_percentile, TRUE, FALSE))) }
    expression_plt + 
      geom_histogram(binwidth = 0.1) +
      labs(title = "Distribution of Normalized Count Values",
           x = "Normalized Counts (VST)",
           y = "Frequency",
           fill = "Greater than or equal\nto input percentile") +
      theme_classic()
  })
  
  output$always_expression_table <- renderDT({
    expression_data()$always_expressed_filtered %>%
      select(-expression_metric) %>% # rm internal col used for filtering
      datatable(options = list(pageLength = 10))
  })
  
  output$search_expression_table <- renderDT({
    always <- expression_data()$always_expressed_df %>%
      select(gene) %>%
      mutate(expression_category = "always")
    sometimes <- expression_data()$sometimes_expressed_df %>%
      select(gene) %>%
      mutate(expression_category = "sometimes") 
    never <- expression_data()$not_expressed_df %>%
      select(gene) %>%
      mutate(expression_category = "never")
    
    bind_rows(always, sometimes, never) %>%
      distinct() %>%
      datatable(options = list(pageLength = 10))
  })

  # include download handlers for tables
  output$download_always <- downloadHandler(
    filename = function() {
      paste0('always_expressed_', input$condition_input, '.csv')
    },
    content = function(file) {
      write_csv(expression_data()$expressed_df, file)
    })
  
  output$download_never <- downloadHandler(
    filename = function() {
      paste0('never_expressed_', input$condition_input, '.csv')
    },
    content = function(file) {
      write_csv(expression_data()$not_expressed_df, file)
    })
  
  output$download_sometimes <- downloadHandler(
    filename = function() {
      paste0('sometimes_expressed_', input$condition_input, '.csv')
    },
    content = function(file) {
      write_csv(expression_data()$sometimes_expressed_df, file)
    })
}


# run server --------------------------------------------------------------

shinyApp(ui = ui, server = server)
