# Load the packages
library(pheatmap)
library(shiny)
library(RColorBrewer)
library(data.table)
library(msigdbr)
library(heatmaply)
library(dplyr)
library(shinyWidgets)
library(rintrojs)
library(readxl)

# Define UI for application that draws a histogram

heatmap_ui <- function(id, label, choices) {
  ns <- NS(id)
  tagList(
    
    useShinyalert(),
    
    add_busy_spinner(
      spin = "cube-grid",
      position = "top-right",
      color = "#01303f",
      margins = c(300, 500),
      height = "60px",
      width = "60px"),
    
    sidebarPanel(
      selectizeInput(inputId =ns("hm_definition_sel"), 
                     multiple=T,
                     label = "Select sample types",
                     choices=NULL, # will be updated dynamically
                     options=list(placeholder = "eg. Primary solid tumor")),
      selectizeInput(inputId = ns("selectors"),
                     label= "Choose how you would like to select genes to be shown",
                     choices = c("Manually enter gene names", "MSigDB gene sets", "Upload a xlsx/xls file"),
                     multiple = FALSE),
      conditionalPanel(
        condition = "input.selectors == 'Manually enter gene names'",
        ns=ns,
        selectizeInput(inputId = ns("genes"),
                       label= "Type gene names",
                       choices = NULL,
                       multiple = TRUE,
                       options=list(placeholder = "eg. TSPAN6, TNMD etc."))
      ),
      conditionalPanel(
        condition= "input.selectors == 'MSigDB gene sets'",
        ns=ns,
        selectizeInput(ns("cat"), "Please select an MSigDB Collection", choices = c("Hallmark gene sets (H)" = "H",
                                                                                    "Positional gene sets (C1)" = "C1",
                                                                                    "Curated gene sets (C2)" = "C2",
                                                                                    "Regulatory target gene sets (C3)" = "C3",
                                                                                    "Computational gene sets (C4)" = "C4",
                                                                                    "Ontology gene sets (C5)" = "C5" ,
                                                                                    "Oncogenic gene sets (C6)" = "C6",
                                                                                    "Immunologic gene sets (C7)" = "C7",
                                                                                    "Cell type signature gene sets (C8)" = "C8")),
        conditionalPanel(condition = "input.cat == 'C2'|input.cat =='C3'|
                                      input.cat =='C4'|
                                      input.cat =='C5'|
                                      input.cat =='C7' ", ns = ns, 
                         selectizeInput(ns("heatmap_subcat"),"Please select a subcategory" ,choices = c(""))
        ),
        selectizeInput(ns("chosen_gse"), #Automatically updates to show subset.
                       "Select a specific gene set", 
                       choices = NULL)
      ),
      conditionalPanel(
        condition = "input.selectors == 'Upload a xlsx/xls file'",
        ns=ns,
        fileInput(ns("heatmap_csv"), 
                  label = tags$span(
                    "Please upload your xlsx/xls file.", 
                    tags$i(
                      class = "glyphicon glyphicon-info-sign", 
                      style = "color:#0072B2;",
                      title = "The xlsx/xls file should contain one unnamed column: the first column should contain the gene names."
                    ), tags$br(),
                    a(href="sample_gene_input.xlsx", "Sample Input File", download=NA, target="_blank")),
                  accept = c(".xls", ".xlsx"),
                  multiple = FALSE)
      ),
      radioButtons(inputId = ns("scale_by_row_column"),
                   label = "Scale by row or column",
                   c("Row" = "row", 
                     "Column" = "column",
                     "None" = "none"),
                   selected = "none"),
      sliderInput(inputId = ns("variable"),
                  label = "Variation filter (Keep top n% variable genes)",
                  min = 0, max = 100, value = 100
      ),
      
      span(style="color:#3382FF",
           
           selectizeInput(ns("annotation"), 
                          "Please select annotation variables (optional)", 
                          choices = NULL,
                          multiple = TRUE,
                          options=list(placeholder = "eg. meta.definition, meta.LYMPHOCYTE.SCORE etc.")),
           selectizeInput(ns("hm_categorized_gene"), 
                          label = tags$span(
                            "Please select a gene/genes to be categorized as high/low (optional)",
                            tags$i(
                              class = "glyphicon glyphicon-info-sign",
                              style = "color:#0072B2;",
                              title = "Also meta variables can be added into this input box, but averaging meta variables and genes might produce meaningless results, try to examine them separately."
                            )),
                          choices = NULL,
                          multiple = TRUE,
                          options=list(placeholder = "eg. TSPAN6, TNMD etc."))
           
      ),
      radioButtons(inputId = ns("hm_gene_categorization_button"),
                   label = "Choose if you would like to categorize genes as high/low separately or by taking the median.",
                   c("Categorize genes after averaging" = "take_median",
                     "Categorize genes separately" = "take_separately"),
                   selected = "take_separately"),
      
      selectizeInput(inputId = ns("clustering_distance_rows"), 
                     multiple=F,
                     label = "Choose the distance calculation method for genes",
                     choices=c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")),
      
      selectizeInput(inputId = ns("clustering_distance_columns"), 
                     multiple=F,
                     label = "Choose the distance calculation method for samples",
                     choices=c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")),
      
      
      selectizeInput(inputId = ns("clustering_method_rows"), 
                     multiple=F,
                     label = "Choose the clustering method for genes",
                     choices=c("ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D")),
      
      
      selectizeInput(inputId = ns("clustering_method_columns"), 
                     multiple=F,
                     label = "Choose the clustering method for samples",
                     choices=c("ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D")),
      
      actionBttn(inputId = ns("heatmap_run"), 
                 label = "Analyze",
                 style = "unite",
                 block = TRUE,
                 color = "primary"),
      br(),
      
      downloadButton(ns("downloadPlot6"), "Download Heatmap Plot", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
      
      #help section UI
      
      introjsUI(),
      actionButton(ns("heatmap_help"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
      
      textOutput(ns("filewarning_6")) ,
      
      width = 3
      
    ),
    
    mainPanel(
      plotlyOutput(outputId = ns("heatmap_plot"), width = "950px", height = "900px")
    )
  )
  
}

# Define server logic required to draw a histogram
heatmap_server <- function(id,Xproj) {
  moduleServer(
    id,
    function(input, output, session) {
      
      output$filewarning_6 <- renderText({
        
        if (!is.null(Xproj$fileInfost())) {
          shinyalert("Warning!", "To perform this analysis using MsigDB gene sets, please ensure that your uploaded data set contains gene symbols rather than Entrez or Ensembl gene IDs. Otherwise you may receive errors.") }
      })
      
      ns <- session$ns
      
      # Heatmap App Tutorial
      
      heatmap_help_tutorial <- reactive({
        
        if(input$selectors == "Manually enter gene names"){
          
          return(
            data.frame(
              
              element = paste0("#", session$ns(c(NA, "hm_definition_sel + .selectize-control", "selectors + .selectize-control ", "genes + .selectize-control", "variable + .selectize-control", "annotation + .selectize-control ", "hm_categorized_gene + .selectize-control",
                                                 "clustering_distance_rows + .selectize-control", "clustering_distance_columns + .selectize-control", "clustering_method_rows + .selectize-control", "clustering_method_columns + .selectize-control"))),
              
              intro = paste(c(
                "This is the heatmap module where you can visualize expression patterns of selected genes. Continue the tutorial to learn how to use this module",
                "You can select sample types to focus the analysis on the specific subsets",
                "Here you can choose how you would like to select genes for the plot",
                "You can manually type genes of interest here (for other input possibilities (and associated tutorial steps), change the selection in the previous box)",
                "Next, you can apply a variance filter to keep only highly variable genes in the plot. 100 (default) means no filtering is applied. If you like to see top 10% variable genes only, set this value to 10. Such filtering can help see more informative genes.",
                "You can select categorical clinical meta data features to show as annotations on top of the heatmap.",
                "You can also create an annotation bar by categorizing the patients based on their gene expression levels. You can specify one or more genes here. When multiple genes are entered, you can categorize features separately and show them in individual annotation bars; or you can categorize after taking the overall average. Categorization is done as 'high' and 'low' at the median gene expression value",
                "You can choose how the distance will be calculated for genes here. You can read <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/dist.html'>R documentation</a> for further details.",
                "You can choose how the distance will be calculated for samples here",
                "You can choose different hierarchical clustering methods for genes here. You can read <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html'>R documentation</a> for further details.",
                "You can choose different hierarchical clustering methods for samples here"
              ))
            )
          )
        } else if (input$selectors == "MSigDB gene sets") {
          return(
            
            data.frame(
              
              element = paste0("#", session$ns(c(NA, "hm_definition_sel + .selectize-control", "selectors + .selectize-control ", "cat + .selectize-control", "chosen_gse + .selectize-control", "variable + .selectize-control", "annotation + .selectize-control ", "hm_categorized_gene + .selectize-control", 
                                                 "clustering_distance_rows + .selectize-control", "clustering_distance_columns + .selectize-control", "clustering_method_rows + .selectize-control", "clustering_method_columns + .selectize-control"))),
              
              intro = paste(c(
                "This is the heatmap module where you can visualize expression patterns of selected genes. Continue the tutorial to learn how to use this module",
                "You can select sample types to focus the analysis on the specific subsets.",
                "Here you can choose how you would like to select genes for the plot",
                "If MSigDB is selected, you can specify the main MSigDB collection you are interested in here. The next selection box will allow you to select a particular gene set",
                "When plotting gene sets from MSigDB, you can specify the particular gene set here.",
                "Next, you can apply a variance filter to keep only highly variable genes in the plot. 100 (default) means no filtering is applied. If you like to see top 10% variable genes only, set this value to 10. Such filtering can help see more informative genes.",
                "You can select categorical clinical meta data features to show as annotations on top of the heatmap.",
                "You can also create an annotation bar by categorizing the patients based on their gene expression levels. You can specify one or more genes here. When multiple genes are entered, you can categorize features separately and show them in individual annotation bars; or you can categorize after taking the overall average. Categorization is done as 'high' and 'low' at the median gene expression value",
                "You can choose how the distance will be calculated for genes here. You can read <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/dist.html'>R documentation</a> for further details.",
                "You can choose how the distance will be calculated for samples here",
                "You can choose different hierarchical clustering methods for genes here. You can read <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html'>R documentation</a> for further details.",
                "You can choose different hierarchical clustering methods for samples here"
              ))
            )
            
          )
        } else if (input$selectors == "Upload a xlsx/xls file"){
          return(
            
            data.frame(
              
              element = paste0("#", session$ns(c(NA, "hm_definition_sel + .selectize-control", "selectors + .selectize-control ", "heatmap_csv + .selectize-control", "variable + .selectize-control", "annotation + .selectize-control ", "hm_categorized_gene + .selectize-control", 
                                                 "clustering_distance_rows + .selectize-control", "clustering_distance_columns + .selectize-control", "clustering_method_rows + .selectize-control", "clustering_method_columns + .selectize-control"))),
              
              intro = paste(c(
                "This is the heatmap module where you can visualize expression patterns of selected genes. Continue the tutorial to learn how to use this module",
                "You can select sample types to focus the analysis on the specific subsets.",
                "Here you can choose how you would like to select genes for the plot",
                "You can upload a xlsx/xls file including your genes of interest to see them on the heatmap",
                "Next, you can apply a variance filter to keep only highly variable genes in the plot. 100 (default) means no filter is applied. If you like to see top 10% variable genes only, set this value to 10. Such filtering can help see more informative genes.",
                "You can select categorical clinical meta data features to show as annotations on top of the heatmap.",
                "You can also create an annotation bar by categorizing the patients based on their gene expression levels. You can specify one or more genes here. When multiple genes are entered, you can categorize features separately and show them in individual annotation bars; or you can categorize after taking the overall average. Categorization is done as 'high' and 'low' at the median gene expression value",
                 "You can choose how the distance will be calculated for genes here. You can read <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/dist.html'>R documentation</a> for further details.",
                "You can choose how the distance will be calculated for samples here",
                "You can choose different hierarchical clustering methods for genes here. You can read <a href='https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html'>R documentation</a> for further details.",
                "You can choose different hierarchical clustering methods for samples here"
              ))
            )
            
          )
        }
        
      })
      
      
      observeEvent(input$heatmap_help, {
        
        introjs(session, options = list(steps = heatmap_help_tutorial() ) )
        
      })
      
      # Reactive selectboxes
      
      observe({
        req(input$cat)
        
        heatmap_subcat = names(heatmap_df_msigdb()[[input$cat]])
        if (length(heatmap_subcat) > 1)  {
          updateSelectizeInput(session,'heatmap_subcat', choices = heatmap_subcat , server = TRUE)
        } 
      })
      
      
      # Remove cols which is including more than 10 levels and numeric
      
      hidden_cols_heatmap_plotvar<-reactive({"meta.definition"})
      
      meta_cols_heatmap_plotvar <- reactive({colnames(Xproj$a())[grep("^meta\\.", colnames(Xproj$a()))]})
      
      sel_cols_heatmap_plotvar <- reactive({meta_cols_heatmap_plotvar()[unlist(lapply(Xproj$a()[, meta_cols_heatmap_plotvar(), with = FALSE], function(x) length(levels(x)))) < 10]})
      
      numeric_cols_heatmap_plotvar <- reactive({unlist(lapply(Xproj$a()[, sel_cols_heatmap_plotvar(), with = FALSE], function(x) is.numeric(x)))})
      
      available_cols_heatmap_plotvar <- reactive({setdiff(sel_cols_heatmap_plotvar(), c(hidden_cols_heatmap_plotvar(), sel_cols_heatmap_plotvar()[numeric_cols_heatmap_plotvar()]))})
      
      
      # observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "selectors", choices = c("Manually enter gene names", "MSigDB gene sets", "Upload a xlsx/xls file"), server = TRUE)})
      
      observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "genes", choices = colnames(Xproj$a()[, lapply(Xproj$a(), is.numeric) == TRUE, with = FALSE]), server = TRUE)})
      observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "annotation", choices = available_cols_heatmap_plotvar(), selected = character(0), server = TRUE)})
      observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "hm_categorized_gene", choices = colnames(Xproj$a()[, lapply(Xproj$a(), is.numeric) == TRUE, with = FALSE]), selected = character(0), server = TRUE)})
      observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "hm_definition_sel", choices = c(Xproj$a()$meta.definition), selected = character(0), server = TRUE)})
      
      heatmap_df_msigdb <- reactive({heatmap_df_msigdb <- readRDS(paste0("genesets/", "msigdb_long_w_immth", ".rds"))})
      
      hm_gene_sets <- reactive({ 
        
        req(heatmap_df_msigdb())
        
        #Preparation of Human MSigDB genesets.
        
        df_msigdb <- heatmap_df_msigdb()
        
        if(input$selectors == 'MSigDB gene sets'){
          
          if(input$cat %in% c("C2","C3","C4","C5","C7")) {
            
            df_msigdb2 = names(df_msigdb[[input$cat]][[input$heatmap_subcat]])
            
          } else {
            df_msigdb2 = names(df_msigdb[[input$cat]][[1]])
          }
          # return(df_msigdb2)
          
        } else {
          
          return({
            df_msigdb3 <- NULL
          })
        }
        
      })
      
      
      observe({
        
        req(heatmap_df_msigdb())
        
        req(hm_gene_sets())
        
        if(input$cat %in% c("C2","C3","C4","C5","C7")) {
          
          
          heatmap_path = hm_gene_sets()
          updateSelectizeInput(session,'chosen_gse', choices = heatmap_path , server = TRUE)
          
        } else {
          
          heatmap_path = hm_gene_sets()
          updateSelectizeInput(session,'chosen_gse', choices = heatmap_path , server = TRUE)
        }
        
      })
      
      pre_data <- eventReactive(input$heatmap_run, {
        
        # Preparing the preliminary data.
        
        if(length(as.vector(input$hm_definition_sel)) > 0) {
          return({
            daf1 <- Xproj$a()
            daf1 <- setDT(daf1, key = 'meta.definition')[J(input$hm_definition_sel)]
            daf1
          })
        } else if (length(as.vector(input$hm_definition_sel)) == 0){
          return({
            daf2 <- Xproj$a()
            daf2
          })
        }
        
      })
      
      mat <- eventReactive(input$heatmap_run, {
        
        # browser()
        
        # Matching the MSigDB gene sets with selected TCGA data.
        
        req(pre_data())
        
        if(input$selectors == "MSigDB gene sets"){
          
          req(hm_gene_sets())
          
          req(input$chosen_gse)
          
          heatmap_df_msigdb <- readRDS(paste0("genesets/", "msigdb_long_w_immth", ".rds"))
          
          if(input$cat %in% c("C2","C3","C4","C5","C7")) {
            
            heatmap_msigdb_genes <- heatmap_df_msigdb[[input$cat]][[input$heatmap_subcat]][[input$chosen_gse]]
            
          } else {
            
            heatmap_msigdb_genes <- heatmap_df_msigdb[[input$cat]][[1]][[input$chosen_gse]]
          }
          
          #Take the subset of chosen Human MsigDB geneset.
          
          final_gene_sets <- heatmap_msigdb_genes
          
          # Drop and the Human MsigDB genes which are not in our data.
          
          daf <- as.data.frame(pre_data())
          
          # rownames(daf) <- daf$meta.barcode
          
          same_hallmarks_names = intersect(final_gene_sets, colnames(daf))
          
          selected_cols <- c(same_hallmarks_names)
          
          daf <- daf[,selected_cols]  
          
        } else if(input$selectors == "Manually enter gene names"){
          
          req(input$genes)
          
          daf <- as.data.frame(pre_data())
          
          rownames(daf)<-NULL
          
          rownames(daf) <- daf$meta.barcode
          
          selected_cols <- c(input$genes)
          
          daf <- daf[,selected_cols]  
          
        } else if(input$selectors == "Upload a xlsx/xls file"){
          
          
          daf <- as.data.frame(pre_data())
          
          # rownames(daf) <- daf$meta.barcode
          
          uploaded_heatmap_csv <- input$heatmap_csv
          
          selected_csv <- as.data.frame(read_excel(uploaded_heatmap_csv$datapath, sheet = 1, col_names = F))
          
          colnames(selected_csv)[1] <- "Genes" 
          
          same_gene_names = intersect(selected_csv[["Genes"]], colnames(daf))
          
          selected_cols <- c(same_gene_names)
          
          daf <- daf[,selected_cols]  
          
        }
        
        
        if (input$variable !=100) {
          
          # daf_variable <- t(daf)
          # 
          # evar <- apply(daf_variable,1,var)
          # 
          # mostVariable <- daf_variable[evar>quantile(evar,as.numeric(input$variable)/100),]
          # mostvariablegenes <- rownames(mostVariable)
          # daf <- daf[,mostvariablegenes]
          
          daf_variation <- apply(daf,2,var, na.rm=T) 
          
          mostVariablegenes <- names(which(daf_variation >= quantile(daf_variation, 1-as.numeric(input$variable)/100)))
          
          daf <- daf[,mostVariablegenes]
          
          
          
        } else if(input$variable == 100){
          daf
        }
        
        # When you have several variables to examine over multiple scales, the scale() function makes more sense. 
        # One variable, for example, is of magnitude 100, whereas another is of magnitude 1000.
        # if(input$scale_by_row == TRUE){
        #   
        # daf <- scale(daf, scale = FALSE)
        ## scaling will be done in heatmaply from now on, scale function is not needed anymore
        
        daf <- as.matrix(daf)
        
        daf <- t(daf)
        
        # } else if (input$scale_by_row == FALSE){
        #   
        #   daf <- as.matrix(daf)
        # 
        # daf <- t(daf)
        # }
      })
      
      meta <- eventReactive(input$heatmap_run, {

        # Preparing the meta data for annotation bar.
        
        req(pre_data())
        
        if(length(as.vector(input$annotation)) > 0 & length(as.vector(input$hm_categorized_gene)) > 0){
          
          return({
            
            # When there are both categorized genes and meta annotation
            
            daf <- as.data.frame(pre_data())
            
            rownames(daf) <- daf$meta.barcode
            
            hm_categorized_gene <- c(input$hm_categorized_gene)
            
            chosen_meta <- c(input$annotation, input$hm_categorized_gene)
            
            meta <- daf[,chosen_meta]
            
            # High-low expression categorization of selected genes
            
            meta <- as.data.table(meta)
            
            if (input$hm_gene_categorization_button == "take_median" & length(as.vector(input$hm_categorized_gene)) > 1){

            
              meta <- meta[, hm_categorized_gene_means := rowMeans(.SD, na.rm = TRUE), .SDcols = c(hm_categorized_gene)]  

              cat_gene_median <-  median(meta$hm_categorized_gene_means, na.rm = TRUE)
              
              meta <- meta[,c(hm_categorized_gene):=NULL]
              
              meta <- as.data.frame(meta)
              
              meta <- meta %>% 
                mutate(categorized_gene_set_expression = case_when(
                  meta$hm_categorized_gene_means > cat_gene_median ~ "high",
                  meta$hm_categorized_gene_means < cat_gene_median ~ "low"
                )) %>% 
                select(-hm_categorized_gene_means)
              
            } else if (input$hm_gene_categorization_button == "take_separately" | length(as.vector(input$hm_categorized_gene)) == 1){
              
              dat <- meta
              
              for (i in input$hm_categorized_gene) {
                
                meta[, (i) := ifelse(dat[[i]] >= median(dat[[i]], na.rm = T), "high", ifelse(dat[[i]] < median(dat[[i]], na.rm = T), "low", "mid"))]
                
              }
              
              return(meta <- as.data.frame(meta))
              
            }
            
          })} else if (length(as.vector(input$annotation)) > 0 & length(as.vector(input$hm_categorized_gene)) == 0){
            
            return({
              
              # When there is no categorized genes entered
              
              daf <- as.data.frame(pre_data())
              
              rownames(daf) <- daf$meta.barcode
              
              chosen_meta <- c(input$annotation)
              
              meta <- daf[,chosen_meta]
              
            })
            
          } else if (length(as.vector(input$annotation)) == 0 & length(as.vector(input$hm_categorized_gene)) > 0){
            
            return({
              
              # When there are categorized genes but no meta annotation
              
              req(input$hm_categorized_gene)
              
              daf <- as.data.frame(pre_data())
              
              rownames(daf) <- daf$meta.barcode
              
              hm_categorized_gene <- c(input$hm_categorized_gene)
              
              meta <- daf[,input$hm_categorized_gene]
              
              # High-low expression categorization of selected genes
              
              meta <- as.data.table(meta)
              
              if (input$hm_gene_categorization_button == "take_median"){
                
                if (length(as.vector(input$hm_categorized_gene)) == 1){
                  
                  setnames(meta, "meta", "hm_categorized_gene_means")
                  
                } else {
                  
                  meta <- meta[, hm_categorized_gene_means := rowMeans(.SD, na.rm = TRUE), .SDcols = c(hm_categorized_gene)]  
                  
                }
                
                
                cat_gene_median <-  median(meta$hm_categorized_gene_means, na.rm = TRUE)
                
                meta <- meta[,c(hm_categorized_gene):=NULL]
                
                meta <- as.data.frame(meta)
                
                meta <- meta %>% 
                  mutate(categorized_gene_set_expression = case_when(
                    meta$hm_categorized_gene_means > cat_gene_median ~ "high",
                    meta$hm_categorized_gene_means < cat_gene_median ~ "low"
                  )) %>% 
                  select(-hm_categorized_gene_means)
                
              } else if (input$hm_gene_categorization_button == "take_separately"){
                
                if (length(as.vector(input$hm_categorized_gene)) == 1){
                  
                  setnames(meta, "meta", "hm_categorized_gene_means")
                  
                  cat_gene_median <-  median(meta$hm_categorized_gene_means, na.rm = TRUE)
                  
                  meta <- meta[,c(hm_categorized_gene):=NULL]
                  
                  meta <- as.data.frame(meta)
                  
                  meta <- meta %>% 
                    mutate(categorized_gene_set_expression = case_when(
                      meta$hm_categorized_gene_means > cat_gene_median ~ "high",
                      meta$hm_categorized_gene_means < cat_gene_median ~ "low"
                    )) %>% 
                    select(-hm_categorized_gene_means)
                  
                } else{
                  
                  dat <- meta
                  
                  for (i in input$hm_categorized_gene) {
                    
                    meta[, (i) := ifelse(dat[[i]] >= median(dat[[i]], na.rm = T), "high", ifelse(dat[[i]] < median(dat[[i]], na.rm = T), "low", "mid"))]
                    
                  }
                  
                  return(meta <- as.data.frame(meta))
                  
                } 
              }
              
              
            })} else if (length(as.vector(input$annotation)) == 0 & length(as.vector(input$hm_categorized_gene)) == 0){
              
              return({
                
                # When there is no categorized genes or meta annotation entered
                meta <- NULL
                
              })
            }
      })
      
      heatmap_object <- eventReactive(input$heatmap_run, {
        
        # Creating the plot.
        
        if(is.null(meta()) == TRUE){
          
          # If no categorical annotations are selected.
          
          validate(need(input$hm_definition_sel, "Please select a sample type"))
          validate(need(input$selectors, "Please select a Human MSigDB Gene set or create a gene set"))
          if(input$selectors == 'Manually enter gene names'){
            validate(need(input$genes, "Create genes to create your gene set"))
          } else if (input$selectors == 'MSigDB gene sets'){
            validate(need(input$cat, "Please select a Human MSigDB Collection"))
            validate(need(input$chosen_gse, "Please select the subset of your chosen Human MSigDB Collection"))
          }
          validate(need(input$heatmap_run, "Please click the 'Analyze' button"))
          
          
          return({
            if(input$clustering_distance_rows == "pearson"| input$clustering_distance_rows == "spearman" | input$clustering_distance_rows == "kendall"){
              distfun_row = function(x) as.dist(1 - cor(t(x), method=input$clustering_distance_rows))
            } else{
              distfun_row = function(x) stats::dist(x, method = input$clustering_distance_rows)}
            
            if(input$clustering_distance_columns == "pearson"| input$clustering_distance_columns == "spearman" | input$clustering_distance_columns == "kendall"){
              distfun_col = function(x) as.dist(1 - cor(t(x), method = input$clustering_distance_columns))
            }else{
              distfun_col =  function(x) stats::dist(x, method = input$clustering_distance_columns)}
            
            mat_data <- mat()
            
            mat_data <- mat_data[ , colSums(is.na(mat_data))==0]
            
            hclustfun_row = function(x) stats::hclust(x, method = input$clustering_method_rows)
            hclustfun_col = function(x) stats::hclust(x, method = input$clustering_method_columns)
            
            heatmap_obj <- heatmaply(mat_data, 
                                     fontsize_row = 9 , 
                                     colors = rev(brewer.pal(n= 10, "RdBu")) , 
                                     showticklabels = c(FALSE, TRUE) ,
                                     distfun_row = distfun_row,
                                     distfun_col = distfun_col,
                                     hclustfun_row = hclustfun_row,
                                     hclustfun_col = hclustfun_col,
                                     scale = input$scale_by_row_column,
                                     plot_method = "plotly")
          })
        } else {
          
          # If categorical annotations are selected to be shown on top.
          
          validate(need(input$hm_definition_sel, "Please select a sample type"))
          validate(need(input$selectors, "Please select a Human MSigDB Gene set or create a gene set"))
          if(input$selectors == 'Manually enter gene names'){
            validate(need(input$genes, "Create genes to create your gene set"))
          } else if (input$selectors == 'MSigDB gene sets'){
            validate(need(input$cat, "Please select a Human MSigDB Collection"))
            validate(need(input$chosen_gse, "Please select the subset of your chosen Human MSigDB Collection"))
          }
          validate(need(input$heatmap_run, "Please click the 'Analyze' button"))
          
          if(input$clustering_distance_rows == "pearson"| input$clustering_distance_rows == "spearman" | input$clustering_distance_rows == "kendall"){
            distfun_row = function(x) as.dist(1 - cor(t(x), method=input$clustering_distance_rows))
          } else{
            distfun_row = function(x) stats::dist(x, method = input$clustering_distance_rows)}
          
          if(input$clustering_distance_columns == "pearson"| input$clustering_distance_columns == "spearman" | input$clustering_distance_columns == "kendall"){
            distfun_col = function(x) as.dist(1 - cor(t(x), method =input$clustering_distance_columns))
          }else{
            distfun_col =  function(x) stats::dist(x, method = input$clustering_distance_columns)}
          
          mat_data <- mat()
          
          # mat_data <- mat_data[ , colSums(is.na(mat_data))==0]
          
          hclustfun_row = function(x) stats::hclust(x, method = input$clustering_method_rows)
          hclustfun_col = function(x) stats::hclust(x, method = input$clustering_method_columns)
          
          heatmap_obj <- heatmaply(mat_data, 
                                   fontsize_row = 9 , 
                                   col_side_colors = meta(),
                                   colors = rev(brewer.pal(n= 10, "RdBu")) , 
                                   showticklabels = c(FALSE, TRUE) ,
                                   distfun_row = distfun_row,
                                   distfun_col = distfun_col,
                                   hclustfun_row = hclustfun_row,
                                   hclustfun_col = hclustfun_col,
                                   scale = input$scale_by_row_column,
                                   plot_method = "plotly")
        }})
      
      output$heatmap_plot <- renderPlotly({
        
        # Output for the heatmap plot
        
        req(heatmap_object())
        
        print(heatmap_object())
        
      })
      
      
      ## Since it is plotly object, it will be downloaded in hml format.
      
      output$downloadPlot6 <- downloadHandler(
        filename = function() {'heatmap.html'},
        content = function(file) {
          
          htmlwidgets::saveWidget(as_widget(heatmap_object()), file)
          
        })
      
    }
  )
}
