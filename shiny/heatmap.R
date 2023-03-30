# Load the packages
library(pheatmap)
library(shiny)
library(RColorBrewer)
library(data.table)
library(msigdbr)
library(ComplexHeatmap)
library(dplyr)
library(shinyWidgets)
library(rintrojs)

# Define UI for application that draws a histogram
heatmap_ui <- function(id, label, choices) {
  ns <- NS(id)
  tagList(
    
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
                     label = "1. Select sample types",
                     choices=NULL, # will be updated dynamically
                     options=list(placeholder = "eg. Primary solid tumor")),
      selectizeInput(inputId = ns("selectors"),
                     label= "2. Choose if you would like to examine your own geneset or Human MSigDB genesetsi or upload an csv file including your genes of interest.",
                     choices = NULL,
                     multiple = FALSE),
      conditionalPanel(
        condition = "input.selectors == 'Create Gene Set'",
        ns=ns,
        selectizeInput(inputId = ns("genes"),
                       label= "3. Choose genes to create your gene set",
                       choices = NULL,
                       multiple = TRUE,
                       options=list(placeholder = "eg. TSPAN6, TNMD etc."))
      ),
      conditionalPanel(
        condition= "input.selectors == 'Human MSigDB Gene sets'",
        ns=ns,
        selectInput(ns("cat"), "3. Please select a Human MSigDB Collection", 
                    choices = list(`Gene Sets` = list("H", "C1", "C6", "C8"),
                                   `C2 Subcategories` = list("CGP", "CP", "CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS"),
                                   `C3 Subcategories` = list("MIR:MIR_Legacy", "MIR:MIRDB", "TFT:GTRD", "TFT:TFT_Legacy"),
                                   `C4 Subcategories` = list("CGN", "CM"),
                                   `C5 Subcategories` = list("GO:BP", "GO:CC", "GO:MF", "HPO"),
                                   `C7 Subcategories` = list("IMMUNESIGDB", "VAX")
                                   
                    )),
        selectizeInput(ns("chosen_gse"), 
                       "3.1 Please select the subset of your chosen Human MSigDB Collection", 
                       choices = NULL)
      ),
      conditionalPanel(
        condition = "input.selectors == 'Upload a csv File'",
        ns=ns,
        fileInput(ns("heatmap_csv"), 
                  label = "3. Please upload your csv file.",
                  multiple = FALSE)
      ),
      
      sliderInput(inputId = ns("variable"),
                               label = "Only Show Top Variable Genes",
                  min = 0, max = 1, value = 0
      ),
  
      selectizeInput(ns("annotation"), 
                     "4. Please select annotation variables (optional)", 
                     choices = NULL,
                     multiple = TRUE,
                     options=list(placeholder = "eg. meta.definition, meta.LYMPHOCYTE.SCORE etc.")),
      selectizeInput(ns("hm_categorized_gene"), 
                     "5. Please select a gene/genes to be categorized as high/low (optional)", 
                     choices = NULL,
                     multiple = TRUE,
                     options=list(placeholder = "eg. TSPAN6, TNMD etc.")),
      p("Click 'Create Heatmap' button in the bottom of the page each time after a parameter is changed."),
      materialSwitch(inputId = ns("cluster_rows"),
                     label = "Cluster the rows",
                     status = "success",
                     inline = TRUE,
                     value = TRUE),
      materialSwitch(inputId = ns("cluster_columns"), 
                     label = "Cluster the columns",
                     status = "success",
                     inline = TRUE,
                     value = TRUE),
      selectizeInput(inputId = ns("clustering_distance_rows"), 
                     multiple=F,
                     label = "choose the distance method for rows",
                     choices=c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")),
      selectizeInput(inputId = ns("clustering_distance_columns"), 
                     multiple=F,
                     label = "Choose the distance method for columns",
                     choices=c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")),
      
      selectizeInput(inputId = ns("clustering_method_rows"), 
                     multiple=F,
                     label = "Choose the clustering method for rows",
                     choices=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")),
      
      selectizeInput(inputId = ns("clustering_method_columns"), 
                     multiple=F,
                     label = "Choose the clustering method for columns",
                     choices=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")),
      actionBttn(inputId = ns("heatmap_run"), 
                 label = "Create Heatmap",
                 style = "unite",
                 block = TRUE,
                 color = "primary"),
      br(),
      
      downloadButton(ns("downloadPlot6"), "Download Heatmap Plot", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
      
      #help section UI
      
      introjsUI(),
      actionButton(ns("heatmap_help"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
      
    ),
    
    mainPanel(
      plotOutput(outputId = ns("heatmap_plot"))
    )
  )
  
}

# Define server logic required to draw a histogram
heatmap_server <- function(id,Xproj) {
  moduleServer(
    id,
    function(input, output, session) {
      
      ns <- session$ns
      
      # Heatmap App Tutorial
      
      heatmap_help_tutorial <- reactive({
        
        if(input$selectors == "Create Gene Set"){
          
          return(
            data.frame(
              
              element = paste0("#", session$ns(c(NA, "hm_definition_sel + .selectize-control", "selectors + .selectize-control ", "genes + .selectize-control", "variable + .selectize-control", "annotation + .selectize-control ", "hm_categorized_gene + .selectize-control", "cluster_rows + .selectize-control", "clustering_distance_rows + .selectize-control", "clustering_distance_columns + .selectize-control", "clustering_method_rows + .selectize-control", "clustering_method_columns + .selectize-control"))),
              
              intro = paste(c(
                "This is the Heatmap Generator App, press the buttons to learn features of the app.",
                "You can choose one or more than one tumor types to create a plot including certain tumor types.",
                "You can either create your own gene set or choose genesets from the Human MSigDB Collections.",
                "Here, you can create your own geneset to observe on the heatmap.",
                "Only Show Top Variable Genes: If this toggle is on, only the top differentially expressed genes will be filtered and shown on the heatmap.",
                "You can add categorical annotation parameters here.",
                "You can choose gene or genes to take its expression value. For the selected gene, patients with lower expression value than the median will be taken as 'low' while patients with higher expression value will be taken as 'high' on the heatmap.",
                "Cluster the rows: 
                Cluster the columns: ",
                "Desired distance methods for rows can be chosen.",
                "Desired distance methods for columns can be chosen.",
                "Clustering methods for rows can be chosen.",
                "Clustering methods for columns can be chosen."
                ))
            )
          )
        } else {
          return(
            
            data.frame(
              
              element = paste0("#", session$ns(c(NA, "hm_definition_sel + .selectize-control", "selectors + .selectize-control ", "cat + .selectize-control", "chosen_gse + .selectize-control", "variable + .selectize-control", "annotation + .selectize-control ", "hm_categorized_gene + .selectize-control", "cluster_rows + .selectize-control", "clustering_distance_rows + .selectize-control", "clustering_distance_columns + .selectize-control", "clustering_method_rows + .selectize-control", "clustering_method_columns + .selectize-control"))),
              
              intro = paste(c(
                "This is the Heatmap Generator App, press the buttons to learn features of the app.",
                "You can choose one or more than one tumor types to create a plot including certain tumor types.",
                "You can either create your own gene set or choose genesets from the Human MSigDB Collections.",
                "Here you can choose between different Human MSigDB Collections and their subcategories.",
                "Each Human MSigDB Collection and subcategory have different sets in it. Here you can filter between those sets. Plot is not created untill these sets are chosen.",
                "Only Show Top Variable Genes: If this toggle is on, only the top differentially expressed genes will be filtered and shown on the heatmap.",
                "You can add categorical annotation parameters here.",
                "You can choose gene or genes to take its expression value. For the selected gene, patients with lower expression value than the median will be taken as 'low' while patients with higher expression value will be taken as 'high' on the heatmap.",
                "Cluster the rows: 
                Cluster the columns: ",
                "Desired distance methods for rows can be chosen.",
                "Desired distance methods for columns can be chosen.",
                "Clustering methods for rows can be chosen.",
                "Clustering methods for columns can be chosen."
              ))
            )
            
          )
        }
        
      })
      
      
      observeEvent(input$heatmap_help, {
        
        introjs(session, options = list(steps = heatmap_help_tutorial() ) )
        
      })
      
      # Reactive selectboxes
      
      observeEvent(input$cat,{
        req(input$cat)
        if(length(as.vector(input$selectors)) == "Create Gene Set"){
          return()
        }
        updateSelectizeInput(session = getDefaultReactiveDomain(),"chosen_gse", choices = gene_sets()$gs_name, selected = character(0) ,server = TRUE)
      })
      
      observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "selectors", choices = c("Create Gene Set", "Human MSigDB Gene sets", "Upload a csv File"), server = TRUE)})
      observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "genes", choices = colnames(Xproj$a()[, lapply(Xproj$a(), is.numeric) == TRUE, with = FALSE]), server = TRUE)})
      observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "annotation", choices = colnames(Xproj$a()[, lapply(Xproj$a(), is.factor) == TRUE, with = FALSE]), selected = character(0), server = TRUE)})
      observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "hm_categorized_gene", choices = colnames(Xproj$a()[, lapply(Xproj$a(), is.numeric) == TRUE, with = FALSE]), selected = character(0), server = TRUE)})
      observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "hm_definition_sel", choices = c(Xproj$a()$meta.definition), selected = character(0), server = TRUE)})
      
      gene_sets <- reactive({ 
        
        # Downloading the MSigDB datasets if selected.
        
        if(input$cat %in% c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")){
          
          return({
            
            genset2 <- msigdbr(species = "human", category = input$cat)
            
            genset2
            
          })
          
          
        } else if(input$cat %in% c("CGP", "CP", "CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS",
                                   "MIR:MIR_Legacy", "MIR:MIRDB", "TFT:GTRD", "TFT:TFT_Legacy",
                                   "CGN", "CM",
                                   "GO:BP", "GO:CC", "GO:MF", "HPO",
                                   "IMMUNESIGDB", "VAX")){
          
          
          return({
            
            genset3 <- msigdbr(species = "human", category = NULL, subcategory = input$cat)
            
            genset3
            
          })}
        
      })
      
      
      pre_data <- eventReactive(input$heatmap_run, {
        
        # Preparing the preliminary data.
        
        if(length(as.vector(input$hm_definition_sel)) > 0) {
          return({
            daf1 <- Xproj$a()
            daf1 <- setDT(daf1, key = 'meta.definition')[J(input$hm_definition_sel)]
            daf1
          })
        } else if (length(as.vector(input$hm_hm_definition_sel)) == 0){
          return({
            daf2 <- Xproj$a()
            daf2
          })
        }
        
      })
      
      mat <- eventReactive(input$heatmap_run, {
        
        # Matching the MSigDB gene sets with selected TCGA data.
        
        req(pre_data())
        
        if(input$selectors == "Human MSigDB Gene sets"){
          
          req(gene_sets())
          
          req(input$chosen_gse)
          
          final_gene_sets <- gene_sets()
          
          final_gene_sets <- setDT(final_gene_sets, key = 'gs_name')[J(input$chosen_gse)]
          
          daf <- as.data.frame(pre_data())
          
          rownames(daf) <- daf$meta.sample
          
          same_hallmarks_names = intersect(final_gene_sets$gene_symbol, colnames(daf))
          
          selected_cols <- c(same_hallmarks_names)
          
          daf <- daf[,selected_cols]  
          
        } else if(input$selectors == "Create Gene Set"){
          
          req(input$genes)
          
          daf <- as.data.frame(pre_data())
          
          rownames(daf) <- daf$meta.sample
          
          selected_cols <- c(input$genes)
          
          daf <- daf[,selected_cols]  
          
        } else if(input$selectors == "Upload a csv File"){

          req(input$heatmap_csv)
          
          daf <- as.data.frame(pre_data())
          
          rownames(daf) <- daf$meta.sample
          
          uploaded_heatmap_csv <- input$heatmap_csv

          selected_csv <- read.csv(uploaded_heatmap_csv$datapath, stringsAsFactors = FALSE, header = FALSE)$V1
          
          same_gene_names = intersect(selected_csv, colnames(daf))
          
          selected_cols <- c(same_gene_names)
          
          daf <- daf[,selected_cols]  
          
        }
        
        if (input$variable !=1) {
          
          daf_variable <- t(daf)
          
          evar <- apply(daf_variable,1,var)
          
          mostVariable <- daf_variable[evar>quantile(evar,as.numeric(input$variable)),]
          mostvariablegenes <- rownames(mostVariable)
          daf <- daf[,mostvariablegenes]
        } else if(input$variable == 0){
          daf
        }
        
        # When you have several variables to examine over multiple scales, the scale() function makes more sense. 
        # One variable, for example, is of magnitude 100, whereas another is of magnitude 1000.
        
        daf <- scale(daf, scale = FALSE)
        
        daf <- as.matrix(daf)
        
        daf <- t(daf)
      })
      
      meta <- eventReactive(input$heatmap_run, {
        
        # Preparing the meta data for annotation bar.
        
        req(pre_data())
        
        if(length(as.vector(input$annotation)) > 0 & length(as.vector(input$hm_categorized_gene)) > 0){
          
          return({
            
            # When there are both categorized genes and meta annotation
            
            daf <- as.data.frame(pre_data())
            
            rownames(daf) <- daf$meta.sample
            
            hm_categorized_gene <- c(input$hm_categorized_gene)
            
            chosen_meta <- c(input$annotation, input$hm_categorized_gene)
            
            meta <- daf[,chosen_meta]
            
            # High-low expression categorization of selected genes
            
            meta <- as.data.table(meta)
            
            meta <- meta[, hm_categorized_gene_means := rowMeans(.SD, na.rm = TRUE), .SDcols = c(hm_categorized_gene)]  
            
            cat_gene_median <-  median(meta$hm_categorized_gene_means, na.rm = TRUE)
            
            meta <- meta[,c(hm_categorized_gene):=NULL]
            
            meta <- as.data.frame(meta)
            
            meta <- meta %>% 
              mutate(hm_categorized_gene_high_low = case_when(
                meta$hm_categorized_gene_means > cat_gene_median ~ "high",
                meta$hm_categorized_gene_means < cat_gene_median ~ "low"
              )) %>% 
              select(-hm_categorized_gene_means)
            
          })} else if (length(as.vector(input$annotation)) > 0 & length(as.vector(input$hm_categorized_gene)) == 0){
            
            return({
              
              # When there is no categorized genes entered
              
              daf <- as.data.frame(pre_data())
              
              rownames(daf) <- daf$meta.sample
              
              chosen_meta <- c(input$annotation)
              
              meta <- daf[,chosen_meta]
              
            })
            
          } else if (length(as.vector(input$annotation)) == 0 & length(as.vector(input$hm_categorized_gene)) > 0){
            
            return({
              
              # When there are categorized genes but no meta annotation
              
              req(input$hm_categorized_gene)
              
              daf <- as.data.frame(pre_data())
              
              rownames(daf) <- daf$meta.sample
              
              hm_categorized_gene <- c(input$hm_categorized_gene)
              
              meta <- daf[,input$hm_categorized_gene]
              
              # High-low expression categorization of selected genes
              
              meta <- as.data.table(meta)
              
              if (length(as.vector(input$hm_categorized_gene)) == 1){
                
                setnames(meta, "meta", "hm_categorized_gene_means")
                
              } else {
                
                meta <- meta[, hm_categorized_gene_means := rowMeans(.SD, na.rm = TRUE), .SDcols = c(hm_categorized_gene)]  
                
              }
              
              
              cat_gene_median <-  median(meta$hm_categorized_gene_means, na.rm = TRUE)
              
              meta <- meta[,c(hm_categorized_gene):=NULL]
              
              meta <- as.data.frame(meta)
              
              meta <- meta %>% 
                mutate(hm_categorized_gene_high_low = case_when(
                  meta$hm_categorized_gene_means > cat_gene_median ~ "high",
                  meta$hm_categorized_gene_means < cat_gene_median ~ "low"
                )) %>% 
                select(-hm_categorized_gene_means)
              
            
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
          if(input$selectors == 'Create Gene Set'){
            validate(need(input$genes, "Create genes to create your gene set"))
          } else if (input$selectors == 'Human MSigDB Gene sets'){
            validate(need(input$cat, "Please select a Human MSigDB Collection"))
            validate(need(input$chosen_gse, "Please select the subset of your chosen Human MSigDB Collection"))
          }
          validate(need(input$heatmap_run, "Please click the create heatmap button"))
          
          
          return({
            heatmap_obj <- Heatmap(mat(),
                                   name = "mat",
                                   cluster_rows = input$cluster_rows,
                                   cluster_columns = input$cluster_columns,
                                   clustering_distance_rows = input$clustering_distance_rows,
                                   clustering_distance_columns = input$clustering_distance_columns,
                                   clustering_method_rows = input$clustering_method_rows,
                                   clustering_method_columns = input$clustering_method_columns,
                                   show_column_names = FALSE)
            
            heatmap_obj <- draw(heatmap_obj)
            heatmap_obj 
          })
        } else {
          
          # If categorical annotations are selected to be shown on top.
          
          validate(need(input$hm_definition_sel, "Please select a sample type"))
          validate(need(input$selectors, "Please select a Human MSigDB Gene set or create a gene set"))
          if(input$selectors == 'Create Gene Set'){
            validate(need(input$genes, "Create genes to create your gene set"))
          } else if (input$selectors == 'Human MSigDB Gene sets'){
            validate(need(input$cat, "Please select a Human MSigDB Collection"))
            validate(need(input$chosen_gse, "Please select the subset of your chosen Human MSigDB Collection"))
          }
          validate(need(input$heatmap_run, "Please click the create heatmap button"))
          
          
          column_ha <- HeatmapAnnotation(df = meta())
          
          heatmap_obj <- Heatmap(mat(),
                                 name = "mat",
                                 cluster_rows = input$cluster_rows,
                                 cluster_columns = input$cluster_columns,
                                 clustering_distance_rows = input$clustering_distance_rows,
                                 clustering_distance_columns = input$clustering_distance_columns,
                                 clustering_method_rows = input$clustering_method_rows,
                                 clustering_method_columns = input$clustering_method_columns,
                                 top_annotation = column_ha,
                                 show_column_names = FALSE)
          
          heatmap_obj <- draw(heatmap_obj)
          heatmap_obj 
        }
      })
      
      output$heatmap_plot <- renderPlot({
        
        # Output for the heatmap plot
        
        req(heatmap_object())
        
        print(heatmap_object())
        
      })
      
      
      output$downloadPlot6 <- downloadHandler(
        filename = function() {
          paste("Heatmap_plot.png")
        },
        content = function(file) {
          
          png(file)
          print(heatmap_object())
          dev.off()
          
        })
      
    }
  )
}
