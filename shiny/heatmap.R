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
                     label= "2. Choose how you would like to select genes to be shown",
                     choices = c("Manually enter gene names", "MSigDB gene sets", "Upload a csv file"),
                     multiple = FALSE),
      conditionalPanel(
        condition = "input.selectors == 'Manually enter gene names'",
        ns=ns,
        selectizeInput(inputId = ns("genes"),
                       label= "3. Type gene names",
                       choices = NULL,
                       multiple = TRUE,
                       options=list(placeholder = "eg. TSPAN6, TNMD etc."))
      ),
      conditionalPanel(
        condition= "input.selectors == 'MSigDB gene sets'",
        ns=ns,
        selectInput(ns("cat"), "3. Select an MSigDB Collection", 
                    choices = list(`H: HALLMARK gene sets` = list("H"), 
                                   `C1: Positional gene sets` = list("C1"),
                                   `C2: Curated gene sets` = list("CGP", "CP", "CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS"),
                                   `C3: Regulatory target gene sets` = list("MIR:MIR_Legacy", "MIR:MIRDB", "TFT:GTRD", "TFT:TFT_Legacy"),
                                   `C4: Computational gene sets` = list("CGN", "CM"),
                                   `C5: Ontology gene sets` = list("GO:BP", "GO:CC", "GO:MF", "HPO"),
                                   `C6: Oncogenic gene sets` = list("C6"),
                                   `C7: Immunologic gene sets` = list("IMMUNESIGDB", "VAX"),
                                   `C8: Cell type signature gene sets` = list("C8")
                                   
                    )),
        selectizeInput(ns("chosen_gse"), 
                       "3.1 Please select a specific gene set", 
                       choices = NULL)
      ),
      conditionalPanel(
        condition = "input.selectors == 'Upload a csv file'",
        ns=ns,
        fileInput(ns("heatmap_csv"), 
                  # label = "3. Please upload your csv file.",
                  label = tags$span(
                    "3. Please upload your csv file.", 
                    tags$i(
                      class = "glyphicon glyphicon-info-sign", 
                      style = "color:#0072B2;",
                      title = "The csv file should contain two unnamed columns: the first column should contain the gene set name, and the second column should contain gene names. Each gene should be associated with a gene set (ie. no missing data), and multiple gene sets can be provided in one file."
                    )),
                  multiple = FALSE)
      ),
      
      sliderInput(inputId = ns("variable"),
                  label = "Variation filter (Keep top n% variable genes)",
                  min = 0, max = 100, value = 0
      ),
      
      span(style="color:#3382FF",
           
           selectizeInput(ns("annotation"), 
                          "4. Please select annotation variables (optional)", 
                          choices = NULL,
                          multiple = TRUE,
                          options=list(placeholder = "eg. meta.definition, meta.LYMPHOCYTE.SCORE etc.")),
           selectizeInput(ns("hm_categorized_gene"), 
                          "5. Please select a gene/genes to be categorized as high/low (optional)", 
                          choices = NULL,
                          multiple = TRUE,
                          options=list(placeholder = "eg. TSPAN6, TNMD etc."))
           
      ),
      
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
      plotlyOutput(outputId = ns("heatmap_plot"), width = "950px", height = "900px")
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
        
        if(input$selectors == "Manually enter gene names"){
          
          return(
            data.frame(
              
              element = paste0("#", session$ns(c(NA, "hm_definition_sel + .selectize-control", "selectors + .selectize-control ", "genes + .selectize-control", "variable + .selectize-control", "annotation + .selectize-control ", "hm_categorized_gene + .selectize-control",
                                                 "clustering_distance_rows + .selectize-control", "clustering_distance_columns + .selectize-control", "clustering_method_rows + .selectize-control", "clustering_method_columns + .selectize-control"))),
              
              intro = paste(c(
                "This is the heatmap module where you can visualize expression patterns of selected genes. Continue the tutorial to learn how to use this module",
                "You can select sample types to focus the analysis on the specific subsets.",
                "Here you can choose how you would like to select genes for the plot",
                "You can manually type genes of interest here (for other input possibilities, go back to the previous selection box)",
                "Next, you can apply a variance filter to keep only highly variable genes in the plot. 100 (default) means no filter is applied. If you like to see top 10% variable genes only, set this value to 10. Such filtering can help see more informative genes.",
                "You can select categorical clinical meta data features to show as annotations on top of the heatmap.",
                "You can also create an annotation bar by categorizing the patients based on their gene expression levels. You can specify one or more genes here. When multiple genes are entered, their average is calculated. Patients are categorized as 'high' and 'low' according to the median gene expression value",
                "You can choose how the distance will be calculated for genes here",
                "You can choose how the distance will be calculated for samples here",
                "You can choose different hierarchical clustering methods for genes here",
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
                "Next, you can apply a variance filter to keep only highly variable genes in the plot. 100 (default) means no filter is applied. If you like to see top 10% variable genes only, set this value to 10. Such filtering can help see more informative genes.",
                "You can select categorical clinical meta data features to show as annotations on top of the heatmap.",
                "You can also create an annotation bar by categorizing the patients based on their gene expression levels. You can specify one or more genes here. When multiple genes are entered, their average is calculated. Patients are categorized as 'high' and 'low' according to the median gene expression value",
                "You can choose how the distance will be calculated for genes here",
                "You can choose how the distance will be calculated for samples here",
                "You can choose different hierarchical clustering methods for genes here",
                "You can choose different hierarchical clustering methods for samples here"
              ))
            )
            
          )
        } else if (input$selectors == "Upload a csv file"){
          return(
            
            data.frame(
              
              element = paste0("#", session$ns(c(NA, "hm_definition_sel + .selectize-control", "selectors + .selectize-control ", "heatmap_csv + .selectize-control", "variable + .selectize-control", "annotation + .selectize-control ", "hm_categorized_gene + .selectize-control", 
                                                 "clustering_distance_rows + .selectize-control", "clustering_distance_columns + .selectize-control", "clustering_method_rows + .selectize-control", "clustering_method_columns + .selectize-control"))),
              
              intro = paste(c(
                "This is the heatmap module where you can visualize expression patterns of selected genes. Continue the tutorial to learn how to use this module",
                "You can select sample types to focus the analysis on the specific subsets.",
                "Here you can choose how you would like to select genes for the plot",
                "You can upload a csv file including your genes of interest to see them on the heatmap",
                "Next, you can apply a variance filter to keep only highly variable genes in the plot. 100 (default) means no filter is applied. If you like to see top 10% variable genes only, set this value to 10. Such filtering can help see more informative genes.",
                "You can select categorical clinical meta data features to show as annotations on top of the heatmap.",
                "You can also create an annotation bar by categorizing the patients based on their gene expression levels. You can specify one or more genes here. When multiple genes are entered, their average is calculated. Patients are categorized as 'high' and 'low' according to the median gene expression value",
                "You can choose how the distance will be calculated for genes here",
                "You can choose how the distance will be calculated for samples here",
                "You can choose different hierarchical clustering methods for genes here",
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
      
      observeEvent(input$cat,{
        req(input$cat)
        if(length(as.vector(input$selectors)) == "Manually enter gene names"){
          return()
        }
        updateSelectizeInput(session = getDefaultReactiveDomain(),"chosen_gse", choices = gene_sets()$gs_name, selected = character(0) ,server = TRUE)
      })
      
      # observe({updateSelectizeInput(session = getDefaultReactiveDomain(), "selectors", choices = c("Manually enter gene names", "MSigDB gene sets", "Upload a csv file"), server = TRUE)})
      
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
        
        if(input$selectors == "MSigDB gene sets"){
          
          req(gene_sets())
          
          req(input$chosen_gse)
          
          final_gene_sets <- gene_sets()
          
          final_gene_sets <- setDT(final_gene_sets, key = 'gs_name')[J(input$chosen_gse)]
          
          daf <- as.data.frame(pre_data())
          
          rownames(daf) <- daf$meta.barcode
          
          same_hallmarks_names = intersect(final_gene_sets$gene_symbol, colnames(daf))
          
          selected_cols <- c(same_hallmarks_names)
          
          daf <- daf[,selected_cols]  
          
        } else if(input$selectors == "Manually enter gene names"){
          
          req(input$genes)
          
          daf <- as.data.frame(pre_data())
          
          rownames(daf) <- daf$meta.barcode
          
          selected_cols <- c(input$genes)
          
          daf <- daf[,selected_cols]  
          
        } else if(input$selectors == "Upload a csv file"){
          
          req(input$heatmap_csv)
          
          daf <- as.data.frame(pre_data())
          
          rownames(daf) <- daf$meta.barcode
          
          uploaded_heatmap_csv <- input$heatmap_csv
          
          selected_csv <- read.csv(uploaded_heatmap_csv$datapath, stringsAsFactors = FALSE, header = FALSE)$V1
          
          same_gene_names = intersect(selected_csv, colnames(daf))
          
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
            
            rownames(daf) <- daf$meta.barcode
            
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
          if(input$selectors == 'Manually enter gene names'){
            validate(need(input$genes, "Create genes to create your gene set"))
          } else if (input$selectors == 'MSigDB gene sets'){
            validate(need(input$cat, "Please select a Human MSigDB Collection"))
            validate(need(input$chosen_gse, "Please select the subset of your chosen Human MSigDB Collection"))
          }
          validate(need(input$heatmap_run, "Please click the create heatmap button"))
          
          
          return({
            distfun_row = function(x) stats::dist(x, method = input$clustering_distance_rows)
            distfun_col =  function(x) stats::dist(x, method = input$clustering_distance_columns)
            
            hclustfun_row = function(x) stats::hclust(x, method = input$clustering_method_rows)
            hclustfun_col = function(x) stats::hclust(x, method = input$clustering_method_columns)
            
            
            heatmap_obj <- heatmaply(mat(), 
                                     fontsize_row = 9 , 
                                     colors = rev(brewer.pal(n= 10, "RdBu")) , 
                                     showticklabels = c(FALSE, TRUE) ,
                                     distfun_row = distfun_row,
                                     distfun_col = distfun_col,
                                     hclustfun_row = hclustfun_row,
                                     hclustfun_col = hclustfun_col,
                                     plot_method = "plotly") 
            
            heatmap_obj
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
          validate(need(input$heatmap_run, "Please click the create heatmap button"))
          
          distfun_row = function(x) stats::dist(x, method = input$clustering_distance_rows)
          distfun_col =  function(x) stats::dist(x, method = input$clustering_distance_columns)
          
          hclustfun_row = function(x) stats::hclust(x, method = input$clustering_method_rows)
          hclustfun_col = function(x) stats::hclust(x, method = input$clustering_method_columns)
          
          
          heatmap_obj <- heatmaply(mat(), 
                                   fontsize_row = 9 , 
                                   col_side_colors = meta(),
                                   colors = rev(brewer.pal(n= 10, "RdBu")) , 
                                   showticklabels = c(FALSE, TRUE) ,
                                   distfun_row = distfun_row,
                                   distfun_col = distfun_col,
                                   hclustfun_row = hclustfun_row,
                                   hclustfun_col = hclustfun_col,
                                   plot_method = "plotly") 
          
          heatmap_obj
        }
      })
      
      output$heatmap_plot <- renderPlotly({
        
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
