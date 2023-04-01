library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(plotly)
library(factoextra)
library(kableExtra)
library(janitor)
library(msigdbr)
library(shinybusy)
library(data.table)
library(rintrojs)
library(plotly)
library(shinyvalidate)
library(shinyWidgets)

pca_ui <- function(id) {
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
      

      
      selectizeInput(ns("genecor_samp_2"), multiple=T,
                     "Please select sample types",
                     choices=NULL,
                     options=list(placeholder = "eg. Primary solid tumor")),
      
      selectizeInput(ns("data"), "*Please select input genes", 
                     choices = NULL,
                     options=list(placeholder = "eg. RNAseq")),
      
      checkboxInput(ns("center"), "Center data", value = TRUE),
      checkboxInput(ns("scale"), "Scale variables", value = TRUE),
      
      
      ## conditional panel for gene sets
      # pathway can be chosen according to gene set input
      
      conditionalPanel(
        
        condition = 'input.data == "Custom gene set"', ns = ns,
        
        
        #'[##########################################################################################]
        #'[##########################################################################################]
        #'[Please show a hover info icon here to tell the user about the file format. See example in the GSEA module]
        #'[##########################################################################################]
        fileInput(inputId = ns("pca_up"),
                  label = "Upload your gene set as .csv file",
                  accept = c("text/csv", "text/comma-separated-values,text/plain",
                             ".csv")), 
        
        #'[##########################################################################################]
        #'[##########################################################################################]
        
        
        
      ),
      
      conditionalPanel(
        
        condition = 'input.data == "MSigDB Gene Sets"', ns = ns,
        
        selectizeInput(ns("cat"), "Please select an MSigDB Collection", 
                       choices = NULL,
                    options=list(placeholder = "eg. CM - Cancer Modules")),
        
        selectizeInput(ns("chosen_gs"), 
                       "*Please select a pathway", 
                       choices = NULL,
                       options=list(placeholder = "eg. MODULE_1"))
        
        
      ),
      
      checkboxInput(ns("variance"), "Apply variance filtering", value = F),
      
      conditionalPanel(
        
        condition = 'input.variance', ns = ns,
        
        sliderInput(inputId = ns("var_gene"),
                    label = "Keep genes with the top n% variation",
                    min = 0, max = 100, value = 10
        )
          
      ),
      
      selectizeInput(ns("gene"), "*Please select feature to annotate", 
                     options=list(placeholder = "eg. meta.gender, gene name"),
                     choices = NULL),

      selectizeInput(ns("palette"), "*Please select a color palette",
                     choices = NULL,
                     options=list(placeholder = "eg. Nature")

      ),

      actionBttn(inputId = ns("did"), 
                 label = "Perform PCA",
                 style = "unite",
                 block = TRUE,
                 color = "primary"),
      br(),
      
      #help section UI
      
      introjsUI(),
      actionButton(ns("intro2"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7")
      
    ),
    
    mainPanel(
      
      plotlyOutput(outputId = ns("ex"))
      
      
    )
  )
  
  
}





pca_server <- function(id,Xproj) {
  moduleServer(
    id,
    function(input, output, session) {
      
      ## help section server
      
      ns <- session$ns
      
      steps <- reactive({
        
        data.frame(
          
          element = paste0("#", session$ns(c(NA, "center", "scale", "genecor_samp_2 + .selectize-control" ,"data + .selectize-control", "variance","gene + .selectize-control", "palette + .selectize-control"))),
          
          intro = paste(c(
            "This is the Principal Compenent Analysis (PCA) module. PCA is a dimensionality reduction method allowing multidimensional datasets to be visualized on a two dimensional graph. Principal components are derived from linear combinations of the original variables and they represent the variation in the data set. <b>NOTE:</b> Creating PCA plots can take some time, please be patient. Continue tutorial to learn more about this module.",
            "By default, gene expression values are centered by subtracting the mean expression value.",
            "A variable that is on a different scale from the others may dominate the variance direction. Scaling (default) gene expression values prevents this effect.",
            "You can select sample types to include in the analysis here.",
            "Here, you can select the genes to be used in PCA. You can select <b>i)<b/> a custom list of genes, <b>ii)<b/> all genes (RNAseq and miRNAseq data), <b>iii)<b/> miRNAs (mature miRNAs from miRNAseq), or <b>iv)<b/> genes annotated in MSigDB gene sets.",
            "Here, you can apply variance filtering to keep most highly variable genes in the analysis. Setting this value to 10, for instance, will select the genes having the top 10% highest variation in the dataset.",
            "You can color code the data points on the graph using gene expression values or clinical meta data. If you select a gene name here, gene expression will be categorized at the median value per sample and points will be annotated. You can also select a clinical meta data (eg. meta.gender) to color points accordingly.",
            "You can change the color palette of the graph here."
          ))
          
        )
        
      })
      
      
      observeEvent(input$intro2, {
        
        introjs(session, options = list(steps = steps() ) )
        
      })
      
      ##validate
      iv <- InputValidator$new()

      iv$add_rule("var_gene", ~ if (input$var_gene == 0) "Your choice cannot be zero")

      iv$enable()
      
      #inputs
      
      observe({updateSelectizeInput(session, "genecor_samp_2",choices = Xproj$a()[["meta.definition"]], server = T)})
      observe(updateSelectizeInput(session, 'data', 
                                   choices = c("All genes", "miRNA", "RNAseq", "MSigDB Gene Sets", "Custom gene set"),
                                   selected = "",
                                   server = TRUE))
      
      observe(updateSelectizeInput(session, 'cat', 
                                   choices = list(`Gene Sets` = list("H - Hallmark Gene Sets" = "H",
                                                                     "C1 - Positional Gene Sets" = "C1",
                                                                     "C6 - Oncogenic Signature Gene Sets" = "C6", 
                                                                     "C8 - Cell Type Signature Gene Sets" = "C8"),
                                                  `C2 - Curated Gene Sets` = list("CGP - Chemical and Genetic Perturbations" = "CGP", 
                                                                                  "CP - Canonical Pathways" = "CP", 
                                                                                  "BioCarta Pathway Subset of CP" = "CP:BIOCARTA", 
                                                                                  "KEGG Pathway Subset of CP" = "CP:KEGG", 
                                                                                  "PID Pathway Subset of CP" = "CP:PID", 
                                                                                  "Reactome Pathway Subset of CP" = "CP:REACTOME", 
                                                                                  "Wikipathways Pathway Subset of CP" = "CP:WIKIPATHWAYS"),
                                                  `C3 - Regulatory Target Gene Sets` = list("MIR_Legacy Subset of microRNA Targets (MIR)" = "MIR:MIR_Legacy",
                                                                                            "miRDB Subset of microRNA Targets (MIR)" = "MIR:MIRDB",
                                                                                            "GTRD Subset of Transcription Factor Targets (TFT)" =  "TFT:GTRD", 
                                                                                            "TFT_Legacy Subset of Transcription Factor Targets (TFT)" =  "TFT:TFT_Legacy"),
                                                  `C4 - Computational Gene Sets` = list("CGN - Cancer Gene Neighborhoods" = "CGN", 
                                                                                        "CM - Cancer Modules" =  "CM"),
                                                  `C6 - Oncogenic Signature Gene Sets` = list("BP - Biological Process Subsets of Gene Ontology (GO)" = "GO:BP", 
                                                                                              "CC - Cellular Component Subsets of Gene Ontology (GO)" = "GO:CC", 
                                                                                              "MF - Molecular Function Subsets of Gene Ontology (GO)" = "GO:MF", 
                                                                                              "HPO - Human Phenotype Ontology subsets of Gene Ontology (GO)" = "HPO"),
                                                  `C7 - Immunological Signature Gene Sets` = list("ImmuneSigDB Subset of C7" = "IMMUNESIGDB", 
                                                                                                  "VAX - Vaccine Response Gene Sets" =  "VAX")
                                                  
                                   ),
                                   selected = "",
                                   server = TRUE))
    
      observe(updateSelectizeInput(session, 'gene', 
                                   choices = colnames(Xproj$a()), 
                                   selected = "",
                                   server = TRUE))

      #'[##########################################################################################]
      #'[##########################################################################################]
      #'[Does this have to be updated dynamically?]

      observe(updateSelectizeInput(session, 'palette',
                                   choices = list(`Pre-made palettes` = list(
                                     "Nature"="npg",
                                     "Science"= "aaas",
                                     "Lancet"= "lancet",
                                     "New Engl J Med" = "nejm",
                                     "J Clin Onc"="jco",
                                     "JAMA" = "jama",
                                     "Int Genomics Viewer (IGV)" = "igv",
                                     "UCSC Genome Browser"="ucscgb"),
                                     `Brewer palettes` = list(
                                       "Reds", "Blues", "Greys",
                                       "Set1", "Set2", "Set3",
                                       "Rainbow" = "Spectral",
                                       "Red-Blue" = "RdBu",
                                       "Red-Yellow-Blue" = "RdYlBu",
                                       "Red-Yellow-Green" = "RdYlGn",
                                       "Red-Grey" = "RdGy",
                                       "Purple-Orange" = "PuOr",
                                       "Purple-Green" = "PRGn",
                                       "Pink-Green" = "PiYG",
                                       "Brown-Green" = "BrBG")),
                                   selected = "J Clin Onc",
                                   server = TRUE))
      
      #'[##########################################################################################]
      #'[##########################################################################################]

      observeEvent(input$cat,{
    
        req(input$cat)
        
        updateSelectizeInput(session, "chosen_gs", choices = gene_sets()[["gs_name"]],selected = "", server = TRUE)
        
      })
      
      ## PATHWAY SELECTION
      # gene sets are downloaded via msigdbr package from the internet 
      gene_sets <- reactive({ 
        
        if(input$cat %in% c("H", "C1", "C6", "C8")){
          
          return({
            
            gene_set_nonsub <- msigdbr(species = "human", category = input$cat)
            
            gene_set_nonsub
            
          })
          
          
        } else if(input$cat %in% c("CGP", "CP", "CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS",
                                   "MIR:MIR_Legacy", "MIR:MIRDB", "TFT:GTRD", "TFT:TFT_Legacy",
                                   "CGN", "CM",
                                   "GO:BP", "GO:CC", "GO:MF", "HPO",
                                   "IMMUNESIGDB", "VAX")){
          
          
          return({
            
            gene_set_sub <- msigdbr(species = "human", category = NULL, subcategory = input$cat)
            
            gene_set_sub
            
          })}
        
      })
      
      
      
      ### DATA SELECTION PART
      
      pca_df <- reactive({
        
        browser()
        
        gene_cols <- Xproj$a() %>%
              select(!starts_with("meta.")) %>%
              select(where(is.numeric))
        gene_names <- colnames(gene_cols)

        pre_d <- Xproj$a()[meta.definition %in% input$genecor_samp_2, ]
        
        sel_cols <- apply(select(pre_d, all_of(gene_names)), 2, sum)
        
        #'[Causes an error in GBM. All rows are dropped. We need to rethink this]
        #'[###############################################################################]
        #'[###############################################################################
        #'[###############################################################################
        # pre_d <- pre_d %>% drop_na(all_of(gene_names))
        #'[###############################################################################
        #'[###############################################################################
        #'[###############################################################################
        
        pre_d
        
        })
      
      all_genes <- reactive({
        
        pc_df <-  pca_df() %>% 
              select(!starts_with("meta.")) %>%
              select(where(is.numeric))
        
        pc_df <- remove_constant(pc_df)
        
        })
      
      cre_gene <- reactive({
        
        req(input$pca_up)
        
        c_gene <- read.table(input$pca_up$datapath, header = FALSE, sep = ";", quote = "\"'",
                             dec = ".")  
        
        colnames(c_gene)[1] <- "genes" 
        
        c_gene
        
      })
      
      int_dat <- reactive({
        
        if(input$data == "All genes" ){
          
          return(all_genes())
          
        }else if(input$data == "miRNA"){
          
          return({
            
            df_mirna <- select(all_genes(), starts_with("hsa."))
            
            
            df_mirna
          })
          
          
        }else if(input$data == "RNAseq"){
          
          return({
            
            df_rna <- select(all_genes(), !starts_with("hsa."))
            
            df_rna
            
          })
          
        }else if(input$data == "MSigDB Gene Sets" & input$did > 0) {
          
          
          pathway <- reactive({filter(gene_sets(), gs_name == input$chosen_gs)})
          df_pathway <- all_genes()[,intersect(pathway()[["human_gene_symbol"]], colnames(all_genes())),with=F]
          
          df_pathway
          
        }else if(input$data == "Custom gene set" & input$did > 0){
          

          
          df_created <-  all_genes()[,intersect(cre_gene()[["genes"]], colnames(all_genes())),with=F]
        
            
        }
        
      })
      
      ## keep highly variant genes
      
      p_dat <- reactive({
        
        # browser()
        
        if(input$variance){
          
          if (input$var_gene != 100) {
            
            # var_dat <- int_dat()
            # 
            # variance_dat <- t(var_dat)
            # 
            # evar <- apply(variance_dat,1,var)
            # 
            # mostVariant <- variance_dat[evar>quantile(evar,as.numeric(input$var_gene)),]
            # mostvariantgenes <- rownames(mostVariant)
            # 
            # pc_d <- int_dat() %>% select(all_of(mostvariantgenes))
            # 
            # pc_d
            
            evar <- apply(int_dat(),2,var)
            mostvariantgenes <- names(evar)[evar > quantile(evar, 1- as.numeric(input$var_gene)/100)]
            return(
              select(int_dat(), all_of(mostvariantgenes))
            )
            
            
          } else {
            
            int_dat()
         
          }
          
        }else {
          
          int_dat()
          
        }
        
        
        
      })
      
      ### extracting gene names 
      
      gene_sel <- reactive({Xproj$a() %>% select(where(is.numeric))})
      
      # observe({return(gene_sel())})
      
      gene_op <- reactive({colnames(gene_sel())})
      
      # observe({return(gene_op())})
      
      
      
      
      ## classification function
      # if input is a gene, high low classification is performed, if input is a clinic, clinic class. is performed
      

      gene_median <- reactive({median(pca_df()[[input$gene]], na.rm = T)})

      # df_class <-  reactive({
      # 
      #   mutate(pca_df(), gene_categorise = case_when(
      #     gene_mean() > pca_df()[[input$gene]] ~ "High",
      #     gene_mean() < pca_df()[[input$gene]]  ~ "Low")
      # 
      #   )
      # 
      # })


      fill <- reactive({


        ifelse(input$gene %in% gene_op(),



               return({

                 # df_class()[["gene_categorise"]]
                 
                 ifelse(pca_df()[[input$gene]] > gene_median(),
                        "High",
                        "Low")
                 
                 #'[##########################################################################################]
                 #'[##########################################################################################]
                 #'[Please ensure that the order of annotation coloring vector and the pca data is the same]
                 #'[We don't want the order to be shuffled somehow, otherwise wrongs plots will be generated]
                 #'[##########################################################################################]
                 #'[##########################################################################################]


               }), return({

                 pca_df()[[input$gene]]

               })

        )

      })

      
      
      
      ### legend title 
      
      tit_legend <- reactive({
        
       
        
        ifelse(input$gene %in% gene_op(), return({paste("Chosen Gene Factor : ", input$gene)}), 
               paste("Chosen Clinical Factor  : ", input$gene) )
        
      })
      
      
      
      ## plot title
      
      tit_plot <- reactive({
        
        
        #'[##########################################################################################]
        #'[##########################################################################################]
        #'[Please ensure these titles look good on different screen sizes (try changing the window size)]
        #'[If it doesn't look good, maybe titles can be omitted altogether]
        #'[##########################################################################################]
        #'[##########################################################################################]
        
        
        if(input$data == "All genes" ){
          
          return("PCA using all genes")
          
        }else if(input$data == "miRNA"){
          
          return("PCA using all miRNAseq data")
          
          
        }else if(input$data == "RNAseq"){
          
          return({
            
            return("PCA using all RNAseq data")
            
          })
          
        }else if(input$data == "MSigDB Gene Sets") {
          
          return(paste("PCA using genes from", input$chosen_gs))
          
        }else if(input$data == "Custom gene set") {
          
          return("PCA using genes from custom gene set")
          
        }
        
      })
      
      #prcomp
      
      pc2<- reactive(prcomp(p_dat(), center = input$center, scale. = input$scale))
      
      #scatting plot 
      
      observeEvent(input$did, {
        
        output$ex <- renderPlotly({
          
          req(input$did)
          input$did
          
          isolate({
            
            req(iv$is_valid()) 
            
            validate(
              need(input$genecor_samp_2, 'Choose at least one sample type'),
              need(input$data, 'Please select input genes'),
              need(input$gene, 'Please select a feature annotate (eg. specific genes or clinical meta data'),
              need(input$palette, 'Please select a color palette'),
            )
            
            if(input$data == "Custom gene set"){
              
              validate(need(input$pca_up, "Please upload your csv file"))
              
            }
            
            if (input$cat < 0  && input$data == "MSigDB Gene Sets" ) {
              validate("Choose a gene set")
            }
            
            if (input$chosen_gs < 0  && input$data == "MSigDB Gene Sets" ) {
              validate("Choose a pathway")
            }
            
            
            # pc2<- prcomp(p_dat(), center = input$center, scale. = input$scale)
            
            pca <- fviz_pca_ind(pc2(), geom.ind = "point", pointshape = 19,
                         fill.ind = fill(),
                         col.ind = "black",
                         palette = input$palette,
                         addEllipses = T,
                         ellipse.level=0.95,
                         label = "var",
                         col.var = "black",
                         repel = TRUE,
                         legend.title = tit_legend(),
                         geom = c("text","point"),
                         stroke = 0.1
                         
                         
            ) +
              ggtitle(tit_plot()) +
              labs(x = "PC1", y = "PC2")+
              theme(plot.title = element_text(hjust = 0.5))
            
            
           

          })
          
           
            
          
          
        })
      })
      
      
    }
    
  )
}