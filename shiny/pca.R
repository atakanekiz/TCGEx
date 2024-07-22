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
library(readxl)
library(shinyBS)

pca_ui <- function(id) {
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
      
      
      
      selectizeInput(ns("genecor_samp_2"), multiple=T,
                     "Please select sample types",
                     choices=NULL,
                     options=list(placeholder = "eg. Primary solid tumor",
                                  plugins = list('restore_on_backspace'))),
      
      selectizeInput(ns("data"), "*Please select input genes", 
                     choices = NULL,
                     options=list(placeholder = "eg. RNAseq",
                                  plugins = list('restore_on_backspace'))),
      
      checkboxInput(ns("center"), "Center data", value = TRUE),
      checkboxInput(ns("scale"), "Scale variables", value = TRUE),
      
      
      ## conditional panel for gene sets
      # pathway can be chosen according to gene set input
      
      conditionalPanel(
        
        condition = 'input.data == "Custom gene set"', ns = ns,
        
        
        
        fileInput(inputId = ns("pca_up"),
                  label = tags$span(
                    "Please upload your xlsx/xls file.",
                    tags$i(
                      class = "glyphicon glyphicon-info-sign",
                      style = "color:#0072B2;",
                      title = "The xlsx/xls file should contain a single unnamed column with human gene names."
                    ), tags$br(),
                    a(href="sample_gene_input.xlsx", "Sample Input File", download=NA, target="_blank")),
                  accept = c(".xls", ".xlsx" # "text/csv", "text/comma-separated-values,text/plain", ".csv")), 
       
        
        
        
      ))),
      
      conditionalPanel(
        
        condition = 'input.data == "MSigDB Gene Sets"', ns = ns,
        
        selectizeInput(ns("pca_cat"), 
                       "Please select an MSigDB Collection",
                       # label = tags$span(
                       #   "Please select an MSigDB Collection",
                       #   tags$i(
                       #     class = "glyphicon glyphicon-info-sign",
                       #     style = "color:#0072B2;",
                       #     title = "For information about gene sets, you can visit",
                       #     tags$a(href='https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp', "here", .noWS = "after")
                       #   )),
                       label = tags$span( "Please select an MSigDB Collection", 
                                          bsButton("pc_msigdb", 
                                                   label = "", 
                                                   icon = icon("info"), 
                                                   style = "info", 
                                                   size = "extra-small")),
                       choices = c("Hallmark gene sets (H)" = "H",
                                   "Positional gene sets (C1)" = "C1",
                                   "Curated gene sets (C2)" = "C2",
                                   "Regulatory target gene sets (C3)" = "C3",
                                   "Computational gene sets (C4)" = "C4",
                                   "Ontology gene sets (C5)" = "C5" ,
                                   "Oncogenic gene sets (C6)" = "C6",
                                   "Immunologic gene sets (C7)" = "C7",
                                   "Cell type signature gene sets (C8)" = "C8"),
                       options = list(plugins = list('restore_on_backspace'))),
        bsPopover(
          id = "pc_msigdb",
          title = "",
          content = paste0(
            "For information about gene sets, you can visit ",
            a("MsigDB", href = "https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp", target="_blank")
          ),
          placement = "right",
          trigger = "focus",
          options = list(container = "body")
        ),
        conditionalPanel(condition = "input.pca_cat == 'C2'|input.msigdb_setnames_response =='C3'|
                                      input.pca_cat =='C4'|
                                      input.pca_cat =='C5'|
                                      input.pca_cat =='C7' ", ns = ns, 
                         selectizeInput(ns("pca_subcat"),
                                        "Please select a subcategory" ,
                                        choices = c(""),
                                        options = list(plugins = list('restore_on_backspace')))
        ),
        
        selectizeInput(ns("pca_pathway"), 
                       "Please select a pathway", 
                       choices = c(""),
                       options = list(plugins = list('restore_on_backspace'))
                       )),
        
      
      checkboxInput(ns("variance"), "Apply variance filtering", value = F),
      
      conditionalPanel(
        
        condition = 'input.variance', ns = ns,
        
        sliderInput(inputId = ns("var_gene"),
                    label = "Keep genes with the top n% variation",
                    min = 0, max = 100, value = 10
        )
        
      ),
      
      selectizeInput(ns("gene"), "*Please select feature to annotate", 
                     options=list(placeholder = "eg. meta.gender, gene name",
                                  plugins = list('restore_on_backspace')),
                     choices = NULL),
      
      selectizeInput(ns("palette"), "*Please select a color palette",
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
                     options=list(placeholder = "eg. Nature",
                                  plugins = list('restore_on_backspace'))
                     
      ),
      
      actionBttn(inputId = ns("did"), 
                 label = "Analyze",
                 style = "unite",
                 block = TRUE,
                 color = "primary"),
      br(),
      
      #help section UI
      
      introjsUI(),
      actionButton(ns("intro2"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
      
      textOutput(ns("filewarning_4")) ,
      
      width = 3
      
    ),
    
    mainPanel(
      
      plotlyOutput(outputId = ns("ex"),
                   width = "100%",
                   height = "600px"
                   ),
      DTOutput(outputId = ns("eigen"))
      
      
    )
  )
  
  
}





pca_server <- function(id,Xproj) {
  moduleServer(
    id,
    function(input, output, session) {
      
      output$filewarning_4 <- renderText({
        
        if (!is.null(Xproj$fileInfost())) {
          shinyalert("Warning!", "To perform this analysis using MsigDB gene sets, please ensure that your uploaded data set contains gene symbols rather than Entrez or Ensembl gene IDs. Otherwise you may receive errors.
                     
                     To perform miRNA-based analysis, remember that the miRNA columns in your data must start with 'hsa.miR.' . Example miRNA column names: 'hsa.miR.155.5p', 'hsa.miR.142.3p', 'hsa.miR.107'.") }
      })
      
      lncrnas_vector <- reactive({readRDS(paste0("genesets/", "filtered_lncrnas", ".rds"))})
      
      ## msigdb_database reading
      
      msigdb_gene_sets =  reactive({readRDS(paste0("genesets/", "msigdb_long_w_immth", ".rds"))})
      
      ## help section server
      
      ns <- session$ns
      
      steps <- reactive({
        
        data.frame(
          
          element = paste0("#", session$ns(c(NA,"genecor_samp_2 + .selectize-control" ,"data + .selectize-control", "center", "scale", "variance","gene + .selectize-control", "palette + .selectize-control"))),
          
          intro = paste(c(
            "This is the Principal Compenent Analysis (PCA) module. PCA is a dimensionality reduction method allowing multidimensional datasets to be visualized on a two dimensional graph. Principal components are derived from linear combinations of the original variables and they represent the variation in the data set. <b>NOTE:</b> Creating PCA plots can take some time, please be patient. Continue tutorial to learn more about this module.",
            "You can select sample types to include in the analysis here.",
            "Here, you can select the genes to be used in PCA. You can select <b>i)</b> a custom list of genes, <b>ii)</b> all genes (RNAseq and miRNAseq data), <b>iii)</b> miRNAs (mature miRNAs from miRNAseq), or <b>iv)</b> genes annotated in MSigDB gene sets.",
            "By default, gene expression values are centered by subtracting the mean expression value.",
            "A variable that is on a different scale from the others may dominate the variance direction. Scaling (default) gene expression values prevents this effect.",
            "Here, you can apply variance filtering to keep most highly variable genes in the analysis. Setting this value to 10, for instance, will select the genes having the top 10% highest variation in the dataset. Such filtering can speed up the analysis.",
            "You can color code the data points on the graph using gene expression values or clinical meta data. If you select a gene name here, gene expression will be categorized at the median value per sample and points will be annotated. You can also select a clinical meta data (eg. meta.gender or a specific gene) to color points accordingly.",
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
      observe({updateSelectizeInput(session, 'data', 
                                   choices = c("All genes", "miRNA", "RNAseq","lncRNA", "MSigDB Gene Sets", "Custom gene set"),
                                   selected = "",
                                   server = TRUE)})
      
      
      observe({
        req(input$pca_cat)
        
        pc_subcat = names(msigdb_gene_sets()[[input$pca_cat]])
        if (length(pc_subcat) > 1)  {
          updateSelectizeInput(session,'pca_subcat', choices = pc_subcat , server = TRUE)
        } 
      })
      

      observe(updateSelectizeInput(session, 'gene', 
                                   choices = colnames(Xproj$a()), 
                                   selected = "",
                                   server = TRUE))
      
      
      
      observe({
        if(input$pca_cat %in% c("C2","C3","C4","C5","C7")) {
          
          pca_path = names(msigdb_gene_sets()[[input$pca_cat]][[input$pca_subcat]])
          updateSelectizeInput(session,'pca_pathway', choices = pca_path , server = TRUE)
          
          
        } else {
          
          pca_path = names(msigdb_gene_sets()[[input$pca_cat]][[1]])
          updateSelectizeInput(session,'pca_pathway', choices = pca_path , server = TRUE)
        }
        
      })
      
      
      
      pca_msigdb_genes <- reactive({
        
        ms = msigdb_gene_sets()
        if(input$pca_cat %in% c("C2","C3","C4","C5","C7")) {
          
          pca_m = ms[[input$pca_cat]][[input$pca_subcat]][[input$pca_pathway]]
          
        } else {
          
          pca_m = ms[[input$pca_cat]][[1]][[input$pca_pathway]]
          
        }
        
        pca_m
        
        })
      
      
      ##NA number
      
      gene_cols <- reactive({
        
        Xproj$a() %>%
          select(!starts_with("meta.")) %>%
          select(where(is.numeric))
        
      })
        

      
      # na_number <- reactive({
      #   
      #   NA_number <- data.frame(colSums(is.na(gene_cols())))
      #   
      #   na_miRNA <- NA_number %>% filter(grepl('hsa.', rownames(NA_number)))
      #   
      #   colnames(na_miRNA)[1] <- "na_numbers"
      #   
      #   if (length(unique(na_miRNA$na_numbers)) == 1) {
      #     
      #     na_num = unique(na_miRNA$na_numbers)
      #     
      #   }else if (length(unique(na_miRNA$na_numbers)) > 1) {
      #     
      #     na_num = max(unique(na_miRNA$na_numbers))
      #     
      #   }
      #   
      #   
      # })
     

      
      
      ### DATA SELECTION PART
      
      pca_df <- reactive({
        
        # browser()
        
        
         pre_d <- Xproj$a()[meta.definition %in% input$genecor_samp_2, ]
        
        # sel_cols <- apply(select(pre_d, all_of(gene_names)), 2, sum)
        
         gene_names <- colnames(gene_cols())
       
        miRNA <- select(gene_cols(),starts_with("hsa.") )
        
        lncRNA <-  gene_cols()[,intersect(colnames(pre_d),lncrnas_vector()),with=F]
        
        RNAseq <- select(gene_cols(),-starts_with("hsa.") )
        
        if(input$data == "All genes" ){
          
            pre_d <- drop_na(pre_d, colnames(miRNA))%>%
              select(!(any_of(colnames(RNAseq))& where(~ any(is.na(.x)))))
          
          
        }else if(input$data == "miRNA"){
          
          pre_d <- drop_na(pre_d, colnames(miRNA))
          
          
        }else if(input$data == "lncRNA"){
          
          pre_d <- drop_na(pre_d,  colnames(lncRNA))
          
          
        }else if(input$data == "RNAseq"){
          
          pre_d <- drop_na(pre_d, colnames(RNAseq))
          
        }
        
        
        pre_d
        
       

        
        
        
      })
      
      all_genes <- reactive({
        
        # browser()
        
        req(pca_df())
        
        pc_df <-  pca_df() %>% 
          select(!starts_with("meta.")) %>%
          select(where(is.numeric))
        
        pc_df <- remove_constant(pc_df)
        
      })
      
      cre_gene <- reactive({
        
        req(input$pca_up)
        
        c_gene <- as.data.frame(read_excel(input$pca_up$datapath, sheet = 1, col_names = F))  
        
        validate(
          need(
            {if(length(colnames(c_gene))== 1) TRUE else FALSE},
            "This file contains more than 1 columns.Please ensure that your file contains a single column with human gene names"))
        
        colnames(c_gene)[1] <- "genes" 
        
        validate(
          need(
            {if(class(c_gene$genes) == "character")TRUE else FALSE}, 
            "The column must contain only character input.Please ensure that the column contains only human gene names"))
        
        
        # gene_var <- Xproj$a() %>%
        #     select(!starts_with("meta.")) %>%
        #     select(where(is.numeric))
        #   
        
        
        
        # NA_number <- data.frame(na_mum = colSums(is.na(gene_cols())))
        # 
        # na_miRNA <- NA_number %>% filter(grepl('hsa.', rownames(NA_number)))
        # 
        # colnames(na_miRNA)[1] <- "na_numbers"
        # 
        # if (length(unique(na_miRNA$na_numbers)) == 1) {
        #   
        #   na_num = unique(na_miRNA$na_numbers)
        #   
        # }else if (length(unique(na_miRNA$na_numbers)) > 1) {
        #   
        #   na_num = max(unique(na_miRNA$na_numbers))
        #   
        # }
        
        
        
          c_gene <- filter(c_gene, !grepl("hsa.", genes))

          showNotification("Genes which starts with 'hsa.' were removed. This tool will allow to use miRNA's coming soon",
                           duration = 7,
                           closeButton = T,
                           type= "message")
        
        
        # if(((NA_number/nrow(gene_var))*100) > 15) {
        #   
        #   c_gene <- filter(c_gene, !grepl("hsa.", genes))
        #   
        #   showNotification("genes which starts with 'hsa.' were removed because there is no enough miRNA data ", 
        #                    duration = 7, 
        #                    closeButton = T, 
        #                    type= "message")
        #   
        # }else {
        #   
        #   c_gene 
        #   
        # }
        
        # ifelse(((NA_number/nrow(gene_var))*100) > 25, return(c_gene <- filter(c_gene, !grepl("hsa.", genes))), c_gene )
        
        # ifelse(((NA_number/nrow(gene_var))*100) > 25,  return({showNotification("genes which starts with 'hsa.' were removed because there is no enough miRNA data ", 
        #                                                                 duration = 7, 
        #                                                                 closeButton = T, 
        #                                                                 type= "message")}), "")
        
        # if((NA_number/(nrow(gene_var)*100)) > 15) {
        # 
        #   showNotification("genes which starts with 'hsa.' were removed because there is no enough miRNA data ",
        #                    duration = 7,
        #                    closeButton = T,
        #                    type= "message")
        # 
        # }

        c_gene
        
      })
      
      int_dat <- reactive({
        
        
        # req(all_genes())
        
        if(input$data == "All genes" ){
          
           df_pc <- all_genes()
           
           
          
        }else if(input$data == "miRNA"){
          
          
            
            df_pc <- select(all_genes(), starts_with("hsa."))
            
            
            
        
          
          
        }else if(input$data == "lncRNA"){
          
          
          
          df_pc <- all_genes()[,intersect(lncrnas_vector(), colnames(all_genes())),with=F]
          
          
          
          
          
          
        }else if(input$data == "RNAseq"){
          
          
            
            df_pc <- select(all_genes(), !starts_with("hsa."))
            
            
            
        
          
        }else if(input$data == "MSigDB Gene Sets"){
          
          df_pc <- all_genes()[,intersect(pca_msigdb_genes(), colnames(all_genes())),with=F]
          
          
          
        }else if(input$data == "Custom gene set"){
          # input$data == "Custom gene set" & input$did > 0
          
          
          df_pc <-  all_genes()[,intersect(cre_gene()[["genes"]], colnames(all_genes())),with=F]
          
          
        }
        
        df_pc
        
        
      })
      
      

      
      
      ## keep highly variant genes
      
      p_dat <- reactive({
        
        # browser()
        
        if(input$variance){
          
          if (input$var_gene != 100) {
            

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
      
      gene_op <- reactive({colnames(gene_sel())})
      
      
      ## classification function
      # if input is a gene, high low classification is performed, if input is a clinic, clinic class. is performed
      
      
      gene_median <- reactive({median(pca_df()[[input$gene]], na.rm = T)})
      
      
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
                 #'[for PCA ,constant columns are removed, some gene columns for classification might be constant, so initial data must be proper because there is no change in row number]
                 #'[##########################################################################################]
                 #'[##########################################################################################]
                 
                 
               }), return({
                 
                 pca_df()[[input$gene]]
                 
               })
               
        )
        
      })
      
      
      
      
      ### legend title 
      
      tit_legend <- reactive({
        
        
        
        ifelse(input$gene %in% gene_op(), return({paste("Chosen Gene : ", input$gene)}), 
               paste("Chosen Clinical Variable  : ", input$gene) )
        
      })
      
      
      
      ## plot title
      
      tit_plot <- reactive({
        
        
        #'[##########################################################################################]
        #'[##########################################################################################]
        #'[Please ensure these titles look good on different screen sizes (try changing the window size)]
        #'[If it doesn't look good, maybe titles can be omitted altogether]
        #'[I updated plot visualization please check]
        #'[##########################################################################################]
        #'[##########################################################################################]
        
        
        if(input$data == "All genes" ){
          
          return("PCA using all genes")
          
        }else if(input$data == "miRNA"){
          
          return("PCA using all miRNAseq data")
          
          
        }else if(input$data == "lncRNA"){
          
          return("PCA using all lncRNA data")
          
          
        }else if(input$data == "RNAseq"){
          
          return({
            
            return("PCA using all RNAseq data")
            
          })
          
        }else if(input$data == "MSigDB Gene Sets") {
          
          return(paste("PCA using genes from", input$pca_pathway))
          
        }else if(input$data == "Custom gene set") {
          
          return("PCA using genes from custom gene set")
          
        }
        
      })
      
      #prcomp
      
      pc2<- reactive({
        
        prcomp(p_dat(), center = input$center, scale. = input$scale)
        
        })
      
      
      eigentable <- reactive({
        
        eigenval = pc2()$rotation
        eigenval = as.data.frame(eigenval)
        eigenval$genes <- row.names(eigenval)
        eigenval_pc1 = select(eigenval, c("PC1", "genes"))
        eigenval_pc2 = select(eigenval, c("PC2","genes"))
        high_pc1 = eigenval_pc1 %>% arrange(desc(PC1)) %>% head(10)
        low_pc1 = eigenval_pc1 %>% arrange(desc(PC1)) %>% tail(10)
        pc1_eigen = rbind(high_pc1,low_pc1)
        high_pc2 = eigenval_pc2 %>% arrange(desc(PC2)) %>% head(10)
        low_pc2 = eigenval_pc2 %>% arrange(desc(PC2)) %>% tail(10)
        pc2_eigen = rbind( high_pc2,low_pc2)
        
        eig = cbind(pc1_eigen,pc2_eigen)
        eig
        
        
      })
      
      
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
              need(input$gene, 'Please select a feature to annotate in the plot (eg. specific genes or clinical metadata)'),
              need(input$palette, 'Please select a color palette'),
            )
            
            
            if(input$data == "Custom gene set"){
              
              validate(need(input$pca_up, "Please upload your xlsx/xls file"))
              
            }
            
            if (input$pca_cat < 0  && input$data == "MSigDB Gene Sets" ) {
              validate("Choose a gene set")
            }
            
            req(int_dat())

            validate(

              need(
                {if(nrow(int_dat()) < 2) FALSE else TRUE},
                "There isn't enough data related to your choice , so please change input genes, sample type or try another project"))
            
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
              theme(plot.title = element_text(size=11, hjust = 0.5))
            
            
            
            
          })
          
          
          
          
          
        })
      })
      
      observeEvent(input$did, {
        
        output$eigen <- renderDT({
          
          req(input$did)
          input$did
          
          isolate({
            
            req(iv$is_valid()) 
            
            validate(
              need(input$genecor_samp_2, ''),
              need(input$data, ''),
              need(input$gene, ''),
              need(input$palette, ''),
            )
            
            
            if(input$data == "Custom gene set"){
              
              validate(need(input$pca_up, ""))
              
            }
            
            if (input$pca_cat < 0  && input$data == "" ) {
              validate("")
            }
            
            req(int_dat())
            
            validate(
              
              need(
                {if(nrow(int_dat()) < 2) FALSE else TRUE},
                ""))
            
              eigentable()
            
            
            
            
          })
          
          
          
          
          
        })
      })
      
      
      
    }
    
  )
}