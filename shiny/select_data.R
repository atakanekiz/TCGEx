library(shiny)
library(shinydashboard)
library(ggplot2)
library(ggpubr)
library(data.table)
library(plotly)
library(dplyr)
# library(zstdlite)
library(shinyWidgets)
# library(shinycustomload)
library(tools)
library(readxl)

#ui

options(shiny.maxRequestSize = 250*1024^2)

select_data_ui <- function(id) {
  
  ns <- NS(id)
  
  fluidPage(
    
    add_busy_spinner(
      spin = "cube-grid",
      position = "top-right",
      color = "#01303f",
      margins = c(300, 500),
      height = "60px",
      width = "60px"),
    
    useShinyjs(),
    
    conditionalPanel(
      condition = "input.select_data_upload_user_data == false",
      ns = ns,
      pickerInput(inputId = ns("proj"), 
                  # "Project List",
                  
                  label = tags$span(
                    "Project List",
                    tags$i(
                      class = "glyphicon glyphicon-info-sign",
                      style = "color:#0072B2;",
                      title = "Please note that all TCGA data have been normalized using the CPM method. However, the normalization method for immunotherapy data is specified in the options. Be cautious when merging and analyzing datasets with different normalization methods, as this may introduce biases into your analysis results."
                    )),
                  
                  choices = list("TCGA Datasets" = 
                                   c("ACC-Adrenocortical carcinoma" = "ACC",
                                     "BLCA-Bladder Urothelial Carcinoma" = "BLCA", 
                                     "BRCA-Breast invasive carcinoma" = "BRCA", 
                                     "CESC-Cervical squamous cell carcinoma and endocervical adenocarcinoma" = "CESC", 
                                     "CHOL-Cholangiocarcinoma" = "CHOL", 
                                     "COAD-Colon adenocarcinoma" = "COAD" ,
                                     "DLBC-Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "DLBC", 
                                     "ESCA-Esophageal carcinoma" = "ESCA", 
                                     "GBM-Glioblastoma multiforme" ="GBM" , 
                                     "HNSC-Head and Neck squamous cell carcinoma" = "HNSC", 
                                     "KICH-Kidney Chromophobe" = "KICH", 
                                     "KIRC-Kidney renal clear cell carcinoma" = "KIRC" ,
                                     "KIRP-Kidney renal papillary cell carcinoma" = "KIRP",
                                     "LAML-Acute Myeloid Leukemia" = "LAML", 
                                     "LGG-Brain Lower Grade Glioma" = "LGG" , 
                                     "LIHC-Liver hepatocellular carcinoma" = "LIHC", 
                                     "LUAD-Lung adenocarcinoma" = "LUAD",
                                     "LUSC-Lung squamous cell carcinoma" = "LUSC", 
                                     "MESO-Mesothelioma" ="MESO",
                                     "OV-Ovarian serous cystadenocarcinoma" ="OV",
                                     "PAAD-	Pancreatic adenocarcinoma" = "PAAD",
                                     "PCPG-Pheochromocytoma and Paraganglioma" = "PCPG",
                                     "PRAD-Prostate adenocarcinoma" = "PRAD",
                                     "READ-Rectum adenocarcinoma" = "READ",
                                     "SARC-Sarcoma" = "SARC",
                                     "SKCM-Skin Cutaneous Melanoma" = "SKCM",
                                     "STAD-Stomach adenocarcinoma" = "STAD", 
                                     "TGCT-Testicular Germ Cell Tumors" = "TGCT", 
                                     "THCA-Thyroid carcinoma" = "THCA", 
                                     "THYM-Thymoma" = "THYM", 
                                     "UCEC-Uterine Corpus Endometrial Carcinoma" = "UCEC", 
                                     "UCS-Uterine Carcinosarcoma" = "UCS", 
                                     "UVM-Uveal Melanoma" = "UVM"),
                                 #'[#### THIS ADDED]
                                 "Immune checkpoint inhibition studies (CRI-iAtlas, quartile normalized)" =
                                   c("SKCM-Gide_Cell_2019 (PD-1, PD-1+CTLA-4)" =  "Gide_Cell_2019", 
                                     "SKCM-Hugo_Cell_2016 (PD-1)" = "HugoLo_IPRES_2016",
                                     "SKCM-Liu_NatMed_2019 (PD-1)" = "Liu_NatMed_2019", 
                                     "SKCM-Riaz_Nivolumab_2017 (PD-1)"  = "Riaz_Nivolumab_2017",
                                     "SKCM-VanAllen_Science_2015 (CTLA-4)" = "VanAllen_antiCTLA4_2015", 
                                     "KIRC-Miao_Science_2018 (PD-1, PD-1+CTLA-4, PD-L1)" = "Miao_Science_2018", 
                                     "KIRC-Choueiri_CCR_2016 (PD-1)" = "Choueiri_CCR_2016", 
                                     "KIRC-McDermott_NatMed_2018 (PD-L1)" = "IMmotion150", 
                                     "BLCA-Balar_Lancet_2017 (PD-L1)" = "IMVigor210", 
                                     "STAD-Kim_NatMed_2018 (PD-1)" =  "Kim_NatMed_2018", 
                                     "GBM-Cloughesy_NatMed_2019 (PD-1)" = "Prins_GBM_2019",
                                     "GBM-Zhao_NatMed_2019 (PD-1)" = "Zhao_NatMed_2019",
                                     "PAN-CANCER-CRI (PD-1, CTLA-4, PD-L1)"="all_icb_tcgex_ICI_TRT"),
                                 #'[#### THIS ADDED]
                                 "Immune checkpoint inhibition studies (cBioPortal)" = 
                                   c("SKCM-Van_Allen_2015 (CTLA-4, FPKM Normalized)"="skcm_dfci_2015_tcgex",
                                     "SKCM-Snyder_2014 (CTLA-4, RPKM Normalized)"="msk_2014_tcgex",
                                     "SKCM-Hugo_2016 (PD-1, RPKM Normalized)"="mel_ucla_2016_tcgex",
                                     "ACRM-Liang_2017 (CTLA-4, PD-1, FPKM Normalized)"="mel_tsam_liang_2017_tcgex",
                                     "SKCM-Liu_2019 (PD-1, TPM Normalized)"="mel_dfci_2019_tcgex",
                                     "SKCM-Riaz_2017 (CTLA-4, FPKM Normalized)"="riaz_2017_tcgex",
                                     "SKCM-Gide_2019 (CTLA-4+PD-1, CPM Normalized)"="dual_gide_tcgex",
                                     "SKCM-Gide_2019 (PD-1, CPM Normalized)"="pd1_gide_tcgex")
                                 
                  ),
                  #selectize = T,
                  # options = list('actions-box' = TRUE), #build buttons for collective selection
                  multiple = TRUE)
    ),
    materialSwitch(inputId = ns("select_data_upload_user_data"),
                   label = "Upload your own dataset",
                   status = "info",
                   value = FALSE),
    conditionalPanel(
      condition = "input.select_data_upload_user_data == true",
      ns = ns,
      fileInput(ns("file"), 
                label = tags$span(
                  "You Can Upload Your RDS or Excel Data Less Than 250 MB", 
                  tags$i(
                    class = "glyphicon glyphicon-info-sign", 
                    style = "color:#0072B2;",
                    title = "Please check sample data for your analysis. The rds/xlsx/xls file should contain gene and category names. Clinical datas should contain 'meta.' before column names, 'meta.gender' etc."
                  ), tags$br(),
                  a(href="sample_data_for_loading.rds", "Sample Input File", download=NA, target="_blank")),
                accept = c(".rds", ".xlsx", ".xls"),
                multiple = FALSE),
      
      materialSwitch(inputId = ns("nrm_dat"),
                     label = tags$span(
                       "Normalize data", 
                       tags$i(
                         class = "glyphicon glyphicon-info-sign", 
                         style = "color:#0072B2;",
                         title = "Data normalization is a process used in data preprocessing to standardize the range of independent variables or features of a dataset. If you have raw data, you can click this button and normalize your data by using log(CPM+1) normalization method"
                       )),
                     status = "info",
                     value = FALSE)
      
      #'[#### THIS COMMENTED OUT]
      # materialSwitch(inputId = ns("flt_dat"),
      #                label = tags$span(
      #                  "Filter  data", 
      #                  tags$i(
      #                    class = "glyphicon glyphicon-info-sign", 
      #                    style = "color:#0072B2;",
      #                    title = "Genes containing a percentage of 'NA' and 'zero count' greater than that specified by the slider will be filtered out and excluded from the data."
      #                  )),
      #                status = "info",
      #                value = FALSE),
      # 
      # 
      # conditionalPanel(
      #   ns=ns,
      #   condition = "input.flt_dat == true",
      #   sliderInput(inputId = ns("filter_percentage"),
      #               label = "Filtering (Keep genes that are expressed in at least n% of the samples)",
      #               min = 0,
      #               max = 100, step = 5,
      #               value = 25))
      #'[#### THIS COMMENTED OUT]
      
    ),
    
    #'[#### THIS MOVED HERE]
    materialSwitch(inputId = ns("flt_dat"),
                   label = tags$span(
                     "Filter  data", 
                     tags$i(
                       class = "glyphicon glyphicon-info-sign", 
                       style = "color:#0072B2;",
                       title = "Genes containing a percentage of 'NA' and 'zero count' greater than that specified by the slider will be filtered out and excluded from the data."
                     )),
                   status = "info",
                   value = FALSE),
    
    
    conditionalPanel(
      ns=ns,
      condition = "input.flt_dat == true",
      sliderInput(inputId = ns("filter_percentage"),
                  label = "Filtering (Keep genes that are expressed in at least n% of the samples)",
                  min = 0,
                  max = 100, step = 5,
                  value = 25)),
    #'[#### THIS MOVED HERE]
    #'
    br(),
    
    actionBttn(inputId = ns("run"), 
               label = "Load Data",
               style = "unite",
               block = FALSE,
               color = "primary"),
    br(),
    br(),
    
    textOutput(ns("fileInfos2")),
    br(),
    textOutput(ns("fileInfos3")),
    br(),
    dataTableOutput(ns("fileInfos4")),
    textOutput(ns("fileInfos5")),
    
    fluidRow(width = 12,
             column(6,
                    # withLoader( plotlyOutput(outputId =   ns("patient_hist")),type = "html",loader = "dnaspin")
                    plotlyOutput(outputId =   ns("patient_hist")), br(),
                    plotlyOutput(outputId =   ns("gender_hist"))
             ),
             column(6,
                    plotlyOutput(outputId =   ns("definition_hist")), br(),
                    plotlyOutput(outputId =   ns("age_hist"))
             ))
  )}

#server

select_data_server <- function(id,Xproj){
  moduleServer(id, function(input, output, session) {
    
    
  #'[#### THIS ADDED]
    all_projects <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM",  "Choueiri_CCR_2016", "Gide_Cell_2019","HugoLo_IPRES_2016","IMmotion150","IMVigor210","Kim_NatMed_2018","Liu_NatMed_2019","Miao_Science_2018","Prins_GBM_2019","Riaz_Nivolumab_2017", "VanAllen_antiCTLA4_2015", "Zhao_NatMed_2019","skcm_dfci_2015_tcgex","msk_2014_tcgex","mel_ucla_2016_tcgex","mel_tsam_liang_2017_tcgex","mel_dfci_2019_tcgex","all_icb_tcgex_ICI_TRT","riaz_2017_tcgex","dual_gide_tcgex","pd1_gide_tcgex")
    #'[#### THIS ADDED]
    
    file_type <- reactive({file_ext(input$file$name)})
    
    Xproj$fileInfost <- reactive({input$file$name})
    
    observeEvent(input$select_data_upload_user_data, {
      
      validate(need(input$file, ""))
      
      shinyjs::reset("proj")
      shinyjs::reset("file")
      shinyjs::reset(Xproj$fileInfost())
      
      shinyjs::reset("fileInfos2")
      shinyjs::reset("fileInfos3")
      shinyjs::reset("fileInfos4")
      shinyjs::reset("fileInfos5")
      
      shinyjs::reset("patient_hist")
      shinyjs::reset("gender_hist")
      shinyjs::reset("definition_hist")
      shinyjs::reset("age_hist")
      
      output$fileInfos2 <- renderText({NULL})
      output$fileInfos3 <- renderText({NULL})
      output$fileInfos4 <- renderDataTable({NULL})
      output$fileInfos5 <- renderText({NULL})
    })
    
    observeEvent(input$select_data_upload_user_data, {
      
      shinyjs::reset("proj")
      shinyjs::reset("file")
      
      shinyjs::reset("patient_hist")
      shinyjs::reset("gender_hist")
      shinyjs::reset("definition_hist")
      shinyjs::reset("age_hist")
      
      shinyjs::reset("fileInfos2")
      shinyjs::reset("fileInfos3")
      shinyjs::reset("fileInfos4")
      shinyjs::reset("fileInfos5")
      
      output$fileInfos2 <- renderText({NULL})
      output$fileInfos3 <- renderText({NULL})
      output$fileInfos4 <- renderDataTable({NULL})
      output$fileInfos5 <- renderText({NULL})
    })
    
    
    observeEvent(input$run,{
      
      validate(need(input$file, ""))
      
      if (!is.null(input$file) && is.null(input$proj) ){
        
        xdata <- input$file
        ydata <- reactive({as.data.table(readRDS(paste0(input$file$datapath)))})
        
        if (file_type() %in% c("rds","RDS","Rds")) {
          
          file_size <- round(xdata$size / (1024^2), 2)
          num_rows <- nrow(Xproj$a())
          num_cols <- ncol(Xproj$a())
          
          file_info <- paste("Size of Uploaded RDS File:", file_size, "MB")
          file_info <- paste(file_info, "Number of Rows:", num_rows)
          file_info <- paste(file_info, "Number of Columns:", num_cols)
          
          output$fileInfos2 <- renderText({
            validate(need(input$run, ""))
            file_info
          })
          
          output$fileInfos3 <- renderText({
            validate(need(input$run, ""))
            "You Can Start Analyzing Your Data by Switching the Tab"
          })
          
          output$fileInfos4 <- renderDataTable({
            validate(need(input$run, ""))
            Xproj$a()[1:10,1:15]
            
          })
          output$fileInfos5 <- renderText({
            validate(need(input$run, ""))
            "The first 10 rows and first 15 columns of your data are shown..."
          })
        }
        
        else if ((file_type() %in% c("xlsx", "xls"))) {
          
          xdata <- input$file
          ydata <- reactive({as.data.table(read_excel(paste0(input$file$datapath)))})
          
          file_size <- round(xdata$size / (1024^2), 2)
          num_rows <- nrow(Xproj$a())
          num_cols <- ncol(Xproj$a())
          
          file_info <- paste("Size of Uploaded EXCEL File:", file_size, "MB")
          file_info <- paste(file_info, "Number of Rows:", num_rows)
          file_info <- paste(file_info, "Number of Columns:", num_cols)
          
          output$fileInfos2 <- renderText({
            validate(need(input$run, ""))
            file_info
          })
          
          output$fileInfos3 <- renderText({
            validate(need(input$run, ""))
            "You Can Start Analyzing Your Data by Switching the Tab"
          })
          
          output$fileInfos4 <- renderDataTable({
            validate(need(input$run, ""))
            Xproj$a()[1:10,1:15]
            
          })
          output$fileInfos5 <- renderText({
            validate(need(input$run, ""))
            "The first 10 rows and first 15 columns of your data are shown..."
          })
        }
        
      }
      
    })
    
    Xproj$cancer_length <- reactive({length(as.vector(input$proj))}) ## a reactive that created for other modules to use the length information for several cancers(Cagatay)
    
    Xproj$a<- eventReactive(input$run, {
      
      req (input$run)
      
      # browser()
      
      validate(need(input$run, ""))
      
      if(!is.null(input$proj))  {
        
        if (Xproj$cancer_length() == 1 ){
          
          req (input$run)
          
          # readRDS(paste0("projects/", input$proj, ".rds"))
          
          #'[ ####### THIS CHANGED]
         preloaded_data <-  readRDS(paste0("projects/", input$proj, ".rds"))
          
         if(input$flt_dat == TRUE){
           
           filter_percentage <- reactive ({1-(input$filter_percentage / 100)})
           
           # filter_percentage <- reactive ({input$filter_percentage / 100})
           
           df_num = preloaded_data %>% select(where(is.numeric))
           
           df_nonnum = preloaded_data %>% select(-where(is.numeric))
           
           # df_nongene = df_num %>% select(starts_with("meta."))
           # 
           # df_gene = df_num %>% select(-starts_with("meta."))
           
           df_nongene = df_num %>% 
             select(starts_with("meta."),starts_with("meta_"))
           
           df_gene <- df_num %>%
             select(-starts_with("meta."), -starts_with("meta_"))
           
           # browser()
           
           na_zero_percent <- apply(df_gene, 2, function(x) mean(is.na(x) | x == 0))
           
           selected_columns <- names(df_gene)[na_zero_percent < filter_percentage()]
           
           df_gene <- df_gene[, ..selected_columns]
           
           preloaded_data = cbind(df_gene, df_nonnum,df_nongene)
           
           preloaded_data
          
         }
         
         preloaded_data
         #'[ ####### THIS CHANGED]
          
          # compressed<-readRDS(paste0("projects/", input$proj, ".rds"))
          #
          # zstd_unserialize(compressed
          
        }
        
        #'[ ##### THIS ADDED]
        # else if (!(Xproj$cancer_length() == 1) && any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj)) {
        
        else if (!(Xproj$cancer_length() == 1) && any(all_projects %in% input$proj)) {
          
          req (input$run)
          
          #'[ ##### THIS ADDED]
          
          validate(need(input$run, ""))
          
          # browser()
          
          dflist <- list()
          
          for (i in 1:length(input$proj)) {
            
            dflist[[i]] <- readRDS(paste0("projects/", input$proj[[i]], ".rds"))
            
            # compressed_data<-readRDS(paste0("projects/", input$proj[[i]], ".rds"))
            #
            # dflist[[i]] <- zstd_unserialize(compressed_data)
          }
          
          same_columns <- Reduce(function(x, y){intersect(x, names(y))}, dflist, init = names(dflist[[1]]))
          
          for (i in 1:length(input$proj)) {
            
            dflist[[i]] <- dflist[[i]][, ..same_columns] ## we added .. this to same columns, as syntax
            
          }
          
          big_data <- do.call(rbind, dflist)
          
          #'[ ####### THIS CHANGED. CONSIDER MOVING THIS INTO PREVIOUS CHUNK TO APPLY FILTERING TO EACH DATASET SEPARATELY]
          
          if(input$flt_dat == TRUE){
            
            filter_percentage <- reactive ({1-(input$filter_percentage / 100)})
            
            # filter_percentage <- reactive ({input$filter_percentage / 100})
            
            df_num = big_data %>% select(where(is.numeric))
            
            df_nonnum = big_data %>% select(-where(is.numeric))
            
            # df_nongene = df_num %>% select(starts_with("meta."))
            # 
            # df_gene = df_num %>% select(-starts_with("meta."))
            
            df_nongene = df_num %>% 
              select(starts_with("meta."),starts_with("meta_"))
            
            df_gene <- df_num %>%
              select(-starts_with("meta."), -starts_with("meta_"))
            
            # browser()
            
            na_zero_percent <- apply(df_gene, 2, function(x) mean(is.na(x) | x == 0))
            
            selected_columns <- names(df_gene)[na_zero_percent < filter_percentage()]
            
            df_gene <- df_gene[, ..selected_columns]
            
            big_data <-  cbind(df_gene, df_nonnum,df_nongene)
            
            }
          #'[ ####### THIS CHANGED]
          
          return(big_data)
          
          
        }}
      
      else if (!is.null(input$file)) {
        
        # filter_per<-reactive({input$filter_percentage})
        
        # browser()
        
        validate(need(input$file, ""))
        
        
        if (file_type() %in% c("rds","RDS","Rds","rDs","rDS","rdS","RDs","RdS")){
          
          
          uploaded_data <- as.data.table(readRDS(paste0(input$file$datapath)))
          
          colnames(uploaded_data)<- gsub("-", ".", colnames(uploaded_data), fixed = TRUE)
          
          
          if (!"meta.definition" %in% names(uploaded_data)) {
            uploaded_data[, meta.definition := "All Samples"]
          }
          
          if(input$nrm_dat == TRUE){
            
            # browser()
            
            df_num = uploaded_data %>% select(where(is.numeric))
            df_nonnum = uploaded_data %>% select(-where(is.numeric))
            
            df_nongene = df_num %>% 
              select(starts_with("meta."),starts_with("meta_"))
            
            df_gene <- df_num %>%
              select(-starts_with("meta."), -starts_with("meta_"))
            
            df_gene <- df_gene[, apply(df_gene,2,var)!=0, with=F]
            
            df_gene = log(cpm(df_gene, prior.count=0)+1, base=10)
            
            uploaded_data = cbind(df_gene, df_nonnum,df_nongene)
            
            uploaded_data
            
          }
          
          
          if(input$flt_dat == TRUE){
            
            filter_percentage <- reactive ({1-(input$filter_percentage / 100)})
            
            # filter_percentage <- reactive ({input$filter_percentage / 100})
            
            df_num = uploaded_data %>% select(where(is.numeric))
            
            df_nonnum = uploaded_data %>% select(-where(is.numeric))
            
            # df_nongene = df_num %>% select(starts_with("meta."))
            # 
            # df_gene = df_num %>% select(-starts_with("meta."))
            
            df_nongene = df_num %>% 
              select(starts_with("meta."),starts_with("meta_"))
            
            df_gene <- df_num %>%
              select(-starts_with("meta."), -starts_with("meta_"))
            
            # browser()
            
            na_zero_percent <- apply(df_gene, 2, function(x) mean(is.na(x) | x == 0))
            
            selected_columns <- names(df_gene)[na_zero_percent < filter_percentage()]
            
            df_gene <- df_gene[, ..selected_columns]
            
            uploaded_data = cbind(df_gene, df_nonnum,df_nongene)
            
            uploaded_data
            
          }
          
          
          return(uploaded_data)
          
          # uploaded_data
          
        }
        
        else if (file_type() %in% c("xlsx", "xls")){
          
          uploaded_data_xl<- as.data.table(read_excel(input$file$datapath))
          
          if (!"meta.definition" %in% names(uploaded_data_xl)) {
            uploaded_data_xl[, meta.definition := "All Samples"]
          }
          
          if(input$nrm_dat == TRUE){
            # browser()
            
            df_num = uploaded_data_xl %>% select(where(is.numeric))
            df_nonnum = uploaded_data_xl %>% select(-where(is.numeric))
            
            df_nongene = df_num %>% 
              select(starts_with("meta."),starts_with("meta_"))
            
            df_gene <- df_num %>%
              select(-starts_with("meta."), -starts_with("meta_"))
            
            df_gene <- df_gene[, apply(df_gene, 2, var) !=0, with=F]
            
            #filter out all-zero columns
            
            df_gene = log(cpm(df_gene, prior.count=0)+1, base=10)
            
            uploaded_data_xl = cbind(df_gene, df_nonnum,df_nongene)
            
            uploaded_data_xl
            
          }
          
          if(input$flt_dat == TRUE){
            
            filter_percentage <- reactive ({input$filter_percentage / 100})
            
            df_num = uploaded_data_xl %>% select(where(is.numeric))
            
            df_nonnum = uploaded_data_xl %>% select(-where(is.numeric))
            
            # df_nongene = df_num %>% select(starts_with("meta."))
            # 
            # df_gene = df_num %>% select(-starts_with("meta."))
            
            df_nongene = df_num %>% 
              select(starts_with("meta."),starts_with("meta_"))
            
            df_gene <- df_num %>%
              select(-starts_with("meta."), -starts_with("meta_"))
            
            na_zero_percent <- apply(df_gene, 2, function(x) mean(is.na(x) | x == 0))
            
            selected_columns <- names(df_gene)[na_zero_percent < filter_percentage()]
            
            df_gene <- df_gene[, ..selected_columns]
            
            uploaded_data_xl = cbind(df_gene, df_nonnum,df_nongene)
            
            uploaded_data_xl
            
          }
          
          return(uploaded_data_xl)
          
          # uploaded_data_xl
        }
        
        
        else if (Xproj$cancer_length() == 0) {
          
          Xproj$a()<- NULL
          
          shinyjs::reset("proj")
          shinyjs::reset("file")
          shinyjs::reset(Xproj$fileInfost())

          shinyjs::reset("fileInfos2")
          shinyjs::reset("fileInfos3")
          shinyjs::reset("fileInfos4")
          shinyjs::reset("fileInfos5")

          shinyjs::reset("patient_hist")
          shinyjs::reset("gender_hist")
          shinyjs::reset("definition_hist")
          shinyjs::reset("age_hist")

          output$fileInfos2 <- renderText({NULL})
          output$fileInfos3 <- renderText({NULL})
          output$fileInfos4 <- renderDataTable({NULL})
          output$fileInfos5 <- renderText({NULL})
          
        }
        
        
      }
      
      
    })
    
    output$gender_hist<- renderPlotly({
      
      req(input$run)
      
      validate(need(input$run, ""))
      
      #'[#### THIS ADDED]
      # if (any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj) ){
        
        if (any(all_projects %in% input$proj) ){
          #'[#### THIS ADDED]
          
        #gghistogram(Xproj$a(), "meta.gender", stat="count", legend="none",
        #font.x=18, font.y=18, font.tickslab = 18,
        #fill="meta.gender", palette = c("skyblue", "coral", "gold"))
        
        fig_gender <- plot_ly(Xproj$a(),  labels = ~meta.gender, type = 'pie',
                              marker = list(color = viridis::viridis_pal(begin = 0.2, end = 0.8)(4)))
        
        fig_gender <- fig_gender %>% layout(title = 'Gender Statistics For Chosen Data',
                                            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        
        fig_gender
      }
      
      else {
        return(NULL)
      }
      
    })
    
    output$patient_hist <- renderPlotly({
      
      req(input$run)
      
      #'[#### THIS CHANGED]
      # if (any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj) ){
        
        if (any(all_projects %in% input$proj) ){
        
        #'[#### THIS CHANGED]
        
        validate(need(input$run, "Load TCGA project to see descriptive statistics"))
        
        fig_patient <- plot_ly(Xproj$a(),  labels = ~meta.project_id, type = 'pie',
                               marker = list(color = viridis::viridis_pal(begin = 0.2, end = 0.8)(4)))
        
        fig_patient <- fig_patient %>% layout(title = 'Total Patient Number For Chosen Data',
                                              xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                              yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        
        fig_patient
      }
      else {
        return(NULL)
      }
    })
    
    
    output$definition_hist <- renderPlotly({
      
      req(input$run)
      
      #'[THIS CHANGED]
      # if (any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj) ){
      
      if (any(all_projects %in% input$proj) ){
        
        #'[THIS CHANGED]
        
        validate(need(input$run, ""))
        
        #gghistogram(Xproj$a(), "meta.definition", stat="count",
        #font.x=18, font.y=18, font.tickslab = 18,
        #fill="meta.definition", legend="none") +
        #rotate_x_text(45)
        
        Xproj$a() %>%
          plot_ly(x = ~meta.definition,
                  textangle = 45,
                  marker = list(color = viridis::viridis_pal(option = "C", direction = -1)(10)))%>%
          add_histogram()
      }
      else {
        return(NULL)
      }
    })
    
    output$age_hist <- renderPlotly({
      
      req(input$run)
      
      #'[THIS CHANGED]
      # if (any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj) ){
      if (any(all_projects %in% input$proj) ){
        
        #'[THIS CHANGED]
        validate(need(input$run, ""))
        
        ##gghistogram(Xproj$a(), "meta.age_at_diagnosis", fill="gold", bins=50,
        #font.x=18, font.y=18,  font.tickslab = 18)
        #Xproj$a() %>%
        # plot_ly(x = ~meta.age_at_diagnosis,
        #        textangle = 45,
        #       marker = list(color = viridis::viridis_pal(option = "C", direction = -1)(4)))%>%
        #add_histogram()
        fig_age_hist <- plot_ly(Xproj$a(), x = ~ meta.age_at_diagnosis, type = 'histogram',
                                marker = list(color = viridis::viridis_pal(option = "C", direction = -1)(40)))
        fig_age_hist
      }
      else {
        return(NULL)
      } 
    })
    
    
  }
  )}

