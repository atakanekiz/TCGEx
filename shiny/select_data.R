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
                  "Project List",
                  choices = c("ACC-Adrenocortical carcinoma" = "ACC","BLCA-Bladder Urothelial Carcinoma" = "BLCA", "BRCA-Breast invasive carcinoma" = "BRCA", "CESC-Cervical squamous cell carcinoma and endocervical adenocarcinoma" = "CESC", "CHOL-Cholangiocarcinoma" = "CHOL", "COAD-Colon adenocarcinoma" = "COAD" ,"DLBC-Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "DLBC", "ESCA-Esophageal carcinoma" = "ESCA", "GBM-Glioblastoma multiforme" ="GBM" , "HNSC-Head and Neck squamous cell carcinoma" = "HNSC", "KICH-Kidney Chromophobe" = "KICH", "KIRC-Kidney renal clear cell carcinoma" = "KIRC" ,"KIRP-Kidney renal papillary cell carcinoma" = "KIRP",
                              "LAML-Acute Myeloid Leukemia" = "LAML", "LGG-Brain Lower Grade Glioma" = "LGG" , "LIHC-Liver hepatocellular carcinoma" = "LIHC", "LUAD-Lung adenocarcinoma" = "LUAD","LUSC-	Lung squamous cell carcinoma" = "LUSC", "MESO-Mesothelioma" ="MESO", "OV-Ovarian serous cystadenocarcinoma" ="OV",   "PAAD-	Pancreatic adenocarcinoma" = "PAAD", "PCPG-Pheochromocytoma and Paraganglioma" = "PCPG", "PRAD-Prostate adenocarcinoma" = "PRAD",
                              "READ-Rectum adenocarcinoma" = "READ", "SARC-Sarcoma" = "SARC", "SKCM-Skin Cutaneous Melanoma" = "SKCM" ,"STAD-Stomach adenocarcinoma" = "STAD", "TGCT-Testicular Germ Cell Tumors" = "TGCT", "THCA-Thyroid carcinoma" = "THCA", "THYM-Thymoma" = "THYM", "UCEC-Uterine Corpus Endometrial Carcinoma" = "UCEC", "UCS-Uterine Carcinosarcoma" = "UCS",  "UVM-Uveal Melanoma" = "UVM"),
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
      # fileInput(inputId = ns("file"),
      #           label = "Also You Can Upload Your RDS or Excel Data Less Than 250 MB",  #I changed this part to below.
      #           accept = c(".rds", ".xlsx", ".xls"),
      #           multiple = FALSE),

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
                     label = "Normalize the data",
                     status = "info",
                     value = FALSE),
  
      ),
    
    # actionButton(inputId = ns("resetBtn"), label = "Clear Uploaded Data"),
    
    # textOutput(ns("fileInfos")),
    
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
   
    file_type <- reactive({file_ext(input$file$name)})
    
    Xproj$fileInfost <- reactive({input$file$name})
    
    # observeEvent(input$file,{
    #   
    #   if (!is.null(input$file)) {
    #   
    #   # Update the pickerInput choices
    #   updatePickerInput(session, "proj", choices = c(paste0(Xproj$fileInfost())))
    #   }
    #   })
    
    # observeEvent(input$resetBtn,{
    # 
    #     # Update the pickerInput choices
    #     updatePickerInput(session, "proj", choices = c("ACC-Adrenocortical carcinoma" = "ACC","BLCA-Bladder Urothelial Carcinoma" = "BLCA", "BRCA-Breast invasive carcinoma" = "BRCA", "CESC-Cervical squamous cell carcinoma and endocervical adenocarcinoma" = "CESC", "CHOL-Cholangiocarcinoma" = "CHOL", "COAD-Colon adenocarcinoma" = "COAD" ,"DLBC-Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "DLBC", "ESCA-Esophageal carcinoma" = "ESCA", "GBM-Glioblastoma multiforme" ="GBM" , "HNSC-Head and Neck squamous cell carcinoma" = "HNSC", "KICH-Kidney Chromophobe" = "KICH", "KIRC-Kidney renal clear cell carcinoma" = "KIRC" ,"KIRP-Kidney renal papillary cell carcinoma" = "KIRP",
    #                                                    "LAML-Acute Myeloid Leukemia" = "LAML", "LGG-Brain Lower Grade Glioma" = "LGG" , "LIHC-Liver hepatocellular carcinoma" = "LIHC", "LUAD-Lung adenocarcinoma" = "LUAD","LUSC-	Lung squamous cell carcinoma" = "LUSC", "MESO-Mesothelioma" ="MESO", "OV-Ovarian serous cystadenocarcinoma" ="OV",   "PAAD-	Pancreatic adenocarcinoma" = "PAAD", "PCPG-Pheochromocytoma and Paraganglioma" = "PCPG", "PRAD-Prostate adenocarcinoma" = "PRAD",
    #                                                    "READ-Rectum adenocarcinoma" = "READ", "SARC-Sarcoma" = "SARC", "SKCM-Skin Cutaneous Melanoma" = "SKCM" ,"STAD-Stomach adenocarcinoma" = "STAD", "TGCT-Testicular Germ Cell Tumors" = "TGCT", "THCA-Thyroid carcinoma" = "THCA", "THYM-Thymoma" = "THYM", "UCEC-Uterine Corpus Endometrial Carcinoma" = "UCEC", "UCS-Uterine Carcinosarcoma" = "UCS",  "UVM-Uveal Melanoma" = "UVM"))
    #   
    # })
    
    # observeEvent(input$file,{ 
    # 
    # output$fileInfos <- renderText({
    #   
    #   if (!is.null(Xproj$fileInfost())) {
    #     paste("Now You Can Select Your Data On Project List and Click Load Data")
    #   } else {
    #     "No file uploaded"
    #   }
    # })
    # })
    
    observeEvent(input$select_data_upload_user_data, {
      
      validate(need(input$file, ""))
      
      shinyjs::reset("proj")
      shinyjs::reset("file")
      shinyjs::reset(Xproj$fileInfost())
      # Xproj$fileInfost()<-NULL,
      
      # shinyjs::reset("fileInfos")
      shinyjs::reset("fileInfos2")
      shinyjs::reset("fileInfos3")
      shinyjs::reset("fileInfos4")
      shinyjs::reset("fileInfos5")
      
      shinyjs::reset("patient_hist")
      shinyjs::reset("gender_hist")
      shinyjs::reset("definition_hist")
      shinyjs::reset("age_hist")

      # output$patient_hist <- renderPlotly({NULL})
      # output$gender_hist <- renderPlotly({NULL})
      # output$definition_hist <- renderPlotly({NULL})
      # output$age_hist <- renderPlotly({NULL})

      # output$fileInfos <- renderText({NULL})
      output$fileInfos2 <- renderText({NULL})
      output$fileInfos3 <- renderText({NULL})
      output$fileInfos4 <- renderDataTable({NULL})
      output$fileInfos5 <- renderText({NULL})
    })
    
    observeEvent(input$select_data_upload_user_data, {
      
      # validate(need(input$select_data_upload_user_data, ""))
      
      shinyjs::reset("proj")
      # shinyjs::reset(Xproj$a())

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
          num_rows <- nrow(ydata())
          num_cols <- ncol(ydata())

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
            Xproj$a()[1:100,1:100]

            })
          output$fileInfos5 <- renderText({
            validate(need(input$run, ""))
            "The first 100 rows and first 100 columns of your data are shown..."
          })
        }
      
      else if ((file_type() %in% c("xlsx", "xls"))) {
        
        xdata <- input$file
        ydata <- reactive({as.data.table(read_excel(paste0(input$file$datapath)))})
          
          file_size <- round(xdata$size / (1024^2), 2)
          num_rows <- nrow(ydata())
          num_cols <- ncol(ydata())
          
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
            ydata()[1:100,1:100]
            
          })
          output$fileInfos5 <- renderText({
            validate(need(input$run, ""))
            "The first 100 rows and first 100 columns of your data are shown..."
          })
        }
        
      }
      
      })
    
    Xproj$cancer_length <- reactive({length(as.vector(input$proj))}) ## a reactive that created for other modules to use the length information for several cancers(Cagatay)
    
    Xproj$a<- eventReactive(list(input$run, input$file, input$proj), {
      
      # browser()

      validate(need(input$run, ""))

      if(!is.null(input$proj))  {

        if (Xproj$cancer_length() == 1 ){
      
      readRDS(paste0("projects/", input$proj, ".rds"))
      
      # compressed<-readRDS(paste0("projects/", input$proj, ".rds"))
      #
      # zstd_unserialize(compressed
      
        }
        
        
        else if (!(Xproj$cancer_length() == 1) && any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj)) {
          
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
          
          big_data = do.call(rbind, dflist)
          
          return(big_data)
          
          
        }}
    
        else if (!is.null(input$file)) {
          
          # browser()
          
          validate(need(input$file, ""))
          
          
          if (file_type() %in% c("rds","RDS","Rds")){
            
            uploaded_data <- as.data.table(readRDS(paste0(input$file$datapath)))
            
            if (!"meta.definition" %in% names(uploaded_data)) {
              uploaded_data[, meta.definition := "All Samples"]
            }
            
            if(input$nrm_dat == TRUE){
              
              # browser()
              
              df_num = uploaded_data %>% select(where(is.numeric))
              df_nonnum = uploaded_data %>% select(-where(is.numeric))
              
              df_nongene = df_num %>% select(where(~length(unique(.))< 10))
              
              df_gene = df_num %>% select(-where(~length(unique(.))< 10))
              
              df_gene = log(df_gene)
              
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
              
              df_nongene = df_num %>% select(where(~length(unique(.))< 10))
              
              df_gene = df_num %>% select(-where(~length(unique(.))< 10))
              
              df_gene = log(df_gene)
              
              uploaded_data_xl = cbind(df_gene, df_nonnum,df_nongene)
              
              uploaded_data_xl
              
            }
            
            return(uploaded_data_xl)
            
            # uploaded_data_xl
          }
          

          
        }

      
    })

    output$gender_hist<- renderPlotly({
      
      validate(need(input$run, ""))
      
     if (any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj) ){
      
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
      
      if (any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj) ){
        
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
      
      if (any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj) ){
        
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
      
      if (any(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") %in% input$proj) ){
        
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

