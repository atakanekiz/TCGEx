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

#ui

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
    
  
        pickerInput(inputId = ns("proj"), 
                    "TCGA projects",
                    choices = c("ACC-Adrenocortical carcinoma" = "ACC","BLCA-Bladder Urothelial Carcinoma" = "BLCA", "BRCA-Breast invasive carcinoma" = "BRCA", "CESC-Cervical squamous cell carcinoma and endocervical adenocarcinoma" = "CESC", "CHOL-Cholangiocarcinoma" = "CHOL", "COAD-Colon adenocarcinoma" = "COAD" ,"DLBC-Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "DLBC", "ESCA-Esophageal carcinoma" = "ESCA", "GBM-Glioblastoma multiforme" ="GBM" , "HNSC-Head and Neck squamous cell carcinoma" = "HNSC", "KICH-Kidney Chromophobe" = "KICH", "KIRC-Kidney renal clear cell carcinoma" = "KIRC" ,"KIRP-Kidney renal papillary cell carcinoma" = "KIRP",
                                "LAML-Acute Myeloid Leukemia" = "LAML", "LGG-Brain Lower Grade Glioma" = "LGG" , "LIHC-Liver hepatocellular carcinoma" = "LIHC", "LUAD-Lung adenocarcinoma" = "LUAD","LUSC-	Lung squamous cell carcinoma" = "LUSC", "MESO-Mesothelioma" ="MESO", "OV-Ovarian serous cystadenocarcinoma" ="OV",   "PAAD-	Pancreatic adenocarcinoma" = "PAAD", "PCPG-Pheochromocytoma and Paraganglioma" = "PCPG", "PRAD-Prostate adenocarcinoma" = "PRAD",
                                "READ-Rectum adenocarcinoma" = "READ", "SARC-Sarcoma" = "SARC", "SKCM-Skin Cutaneous Melanoma" = "SKCM" ,"STAD-Stomach adenocarcinoma" = "STAD", "TGCT-Testicular Germ Cell Tumors" = "TGCT", "THCA-Thyroid carcinoma" = "THCA", "THYM-Thymoma" = "THYM", "UCEC-Uterine Corpus Endometrial Carcinoma" = "UCEC", "UCS-Uterine Carcinosarcoma" = "UCS",  "UVM-Uveal Melanoma" = "UVM"),
                    #selectize = T,
                    # options = list('actions-box' = TRUE), #build buttons for collective selection
                    multiple = TRUE),
    actionBttn(inputId = ns("run"), 
               label = "Load Data",
               style = "unite",
               block = FALSE,
               color = "primary"),
    br(),
    br(),
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

select_data_server<- function(id,Xproj){
  moduleServer(id, function(input, output, session) {
   
    Xproj$cancer_length <- reactive({length(as.vector(input$proj))}) ## a reactive that created for other modules to use the length information for several cancers(Cagatay)
    Xproj$a <- eventReactive(input$run, {
      
     
      
      if(Xproj$cancer_length() == 1){
        
        readRDS(paste0("projects/", input$proj, ".rds"))
        
        # compressed<-readRDS(paste0("projects/", input$proj, ".rds"))
        # 
        # zstd_unserialize(compressed)

          
      } else {
        
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

        
      }
        
      
    }
    
    
    
    )
    
    output$gender_hist<- renderPlotly({
      
      validate(need(input$run, "Load TCGA project to see descriptive statistics"))
      
      #gghistogram(Xproj$a(), "meta.gender", stat="count", legend="none",
                  #font.x=18, font.y=18, font.tickslab = 18,
                  #fill="meta.gender", palette = c("skyblue", "coral", "gold"))
      
      fig_gender <- plot_ly(Xproj$a(),  labels = ~meta.gender, type = 'pie',
                            marker = list(color = viridis::viridis_pal(begin = 0.2, end = 0.8)(4)))
      
      fig_gender <- fig_gender %>% layout(title = 'Gender Statistics For Chosen Data',
                            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        
      fig_gender
      
    })
    
    output$patient_hist <- renderPlotly({
      
      validate(need(input$run, "Load TCGA project to see descriptive statistics"))
      
      fig_patient <- plot_ly(Xproj$a(),  labels = ~meta.project_id, type = 'pie',
                             marker = list(color = viridis::viridis_pal(begin = 0.2, end = 0.8)(4)))
      
      fig_patient <- fig_patient %>% layout(title = 'Total Patient Number For Chosen Data',
                                          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      
      fig_patient

    })
    
    
    output$definition_hist <- renderPlotly({
      
      #gghistogram(Xproj$a(), "meta.definition", stat="count",
                  #font.x=18, font.y=18, font.tickslab = 18,
                  #fill="meta.definition", legend="none") +
        #rotate_x_text(45)
      
      Xproj$a() %>%
        plot_ly(x = ~meta.definition,
                textangle = 45,
                marker = list(color = viridis::viridis_pal(option = "C", direction = -1)(10)))%>%
        add_histogram()
      
    })
    
    output$age_hist <- renderPlotly({
      
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
    })
  
    
  }
  )}

