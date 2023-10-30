###########
library(dplyr)
library(edgeR)
library(tidyverse)
library(glmnet)
library(msigdbr)
library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(plotly)
library(shinyjs)
library(readxl)
library(shinydashboardPlus) 
library(fresh) 
library(RColorBrewer) 
library(rintrojs)  
library(shinybusy) 
library(shinyWidgets)
# library(zstdlite)
############
#Functions
zero_adjuster <-  function(mdata, max_zero_percent, subgroup = "") {
  
  if (subgroup == "alltranscripts") {
    mdata <- select(mdata,-c(starts_with("meta.")))
    zerofreq_list <- sapply(mdata, function(x){ length(which(x==0))/length(x) *100 })
    available <- names(zerofreq_list[zerofreq_list[] < max_zero_percent])
  } else if (subgroup == "cibersort") {
    
    
    mdata <- mdata %>% 
      
      #'[#########################################################################################################################]   
      #'[would be better to explicitly select columns instead of position based selection]
      # select((which(colnames(mdata) == "meta.immune.subtype") + 1):ncol(mdata))
      select(matches("cells$|cells.(?!\\.)|lymphocytes|macrophages|neutrophils|monocytes|eosinophils", perl = T))
    #'[#########################################################################################################################]
    
    zerofreq_list <- sapply(mdata, function(x){ length(which(x==0))/length(x) *100 })
    available <- names(zerofreq_list[zerofreq_list[] < max_zero_percent])
  } 
  return(as.data.frame(available))
}
ggplot_ml_coefvlambda = function(cvfit_glm, all ,predictor_name) {
  fit <- cvfit_glm$glmnet.fit
  l1se <- cvfit_glm$lambda.1se
  lmin <- cvfit_glm$lambda.min
  
  lambda_df <- data.frame(lambda = fit$lambda, penalty = names(fit$a0))
  if (all == TRUE) {
    
    results <- cbind(
      Transcript = rownames(fit$beta),
      as.data.frame(as.matrix(fit$beta))
    ) %>%
      tidyr::gather(penalty, coefficients, -Transcript) %>%
      inner_join(lambda_df) 
    
  } else {
    results <- cbind(
      Transcript = rownames(fit$beta),
      as.data.frame(as.matrix(fit$beta))
    ) %>%
      tidyr::gather(penalty, coefficients, -Transcript) %>%
      inner_join(lambda_df) %>% filter(Transcript %in% c(predictor_name))
    
  }
  num_groups <- length(unique(results$Transcript))
  
  color_palette <- colorRampPalette(brewer.pal(8, "Paired"))(num_groups)
  g <- ggplot() +
    geom_line(data = results, aes(log(lambda), coefficients, color = Transcript), show.legend = TRUE) +
    geom_vline(xintercept = log(lmin)) +
    geom_vline(xintercept = log(l1se)) +
    theme(legend.title = element_blank(), legend.position = 'none') +
    scale_color_manual(values = color_palette)
  
  return(list('plot' = g))
}
ggplot_boxplot_variables <- function(transcript_list, pre_data) {
  transcript_list = c(transcript_list, "response")
  predata_melt <- pre_data %>%
    select(all_of(transcript_list)) %>%
    gather(key = "Gene", value = "Expression")
  num_groups <- length(unique(predata_melt$Gene))
  color_palette <- colorRampPalette(brewer.pal(8, "Paired"))(num_groups)
  
  p <- ggplot(predata_melt, aes(x = Gene, y = Expression, fill = Gene)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "text", label = "*", size = 6, color = "black") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_manual(values = color_palette)
  return(p)
}

#UI 
dataprepInputControl_UI <- function(id) {
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
    tags$style(type = 'text/css',
               ".selectize-input { word-wrap : break-word;} .selectize-input { word-break: break-word;} .selectize-dropdown {word-wrap : break-word;} "
    ),
    h3("Set up variables for analysis"),
    
    
    selectizeInput(NS(id,"sample_type_reg"), "Select Sample Type", choices = NULL, multiple = TRUE,
                   options=list(placeholder = "eg. Primary solid tumor")),
    sliderInput(NS(id,"max_zero_percent"),"Maximum allowed percentage of zero expressors",
                min = 0, max = 100, value = 100, step = 1),
    
    h3("Response Variable"),
    radioButtons(NS(id,"response_prep_method"), "Select an input method",
                 c("Manually select gene(s) as response variable" = "gene_list",
                   "Create a response variable from MSigDB gene sets" = "msigdb",
                   "Use CIBERSORT immune cell signatures as response variable" = "cibersort")
    ),
    conditionalPanel(condition = "input.response_prep_method == 'msigdb' ", ns = ns,
                     selectizeInput(NS(id,"msigdb_setnames_response"), "MSigDB Human Collections", choices = c("Hallmark gene sets (H)" = "H",
                                                                                                               "Positional gene sets (C1)" = "C1",
                                                                                                               "Curated gene sets (C2)" = "C2",
                                                                                                               "Regulatory target gene sets (C3)" = "C3",
                                                                                                               "Computational gene sets (C4)" = "C4",
                                                                                                               "Ontology gene sets (C5)" = "C5" ,
                                                                                                               "Oncogenic gene sets (C6)" = "C6",
                                                                                                               "Immunologic gene sets (C7)" = "C7",
                                                                                                               "Cell type signature gene sets (C8)" = "C8"
                     )),
                     conditionalPanel(condition = "input.msigdb_setnames_response == 'C2'|input.msigdb_setnames_response =='C3'|
                                      input.msigdb_setnames_response =='C4'|
                                      input.msigdb_setnames_response =='C5'|
                                      input.msigdb_setnames_response =='C7' ", ns = ns, 
                                      selectizeInput(NS(id,"msigdb_subc_response"),"Select subcategory" ,choices = c(""))
                     ),
                     
                     
                     
                     selectizeInput(NS(id,"msigdb_gene_set_response"), "Select a gene set to create an averaged response variable",
                                    choices = c(""))
    ),
    conditionalPanel(condition = "input.response_prep_method == 'gene_list'", ns = ns,
                     radioButtons(NS(id,"obtain_response"), "Select how to specify genes",c("Type gene names" = "gene_enter", "Upload xlsx/xls file" = "txt_upload")),
                     conditionalPanel(condition = "input.obtain_response == 'gene_enter' ", ns = ns,
                                      selectizeInput(NS(id,"gene_list_response"),"Select gene(s)", choices = NULL, multiple = TRUE, options=list(placeholder = "CD8A, hsa.miR.935, etc. "))
                     ),
                     conditionalPanel(condition = "input.obtain_response == 'txt_upload' ", ns = ns,
                                      fileInput(NS(id,"response_set_file"), label = tags$span(
                                        "Upload xlsx/xls File",
                                        tags$i(
                                          class = "glyphicon glyphicon-info-sign", 
                                          style = "color:#0072B2;",
                                          title = "The xlsx/xls file should contain a list of gene symbols. Neither matrix/dataframe structure nor named columns are required. If a matrix/dataframe provided, first column will be extracted as a gene list."
                                        ), tags$br(),
                                        a(href="sample_gene_input.xlsx", "Sample Input File", download=NA, target="_blank")), 
                                        accept =  c(".xls", ".xlsx") 
                                      )
                     )
                     
                     
    ),
    conditionalPanel(condition = "input.response_prep_method == 'cibersort' ", ns = ns,
                     selectizeInput(NS(id,"cibersort_response_var"), "Select immune cell signature", choices = NULL, multiple = FALSE)
    ),
    h3("Predictor Variables"),
    radioButtons(NS(id,"predictor_prep_method"), "Select an input method",
                 c("Manually select genes as predictor variables" = "gene_list",
                   "Select predictor variables from MSigDB gene sets" = "msigdb"
                 )
    ),
    conditionalPanel(condition = "input.predictor_prep_method == 'msigdb' ", ns = ns,
                     selectizeInput(NS(id,"msigdb_setnames_predictor"), "MSigDB Human Collections", choices = c("Hallmark gene sets (H1)" = "H",
                                                                                                                "Positional gene sets (C1)" = "C1",
                                                                                                                "Curated gene sets (C2)" = "C2",
                                                                                                                "Regulatory target gene sets (C3)" = "C3",
                                                                                                                "Computational gene sets (C4)" = "C4",
                                                                                                                "Ontology gene sets (C5)" = "C5" ,
                                                                                                                "Oncogenic gene sets (C6)" = "C6",
                                                                                                                "Immunologic gene sets (C7)" = "C7",
                                                                                                                "Cell type signature gene sets (C8)" = "C8"
                     )),
                     conditionalPanel(condition = "input.msigdb_setnames_predictor == 'C2'|input.msigdb_setnames_predictor =='C3'|
                                      input.msigdb_setnames_predictor =='C4'|
                                      input.msigdb_setnames_predictor =='C5'|
                                      input.msigdb_setnames_predictor =='C7' ", ns = ns, 
                                      selectizeInput(NS(id,"msigdb_subc_predictor"),"Select subcategory" ,choices = c(""))
                     ),
                     
                     
                     
                     
                     selectizeInput(NS(id,"msigdb_gene_set_predictor"), "Select response variables from a gene set",
                                    choices = c(""))
    ),
    conditionalPanel(condition = "input.predictor_prep_method == 'gene_list' ", ns = ns,
                     radioButtons(NS(id,"obtain_predictor"), "Select how to specify genes",c("Type gene names" = "gene_enter",
                                                                                             "Upload xlsx/xls file" = "txt_upload",
                                                                                             "Use all miRNAs available on the selected data" = "allmiRNA_aspredictor")),
                     conditionalPanel(condition = "input.obtain_predictor == 'gene_enter' ", ns = ns,
                                      selectizeInput(NS(id,"gene_list_predictor"),"Select gene(s)", choices = c(""), multiple = TRUE, options=list(placeholder = "FGR, hsa.miR.107, etc. "))
                     ),
                     conditionalPanel(condition = "input.obtain_predictor == 'txt_upload' ", ns = ns,
                                      fileInput(NS(id,"predictor_set_file"), label = tags$span(
                                        "Upload xlsx/xls File",
                                        tags$i(
                                          class = "glyphicon glyphicon-info-sign", 
                                          style = "color:#0072B2;",
                                          title = "The xlsx/xls file should contain a list of gene symbols. Neither matrix/dataframe structure nor named columns are required. If a matrix/dataframe provided, first column will be extracted as a gene list." 
                                        ),tags$br(),
                                        a(href="sample_gene_input.xlsx", "Sample Input File", download=NA, target="_blank")
                                      ), accept =  c(".xls", ".xlsx") 
                                      )
                                      
                                      
                                      
                     )
    ),
    
    # hr(),
    p("Please move on to the 'Ridge/Elastic Net/Lasso Regression' tab after variable selection"),
    hr()
    
  )
}

regression_sidecontrols <- function(id) {  
  ns <- NS(id)
  tagList(
    h3("Workflow"),
    radioButtons(NS(id,"regression_workflow"), "Select a workflow type to train the model",
                 c("Use 100% of the data to train a model" = "create_model",
                   "Split data into train and test sets and evaluate model accuracy" = "ctest_model")
    ),
    conditionalPanel(
      condition = "input.regression_workflow == 'ctest_model' ", ns = ns, 
      sliderInput(NS(id,"train_percentage"),
                  "Set the training set percentage (If 30% selected, 30% of the data will be used as training set)",
                  min = 10, max = 90, value = 10, step = 10)
      
    ),
    h3("Regression: Ridge/ElasticNet/Lasso"),
    fluidRow(
      column(12,
             sliderInput(NS(id,"user_alpha"),
                         
                         label = tags$span(
                           "Select regression technique by setting shrinkage penalty term (alpha)",
                           tags$i(
                             class = "glyphicon glyphicon-info-sign", 
                             style = "color:#0072B2;",
                             title = "Value of 1 corresponds to LASSO regression where some coefficients will be shrunken (ie. penalized) all the way to zero. Value of 0 corresponds to Ridge regression where some coefficients will converge to (but not reach) zero. Other values correspond to elastic net regression where the penalty is a mixture of the previous approaches."
                             
                           )), min = 0, max = 1, value = 1, step = 0.2
                         
                         
                         
                         
                         
             )
      )
    ),
    conditionalPanel(
      condition = "input.regression_workflow == 'create_model' ", ns = ns,
      actionButton(ns("run_train"), "Train the model")
    ),
    conditionalPanel(
      condition = "input.regression_workflow == 'ctest_model' ", ns = ns,
      actionButton(ns("run_traintest"), "Train and test the model")
    ),
    
    downloadButton(NS(id,"download_coef"),"Download Coefficients")
    
  )
  
}

ml_ui <- function(id) {
  ns <- NS(id)
  navbarPage(
    HTML("<i class='fas fa-atom'></i> TCGEx ML"),
    tags$style(HTML("
  .navbar-default {
    background-color: #80a1c0; /* Custom blue background color for navbar */
  }
  .navbar-default .navbar-brand {
    color: white; /* White text color for navbarPage title */
    pointer-events: none; /* Disable all cursor interactions with the title element */
    margin-right: 24px; /* Add some space between the title and the icon */
  }
  .navbar-default .navbar-brand i {
    float: right; /* Float the icon to the right */
    margin-left: 10px; /* Add some space between the icon and the title */
  }
  .nav.navbar-nav > li > a {
    color: white; /* White text color for navbar items */
  }
")),
    tabPanel(
      "Variable selection",
      add_busy_spinner(
        spin = "cube-grid",
        position = "top-right",
        color = "#01303f",
        margins = c(300, 500),
        height = "60px",
        width = "60px"),
      sidebarPanel(
        dataprepInputControl_UI("ml"),
        introjsUI(),
        actionButton(ns("mldata_help"), "Tutorial"),
        
        textOutput(ns("filewarning_7")) ,
        
        width = 3
        
      ),
      mainPanel(
        fluidRow(
          column(6,h3("Response variable(s)")),
          column(6,h3("Predictor variables"))
        ),
        fluidRow(
          column(6,verbatimTextOutput(NS(id,"response_set"))),
          column(6,verbatimTextOutput(NS(id,"predictor_set")))
        ),
        
        fluidRow(
          column(6,textOutput(NS(id,"validation_message_response"))),
          column(6,textOutput(NS(id,"validation_message_predictor")))
        )
        
        
        
      )
    ),
    tabPanel(
      "Ridge/Elastic Net/LASSO Regression",
      add_busy_spinner(
        spin = "cube-grid",
        position = "top-right",
        color = "#01303f",
        margins = c(300, 500),
        height = "60px",
        width = "60px"),
      
      sidebarPanel(
        regression_sidecontrols("ml"),          
        introjsUI(),
        actionButton(ns("mlregression_help"), "Tutorial"),
        width = 3
      ),
      mainPanel(
        fluidRow(
          column(3,
                 selectizeInput(NS(id,"lambda_for_coef"), "Select lambda value at which model coefficients are being displayed: ",
                                choices = c("Lambda + 1 standard error (most regularized model)" = "lambda.1se",
                                            "Minimum lambda (minimum cross-validation error)" = "lambda.min"
                                )),
                 textOutput(NS(id,"lambda_value_min")),
                 textOutput(NS(id,"lambda_value_1se")),
                 conditionalPanel(
                   condition = "input.regression_workflow == 'ctest_model' ", ns = ns, 
                   selectizeInput(NS(id,"lambda_for_prediction"), "Select lambda value for prediction: ", 
                                  choices = c("lambda.min","lambda.1se")),
                   actionButton(NS(id,"run_prediction"), "Run prediction"),
                   textOutput(NS(id,"model_error"))
                   
                 )
          ),
          column(3,
                 DTOutput(NS(id,"coef_data")),
                 hr(),
                 materialSwitch(
                   NS(id,"plot_setting"),
                   label = "Plot selected coefficient(s)",
                   value = FALSE,
                   status = "primary"
                 )
          ),
          column(6,
                 conditionalPanel(
                   condition = "input.plot_setting == false ", ns = ns, 
                   plotlyOutput(NS(id,"coef_lambda_plot"))
                   
                 ),
                 conditionalPanel(
                   condition = "input.plot_setting == true ", ns = ns,
                   plotlyOutput(NS(id,"coef_lambda_plot_some"))
                   
                 )
                 
                 
                 
                 
                 
          )
          
          
          
        ),
        fluidRow(
          column(6,
                 conditionalPanel(
                   condition = "input.plot_setting == true ", ns = ns,
                   plotlyOutput(NS(id,"coef_box_plot"))
                   
                 )),
          column(6, 
                 plotOutput(NS(id,"lambda_error_plot")))
        )
        
        
      )
    )
    
  )
}

#Server Modules

data_prep_ml_server <- function(id,Xproj) {
  moduleServer(id,function(input,output,session){
    
    output$filewarning_7 <- renderText({
      
      if (!is.null(Xproj$fileInfost())) {
        shinyalert("Warning!", "To perform analysis using MsigDB gensets, please make sure that your gene names are in MsigDB format, otherwise you may receive errors.") }
    })
    
    msigdb_gene_sets =  reactive({readRDS(paste0("genesets/", "msigdb_collections", ".rds"))})
    
    
    #Help
    help_dataprep = reactive({
      
      data.frame(
        
        element = paste0("#", session$ns(c(NA,NA, "response_prep_method", "predictor_prep_method")) ),
        intro = paste(c(
          "Welcome to the TCGEx machine learning module. This module is based on generalized linear models and it allows you to perform regularized regression analysis using custom response and predictor variables. Please continue with the tutorial to learn how to use this module.",
          
          "This slider allows you to eliminate genes that are expressed at low levels. The selection here specifies the maximum allowed percentage of zero expression in a given gene. For instance, if this number is set to 50, genes that are not expressed in 50% or more of the samples in the analysis. The default value of 100 indicates that there is no filtering applied.",  
          
          "In this panel, you can specify response variables either by <b>i)</b> entering gene names manually, <b>ii)</b>using genes from MSigDB gene sets, or <b>iii)</b> selecting one of the previously calculated immune cell signatures. When multiple genes are entered or gene sets are selected, a single response variable is calculated by averaging the expression values. For manual gene selection, you can type gene names in the box, or upload a xlsx/xls file containing gene names.",
          
          "In this panel, you can enter predictor variables either by <b>i)</b> entering them manually (you can type or upload a file), or <b>ii)</b> using genes from MSigDB gene sets. The relationship of these predictor variables and the previously specified response variable will be examined in regularized regression models. After making your selections here, please continue to the regression tab."
          
        ))
      )
      
      
      
    })
    
    observeEvent(input$mldata_help, {
      introjs(session, options = list(steps = help_dataprep() ) )
    })
    
    
    #########Gene names from TCGA data
    available_genelist <- reactive({zero_adjuster(mdata = Xproj$a(),max_zero_percent = input$max_zero_percent, subgroup = "alltranscripts")})
    available_cibersort <- reactive({
      
      zero_adjuster(mdata = Xproj$a(),max_zero_percent = input$max_zero_percent, subgroup = "cibersort")})
    #Observers
    observe({
      subcat_response = names(msigdb_gene_sets()[[input$msigdb_setnames_response]])
      if (length(subcat_response) > 1)  {
        updateSelectizeInput(session,'msigdb_subc_response', choices = subcat_response , server = TRUE)
      } 
      
    })
    observe({
      subcat_predictor = names(msigdb_gene_sets()[[input$msigdb_setnames_predictor]])
      if (length(subcat_predictor) > 1)  {
        updateSelectizeInput(session,'msigdb_subc_predictor', choices = subcat_predictor , server = TRUE)
      }
    })
    observe({
      if (!is.null(Xproj$a())) {
        
        updateSelectizeInput(session,'sample_type_reg', choices = names(table(Xproj$a()$meta.sample_type)),selected = "", server = TRUE)
        cibersort_metrics <- available_cibersort()$available
        genelist <- available_genelist()$available
        updateSelectizeInput(session,'gene_list_response', choices = genelist, selected = "", server = TRUE)
        updateSelectizeInput(session, 'gene_list_predictor', choices = genelist,selected = "",server = TRUE)
        updateSelectizeInput(session,'cibersort_response_var', choices = cibersort_metrics, server = TRUE)
      }
    })
    
    #msigdb
    
    response_msigdb_handler <- reactive({
      ms = msigdb_gene_sets()
      if(input$msigdb_setnames_response == "C2"|input$msigdb_setnames_response =="C3"|
         input$msigdb_setnames_response =="C4"|input$msigdb_setnames_response =="C5"|input$msigdb_setnames_response =="C7") {
        
        response_m = ms[[input$msigdb_setnames_response]][[input$msigdb_subc_response]]
        gene_sets = names(table(response_m$gs_name))
        
        updateSelectizeInput(session,'msigdb_gene_set_response', choices = gene_sets , server = TRUE)
        response_m
        
      } else {
        response_m = ms[[input$msigdb_setnames_response]]
        gene_sets = names(table(response_m$gs_name))
        updateSelectizeInput(session,'msigdb_gene_set_response', choices = gene_sets , server = TRUE)
      }
      response_m
    })
    
    response_msigdb_genes <- reactive({
      if (!is.null(response_msigdb_handler())) {
        filter(response_msigdb_handler(), gs_name == input$msigdb_gene_set_response ) %>% select(gene_symbol)
      }
    })
    
    predictor_msigdb_handler <- reactive({
      ms = msigdb_gene_sets()
      if(input$msigdb_setnames_predictor == "C2"|input$msigdb_setnames_predictor =="C3"|
         input$msigdb_setnames_predictor =="C4"|input$msigdb_setnames_predictor =="C5"|input$msigdb_setnames_predictor =="C7") {
        
        predictor_m = ms[[input$msigdb_setnames_predictor]][[input$msigdb_subc_predictor]]
        gene_sets = names(table(predictor_m$gs_name))
        
        updateSelectizeInput(session,'msigdb_gene_set_predictor', choices = gene_sets , server = TRUE)
        predictor_m
        
      } else {
        predictor_m = ms[[input$msigdb_setnames_predictor]]
        gene_sets = names(table(predictor_m$gs_name))
        updateSelectizeInput(session,'msigdb_gene_set_predictor', choices = gene_sets , server = TRUE)
      }
      predictor_m
    })
    
    predictor_msigdb_genes <- reactive({
      if (!is.null(predictor_msigdb_handler())) {
        filter(predictor_msigdb_handler(), gs_name == input$msigdb_gene_set_predictor ) %>% select(gene_symbol)
      }
    })
    
    
    gene_list_selected_df_response <- reactive({as.data.frame(input$gene_list_response)})
    cibersort_list <- reactive({as.data.frame(input$cibersort_response_var)})
    gene_list_selected_df_predictor <- reactive({as.data.frame(input$gene_list_predictor)})
    
    #
    response_det <- reactive({
      if (input$response_prep_method == "msigdb") {
        list_r = response_msigdb_genes()
      } else if (input$response_prep_method == "gene_list") {
        if (input$obtain_response == "gene_enter") {
          list_r = gene_list_selected_df_response()
        } else {
          file <- input$response_set_file
          ext <- tools::file_ext(file$datapath)
          req(file)
          validate(need(ext == c("xls","xlsx"), "Please upload a xlsx/xls file"))
          list_r = as.data.frame(read_excel(file$datapath, sheet = 1, col_names = F))
        }
        
      } else if (input$response_prep_method == "cibersort") {
        list_r = cibersort_list()
      }
      if (length(list_r) > 0) {
        colnames(list_r) = "gene_symbol"
        
      }
      
      list_r = unique(list_r)
      ########################################
      
      genelist <- list_r$gene_symbol
      is.exist <- genelist %in% colnames(Xproj$a())
      response_missing <- list_r[which(is.exist == "FALSE"),]
      
      if (input$response_prep_method == "cibersort") {
        clean_response_set <- list_r
      } else {
        missing_genes <-  as.vector(response_missing)
        clean_response_set <- list_r
        for (i in missing_genes) {
          clean_response_set <- filter(clean_response_set, gene_symbol != i)
        }
      }
      return(list(c = clean_response_set, r = response_missing))
    })
    
    output$response_set <- renderPrint({
      validate(need(response_det()$c$gene_symbol, "Select response set."))
      
      clean_response_set = response_det()$c  
      return(clean_response_set$gene_symbol)
      
      # if(is.null(clean_response_set()$gene_symbol)) {
      #   hide("output")
      # } else {
      #   show("output")
      #   return(clean_response_set()$gene_symbol)
      # }
    }) 
    
    output$validation_message_response <- renderText({
      response_missing = response_det()$r
      if (length(response_missing) > 0) {
        paste("Listed genes are not present in the data and they are removed from the response list:\n",
              paste(response_missing, collapse = ","))
      } else {
        hide("output")
      }
    })
    
    #
    predictor_det <- reactive({
      
      if (input$predictor_prep_method == "msigdb") {
        list_p <- predictor_msigdb_genes()
      } else if (input$predictor_prep_method == "gene_list") {
        if (input$obtain_predictor == "gene_enter") {
          list_p <- gene_list_selected_df_predictor()
        } else if (input$obtain_predictor == "allmiRNA_aspredictor") {
          
          mir_list <-  available_genelist() %>% filter(startsWith(available,"hsa."))
          list_p <- as.data.frame(mir_list)
          
        } else {
          file <- input$predictor_set_file
          ext <- tools::file_ext(file$datapath)
          req(file)
          validate(need(ext == c("xls","xlsx"), "Please upload a xlsx/xls file"))
          list_p <- as.data.frame(read_excel(file$datapath,  sheet = 1, col_names = F))
          
        }
      }
      if (length(list_p) > 0) {
        colnames(list_p) <- "gene_symbol"
      }
      list_p = unique(list_p)
      
      genelist <- list_p$gene_symbol
      is.exist <- genelist %in% colnames(Xproj$a())
      predictor_missing <- list_p[which(is.exist == "FALSE"),]
      
      if (input$predictor_prep_method == "cibersort") {
        clean_predictor_set <- list_p
      } else {
        missing_genes <-  as.vector(predictor_missing)
        clean_predictor_set <- list_p
        for (i in missing_genes) {
          clean_predictor_set <- filter(clean_predictor_set, gene_symbol != i)
        }
      }
      return(list(c = clean_predictor_set, r = predictor_missing))
    })
    
    output$predictor_set <- renderPrint({
      validate(need(predictor_det()$c$gene_symbol, "Select predictor set."))
      
      clean_predictor_set = predictor_det()$c  
      return(clean_predictor_set$gene_symbol)
      
      # if(is.null(clean_response_set()$gene_symbol)) {
      #   hide("output")
      # } else {
      #   show("output")
      #   return(clean_response_set()$gene_symbol)
      # }
    }) 
    
    output$validation_message_predictor <- renderText({
      predictor_missing = predictor_det()$r
      if (length(predictor_missing) > 0) {
        paste("Listed genes are not present in the data and they are removed from the predictor list:\n",
              paste(predictor_missing, collapse = ","))
      } else {
        hide("output")
      }
    })
    
    
    
    
    
    #
    reg_data <- reactive({
      
      clean_response_set = response_det()$c
      clean_predictor_set = predictor_det()$c
      
      filtered_data = Xproj$a() %>% filter(meta.sample_type %in% input$sample_type_reg)
      
      select(filtered_data, clean_response_set[,1], clean_predictor_set[,1]) %>%
        
        mutate(response = select(., clean_response_set[,1]) %>% rowMeans()) %>%
        select(response, clean_predictor_set[,1]) %>% 
        na.omit()
    })
    
    
    return(reg_data)
  } )
}

ml_main_server <- function(id,regress_data,Xproj) {
  moduleServer(id,function(input,output,session) {
    
    
    #################################################HELP(INTRO.)
    help_regression = reactive({
      
      if (input$regression_workflow == "create_model") {
        
        return(
          data.frame(
            element = paste0("#", session$ns(c(NA, "regression_workflow","user_alpha", "lambda_for_coef + .selectize-control", NA)) ),
            intro = paste(c(
              "After specifying response and predictor variables in the previous tab, you can determine necessary parameters for regularized regression analysis here.",
              "You can perform regression on the whole data or split the data set into 'training' and 'test' subsets. Splitting allows examining the model accuracy through the mean-squared-error.",
              "You can choose the method of regularized regression using this slider. alpha = 1 corresponds to LASSO regression where some coefficients will be shrunk (ie. penalized) all the way to zero. alpha = 0 corresponds to Ridge regression where some coefficients will converge to (but not reach) zero. 0 < alpha < 1 corresponds to Elastic-Net regression where the penalty is a mixture of both.",
              "You can choose the lambda value at which the variable coefficients of the model are displayed. Lambda is the regularization parameter in the model. Minimum lambda is the value that gives the minimum cross-validation error in the regression. Lambda + 1 se is the value of lambda that gives the most regularized (ie. more penalized and simpler) model where the cross-validation error is within the one strandard error of the minimum. 
            <br>
            <br>
            For further details, please see Friedman J, Hastie T, and Tibshirani R. “Regularization Paths for Generalized Linear Models via Coordinate Descent.” Journal of Statistical Software, Articles 33 (1): 1–22. (2010) https://doi.org/10.18637/jss.v033.i01",
              
              "Graphs that will appear in the main panel will show how the increasing levels of model penalization affects predictor coefficient shrinkage and the overall mean-squared-error."
              
            ))
          ))
        
        
      } else if (input$regression_workflow == "ctest_model") {
        
        return(
          data.frame(
            element = paste0("#", session$ns(c(NA,"regression_workflow","train_percentage","lambda_for_prediction + .selectize-control","user_alpha","lambda_for_coef + .selectize-control")) ),
            intro = paste(c(
              "In this tab, you can determine necessary parameters to conduct regularized regression with the variables that you have already specified.",
              "You can either split the data into train and test sets or use all of it to create the model. If you choose splitting option, you will be able to
            make a prediction with your model and examine mean-squared-error.",
              "You can determine how much of the data will be used as train set and how much of it will be reserved for the test set by using the data spliting toggle.",
              "You can choose the lambda value for coefficients at which prediction of the test set will be performed.",
              "You can choose the method of regularized regression from the alpha parameter slider control. ɑ=0 for  Ridge, ɑ=1 for Lasso and 0<ɑ<1 for Elastic Net regression.",
              "You can choose the lambda value at which the variable coefficients of the model are displayed."
              
            ))
          ))
        
      }
    })
    
    observeEvent(input$mlregression_help, {
      introjs(session, options = list(steps = help_regression() ) )
    })
    
    ################################################################REGRESSION
    trows <- eventReactive(input$run_traintest, {
      if (input$regression_workflow == "ctest_model") {
        nrowd <- nrow(regress_data())
        set.seed(30)
        trowsx <- sample(1:nrowd, (input$train_percentage*nrowd)/100 ,replace = FALSE)
        trowsx
      } else {
        NULL
        
      }
    })
    
    
    cvfit <- eventReactive(c(input$run_train, input$run_traintest), {
      validate(need(input$sample_type_reg, "Choose at least one sample type."))
      if (input$run_train == 0 && input$run_traintest == 0) {
        return(NULL) 
      }
      if (input$regression_workflow == "create_model") {
        response_var <- as.matrix(select(regress_data(), response)) 
        predictor_var <-  as.matrix(select(regress_data(), -c("response")))
        set.seed(7)
        foldid <- sample(1:10,size = length(response_var), replace = TRUE)
        fitty <- cv.glmnet(predictor_var, response_var,weights = NULL, foldid = foldid, alpha = input$user_alpha)
      } else if (input$regression_workflow == "ctest_model") {
        train_regress_data <- regress_data()[trows(),]
        response_var_train <- as.matrix(select(train_regress_data, response))
        predictor_var_train <-  as.matrix(select(train_regress_data, -c("response")))
        set.seed(7)
        foldid <- sample(1:10,size = length(response_var_train), replace = TRUE)
        fitty <- cv.glmnet(predictor_var_train, response_var_train, foldid = foldid, alpha = input$user_alpha)
      } 
      fitty
    })
    
    
    
    
    prediction_error <- eventReactive(input$run_prediction, {
      if(input$regression_workflow == "ctest_model" & !is.null(trows())){
        
        test_regress_data <- regress_data()[-trows(),]
        response_var_test <- as.matrix(select(test_regress_data, response)) 
        predictor_var_test <-  as.matrix(select(test_regress_data, -c("response"))) 
        prediction <- predict(cvfit(), s = input$lambda_for_prediction, newx = predictor_var_test)
        
        error <- mean((response_var_test - prediction)^2)
        error
      } else  {
        
        
      }
      
    })
    
    coef_data <- reactive({
      
      if(!is.null(cvfit())) {
        c <- as.matrix(coef(cvfit(), s = input$lambda_for_coef))
        c <-  c[2:length(rownames(c)), ]
        c = c[order(abs(c), decreasing = TRUE)]
        as.data.frame(c)
      } else {
        NULL
      }
    }) 
    
    
    ######################################################OUTPUTS
    
    output$coef_data <- renderDT({
      if(is.null(coef_data())) {
        hide("output")
      } else {
        show("output")
        return(
          datatable(coef_data())
          
        )
      }
      
    })
    coefholder <- reactive({
      indexes = input$coef_data_rows_selected
      pl = rownames(coef_data())
      coefsel = pl[c(indexes)]
      coefsel
      
    }) 
    
    
    coef_lambda_all = reactive({
      ggplot_ml_coefvlambda(cvfit_glm = cvfit(),all = TRUE)$plot
    })
    coef_lambda_some = reactive({
      
      if (length(coefholder() >0)) {
        ggplot_ml_coefvlambda(cvfit_glm = cvfit(),all = FALSE,predictor_name = coefholder())$plot
      }
    })
    
    output$coef_lambda_plot_some <- renderPlotly({
      validate(
        need(input$sample_type_reg, "Choose at least one sample type."),
        need(cvfit(), "Run the analysis.")
        
      )
      coef_lambda_some()
    })
    output$coef_lambda_plot <- renderPlotly({
      validate(
        need(input$sample_type_reg, "Choose at least one sample type."),
        need(cvfit(), "Run the analysis.")
        
      )
      coef_lambda_all()
      
    })
    
    output$coef_box_plot = renderPlotly({
      ggplot_boxplot_variables(coefholder(),regress_data())
    })
    
    
    ############################################################
    output$lambda_error_plot <- renderPlot({
      validate(
        need(input$sample_type_reg, "Choose at least one sample type.")
      )
      if(!is.null(cvfit())) {
        plot(cvfit(),xvar = "lambda", label = TRUE)
      }
      
    })
    
    output$lambda_value_min <- renderText({
      paste("lambda.min: ", cvfit()$lambda.min)
    })
    
    output$lambda_value_1se <- renderText({
      paste("lambda.1se: ", cvfit()$lambda.1se)})
    
    output$model_error <- renderText({
      if(is.null(prediction_error())) {
        hide("output")
      } else {
        show("output")
        paste(" MSE Error: ", prediction_error() )
      }
      
    })
    
    output$download_coef <- downloadHandler(
      filename <-  function() {
        paste(input$lambda_for_coef, ".xlsx", sep = "")
      },
      content <-  function(file) {
        write.xlsx(coef_data(), file, row.names = TRUE)
      }
    )
    ############################################################
  })
}
