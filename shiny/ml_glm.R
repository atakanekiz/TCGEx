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

#Functions

zero_adjuster <-  function(mdata, max_zero_percent, subgroup = "") {
  
  if (subgroup == "alltranscripts") {
    mdata <- select(mdata,-c(starts_with("meta.")))
    zerofreq_list <- sapply(mdata, function(x){ length(which(x==0))/length(x) *100 })
    available <- names(zerofreq_list[zerofreq_list[] < max_zero_percent])
  } else if (subgroup == "cibersort") {
    mdata <- mdata %>% 
      select((which(colnames(mdata) == "meta.immune.subtype") + 1):ncol(mdata))
    zerofreq_list <- sapply(mdata, function(x){ length(which(x==0))/length(x) *100 })
    available <- names(zerofreq_list[zerofreq_list[] < max_zero_percent])
  } 
  return(as.data.frame(available))
}

#UI
dataprepInputControl_UI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Metadata"),
    sliderInput(NS(id,"max_zero_percent"),"Set the maximum percentage of acceptable zero count (for manual selections)",
                min = 0, max = 100, value = 0, step = 1),
    h3("Response Variable"),
    radioButtons(NS(id,"response_prep_method"), "Data preparation pipeline: ",
                 c("Create response set from gene sets available on MSigDB" = "msigdb",
                   "Create response set manually" = "gene_list",
                   "Use CIBERSORT metrics as a response variable" = "cibersort")
    ),
    conditionalPanel(condition = "input.response_prep_method == 'msigdb' ", ns = ns,
                     selectizeInput(NS(id,"msigdb_setnames_response"), "MSigDB Human Collections", choices = c("hallmark gene sets" = "h","ontology gene sets" = "c5" ,"oncogenic gene sets" = "c6",
                                                                                                               "immunologic gene sets" = "c7")),
                     conditionalPanel(condition = "input.msigdb_setnames_response == 'c5'", ns = ns, 
                                      radioButtons(NS(id,"go_term_response"),"Select go term" ,choices = c("GO:BP" = "go_bp","GO:CC" = "go_cc","GO:MF" = "go_mf"))
                     ),
                     conditionalPanel(condition = "input.msigdb_setnames_response == 'c7'", ns = ns,
                                      radioButtons(NS(id,"immuno_response"), "Select subcategory", choices = c("IMMUNESIGDB"= "immunesigdb","VAX" = "vax"))
                     ),
                     selectizeInput(NS(id,"msigdb_gene_set_response"), "Select a gene set as a response set",
                                    choices = c(""))
    ),
    conditionalPanel(condition = "input.response_prep_method == 'gene_list' ", ns = ns,
                     radioButtons(NS(id,"obtain_response"), "Create response set by",c("Upload TXT file" = "txt_upload", "Enter gene" = "gene_enter" )),
                     conditionalPanel(condition = "input.obtain_response == 'gene_enter' ", ns = ns,
                                      selectizeInput(NS(id,"gene_list_response"),"Select gene(s)", choices = NULL, multiple = TRUE)
                     ),
                     conditionalPanel(condition = "input.obtain_response == 'txt_upload' ", ns = ns,
                                      fileInput(NS(id,"response_set_file"), "Upload File",
                                                accept =  c(".txt",".csv"))
                     )
    ),
    conditionalPanel(condition = "input.response_prep_method == 'cibersort' ", ns = ns,
                     selectizeInput(NS(id,"cibersort_response_var"), "Select CIBERSORT metric", choices = NULL, multiple = FALSE)
    ),
    h3("Predictor Variables"),
    radioButtons(NS(id,"predictor_prep_method"), "Data preparation pipeline: ",
                 c("Create predictor set from gene sets available on MSigDB" = "msigdb",
                   "Create predictor set manually" = "gene_list")
    ),
    conditionalPanel(condition = "input.predictor_prep_method == 'msigdb' ", ns = ns,
                     selectizeInput(NS(id,"msigdb_setnames_predictor"), "MSigDB Human Collections", choices = c("hallmark gene sets" = "h","ontology gene sets" = "c5" ,
                                                                                                                "oncogenic gene sets" = "c6","immunologic gene sets" = "c7")),
                     conditionalPanel(condition = "input.msigdb_setnames_predictor == 'c5'", ns = ns,
                                      radioButtons(NS(id,"go_term_predictor"),"Select go term" ,choices = c("GO:BP" = "go_bp","GO:CC"= "go_cc","GO:MF"="go_mf"))
                     ),
                     conditionalPanel(condition = "input.msigdb_setnames_predictor == 'c7'", ns = ns,
                                      radioButtons(NS(id,"immuno_predictor"), "Select subcategory", choices = c("IMMUNESIGDB" = "immunesigdb","VAX" = "vax"))
                     ),
                     selectizeInput(NS(id,"msigdb_gene_set_predictor"), "Select a gene set as a response set",
                                    choices = c(""))
    ),
    conditionalPanel(condition = "input.predictor_prep_method == 'gene_list' ", ns = ns,
                     radioButtons(NS(id,"obtain_predictor"), "Create predictor set by",c("Upload TXT file" = "txt_upload", "Enter gene" = "gene_enter",
                                                                                         "Use all mRNAs available on the selected data(!)" = "allmRNA_aspredictor",
                                                                                         "Use all miRNAs available on the selected data" = "allmiRNA_aspredictor")),
                     conditionalPanel(condition = "input.obtain_predictor == 'gene_enter' ", ns = ns,
                                      selectizeInput(NS(id,"gene_list_predictor"),"Select gene(s)", choices = c(""), multiple = TRUE)
                     ),
                     conditionalPanel(condition = "input.obtain_predictor == 'txt_upload' ", ns = ns,
                                      fileInput(NS(id,"predictor_set_file"), "Upload File",
                                                accept =  c(".txt","csv"))
                     )
    )
    
    
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
                  min = 10, max = 90, value = 10, step = 10),
      selectizeInput(NS(id,"lambda_for_prediction"), "Select lambda value for prediction: ", 
                     choices = c("lambda.min","lambda.1se"))
      
    ),
    h3("Regression: Ridge/ElasticNet/Lasso"),
    fluidRow(
      column(12,
             sliderInput(NS(id,"user_alpha"),"Enter an alpha parameter (shrinkage penalty term) to set the regression technique", min = 0, max = 1, value = 1, step = 0.2)
      )
    ),
    
    actionButton(ns("run"), "Run"),
    downloadButton(NS(id,"download_coef"),"Download Coefficients")
          
  )
  
}
ml_ui <- function(id) {
  ns <- NS(id)
  navbarPage(
    "TCGExplorer ML",
    tabPanel(
      "Data&Prep",
      sidebarPanel(
        dataprepInputControl_UI("ml"),
        introjsUI(),
        actionButton(ns("mldata_help"), "Tutorial")
        
        
      ),
      mainPanel(
        fluidRow(
          column(6,h3("Response")),
          column(6,h3("Predictor"))
        ),
        fluidRow(
          column(6,verbatimTextOutput(NS(id,"response_set"))),
          column(6,verbatimTextOutput(NS(id,"predictor_set")))
        ),
        fluidRow(
          column(6,h6("Listed genes are not present in the data and they are removed from the response list.")),
          column(6,h6("Listed genes are not present in the data and they are removed from the predictor list."))
        ),
        fluidRow(
          column(6,verbatimTextOutput(NS(id,"validation_message_response"))),
          column(6,verbatimTextOutput(NS(id,"validation_message_predictor")))
        )
        
        
        
      )
    ),
    tabPanel(
      "Regression R/EN/L",
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
                                choices = c("lambda.min","lambda.1se")),
                 textOutput(NS(id,"lambda_value_min")),
                 textOutput(NS(id,"lambda_value_1se")),
                 textOutput(NS(id,"model_error"))
                 
                 ),
          column(3,
                 dataTableOutput(NS(id,"coef_data"))
                 
                 ),
          column(6,
                 plotlyOutput(NS(id,"coef_lambda_plot")),
                 
                 plotOutput(NS(id,"lambda_error_plot"))
                 
                 
                 )
          
          
          
        ),
        width = 9
        
      )
    )
    
  )
}


#Server Modules
data_prep_ml_server <- function(id,Xproj) {
  moduleServer(id,function(input,output,session){
    
    #Response, obtain from msigdb
    
    # msigdb_gene_sets =  reactive({readRDS(paste0("projects/", "msigdb_gene_sets", ".rds"))})
    
    msigdb_gene_sets =  reactive({zstd_unserialize(readRDS(paste0("genesets/", "compressed_msigdb_gene_sets", ".rds")))})
    
    
    
    help_dataprep = reactive({
      
      data.frame(
        
        element = paste0("#", session$ns(c(NA,NA, "response_prep_method", "predictor_prep_method")) ),
        intro = paste(c(
          "Welcome to the TCGEx ML module. In this tab you can construct and choose dependent and independent variables (response & predictors)
          which will be used to construct a model via regularized regression.",
          
          "You can filter out transcriptomes by specifying maximum zero-count percantage toggle prior to the variable selection procedure.
          This way you can eliminate the transcripts that have above-threshold zero-count.",  
          
          "From this panel, you can enter response variable determinants either by searching for available MsigDB collections or entering them manually.
          You can also upload your set of interest from the manual selection menu.",
          
          "From this panel, you can enter predictor variables either by searching for available MsigDB collections or entering them manually.
          You can also upload your set of interest or use all available transcripts from the manual selection menu."
          
        ))
      )
      
      
      
    })
    
    observeEvent(input$mldata_help, {
      introjs(session, options = list(steps = help_dataprep() ) )
    })
    
    
    response_msigdb_data <- reactive({
      if (input$msigdb_setnames_response == "c5") {
        msigdb_response <- msigdb_gene_sets()[["c5"]]
        msigdb_response <- msigdb_response[[input$go_term_response]]
        updateSelectizeInput(session, 'msigdb_gene_set_response', choices = msigdb_response$names$gs_name,server = TRUE)
        
        response_m <- as.data.frame(msigdb_response$data)
      } else if (input$msigdb_setnames_response == "c7") {
        
        msigdb_response <- msigdb_gene_sets()[["c7"]]
        msigdb_response <- msigdb_response[[input$immuno_response]]
        
        updateSelectizeInput(session, 'msigdb_gene_set_response', choices = msigdb_response$names$gs_name,server = TRUE)
        
        response_m <- as.data.frame(msigdb_response$data)
      } else {
        
        msigdb_response <- msigdb_gene_sets()[[input$msigdb_setnames_response]]
        updateSelectizeInput(session, 'msigdb_gene_set_response', choices = msigdb_response$names$gs_name,server = TRUE)
        
        response_m <- as.data.frame(msigdb_response$data)
      }
      response_m
    })
    
    response_msigdb_genes <- reactive({
      filter(response_msigdb_data(), gs_name == input$msigdb_gene_set_response) %>% select(gene_symbol)
    })
    
    
    
    #Predictor, obtain from msigdb.
    predictor_msigdb_data <- reactive({
      if (input$msigdb_setnames_predictor == "c5") {
        
        msigdb_predictor <- msigdb_gene_sets()[["c5"]]
        msigdb_predictor <- msigdb_predictor[[input$go_term_predictor]]
        
        updateSelectizeInput(session, 'msigdb_gene_set_predictor', choices = msigdb_predictor$names$gs_name,server = TRUE)
        
        predictor_m <- as.data.frame(msigdb_predictor$data)
        
      } else if (input$msigdb_setnames_predictor == "c7") {
        
        msigdb_predictor <- msigdb_gene_sets()[["c7"]]
        msigdb_predictor <- msigdb_predictor[[input$immuno_predictor]]
        
        updateSelectizeInput(session,'msigdb_gene_set_predictor', choices = msigdb_predictor$names$gs_name,server = TRUE)
        
        predictor_m <- as.data.frame(msigdb_predictor$data)
      } else {
        
        msigdb_predictor <- msigdb_gene_sets()[[input$msigdb_setnames_predictor]]
        updateSelectizeInput(session,'msigdb_gene_set_predictor', choices = msigdb_predictor$names$gs_name,server = TRUE)
        
        predictor_m <- as.data.frame(msigdb_predictor$data)
      }
      predictor_m
    })
    
    predictor_msigdb_genes <- reactive({
      filter(predictor_msigdb_data(), gs_name == input$msigdb_gene_set_predictor) %>% select(gene_symbol)
    })
    
    
    # create a list.
    
    available_genelist <- reactive({zero_adjuster(mdata = Xproj$a(),max_zero_percent = input$max_zero_percent, subgroup = "alltranscripts")})
    available_cibersort <- reactive({zero_adjuster(mdata = Xproj$a(),max_zero_percent = input$max_zero_percent, subgroup = "cibersort")})
    
    
    observeEvent(input$max_zero_percent, {
      cibersort_metrics <- available_cibersort()$available
      genelist <- available_genelist()$available
      updateSelectizeInput(session,'gene_list_response', choices = genelist, server = TRUE)
      updateSelectizeInput(session, 'gene_list_predictor', choices = genelist, server = TRUE)
      updateSelectizeInput(session,'cibersort_response_var', choices = cibersort_metrics, server = TRUE)
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
          validate(need(ext == "txt", "Please upload a txt file"))
          list_r = read.table(file$datapath, header = TRUE, sep = "", dec = ".")
        }
        
      } else if (input$response_prep_method == "cibersort") {
        list_r = cibersort_list()
      }
      colnames(list_r) = "gene_symbol"
      unique(list_r)
    })
    
    
    predictor_det <- reactive({
      if (input$predictor_prep_method == "msigdb") {
        list_p <- predictor_msigdb_genes()
      } else if (input$predictor_prep_method == "gene_list") {
        if (input$obtain_predictor == "gene_enter") {
          list_p <- gene_list_selected_df_predictor()
        } else if (input$obtain_predictor == "allmiRNA_aspredictor") {
          mir_list <- colnames(select(Xproj$a(),starts_with("hsa")))
          list_p <- as.data.frame(mir_list)
          
        } else {
          file <- input$predictor_set_file
          ext <- tools::file_ext(file$datapath)
          req(file)
          validate(need(ext == "txt", "Please upload a txt file"))
          list_p <- read.table(file$datapath, header = FALSE, sep = "", dec = ".")
          
        }
      }
      colnames(list_p) <- "gene_symbol"
      unique(list_p)
    })
    
    predictor_missing <- reactive({
      genelist <- predictor_det()$gene_symbol
      is.exist <- genelist %in% colnames(Xproj$a())
      predictor_missing <- predictor_det()[which(is.exist == "FALSE"),]
      predictor_missing
    })
    
    response_missing <- reactive({
      genelist <- response_det()$gene_symbol
      is.exist <- genelist %in% colnames(Xproj$a())
      response_missing <- response_det()[which(is.exist == "FALSE"),]
      response_missing
    })
    output$validation_message_predictor <- renderPrint({
      predictor_missing()
    })
    
    output$validation_message_response <- renderPrint({
      response_missing()
    })
    
    
    clean_response_set = reactive({
      if (input$response_prep_method == "cibersort") {
        response_det <- response_det()
      } else {
        missing_genes <-  as.vector(response_missing())
        response_det <- response_det()
        for (i in missing_genes) {
          response_det <- filter(response_det, gene_symbol != i)
        }
      }
      response_det
    })
    
    clean_predictor_set = reactive({
      missing_genes <- as.vector(predictor_missing())
      predictor_det <- predictor_det()
      for (i in missing_genes) {
        predictor_det <- filter(predictor_det, gene_symbol != i)
      }
      predictor_det
    })
    
    output$response_set <- renderPrint({clean_response_set()$gene_symbol}) #response genes 
    output$predictor_set <- renderPrint({clean_predictor_set()$gene_symbol }) #predictor genes
  
    reg_data <- reactive({
      select(Xproj$a(), clean_response_set()[,1], clean_predictor_set()[,1]) %>%
        mutate(response = select(., clean_response_set()[,1]) %>% rowMeans()) %>%
        select(response, clean_predictor_set()[,1]) %>% 
        na.omit()
    })
    
    
    return(reg_data)
  } )
}
ml_main_server <- function(id,regress_data,Xproj) {
  moduleServer(id,function(input,output,session) {
    
    help_regression = reactive({
      
      if (input$regression_workflow == "create_model") {
        
        return(
          data.frame(
          element = paste0("#", session$ns(c(NA, "regression_workflow","user_alpha", "lambda_for_coef + .selectize-control")) ),
          intro = paste(c(
            "In this tab, you can determine necessary parameters to conduct regularized regression with the variables that you have already specified.",
            "You can either split the data into train and test sets or use all of it to create the model. If you choose splitting option, you will be able to
            make a prediction with your model and examine mean-squared-error.",
            "You can choose the method of regularized regression from the alpha parameter slider control. ɑ=0 for  Ridge, ɑ=1 for Lasso and 0<ɑ<1 for Elastic-Net regression.",
            "You can choose the lambda value at which the variable coefficients of the model are displayed."
            
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
    
    trows <- reactive( { 
      nrowd <- nrow(regress_data())
      set.seed(30)
      trowsx <- sample(1:nrowd, (input$train_percentage*nrowd)/100 ,replace = FALSE)
      trowsx
    })
    
    
    cvfit <- eventReactive(input$run, {
      if (input$regression_workflow == "create_model") {
        response_var <- as.matrix(select(regress_data(), response)) 
        predictor_var <-  as.matrix(select(regress_data(), -c("response"))) 
        foldid <- sample(1:10,size = length(response_var), replace = TRUE)
        fitty <- cv.glmnet(predictor_var, response_var, foldid = foldid, alpha = input$user_alpha)
      } else if (input$regression_workflow == "ctest_model") {
        train_regress_data <- regress_data()[trows(),]
        response_var_train <- as.matrix(select(train_regress_data, response)) 
        predictor_var_train <-  as.matrix(select(train_regress_data, -c("response"))) 
        foldid <- sample(1:10,size = length(response_var_train), replace = TRUE)
        fitty <- cv.glmnet(predictor_var_train, response_var_train, foldid = foldid, alpha = input$user_alpha)
      }
      fitty
    })
    
    
    prediction_error <- eventReactive(input$lambda_for_prediction, {
      if(input$regression_workflow == "ctest_model"){
        
        test_regress_data <- regress_data()[-trows(),]
        response_var_test <- as.matrix(select(test_regress_data, response)) 
        predictor_var_test <-  as.matrix(select(test_regress_data, -c("response"))) 
        prediction <- predict(cvfit(), s = input$lambda_for_prediction, newx = predictor_var_test)
        
        error <- mean((response_var_test - prediction)^2)
        
      }
      error
    })
    
    coef_data <- reactive({
      c <- as.matrix(coef(cvfit(), s = input$lambda_for_coef))
      c <-  c[2:length(rownames(c)), ]
      as.data.frame(c)})
    
    
    output$lambda_error_plot <- renderPlot({plot(cvfit(),xvar = "lambda", label = TRUE)})
    
    output$coef_lambda_plot <- renderPlotly({
      
      fit <- cvfit()$glmnet.fit
      l1se <- cvfit()$lambda.1se
      lmin <- cvfit()$lambda.min
      
      lam <- fit$lambda %>%
        as.data.frame() %>%
        mutate(penalty = fit$a0 %>% names()) %>%
        dplyr:: rename(lambda = ".")
      results <- fit$beta %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        gather(penalty, coefficients, -rowname) %>%
        left_join(lam)
      
      result_labels <- results %>%
        filter(lambda == min(lambda))
      g=ggplot() +
        geom_line(data = results, aes(lambda, coefficients, group = rowname, color = rowname), show.legend = FALSE) +
        scale_x_log10() +
        geom_text(data = result_labels, aes(lambda, coefficients, label = rowname, color = rowname), nudge_x = .15,nudge_y = 0.0035,
                  show.legend = FALSE)+
        geom_vline(xintercept = lmin) +
        geom_vline(xintercept = l1se) 

      g <- g + theme(legend.title = element_blank())
      g <- g + theme(legend.position='none')
      
      
      ggplotly(g)
      
    })
    
    output$coef_data <- renderDataTable({
      coef_data()
    })
    
    output$lambda_value_min <- renderText({
      paste("lambda.min: ", cvfit()$lambda.min) })
    
    output$lambda_value_1se <- renderText({
      paste("lambda.1se: ", cvfit()$lambda.1se)})
    
    output$model_error <- renderText({
      paste(" MSE Error: ", prediction_error() ) 
    })
    
    output$download_coef <- downloadHandler(
      filename <-  function() {
        paste(input$lambda_for_coef, ".csv", sep = "")
      },
      content <-  function(file) {
        write.csv(coef_data(), file, row.names = FALSE)
      }
    )
    
  })
}