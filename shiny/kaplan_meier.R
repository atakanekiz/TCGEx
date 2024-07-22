library(survival)
library(survminer)
library(ggpubr)
library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)

##KM survival modularization

km_ui <- function(id, label, choices) {
  
  ns <- NS(id)
  
  tagList(
    
    useShinyalert(),
    
    ui <- fluidPage(
      
      sidebarPanel(
      
        selectizeInput(inputId = ns("km_samptyp"), 
                       multiple=T,
                       label = "1. Select sample types",
                       choices=NULL, # will be updated dynamically
                       options=list(placeholder = "eg. Primary solid tumor")),
        
        hr(),
        
        
        selectizeInput(inputId = ns("km_feat"),
                       label = "2. Select feature",
                       choices = NULL, # will be updated dynamically
                       options=list(placeholder = "eg. IFNG",
                                    plugins = list('restore_on_backspace'))),
        
        conditionalPanel(
          
          condition = "output.km_var_status",
          ns=ns, 
          numericInput(inputId = ns("lo_cutoff"), "Low cutoff percent", 50, 1, 100),
          
          numericInput(inputId = ns("hi_cutoff"), "High cutoff percent", 50, 1, 100),
          
          checkboxInput(inputId = ns("keep_mid"), "Keep samples in the middle?", value = F)
          
        ),
        
        conditionalPanel(
          
          condition = "output.km_var_status2",
          ns=ns, 
          selectizeInput(inputId = ns ("sel_feat_meta_groups"), "Select feature subsets", multiple=T,
                         choices = NULL, # will be updated dynamically
                         options = list(placeholder = "eg. male, female")),
          
          
        ),
        
        hr(),
        
        conditionalPanel(
          
          condition = " output.proj_length_KM > '1' ",   ## Conditions of conditional panel does not working with the reactives from other modules, so I created Xproj$b as reactive shows the length of cancers, Since conditions needs an input or output result, I made a output that takes info from selectdata Xproj$b reactive(Cagatay)
          ns=ns,
          radioButtons(inputId = ns ("choose_KM"), "Categorization rule",
                       selected = "Aggregated",
                       c("Aggregated", "Per cancer type")),
          
          
          
        ),
        
        span(style="color:#3382FF",
             
             selectizeInput(inputId = ns ("km_covar") , "3. Select covariate (optional)",
                            choices=NULL, # will be updated dynamically
                            options=list(placeholder = "eg. meta.gender",
                                         plugins = list('restore_on_backspace'))
             ),
             
             conditionalPanel(
               
               condition="output.km_covar_status",
               ns=ns, 
               numericInput(inputId = ns("lo_cutoff_covar"), "Low cutoff percent", 50, 1, 100),
               
               numericInput(inputId = ns("hi_cutoff_covar"), "High cutoff percent", 50, 1, 100),
               
               checkboxInput(inputId = ns("keep_mid_covar"), "Keep samples in the middle?", value = F)
               
               
             ),
             
             conditionalPanel(
               
               condition = "output.km_covar_status2",
               ns=ns,
               selectizeInput(inputId = ns("sel_covar_meta_groups"), "Select feature subsets", multiple=T,
                              choices = NULL, # will be updated dynamically
                              options = list(placeholder = "eg. male, female")),
               
               
             )
             
        ),
        
        actionBttn(inputId = ns("km_run"), 
                   label = "Analyze",
                   style = "unite",
                   block = TRUE,
                   color = "primary"),
        br(),
        
        downloadButton(ns("downloadPlot3"), "Download Kaplan Meier Plot", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
        br(),
        br(),
        #help section UI
        
        introjsUI(),
        actionButton(ns("KM_help"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
        
        textOutput(ns("filewarning_two")) ,     
        
        width = 3
        
    ),
      
      mainPanel(
       
        
        plotOutput(outputId = ns("km_plot")),
        
        
        
        conditionalPanel(
          
          h3("Graphing options"),
          
          condition = "output.km_plot" ,
          ns=ns, 
          column(4,
                 
                 checkboxInput(inputId = ns("km_risk"), label="Show risk table", F),
                 checkboxInput(inputId = ns("km_pval"), label = "Show p-value", T),
                 checkboxInput(inputId = ns("km_confint"), label = "Show conf. intervals", F),
          ),
          
          column(4,
                 selectInput(inputId = ns ("km_medline"), label = "Mark median survival", 
                             selected = "hv",
                             choices = c("None" = "none",
                                         "Horizontal-Vertical" = "hv",
                                         "Vertical" = "v",
                                         "Horizontal" = "h")),
                 
                 selectInput(inputId = ns ("km_pal"), "Select colors", multiple = F, 
                             selectize = T, 
                             selected = "npg",
                             choices = list(`Pre-made palettes` = list(
                               "Nature"="npg",
                               "Science"= "aaas",
                               "Lancet"= "lancet",
                               "New Engl J Med" = "nejm",
                               "J Clin Onc"="jco",
                               "JAMA" = "jama",
                               "Int Genomics Viewer (IGV)" = "igv",
                               "UCSC Genome Browser"="ucscgb",
                               "Default" = NULL),
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
                                 "Brown-Green" = "BrBG"
                               )) 
                 ),
                 actionButton(inputId = ns("km_update"), "update graph", icon("refresh")),      
          ),
          
          column(4,
                 
                 sliderInput(inputId = ns("km_xlim"), label="Define endpoint (days)", value = 6000, min = 10, max=10000),
                 sliderInput(inputId = ns("km_breaktime"), label="Break time invervals", value=500, min=100, max = 1000, step = 100)),
          h3("Survival fit summary")
          
        ),
        
        
        
        verbatimTextOutput(outputId = ns("km_text"))
        
      )
    )
  )
} 





km_server <- function(id,Xproj) {
  moduleServer(id,function(input, output, session){
    
    ns <- session$ns
    
    # v <- eventReactive(Xproj$fileInfost(), {
    #   
    # 
    #   shinyalert("Warning!", "To perform Kaplan-Meier Survival analysis, the data you upload must contain columns containing survival information such as 'vital_status' and 'days _to_event'.",html = FALSE,imageUrl = "",closeOnEsc = TRUE,
    #                          closeOnClickOutside = TRUE)
    #   
    # })
    # 
    # output$filewarning_two <- renderText({
    #   v()
    # })
    # 

    
    observeEvent(Xproj$fileInfost(),{

    output$filewarning_two <- renderText({

      if (!is.null(Xproj$fileInfost())) {
        shinyalert("Warning!", "To perform Kaplan-Meier Survival analysis, the data you upload must contain columns containing survival information such as 'vital_status' and 'days _to_event'.",html = FALSE,imageUrl = "",closeOnEsc = TRUE,
                   closeOnClickOutside = TRUE)

        # tags$style(HTML(".km-shinyalert-12a326a245e24d08919ea6ee5a3110e7 { display: none !important; }"))
      }
    }) })
    
    KM_steps <- reactive({
      
      if(Xproj$cancer_length() ==1 & input$km_feat != "" | !is.null(Xproj$fileInfost())) {
        
        return(
          
          data.frame(
            
            element = paste0("#", session$ns(c(NA, "km_samptyp + .selectize-control", "km_feat + .selectize-control","lo_cutoff", "hi_cutoff", "keep_mid", "km_covar + .selectize-control", "km_risk", "km_pval", "km_confint", "km_medline + .selectize-control", "km_pal + .selectize-control", "km_xlim", "km_breaktime"))),
            
            intro = paste(c(
              "This is Kaplan-meier (KM) survival analysis module. Here, you can examine how different data subsets differ in terms of survival. You can define data subsets by categorizing gene expression at desired cutoffs and/or use metadata features that are already categorical. The log-rank test p-value is reported on the graph and the survival model fit is shown to provide further details about the analysis. Continue tutorial to learn how to use the module.",
              "You can select the sample types (eg. primary and/or metastatic) to tailor the analysis to your needs.",
              "KM analysis is performed between groups of data. You can select genes, miRNAs, or clinical meta data features here. If your selection is a categorical data type (eg. patient gender, tumor subtype), you will be asked to select which subsets to be included in the analysis. If your selection is a numerical data type (eg. gene expression), you will be asked to define quantile cutoffs to categorize gene expression as 'high' and 'low'",
              
              "Define the quantile cutoff for the low expression group. 50 (default) means that the samples expressing the gene of interest at lower levels than the median value will be categorized as 'low'. Setting this value to 25, for instance, will categorize the bottom 25% of the data as the low group.",
              
              "Define the quantile cutoff for the high expression group. 50 (default) means that the samples expressing the gene of interest at the median value or higher will be categorized as 'high'. Setting this value to 25, for instance, will categorize the top 25% of the data as the high group.",
              
              "If your numeric categorization results in three groups (ie low, middle, high), you can hide (default) or show the middle group in the graph",
              "You can add a covariate into the analysis and define data subsets as described before",
              "You can show the risk table by clicking this box. A table will be added below the KM curves showing the number of surviving patients at different timepoints. <i>(This option will be visible after plotting the graph.)</i>",
              "You can show the logrank p-value on the graph by clicking this box. If there are more than two groups in the analysis, the p-value is calculated by testing the null hypothesis that all the samples come from populations with identical survival. You can select two specific data subsets to show pair-wise p-values. <i>This option will be visible after plotting the graph.)</i>",
              "Confidence interval bands can be added to the graph by clicking this box. <i>This option will be visible after plotting the graph.)</i>",
              "You can plot dashed lines to highlight median survival in data subsets. <i>This option will be visible after plotting the graph.)</i> ",
              "You can change the color palette of the graph here. <i>This option will be visible after plotting the graph.)</i>",
              "You can manually change the plotted time interval here. This does not affect the results of survival analysis. <i>This option will be visible after plotting the graph.)</i>",
              "The breaks on the x-axis can be changed here. <i>This option will be visible after plotting the graph.)</i>"
            ))
            
          )
          
        )
      } else if (Xproj$cancer_length() > 1 & input$km_feat != "") {
        return(
          
          data.frame(
            
            element = paste0("#", session$ns(c(NA, "km_samptyp + .selectize-control", "km_feat + .selectize-control", "choose_KM", "lo_cutoff", "hi_cutoff", "keep_mid", "km_covar + .selectize-control", "km_risk", "km_pval", "km_confint", "km_medline + .selectize-control", "km_pal + .selectize-control", "km_xlim", "km_breaktime"))),
            
            intro = paste(c(
              "This is Kaplan-meier (KM) survival analysis module. Here, you can examine how different data subsets differ in terms of survival. You can define data subsets by categorizing gene expression at desired cutoffs and/or use metadata features that are already categorical. The log-rank test p-value is reported on the graph and the survival model fit is shown to provide further details about the analysis. Continue tutorial to learn how to use the module.",
              "You can select the sample types (eg. primary and/or metastatic) to tailor the analysis to your needs.",
              "When analyzing multiple cancer types together, you can categorize gene expression across the aggregated dataset as a whole or separately for each cancer types.",
              "KM analysis is performed between groups of data. You can select genes, miRNAs, or clinical meta data features here. If your selection is a categorical data type (eg. patient gender, tumor subtype), you will be asked to select which subsets to be included in the analysis. If your selection is a numerical data type (eg. gene expression), you will be asked to define quantile cutoffs to categorize gene expression as 'high' and 'low'",
              
              "Define the quantile cutoff for the low expression group. 50 (default) means that the samples expressing the gene of interest at lower levels than the median value will be categorized as 'low'. Setting this value to 25, for instance, will categorize the bottom 25% of the data as the low group.",
              
              "Define the quantile cutoff for the high expression group. 50 (default) means that the samples expressing the gene of interest at the median value or higher will be categorized as 'high'. Setting this value to 25, for instance, will categorize the top 25% of the data as the high group.",
              
              "If your numeric categorization results in three groups (ie low, middle, high), you can hide (default) or show the middle group in the graph",
              "You can add a covariate into the analysis and define data subsets as described before",
              "<i>(This option will be visible after plotting)</i> You can show the risk table by clicking this box. A table will be added below the KM curves showing the number of surviving patients at different timepoints. ",
              "<i>(This option will be visible after plotting)</i> You can show the logrank p-value on the graph by clicking this box. If there are more than two groups in the analysis, the p-value is calculated by testing the null hypothesis that all the samples come from populations with identical survival. You can select two specific data subsets to show pair-wise p-values. ",
              "<i>(This option will be visible after plotting)</i> Confidence interval bands can be added to the graph by clicking this box.",
              "<i>(This option will be visible after plotting)</i> You can plot dashed lines to highlight median survival in data subsets.",
              "<i>(This option will be visible after plotting)</i>You can change the color palette of the graph here.",
              "<i>(This option will be visible after plotting)</i> You can manually change the plotted time interval here. This does not affect the results of survival analysis.",
              "<i>(This option will be visible after plotting)</i> The breaks on the x-axis can be changed here."
            ))
            
          )
          
        )
        
      } else if (input$km_feat == "") {
        return(
          data.frame(
            
            element = paste0("#", session$ns(c(NA, "km_samptyp + .selectize-control", "km_feat + .selectize-control","km_covar + .selectize-control"))),
            
            intro = paste(c(
              "This is Kaplan-meier (KM) survival analysis module. Here, you can examine how different data subsets differ in terms of survival. You can define data subsets by categorizing gene expression at desired cutoffs and/or use metadata features that are already categorical. The log-rank test p-value is reported on the graph and the survival model fit is shown to provide further details about the analysis. Continue tutorial to learn how to use the module. <b>Note:</b> Inputs in this interface are updated according to the user selection and new tutorial steps will be available as you perform analysis.",
              "You can select the sample types (eg. primary and/or metastatic) to tailor the analysis to your needs.",
              "KM analysis is performed between groups of data. You can select genes, miRNAs, or clinical meta data features here. If your selection is a categorical data type (eg. patient gender, tumor subtype), you will be asked to select which subsets to be included in the analysis. If your selection is a numerical data type (eg. gene expression), you will be asked to define quantile cutoffs to categorize gene expression as 'high' and 'low'",
              "You can select a covariate group in order to compare it with first feature. This is optional and not needed for single feature analysis. <i>Note: Other tutorial steps will be available when you start the analysis.</i>"
            ))
          )
        )
      }
      
      
    })
    
    
    
    
    observeEvent(input$KM_help, {
      
      introjs(session, options = list(steps = KM_steps()) )
      
    })
    
    
    
    # shiny invalidate
    observe(req(input$km_feat))
    
    panel_iv <- InputValidator$new()
    
    numeric_iv <- InputValidator$new()
    
    # Does not work in conditional panel
    numeric_iv_covar <- InputValidator$new()
    
    numeric_iv$add_rule("hi_cutoff", ~ if (input$hi_cutoff + input$lo_cutoff > 100 & !anyNA(input$lo_cutoff) & !anyNA(input$hi_cutoff)) "For proper categorization, high and low cutoffs can't exceed 100 when added together")
    numeric_iv$add_rule("lo_cutoff", ~ if (input$hi_cutoff + input$lo_cutoff > 100 & !anyNA(input$lo_cutoff) & !anyNA(input$hi_cutoff)) "For proper categorization, high and low cutoffs can't exceed 100 when added together")
    numeric_iv_covar$add_rule("hi_cutoff_covar", ~ if (input$hi_cutoff_covar + input$lo_cutoff_covar > 100 & !anyNA(input$lo_cutoff_covar) & !anyNA(input$hi_cutoff_covar)) "For proper categorization, high and low cutoffs can't exceed 100 when added together")
    numeric_iv_covar$add_rule("lo_cutoff_covar", ~ if (input$hi_cutoff_covar + input$lo_cutoff_covar > 100 & !anyNA(input$lo_cutoff_covar) & !anyNA(input$hi_cutoff_covar)) "For proper categorization, high and low cutoffs can't exceed 100 when added together")
    
    
    iv <- InputValidator$new()
    iv$add_validator(panel_iv)
    iv$add_validator(numeric_iv)
    iv$add_validator(numeric_iv_covar)
    iv$enable()
    
    
    
    
    
    
    observe({updateSelectizeInput(session, 
                                  "km_samptyp",
                                  choices = Xproj$a()$meta.definition, 
                                  server = T)})
    
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session, 
                                  "km_feat", selected="",
                                  choices = colnames(Xproj$a()), 
                                  server = T)})
    
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session,
                                  "sel_feat_meta_groups",
                                  choices = levels(as.factor(Xproj$a()[[input$km_feat]])),
                                  server = T)})
    
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session, 
                                  "km_covar", selected = "",
                                  choices = c(colnames(Xproj$a())), #projcolnames()[grepl("meta.", projcolnames())], 
                                  server = T)})
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session, 
                                  "sel_covar_meta_groups",
                                  choices = levels(as.factor(Xproj$a()[[input$km_covar]])), 
                                  server = T)})
    
    
    # Prompt conditional panel if variable is a gene
    output$km_var_status <- reactive({
      
      is.numeric(Xproj$a()[[input$km_feat]])
      
    })
    
    outputOptions(output, "km_var_status", suspendWhenHidden = FALSE)
    
    
    # Prompt conditional panel if variable is a categorical variable
    output$km_var_status2 <- reactive({
      
      is.character(Xproj$a()[[input$km_feat]]) | is.factor(Xproj$a()[[input$km_feat]])
      
    })
    
    outputOptions(output, "km_var_status2", suspendWhenHidden = FALSE)
    
    
    
    # Prompt conditional panel if any covariable is a gene
    output$km_covar_status <- reactive({
      
      is.numeric(Xproj$a()[[input$km_covar]])
      
    })
    
    outputOptions(output, "km_covar_status", suspendWhenHidden = FALSE)
    
    
    # Prompt conditional panel if variable is a meta data column
    output$km_covar_status2 <- reactive({
      
      is.character(Xproj$a()[[input$km_covar]]) | is.factor(Xproj$a()[[input$km_covar]])
      
    })
    
    outputOptions(output, "km_covar_status2", suspendWhenHidden = FALSE)
    
    output$proj_length_KM <- reactive({      #an output takes info form selectdata's Xproj$length() reactive
      
      
      Xproj$cancer_length()
      
    })
    
    outputOptions(output, "proj_length_KM", suspendWhenHidden = FALSE)  
    
    km_feat_zero_count <- reactiveValues(value = NULL)
    km_feat_na_count <- reactiveValues(value = NULL)
    km_covar_zero_count <- reactiveValues(value = NULL)
    km_covar_na_count <- reactiveValues(value = NULL)
    
    # Prepare trimmed km_dat object
    
    km_dat <- reactive({
      
      validate(
        need(input$km_samptyp, "Select data subsets"),
        need(input$km_feat, "Select feature"),
        need(input$km_run, "Please click 'Analyze' button to display KM curves"),
        if(input$km_covar != ""){
          if(class(Xproj$a()[[input$km_covar]]) %in% c("character", "factor")){
            need(input$sel_covar_meta_groups, "Select at least one covariate group")}} 
      )
      
      if(Xproj$cancer_length() ==1 | !is.null(Xproj$fileInfost())) {
        
        
        sel_cols <- c(input$km_feat, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient")
        
        # browser()
        
        dat <- Xproj$a()[, ..sel_cols]
        
        dat <- dat[meta.definition %in% input$km_samptyp, ]
        
        # De-duplicate patients. Expand on this to drop certain sample groups when duplicates are present
        dat <- dat[!duplicated(meta.patient), ]
        
        validate(
          if(!is.numeric(dat[[input$km_feat]]))  need(input$sel_feat_meta_groups, "Select at least one feature group"))
        
        if(is.numeric(dat[[input$km_feat]])){
          
          km_feat_zero_count$value <- sum(dat[[input$km_feat]] == 0, na.rm = TRUE)
          
          km_feat_na_count$value <- sum(is.na(dat[[input$km_feat]]))
  
          mid_value_feat <- ifelse(input$keep_mid, "mid", NA)
          
          dat[, (input$km_feat) := ifelse(dat[[input$km_feat]] >= quantile(dat[[input$km_feat]], (100-input$hi_cutoff)/100, na.rm = T), "high", ifelse(dat[[input$km_feat]] < quantile(dat[[input$km_feat]], input$lo_cutoff/100, na.rm = T), "low", mid_value_feat))]
          
        } else if(is.character(dat[[input$km_feat]]) | is.factor(dat[[input$km_feat]])){
          
          
          dat <- dat[get(input$km_feat) %in% input$sel_feat_meta_groups, ]
          
        }
      }
      
      
      else if (Xproj$cancer_length() > 1){
        
        validate(
          need(input$choose_KM, "Please select categorization rule"))
        
        if(input$choose_KM == "Aggregated") {
          
          sel_cols <- c(input$km_feat, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient")
          
          dat <- Xproj$a()[, ..sel_cols]
          
          dat <- dat[meta.definition %in% input$km_samptyp, ]
          
          # De-duplicate patients. Expand on this to drop certain sample groups when duplicates are present
          dat <- dat[!duplicated(meta.patient), ]
          
          validate(
            if(!is.numeric(dat[[input$km_feat]]))  need(input$sel_feat_meta_groups, "Select at least one feature group"))
          
          if(is.numeric(dat[[input$km_feat]])){
            
            km_feat_zero_count$value <- sum(dat[[input$km_feat]] == 0, na.rm = TRUE)
            km_feat_na_count$value <- sum(is.na(dat[[input$km_feat]]))

            
            mid_value_feat <- ifelse(input$keep_mid, "mid", NA)
            
            dat[, (input$km_feat) := ifelse(dat[[input$km_feat]] >= quantile(dat[[input$km_feat]], (100-input$hi_cutoff)/100, na.rm = T), "high", ifelse(dat[[input$km_feat]] < quantile(dat[[input$km_feat]], input$lo_cutoff/100, na.rm = T), "low", mid_value_feat))]
            
          } else if(is.character(dat[[input$km_feat]]) | is.factor(dat[[input$km_feat]])){
            
            
            dat <- dat[get(input$km_feat) %in% input$sel_feat_meta_groups, ]   
            
          }
        } else if (input$choose_KM == "Per cancer type") {
          
          sel_cols_sep <- c(input$km_feat, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient", "meta.project_id")
          
          dat <- Xproj$a()[, ..sel_cols_sep]
          
          dat <- dat[meta.definition %in% input$km_samptyp, ]
          
          # De-duplicate patients. Expand on this to drop certain sample groups when duplicates are present
          dat <- dat[!duplicated(meta.patient), ]
          
          validate(
            if(!is.numeric(dat[[input$km_feat]]))  need(input$sel_feat_meta_groups, "Select at least one feature group"))
          
          if(is.numeric(dat[[input$km_feat]])){
            
            km_feat_zero_count$value <- sum(dat[[input$km_feat]] == 0, na.rm = TRUE)
            km_covar_na_count$value <- sum(is.na(dat[[input$km_feat]]))

            
            mid_value_feat <- ifelse(input$keep_mid, "mid", NA)
            
            dat[, high_sep := quantile(get(input$km_feat), (100-input$hi_cutoff)/100, na.rm = T), by=c("meta.project_id")]  ##GET  ##CHANGE THE PLACE TO THE MIDDLE
            dat[, low_sep := quantile(get(input$km_feat), (input$lo_cutoff)/100, na.rm = T), by=c("meta.project_id")]
            dat[, (input$km_feat) := ifelse(dat[[input$km_feat]] >= high_sep, "high", ifelse(dat[[input$km_feat]] < low_sep, "low", mid_value_feat))]
            
          } else if(is.character(dat[[input$km_feat]]) | is.factor(dat[[input$km_feat]])){
            
            
            dat <- dat[get(input$km_feat) %in% input$sel_feat_meta_groups, ]
            
          } 
          
        }
      }
      
      if(input$km_covar != "") {
        
        
        if(Xproj$cancer_length() ==1 | !is.null(Xproj$fileInfost())) {
          
          # browser()
          
          sel_cols2 <- c(input$km_covar, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient")
          
          covar_dat <- Xproj$a()[, ..sel_cols2][meta.definition %in% input$km_samptyp,][!duplicated(meta.patient), ]
          
          dat <- left_join(dat, covar_dat, by=c("meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient"))
          

          km_covar_zero_count$value <- sum(dat[[input$km_covar]] == 0, na.rm=T)
          km_covar_na_count$value <- sum(is.na(dat[[input$km_covar]]))
          
          if(is.numeric(dat[[input$km_covar]])){
            
            km_covar_zero_count$value <- sum(dat[[input$km_covar]] == 0, na.rm=T)
            km_covar_na_count$value <- sum(is.na(dat[[input$km_covar]]))

            mid_value_covar <- ifelse(input$keep_mid_covar, "mid", NA)
            
            dat[, (input$km_covar) := ifelse(dat[[input$km_covar]] >= quantile(dat[[input$km_covar]], (100-input$hi_cutoff_covar)/100, na.rm = T), "high", ifelse(dat[[input$km_covar]] < quantile(dat[[input$km_covar]], input$lo_cutoff_covar/100, na.rm = T), "low", mid_value_covar))]
            
          } else if(is.character(dat[[input$km_covar]]) | is.factor(dat[[input$km_covar]])){
            
            
            dat <- dat[get(input$km_covar) %in% input$sel_covar_meta_groups, ]
          
            
          }    
        } else if (Xproj$cancer_length() > 1 ) { 
          
          
          if(input$choose_KM == "Aggregated") {
            sel_cols2 <- c(input$km_covar, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient")
            
            covar_dat <- Xproj$a()[, ..sel_cols2][meta.definition %in% input$km_samptyp,][!duplicated(meta.patient), ]
            
            dat <- left_join(dat, covar_dat, by=c("meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient"))
            
            km_covar_zero_count$value <- sum(dat[[input$km_covar]] == 0, na.rm=T)
            km_covar_na_count$value <- sum(is.na(dat[[input$km_covar]]))
            
            if(is.numeric(dat[[input$km_covar]])){
              
              km_covar_zero_count$value <- sum(dat[[input$km_covar]] == 0, na.rm=T)
              km_covar_na_count$value <- sum(is.na(dat[[input$km_covar]]))

              
              mid_value_covar <- ifelse(input$keep_mid_covar, "mid", NA)
              
              dat[, (input$km_covar) := ifelse(dat[[input$km_covar]] >= quantile(dat[[input$km_covar]], (100-input$hi_cutoff_covar)/100, na.rm = T), "high", ifelse(dat[[input$km_covar]] < quantile(dat[[input$km_covar]], input$lo_cutoff_covar/100, na.rm = T), "low", mid_value_covar))]
              
            } else if(is.character(dat[[input$km_covar]]) | is.factor(dat[[input$km_covar]])){
              
              
              dat <- dat[get(input$km_covar) %in% input$sel_covar_meta_groups, ]
              
            }
            
          } else if (input$choose_KM == "Per cancer type") {
            sel_cols2 <- c(input$km_covar, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient", "meta.project_id")
            
            covar_dat <- Xproj$a()[, ..sel_cols2][meta.definition %in% input$km_samptyp,][!duplicated(meta.patient), ]
            
            dat <- left_join(dat, covar_dat, by=c("meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient", "meta.project_id"))
            

            km_covar_zero_count$value <- sum(dat[[input$km_covar]] == 0, na.rm=T)
            km_covar_na_count$value <- sum(is.na(dat[[input$km_covar]]))
            
            if(is.numeric(dat[[input$km_covar]])){
              
              km_covar_zero_count$value <- sum(dat[[input$km_covar]] == 0, na.rm=T)
              km_covar_na_count$value <- sum(is.na(dat[[input$km_covar]]))

              mid_value_covar <- ifelse(input$keep_mid_covar, "mid", NA)
              
              #dat[, (input$km_covar) := ifelse(dat[[input$km_covar]] >= quantile(dat[[input$km_covar]], (100-input$hi_cutoff_covar)/100, na.rm = T), "high", ifelse(dat[[input$km_covar]] < quantile(dat[[input$km_covar]], input$lo_cutoff_covar/100, na.rm = T), "low", mid_value_covar))]
              
              dat[, high_sep_covar := quantile(get(input$km_covar), (100-input$hi_cutoff_covar)/100, na.rm = T), by=c("meta.project_id")]  ##GET  ##CHANGE THE PLACE TO THE MIDDLE
              dat[, low_sep_covar := quantile(get(input$km_covar), (input$lo_cutoff_covar)/100, na.rm = T), by=c("meta.project_id")]
              dat[, (input$km_covar) := ifelse(dat[[input$km_covar]] >= high_sep_covar, "high", ifelse(dat[[input$km_covar]] < low_sep_covar, "low", mid_value_covar))]
              
            } else if(is.character(dat[[input$km_covar]]) | is.factor(dat[[input$km_covar]])){
              
              
              dat <- dat[get(input$km_covar) %in% input$sel_covar_meta_groups, ]
              
              
            }
            
          }
        }
      }
      
      
      dat
    })
    
    
    
    
    
    
    # KM analysis engine (fit object)
    km_form <- reactive({
      
      survform <- paste("Surv(meta.days_to_event, meta.vital_status) ~",
                        input$km_feat)
      
      if(input$km_covar != "") survform <- paste(survform, "+", input$km_covar)
      
      survform <- as.formula(survform)
      
      survform
      
      
    })
    
    val2 <- reactiveValues()
    
    # KM analysis output (text)
    
    output$km_plot <- renderPlot({
      
      
      km_results <- surv_fit(km_form(), dat=km_dat())
      
      
      pl2 <- survminer::ggsurvplot(km_results, data = km_dat(),
                                   pval = input$km_pval,
                                   pval.method = input$km_pval,
                                   tables.height= 0.35,
                                   risk.table = input$km_risk, 
                                   conf.int = input$km_confint,
                                   surv.median.line = input$km_medline,
                                   break.time.by = input$km_breaktime,  
                                   legend="bottom",
                                   xlim=c(0, input$km_xlim),
                                   palette = input$km_pal,
                                   title = input$km_feat,
                                   font.x=18, font.y=18, font.tickslab = 18)
      
      
      val2$pl2 <- pl2
      
      if(length(levels(as.factor(km_results$strata))) > 4) {  pl2$plot <- pl2$plot + guides(color = guide_legend(nrow = 3))
      } else { pl2$plot <- pl2$plot + guides(color = guide_legend(nrow = 2))}
      
      print(pl2)
      
    })
    
    # close observe event
    
    
    output$km_text <- renderPrint({
      
      req(input$km_feat)
      req(input$km_samptyp)
      
      validate(
        if(!is.numeric(Xproj$a()[[input$km_feat]]))  need(input$sel_feat_meta_groups, "Select at least one feature group"),            
        if(input$km_covar != ""){
          if(class(Xproj$a()[[input$km_covar]]) %in% c("character", "factor")){
            need(input$sel_covar_meta_groups, "Select at least one covariate group")}} )
      
      
      cat(paste("Survival formula:\n"))
      print(km_form())
      
      
      fit_txt <- survfit(km_form(), data=km_dat())
      survdiff_txt <- survdiff(km_form(), data = km_dat())
      
      
      print(fit_txt)
      
      writeLines("\n\n\nTest survival curve differences\n")
      
      print(survdiff_txt)
      
      writeLines("\n\n\n")
      
      print(paste("Number of samples with zero counts of the numeric feature:", km_feat_zero_count$value))
      
      writeLines("\n")
      
      print(paste("Number of samples lacking information on selected feature:", km_feat_na_count$value))
      
      writeLines("\n")
      
      if(input$km_covar != ""){
        
        print(paste("Number of samples with zero counts of the numeric covariate:", km_covar_zero_count$value))
        
        writeLines("\n")
        
        print(paste("Number of samples lacking information on selected covariate:", km_covar_na_count$value))
        
      }
    })
    output$downloadPlot3 <- downloadHandler(
      filename = function() {
        paste("Kaplan_meier.png")
      },
      content = function(file) {
        
        png(file, width = 800, height = 600)
        print(val2$pl2)
        dev.off()
        
      })
  })
}
