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
    
    ui <- fluidPage(
      sidebarPanel(
        
        
        selectizeInput(inputId = ns ("km_samptyp"), 
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
          
          condition = " output.proj_lenght_KM > '1' ",   ## Conditions of conditional panel does not working with the reactives from other modules, so I created Xproj$b as reactive shows the lenght of cancers, Since conditions needs an input or output result, I made a output that takes info from selectdata Xproj$b reactive(Cagatay)
          ns=ns,
          radioButtons(inputId = ns ("choose_KM"), "Distribution type:",
                       selected = character(0),
                       c("All", "Separate")),
          
          
          
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
                   label = "Plot KM curves",
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
                 sliderInput(inputId = ns("km_breaktime"), label="Break time invervals", value=1000, min=100, max = 1000, step = 200)),
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
    
    KM_steps <- reactive({
      
      if(Xproj$cancer_lenght() ==1 & input$km_feat != "") {
        
        return(
      
      data.frame(
        
        element = paste0("#", session$ns(c(NA, "km_samptyp + .selectize-control", "km_feat + .selectize-control","lo_cutoff", "hi_cutoff", "keep_mid", "km_covar + .selectize-control", "km_risk", "km_pval", "km_confint", "km_medline + .selectize-control", "km_pal + .selectize-control", "km_xlim", "km_breaktime"))),
        
        intro = paste(c(
          "This is Kaplan-meier (KM) Survival Analysis app.Press the buttons to learn features of the app.",
          "You can select the sample types (primary solid tumors etc.) in order to target the data subsets.",
          "You can select the features (thousands of genes, miRNAs and clinical metadata). Numerical values such as gene expression values will be shown as low/high and the mean will be calculated eith the selected low/high percentages. Categorical choices will be shown as levels.",
          "A probability distribution or sorted data can be divided into equal pieces using quantile values. A q-quantile divides sorted data into q components in general.",
          "A probability distribution or sorted data can be divided into equal pieces using quantile values. A q-quantile divides sorted data into q components in general.",
          "Normally, samples in the middle are removed from the data. However if you check the box, they will be kept.",
          "You can select a covariate group in order to compare it with first feature.",
          "Risk table shows the samples that had not yet experience the event of interest, in that case, death. You can see the values both for high/low options.",
          "The likelihood that your data would have occurred under the null hypothesis of your statistical test is expressed as a number known as a p-value, or probability value. This button shows the p value of the analysis on the graph.",
          "The confidence interval is the range of values that, if you repeated your experiment or resampled the population in the same manner, you would anticipate your estimate to fall within a specific proportion of the time. You can observe the confidence interval for both high-low stratas.",
          "With the mark median survival option, you can select the way of median mark on the graph, vertical, horizontal, both or none.",
          "With the select colors option you can select the colors you desire. Color palettes for journals etc.",
          "With the Define endpoint (days), you can select the endpoint from 10 to 10.000 days.",
          "With the breaktime intervals (days), you can select the breaktime intervals, for example, if you select 500, it makes the breaktime intervals as 500-1000-..."
        ))
        
      )
      
    )
      } else if (Xproj$cancer_lenght() > 1 & input$km_feat != "") {
        return(
          
          data.frame(
            
            element = paste0("#", session$ns(c(NA, "km_samptyp + .selectize-control", "km_feat + .selectize-control", "choose_KM", "lo_cutoff", "hi_cutoff", "keep_mid", "km_covar + .selectize-control", "km_risk", "km_pval", "km_confint", "km_medline + .selectize-control", "km_pal + .selectize-control", "km_xlim", "km_breaktime"))),
            
            intro = paste(c(
              "This is Kaplan-meier (KM) Survival Analysis app.Press the buttons to learn features of the app.",
              "You can select the sample types (primary solid tumors etc.) in order to target the data subsets.",
              "You can select the features (thousands of genes, miRNAs and clinical metadata). Numerical values such as gene expression values will be shown as low/high and the mean will be calculated eith the selected low/high percentages. Categorical choices will be shown as levels.",
              "In the case of several cancer types, if you select the 'All' calculations of high/low separation will be performed without separating the patients with their cancer types. However, if you select the 'Separate', calculations will be performed cancer type specific. For example, there will be different values for high/low separations for different cancers. ", 
              "A probability distribution or sorted data can be divided into equal pieces using quantile values. A q-quantile divides sorted data into q components in general. Bottom percentile will be analyzed here.",
              "Top percentile will be performed here for the high.",
              "Normally, samples in the middle are removed from the data. However if you check the box, they will be kept.",
              "You can select a covariate group in order to compare it with first feature.",
              "Risk table shows the samples that had not yet experience the event of interest, in that case, death. You can see the values both for high/low options.",
              "The likelihood that your data would have occurred under the null hypothesis of your statistical test is expressed as a number known as a p-value, or probability value. This button shows the p value of the analysis on the graph.",
              "The confidence interval is the range of values that, if you repeated your experiment or resampled the population in the same manner, you would anticipate your estimate to fall within a specific proportion of the time. You can observe the confidence interval for both high-low stratas.",
              "With the mark median survival option, you can select the way of median mark on the graph, vertical, horizontal, both or none.",
              "With the select colors option you can select the colors you desire. Color palettes for journals etc.",
              "With the Define endpoint (days), you can select the endpoint from 10 to 10.000 days.",
              "With the breaktime intervals (days), you can select the breaktime intervals, for example, if you select 500, it makes the breaktime intervals as 500-1000-..."
            ))
            
          )
          
        )
        
      } else if (input$km_feat == "") {
        return(
          data.frame(
            
            element = paste0("#", session$ns(c(NA, "km_samptyp + .selectize-control", "km_feat + .selectize-control","km_covar + .selectize-control"))),
            
            intro = paste(c(
            "This is Kaplan-meier (KM) Survival Analysis app.Press the buttons to learn features of the app.",
            "You can select the sample types (primary solid tumors etc.) in order to target the data subsets.",
            "You can select the features (thousands of genes, miRNAs and clinical metadata). Numerical values such as gene expression values will be shown as low/high and the mean will be calculated eith the selected low/high percentages. Categorical choices will be shown as levels.",
            "You can select a covariate group in order to compare it with first feature."
          ))
    )
  )
      }
        
        
        })
    
    
    observeEvent(input$KM_help, {
      
      introjs(session, options = list(steps = KM_steps()) )
      
    })
    
    
    
    
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

    output$proj_lenght_KM <- reactive({      #an output takes info form selectdatas Xproj$b() reactive
      
      
      Xproj$cancer_lenght()
      
    })
    
    outputOptions(output, "proj_lenght_KM", suspendWhenHidden = FALSE)  ##not so sure about this part
    
    
    
    # Prepare trimmed km_dat object
    
    km_dat <- reactive({
      
      validate(
        need(input$km_samptyp, "Select data subsets"),
        need(input$km_feat, "Select feature"),
        need(input$km_run, "Please click the Plot KM curves button to display KM curves"),
        if(input$km_covar != ""){
          if(class(Xproj$a()[[input$km_covar]]) %in% c("character", "factor")){
            need(input$sel_covar_meta_groups, "Select at least one covariate group")}} 
      )
      
      if(Xproj$cancer_lenght() ==1) {
        
        sel_cols <- c(input$km_feat, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient")
        
        dat <- Xproj$a()[, ..sel_cols]
        
        dat <- dat[meta.definition %in% input$km_samptyp, ]
        
        # De-duplicate patients. Expand on this to drop certain sample groups when duplicates are present
        dat <- dat[!duplicated(meta.patient), ]
        
        validate(
          if(!is.numeric(dat[[input$km_feat]]))  need(input$sel_feat_meta_groups, "Select at least one feature group"))
        
        if(is.numeric(dat[[input$km_feat]])){
          
          km_feat_zero_count <<- sum(dat[[input$km_feat]] == 0, na.rm=T)
          km_feat_na_count <<- sum(is.na(dat[[input$km_feat]]))
          
          mid_value_feat <- ifelse(input$keep_mid, "mid", NA)
          
          dat[, (input$km_feat) := ifelse(dat[[input$km_feat]] >= quantile(dat[[input$km_feat]], (100-input$hi_cutoff)/100, na.rm = T), "high", ifelse(dat[[input$km_feat]] <= quantile(dat[[input$km_feat]], input$lo_cutoff/100, na.rm = T), "low", mid_value_feat))]
          
        } else if(is.character(dat[[input$km_feat]]) | is.factor(dat[[input$km_feat]])){
          
          
          dat <- dat[get(input$km_feat) %in% input$sel_feat_meta_groups, ]
          
        }
      }
      
      
      else if (Xproj$cancer_lenght() > 1) {
        
        validate(
          need(input$choose_KM, "Select All or separate"))
      
      if(input$choose_KM == "All") {
        
        sel_cols <- c(input$km_feat, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient")
        
        dat <- Xproj$a()[, ..sel_cols]
        
        dat <- dat[meta.definition %in% input$km_samptyp, ]
        
        # De-duplicate patients. Expand on this to drop certain sample groups when duplicates are present
        dat <- dat[!duplicated(meta.patient), ]
        
        validate(
          if(!is.numeric(dat[[input$km_feat]]))  need(input$sel_feat_meta_groups, "Select at least one feature group"))
        
        if(is.numeric(dat[[input$km_feat]])){
          
          km_feat_zero_count <<- sum(dat[[input$km_feat]] == 0, na.rm=T)
          km_feat_na_count <<- sum(is.na(dat[[input$km_feat]]))
          
          mid_value_feat <- ifelse(input$keep_mid, "mid", NA)
          
          dat[, (input$km_feat) := ifelse(dat[[input$km_feat]] >= quantile(dat[[input$km_feat]], (100-input$hi_cutoff)/100, na.rm = T), "high", ifelse(dat[[input$km_feat]] <= quantile(dat[[input$km_feat]], input$lo_cutoff/100, na.rm = T), "low", mid_value_feat))]
          
        } else if(is.character(dat[[input$km_feat]]) | is.factor(dat[[input$km_feat]])){
          
          
          dat <- dat[get(input$km_feat) %in% input$sel_feat_meta_groups, ]
          
        }
      } else if (input$choose_KM == "Separate") {
        
        sel_cols_sep <- c(input$km_feat, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient", "meta.project_id")
        
        dat <- Xproj$a()[, ..sel_cols_sep]
        
        dat <- dat[meta.definition %in% input$km_samptyp, ]
        
        # De-duplicate patients. Expand on this to drop certain sample groups when duplicates are present
        dat <- dat[!duplicated(meta.patient), ]
        
        validate(
          if(!is.numeric(dat[[input$km_feat]]))  need(input$sel_feat_meta_groups, "Select at least one feature group"))
        
        if(is.numeric(dat[[input$km_feat]])){
          
          km_feat_zero_count <<- sum(dat[[input$km_feat]] == 0, na.rm=T)
          km_feat_na_count <<- sum(is.na(dat[[input$km_feat]]))
          
          mid_value_feat <- ifelse(input$keep_mid, "mid", NA)
          
          dat[, high_sep := quantile(get(input$km_feat), (100-input$hi_cutoff)/100, na.rm = T), by=c("meta.project_id")]  ##GET  ##CHANGE THE PLACE TO THE MIDDLE
          dat[, low_sep := quantile(get(input$km_feat), (input$lo_cutoff)/100, na.rm = T), by=c("meta.project_id")]
          dat[, (input$km_feat) := ifelse(dat[[input$km_feat]] >= high_sep, "high", ifelse(dat[[input$km_feat]] <= low_sep, "low", mid_value_feat))]
          
        } else if(is.character(dat[[input$km_feat]]) | is.factor(dat[[input$km_feat]])){
          
          
          dat <- dat[get(input$km_feat) %in% input$sel_feat_meta_groups, ]
          
        } 
        
      }
      }
      
      if(input$km_covar != "") {
        
        
        sel_cols2 <- c(input$km_covar, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient")
        
        covar_dat <- Xproj$a()[, ..sel_cols2][meta.definition %in% input$km_samptyp,][!duplicated(meta.patient), ]
        
        dat[, (input$km_covar) := covar_dat[, input$km_covar, with=F]]
        
        km_covar_zero_count <<- sum(dat[[input$km_covar]] == 0, na.rm=T)
        km_covar_na_count <<- sum(is.na(dat[[input$km_covar]]))
        
        if(is.numeric(dat[[input$km_covar]])){
          
          km_covar_zero_count <<- sum(dat[[input$km_covar]] == 0, na.rm=T)
          km_covar_na_count <<- sum(is.na(dat[[input$km_covar]]))
          
          mid_value_covar <- ifelse(input$keep_mid_covar, "mid", NA)
          
          dat[, (input$km_covar) := ifelse(dat[[input$km_covar]] >= quantile(dat[[input$km_covar]], (100-input$hi_cutoff_covar)/100, na.rm = T), "high", ifelse(dat[[input$km_covar]] <= quantile(dat[[input$km_covar]], input$lo_cutoff_covar/100, na.rm = T), "low", mid_value_covar))]
          
        } else if(is.character(dat[[input$km_covar]]) | is.factor(dat[[input$km_covar]])){
          
          
          dat <- dat[get(input$km_covar) %in% input$sel_covar_meta_groups, ]
          
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
      
      print(paste("Number of samples with zero counts of the numeric feature:", km_feat_zero_count))
      
      writeLines("\n")
      
      print(paste("Number of samples lacking information on selected feature:", km_feat_na_count))
      
      writeLines("\n")
      
      if(input$km_covar != ""){
        
        print(paste("Number of samples with zero counts of the numeric covariate:", km_covar_zero_count))
        
        writeLines("\n")
        
        print(paste("Number of samples lacking information on selected covariate:", km_covar_na_count))
        
      }
    })
    output$downloadPlot3 <- downloadHandler(
      filename = function() {
        paste("Kaplan_meier.png")
      },
      content = function(file) {
        
        png(file)
        print(val2$pl2)
        dev.off()
        
      })
  })
}
