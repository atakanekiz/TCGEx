library(survminer)
library(survival)
library(data.table)
library(shinyWidgets)
library(shinyalert)

cox_ui <- function(id, label, choices) {
  
  ns <- NS(id)
  
  tagList(
    
    useShinyalert(),
    
    sidebarPanel(
      
      selectizeInput(inputId = ns("cox_samptyp"), 
                     label = "Select sample types",
                     choices = NULL,  # will be updated dynamically,
                     options=list(placeholder = "eg. Primary solid tumor"),
                     multiple=T),
      
      hr(),
      
      selectizeInput(inputId = ns("cox_feat"), 
                     label = "Select features",
                     choices = NULL,
                     multiple=T, # will be updated dynamically,
                     options=list(placeholder = "eg. LAG3",
                                  plugins = list('restore_on_backspace'))),
      hr(),
      
      conditionalPanel(
        
        condition = " output.proj_length_cox > '1' ",    
        ns=ns,                                           
        radioButtons(inputId = ns("choose_cox"), "Select multi-cancer analysis approach:",
                     selected = "Include cancer type as covariate",
                     c("Include cancer type as covariate",
                       "Aggregate cancer types")),
        
        
        
      ),
      
      checkboxInput(inputId = ns("cox_interaction"), "Model interactions?", F), 
      
      conditionalPanel(
        
        condition = "input.cox_interaction",
        ns=ns,  
        
        textInput(inputId = ns("cox_int_feat"), 
                  "Enter feature interactions with asterisks (exact names are needed and capitalization matters). Multiple comma-separated interactions can be provided", 
                  placeholder = "eg. FEAT1*FEAT2, FEAT1*FEAT3")
      ),
      actionBttn(inputId = ns("cox_run"), 
                 label = "Analyze",
                 style = "unite",
                 block = TRUE,
                 color = "primary"),
      br(),
      
      downloadButton(ns("downloadPlot2"), "Download Cox_ph Plot", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
      
      #help section UI
      
      introjsUI(),
      actionButton(ns("COX_help"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
      

      textOutput(ns("filewarning")) ,     
      
      
      width = 3
      
    ),
    
    
    mainPanel(
      plotOutput(outputId = ns("cox_plot")),
      verbatimTextOutput( outputId = ns("cox_text")),
    
    )
  )
  
}


cox_server <- function(id,Xproj) {
  moduleServer(id,function(input, output, session) {
    
      output$filewarning <- renderText({
        
        if (!is.null(Xproj$fileInfost())) {
          shinyalert("Warning!", "To perform Cox Proportional-Hazard Survival analysis, the data you upload must contain columns containing survival information such as 'vital_status' and 'days _to_event'.") }
      })
    
    ns <- session$ns
    
    COX_steps <- reactive({
      
      if(Xproj$cancer_length() ==1) {
        
        return(
          
          data.frame(
            
            element = paste0("#", session$ns(c(NA, "cox_samptyp + .selectize-control", "cox_feat + .selectize-control", "cox_interaction"))),
            
            intro = paste(c(
              "This is  Cox Proportional Hazards (CoxPH) survival analysis module. The CoxPH is used examine whether a feature has an impact on the differential survival outcome. Features can be numeric (eg. expression of a gene) or categorical (eg. mutation subtype), and multiple features can be simultaneously modeled. In CoxPH analysis, hazard ratio (HR) >1 is associated with an increased risk, while HR<1 is associated with reduced risk. Continue tutorial to understand how the module works.",
              "You can select the sample types (eg. primary and/or metastatic samples) in order to target specific data subsets in the analysis.",
              "You can select one or more features you would like to analyze (eg. genes or clinical metadata). If a single feature is selected, univariate CoxPH analysis is performed. When two or more features are selected, multivariate CoxPH analysis is performed where the effects of individual features are reported along with the overall effects.",
              "You can also model interactions between features (optional). This more complex modelling allows examining whether covariates have an impact on each other's effect. The user can, for instance, investigate whether Gene_A has a different survival impact in males and females by specifying Gene_A*meta.gender. Interaction terms must match exactly to the selected features and capitalization matters here."
            ))
            
          )
          
        )
      } else if (Xproj$cancer_length() > 1){
        return(
          data.frame(
            
            element = paste0("#", session$ns(c(NA, "cox_samptyp + .selectize-control", "cox_feat + .selectize-control","choose_cox", "cox_interaction"))),
            
            intro = paste(c(
              "This is  Cox Proportional Hazards (CoxPH) survival analysis module. The CoxPH is used examine whether a feature has an impact on the differential survival outcome. Features can be numeric (eg. expression of a gene) or categorical (eg. mutation subtype), and multiple features can be simultaneously modeled. In CoxPH analysis, hazard ratio (HR) >1 is associated with an increased risk, while HR<1 is associated with reduced risk. Continue tutorial to understand how the module works.",
              "You can select the sample types (eg. primary and/or metastatic samples) in order to target specific data subsets in the analysis.",
              "You can select one or more features you would like to analyze (eg. genes or clinical metadata). If a single feature is selected, univariate CoxPH analysis is performed. When two or more features are selected, multivariate CoxPH analysis is performed where the effects of individual features are reported along with the overall effects.",
              "If you selected more than one TCGA dataset, you have the option to include 'cancer type' as a covariate in survival modeling (default behavior). Alternatively, you can perform analysis agnostic to the cancer type by selecting 'Aggregate'. In this case, cancer type will not be included in the model as a covariate. The latter can be beneficial for examining closely related cancers with similar life expectancies.",
              "You can also model interactions between features (optional). This more complex modelling allows examining whether covariates have an impact on each other's effect. The user can, for instance, investigate whether Gene_A has a different survival impact in males and females by specifying Gene_A*meta.gender. Interaction terms must match exactly to the selected features and capitalization matters here."
            ))
          )
        )
      }
      
      
    })
    
    
    observeEvent(input$COX_help, {
      
      introjs(session, options = list(steps = COX_steps()) )
      
    })
    
    
    observe({updateSelectizeInput(session,
                                  "cox_samptyp", 
                                  choices = Xproj$a()$meta.definition, 
                                  server = T)})
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session, 
                                  "cox_feat", selected="",
                                  choices = colnames(Xproj$a()), 
                                  server = T)})
    output$proj_length_cox <- reactive({                                      
      
      
      
      Xproj$cancer_length()                          ## this is a reactive that takes information form selecdata module's Xproj$b() 
      
    })
    
    outputOptions(output, "proj_length_cox", suspendWhenHidden = FALSE)  
    
    cox_dat <- eventReactive(input$cox_run,{
      
      
      validate(need(input$cox_samptyp, "Select sample type"),
               need(input$cox_feat, "Select features"))
      
      sel_cols <- c(input$cox_feat, "meta.vital_status", "meta.days_to_event", "meta.definition", "meta.patient", "meta.project_id")
      
      df_data <- Xproj$a()[, ..sel_cols]      # df_data is a dataframe that contains needed data for Cox
      
      df_data <- df_data[meta.definition %in% input$cox_samptyp, ]
      
      # De-duplicate patients. Revise this to drop certain subgroups later.
      df_data <- df_data[!duplicated(meta.patient), ]
      
      feats <- input$cox_feat
      
      numeric_cols <- sapply(df_data[, ..feats], is.numeric)
      
      numeric_cols <- names(numeric_cols[numeric_cols == TRUE])
      
      df_data[, (numeric_cols) := lapply(.SD, scale), .SDcols = numeric_cols]
      
      df_data
      
      
    })
    
    
    
    cox_form <- eventReactive(input$cox_run, {
      
      validate(need(cox_dat(), ""))
      
      if(Xproj$cancer_length() ==1) {                     ## in order to avoid zero length warning when there is only one cancer
        
        predvars <- paste(input$cox_feat, collapse="+") 
        
      } else if (Xproj$cancer_length() > 1) {
        
        validate(
          need(input$choose_cox, "Select 'Include cancer type as covariate' or 'Aggregate cancer types'"))
        
        if( input$choose_cox == "Aggregate cancer types") {      ### If user chooses all, there will be no project type covariate
          
          predvars <- paste(input$cox_feat, collapse="+")
          
        } else if (input$choose_cox == "Include cancer type as covariate")  {   
          
          project_type <- "meta.project_id"
          predvars_sep <- paste(input$cox_feat, collapse="+")
          predvars <-  paste(predvars_sep, "+", project_type)
        }
      }
      
      if(input$cox_interaction){
        
        interactvar <- input$cox_int_feat
        
        interactvar <- gsub(" ", "", interactvar)
        
        interactvar <- unlist(strsplit(interactvar, ","))
        
        interactvar <- paste(interactvar,collapse = "+")
        
        predvars <- paste(predvars, "+", interactvar)
        
      }
      
      
      survform <- paste("Surv(meta.days_to_event, meta.vital_status) ~",
                        predvars)
      
      
      survform <- as.formula(survform)
      
      survform
      
    })
    
    fit <- reactive(
      coxph(cox_form(), data=cox_dat())
      
    )
    
    
    val <- reactiveValues()
    
    observeEvent(input$cox_run, {
      
      
      output$cox_plot <- renderPlot({
        
        
        if(input$cox_interaction) validate(need(input$cox_run, "Enter interactions and click analyze"))
        
        pl <- ggforest(model=fit(), fontsize = 1.3, 
                       data = cox_dat())
        
        val$pl <- pl
        
        print(pl)
        
      })
      
      
      output$cox_text <- renderPrint({
        
        if(input$cox_interaction) validate(need(session, input$cox_run, "Enter interactions and click analyze"))
        
        
        cat("\n########################################################################\n")
        
        cat("Survival formula\n")
        print(cox_form())
        
        cat("########################################################################\n\n\n")
        
        
        print(summary(fit()))
        
        
        
        cat("\n########################################################################\n")
        
        cat("Testing proportional hazards assumption.\np-values less than 0.05 indicates the proportionality assumption is violated.\nConsider using time-dependent modeling\n")
        
        
        print(cox.zph(fit()))
        
        cat("########################################################################\n")
        
      })
      
      
      output$downloadPlot2 <- downloadHandler(
        filename = function() {
          paste("Cox_ph_plot.png")
        },
        content = function(file) {
          
          png(file, width = 480, height = 480)
          print(val$pl)
          dev.off()
          
        })
      
      
    })
    
  }
  )
  
  
}


