library(survminer)
library(survival)
library(data.table)
library(shinyWidgets)

cox_ui <- function(id, label, choices) {
  
  ns <- NS(id)
  
  tagList(
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
        
        condition = " output.proj_lenght_cox > '1' ",   ## Conditions does not working with direct ractive values from selectdata which is Xproj$b, conditions are working with 
        ns=ns,                                           ## input or output result, that is wht I created an output that takes information form selecdata Xproj$b reavtive in server(Cagatay)
        radioButtons(inputId = ns ("choose_cox"), "Distribution type:",
                     selected = character(0),
                     c("All", "Separate")),                           ##We need to change All, separate names and make it more informative(Cagatay)
        
        
        
      ),
      
      checkboxInput(inputId = ns("cox_interaction"), "Model interactions?", F), 
      
      conditionalPanel(
        
        condition = "input.cox_interaction",
        ns=ns,  
        
        textInput(inputId = ns("cox_int_feat"), 
                  "Enter feature interactions with asterisks (exact names and capitalization are required). Multiple comma-separated interactions can be provided", 
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
      
    ),
    
    
    mainPanel(
      plotOutput(outputId = ns("cox_plot")),
      
      verbatimTextOutput( outputId = ns("cox_text"))
    )
  )
  
}


cox_server <- function(id,Xproj) {
  moduleServer(id,function(input, output, session) {
    
    ns <- session$ns
    
    COX_steps <- reactive({
      
      if(Xproj$cancer_lenght() ==1) {
        
        return(
          
          data.frame(
            
            element = paste0("#", session$ns(c(NA, "cox_samptyp + .selectize-control", "cox_feat + .selectize-control", "cox_interaction"))),
            
            intro = paste(c(
              "This is  Cox Proportional-Hazards Model app. The Cox proportional-hazards model,is essentially a regression model frequently used in medical research to examine the relationship between patient survival time and one or more predictor variables.Press the buttons to learn features of the app.",
              "You can select the sample types (primary solid tumors etc.) in order to target the data subsets.",
              "You can select the features (thousands of genes, miRNAs and clinical metadata).",
              "When there are several factors that have the potential to interact, multivariate analysis utilizing the method of Cox regression is used. For examplei if you choose A and B genes from select the features part, when you click the Model interactions? button, you can analyze A, B and A and B(A:B). Please write the exact names of the genes with asterisks."
            ))
            
          )
          
        )
      } else if (Xproj$cancer_lenght() > 1){
        return(
          data.frame(
            
            element = paste0("#", session$ns(c(NA, "cox_samptyp + .selectize-control", "cox_feat + .selectize-control","choose_cox", "cox_interaction"))),
            
            intro = paste(c(
              "This is  Cox Proportional-Hazards Model app. The Cox proportional-hazards model,is essentially a regression model frequently used in medical research to examine the relationship between patient survival time and one or more predictor variables.Press the buttons to learn features of the app.",
              "You can select the sample types (primary solid tumors etc.) in order to target the data subsets.",
              "You can select the features (thousands of genes, miRNAs and clinical metadata).",
              "In the case of several cancers, when ou select the 'All' option, calculations will be done without separating the data of different cancers. However, with the 'Separate' option, calculations will be done separately for the spesific cancers.",
              "When there are several factors that have the potential to interact, multivariate analysis utilizing the method of Cox regression is used. For examplei if you choose A and B genes from select the features part, when you click the Model interactions? button, you can analyze A, B and A and B(A:B). Please write the exact names of the genes with asterisks."
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
    output$proj_lenght_cox <- reactive({                                      
      
      
      
      Xproj$cancer_lenght()                          ## this is a reactive that takes information form selecdata module's Xproj$b() 
      
    })
    
    outputOptions(output, "proj_lenght_cox", suspendWhenHidden = FALSE)  
    
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
      
      if(Xproj$cancer_lenght() ==1) {                     ## in order to avoid zero lenght warning when there is omly one cancer
        
        predvars <- paste(input$cox_feat, collapse="+") 
        
      } else if (Xproj$cancer_lenght() > 1) {
        
        validate(
          need(input$choose_cox, "Select All or separate"))
      
      if( input$choose_cox == "All") {      ###İf user chooces all, there will be no project type addition
        
        predvars <- paste(input$cox_feat, collapse="+")
        
      } else if (input$choose_cox == "Separate")  {   
        
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


