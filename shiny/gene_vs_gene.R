library(plotly)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(shiny)
library(shinydashboard)
library(ggiraph)
library(shinybusy)
library(extrafont)
library(rintrojs)
library(shinyvalidate)
library(shinyWidgets)

gene_vs_gene_ui <- function(id) {
  ns <- NS(id)
  tagList(
    
    add_busy_spinner(
      spin = "cube-grid",
      position = "top-right",
      color = "#01303f",
      margins = c(300, 500),
      height = "60px",
      width = "60px"),
    
    
    sidebarPanel(
      
      
      selectizeInput(ns("genecor_samp"), multiple=T,
                     "*Please select sample types",
                     choices=NULL,
                     options=list(placeholder = "eg. Primary solid tumor")),
      
      selectizeInput(
        inputId = ns("Gene1"),
        label = "*Please select Gene-1",
        choices = NULL,
        options=list(placeholder = "eg. TSPAN6")
        
      ),
      
      selectizeInput(
        inputId = ns("Gene2"),
        label = "*Please select Gene-2",
        choices = NULL,
        options=list(placeholder = "eg. TOX"),
        
      ),
      
      conditionalPanel(
        
        condition = "input.facet", ns=ns,
        
        selectizeInput(
          inputId = ns("Facet"),
          label = "Please select faceting variable",
          choices = NULL,
          options=list(placeholder = "eg. meta.gender"))
        
      ),
      
      checkboxInput(ns("notification"), "Show patient information", value = T),
      checkboxInput(ns("genecor_regline"), "Show regression line", value = F),
      checkboxInput(ns("facet"), "Add faceting variable", value = F),
      checkboxInput(ns("formula"), "Show statistics", value = F),
      actionBttn(inputId = ns("do"), 
                 label = "Generate Correlation Plot",
                 style = "unite",
                 block = TRUE,
                 color = "primary"),
      br(),
      #help section UI
      
      introjsUI(),
      actionButton(ns("intro"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
      
      
      
    ),
    
    
    
    mainPanel(
      girafeOutput(ns("plot"))
      
    )
    
    
    
  )
  
}

gene_vs_gene_server <- function(id,Xproj) {
  moduleServer(
    id,
    function(input, output, session) {
      
      
      
      
      
      ## help section server
      
      ns <- session$ns
      
      intro <- reactive({
        
        data.frame(
          
          element = paste0("#", session$ns(c(NA, "genecor_samp + .selectize-control", "Gene1 + .selectize-control ", "Gene2+ .selectize-control ", "notification","genecor_regline", "facet", "formula"))),
          
          intro = paste(c(
            "This is the gene-to-gene visualization module. You can calculate the correlation between two genes and generate scatter plots. Continue with the tutorial to learn features of the module.",
            "You can select the sample type here (eg. primary and/or metastatic). Only the selected subset(s) will be used in the analysis",
            "Choose a gene for the x-axis of the plot.",
            "Choose a gene for the y-axis of the plot.",
            "Here, you can select whether or not to show associated metadata when you hover the mouse cursor over data points",
            "You can show or hide the best-fitting line to the data points",
            "Enter a covariate here to plot the correlations in different data subsets (eg. male and female patients).",
            "You can also show the correlation coefficient and the p-value of the linear regression."
          ))
          
        )
        
      })
      
      
      observeEvent(input$intro, {
        
        introjs(session, options = list(steps = intro() ) )
        
      })
      
      
      
      observe({updateSelectizeInput(session, "genecor_samp",choices = Xproj$a()[["meta.definition"]], server = T)})
      observe({updateSelectizeInput(session, 'Gene1', choices = colnames(gene_choice()), server = TRUE, selected = "")})
      observe({updateSelectizeInput(session, 'Gene2', choices = colnames(gene_choice()), server = TRUE, selected = "" )})
      observe({updateSelectizeInput(session, 'Facet', choices = colnames(facet_choice()), server = TRUE, selected = "")})
      
      
      ## input choices for facet category selection 
      
      facet_choice <- reactive({
        
        Xproj$a() %>% 
          select(starts_with("meta.")) %>% 
          select(where(is.factor))
        
      })
      
      
      ## input choices for gene selection
      
      gene_choice <- reactive({
        
        Xproj$a() %>%
          select(!starts_with("meta.")) %>% 
          select(where(is.numeric))
        
      })
      
      
      
      df <- reactive({Xproj$a()%>% 
          filter( meta.definition %in% input$genecor_samp) %>% 
          filter(!is.na(input$Gene1)) %>% 
          filter(!is.na(input$Gene2))
      })
      
      
      
      ## to make reactive facet_grid function 
      
      facet_cat <- reactive({paste("~", input$Facet)})
      
      
      ## ggiraph package is used to form interactive graph from ggplot 2 
      
      
      
      observeEvent(input$do, {
        
        
        
        
        output$plot <- renderGirafe({
          
          
          req(input$do)
          input$do
          
          isolate({
            
            validate(
              need(input$genecor_samp, "Choose at least one sample type"),
              need(input$Gene1, "Choose the first gene"),
              need(input$Gene2, "Choose the second gene"),
            )
            
            if (input$Facet < 0  && input$facet == T ) {
              validate("Choose a faceting parameter")
            }
            
            dataset = df()
            
            p = ggplot(data = dataset) 
            
            {if(input$notification == T) p <- p + geom_point_interactive(aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]]),
                                                                         tooltip = paste0( "Gender: ",dataset[["meta.gender"]] , "<br/>",
                                                                                           "Patient ID: ", dataset[["meta.patient"]], "<br/>",
                                                                                           "Race: ", dataset[["meta.race"]]))}
            
            
            
            {if(input$notification == F) p <- p + geom_point_interactive(aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]]))}       
            
            
            {if(input$formula) p <- p + stat_cor(mapping = aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]]), family = "Arial", size = 7, color = "black", geom = "label")}        
            
            {if(input$genecor_regline) p <- p + stat_smooth(mapping = aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]]), method = "lm", color = "#ff5317")}   
            
            {if(input$facet) p <- p + facet_wrap(facet_cat(), labeller = as_labeller(dataset, default=label_wrap_gen(16)), scales = "free_y")+ 
                
                theme(aspect.ratio = 1)}   
            
            p <- p + theme_bw(base_size = 16)
            
            p <- p + theme(text=element_text(family="Arial", face="bold", size=30),
                           plot.title = element_text(color="black", face="bold.italic", hjust = 0.5),
                           strip.text.x = element_text(size = 20))
            
            girafe(ggobj = p, 
                   width_svg = 15 , height_svg = 9,
                   options = list(
                     opts_sizing(rescale = TRUE, width = .15),
                     opts_hover_inv(css = "opacity:0.1;"),
                     
                     opts_zoom(max = 5)
                   )
                   
            )
            
          })
          
        })
        
      })
      
    }
  )
}