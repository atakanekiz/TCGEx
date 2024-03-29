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
library(colourpicker)

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
                     "Please select sample types",
                     choices=NULL,
                     options=list(placeholder = "eg. Primary solid tumor")),
      
     radioButtons(ns("gvar_cat_x"), "Please select the type of x-axis variable", choices = c("Gene", "Meta"), selected = "Gene"),
      
       selectizeInput(
        inputId = ns("Gene1"),
        label = "Please select the x-axis variable",
        choices = NULL,
        # options=list(placeholder = "eg. TSPAN6 or meta.gender",
        #              plugins = list('restore_on_backspace'))
        
      ),
      
     radioButtons(ns("gvar_cat_y"), "Please select the type of y-axis variable", choices = c("Gene", "Meta"), selected = "Gene"),
      
     
     selectizeInput(
        inputId = ns("Gene2"),
        label = "Please select the y-axis variable",
        choices = NULL,
        # options=list(placeholder = "eg. TOX or meta.Histology",
        #              plugins = list('restore_on_backspace'))
        
      ),
     
     span(style="color:#3382FF",
          
          selectizeInput(
            inputId = ns("Gene3"),
            label = "Please select a gene to change point size (optional)",
            choices = NULL,
            options=list(placeholder = "eg. WNT",
                         plugins = list('restore_on_backspace'))
            
                 )
          
          ), 
     
     
     span(style="color:#3382FF",
          
          selectizeInput(
            inputId = ns("Gene4"),
            label = "Please select a gene to change point transparency (optional)",
            choices = NULL,
            options=list(placeholder = "eg. BRCA2",
                         plugins = list('restore_on_backspace'))
            
               )
          
          ), 
     
     
     span(style="color:#3382FF",
          
          selectizeInput(
            inputId = ns("Gene5"),
            label = "Please select a gene to change point color (optional)",
            choices = NULL,
            options=list(placeholder = "eg. MYC",
                         plugins = list('restore_on_backspace'))
            
             )
          
          ),
     
     numericInput(inputId = ns("gene_width"), "Choose the width of the plot", min = 1, step = 1, value = 18),
     numericInput(inputId = ns("gene_height"), "Choose the height of the plot", min = 1, step = 1, value = 9),
     numericInput(inputId = ns("text_size"), "Choose the font size of the plot", min = 1, step = 1, value = 25),
   
      
      checkboxInput(ns("notification"), "Show patient information", value = T),
      checkboxInput(ns("genecor_regline"), "Show best fitting line or curve", value = F),
      conditionalPanel(
       
       condition = "input.genecor_regline", ns=ns,
       
       radioButtons(ns("regline_type"), "Please select a fitting type", choices = c("linear", "loess")),
       
       conditionalPanel(
         
         condition = "input.regline_type == 'linear'", ns=ns,
         
         checkboxInput(ns("formula"), "Show statistics", value = T)
         
       ),
       
       conditionalPanel(
         condition = "input.regline_type == 'loess'", ns=ns,
         
         numericInput(ns("span"),
                      "Please pick a number for span",
                      min = 0.01,
                      step = 0.01,
                      value = 0.75)
       ),
       
       numericInput(ns("conf_lev"),
                    "Please select the confidance interval",
                    min = 0.01,
                    step = 0.01,
                    value = 0.95),
       
       colourInput(
         ns("col2"), "Please select a color for the regression line", "darkorange1",
         returnName = TRUE, 
         palette = "limited",
         closeOnClick = TRUE)
       
     ),
      checkboxInput(ns("facet"), "Add faceting variable", value = F),
      conditionalPanel(
       
       condition = "input.facet", ns=ns,
       
       selectizeInput(
         inputId = ns("Facet"),
         label = "Please select faceting variable",
         choices = NULL,
         options=list(placeholder = "eg. meta.gender",
                      plugins = list('restore_on_backspace')))
       
     ),
     
      # checkboxInput(ns("formula"), "Show statistics", value = T),
      actionBttn(inputId = ns("do"), 
                 label = "Analyze",
                 style = "unite",
                 block = TRUE,
                 color = "primary"),
      br(),
     
      #help section UI
      
      introjsUI(),
      actionButton(ns("intro"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
      
      width = 3
      
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
          
          element = paste0("#", session$ns(c(NA, "genecor_samp + .selectize-control", "gvar_cat_x","Gene1 + .selectize-control ", "gvar_cat_y","Gene2+ .selectize-control ","Gene3+ .selectize-control ","Gene4+ .selectize-control ","Gene5+ .selectize-control ","gene_width", "gene_height","text_size","notification","genecor_regline", "facet"))),
          
          intro = paste(c(
            "This is the gene-to-gene visualization module. You can calculate the correlation between two genes and generate scatter plots. Continue with the tutorial to learn features of the module.",
            "You can select the sample type here (eg. primary and/or metastatic). Only the selected subset(s) will be used in the analysis",
            "You can select gene expression or continuous metadata for the x-axis of the plot",
            "Choose a variable for the x-axis of the plot.",
            "You can select gene expression or continuous metadata for the y-axis of the plot ",
            "Choose a variable for the y-axis of the plot.",
            "Size of the data points can be changed according to the levels of the selected variable (optional).",
            "Transparency of the data points can be changed according to the levels of the selected variable (optional).",
            "Color of the data points can be changed according to the levels of the selected variable (optional).",
            "You can change the plot width (in inches)",
            "You can change the plot height (in inches)",
            "You can change the font size here",
            "Here, you can select whether or not to show patient information when you hover the cursor over data points",
            "You can show or hide the best-fitting line or curve to the data points. If you choose to show the best-fitting curve, you can also show the correlation coefficient and the p-value of the linear regression. If you choose to show the best-fitting curve, you can change the span parameter which controls the amount of smoothing (smaller values lead to wigglier lines, larger values lead to smoother lines). In addition, you can also show confidence intervals.",
            "You can select a categorical variable here to plot correlations in different data subsets (eg. male and female patients)."
          ))
          
        )
        
      })
      
      
      observeEvent(input$intro, {
        
        introjs(session, options = list(steps = intro() ) )
        
      })
      
      
      iv <- InputValidator$new()
      iv$add_rule("gene_width", ~ if (input$gene_width == 0 & !anyNA(input$gene_width)) "The number must be greater than 0")
      iv$add_rule("gene_height", ~ if (input$gene_height == 0 & !anyNA(input$gene_height)) "The number must be greater than 0")
      iv$add_rule("text_size", ~ if (input$text_size == 0 & !anyNA(input$text_size)) "The number must be greater than 0")
      iv$add_rule("span", ~ if ((input$span == 0 | input$span > 1) & !anyNA(input$span)) "The number must be in a range between 0 and 1")
      iv$add_rule("conf_lev", ~ if ((input$conf_lev == 0 | input$conf_lev > 1) & !anyNA(input$conf_lev)) "The number must be in a range between 0 and 1")
      iv$enable()
      
      
      observe({updateSelectizeInput(session, "genecor_samp",choices = Xproj$a()[["meta.definition"]], server = T)})
      # observe({updateSelectizeInput(session, 'Gene1', choices = gene_choice(), server = TRUE, selected = "")})
      observe({
        
        if(gvar_x() == "Gene") {
          
          updateSelectizeInput(session, 'Gene1', 
                               choices = gene_choice(), 
                               server = TRUE, 
                               selected = "",
                               options=list(placeholder = "eg. TSPAN6",
                                            plugins = list('restore_on_backspace')))
          
        }else if(gvar_x() == "Meta"){
          
          updateSelectizeInput(session, 'Gene1', 
                               choices = gene_choice(), 
                               server = TRUE, 
                               selected = "",
                               options=list(placeholder = "eg. Meta.t.cells.cd8",
                                            plugins = list('restore_on_backspace')))
          
        }
        
      })
      #observe({updateSelectizeInput(session, 'Gene2', choices = gene_choice_2(), server = TRUE, selected = "" )})
      
      observe({
        
        if(gvar_y() == "Gene") {
          
          updateSelectizeInput(session, 
                               'Gene2', 
                               choices = gene_choice_2(), 
                               server = TRUE, 
                               selected = "",
                               options=list(placeholder = "eg. TOX",
                                            plugins = list('restore_on_backspace')))
          
        }else if(gvar_y() == "Meta"){
          
          updateSelectizeInput(session, 
                               'Gene2', 
                               choices = gene_choice_2(), 
                               server = TRUE, 
                               selected = "",
                               options=list(placeholder = "eg. meta.Histology",
                                            plugins = list('restore_on_backspace')))
          
        }
        
      })
      observe({updateSelectizeInput(session, 'Gene3', choices = colnames(gen_dat()), server = TRUE, selected = "")})
      observe({updateSelectizeInput(session, 'Gene4', choices = colnames(gen_dat()), server = TRUE, selected = "" )})
      observe({updateSelectizeInput(session, 'Gene5', choices = colnames(gen_dat()), server = TRUE, selected = "")})
      observe({updateSelectizeInput(session, 'Facet', choices = available_cols_gplot(), server = TRUE, selected = "")})
      
      ## input choices for facet variable selection 
      
      # Remove cols which is including more than 10 levels and numeric
      
      hidden_cols_gplot<-reactive({"meta.definition"})

      meta_cols_gplot <- reactive({colnames(Xproj$a())[grep("^meta\\.", colnames(Xproj$a()))]})

      sel_cols_gplot <- reactive({meta_cols_gplot()[unlist(lapply(Xproj$a()[, meta_cols_gplot(), with = FALSE], function(x) length(levels(x)))) < 10]})

      numeric_cols_gplot <- reactive({unlist(lapply(Xproj$a()[, sel_cols_gplot(), with = FALSE], function(x) is.numeric(x)))})

      available_cols_gplot <- reactive({setdiff(sel_cols_gplot(), c(hidden_cols_gplot(), sel_cols_gplot()[numeric_cols_gplot()]))})
      
      
      
      ## input choices for gene selection
      
      
  gen_dat <- reactive({
      
      Xproj$a() %>%
        select(!starts_with("meta.")) %>%
        select(where(is.numeric))
      
    })
      
    
    meta_dat <- reactive({
      
      Xproj$a() %>%
        select(starts_with("meta.")) %>%
        select(where(is.numeric))
      
    })
      
     gvar_x <- reactive(input$gvar_cat_x)
     gvar_y <- reactive(input$gvar_cat_y)
      
      gene_choice <- reactive({
        
        if(gvar_x() == "Gene") {
          
          DF <- gen_dat()
     
       }else if(gvar_x() == "Meta"){
         
         DF <- meta_dat()
          
       }
        
          colnames(DF)
        
      })
      
      
      
      gene_choice_2 <- reactive({
        
        if(gvar_y() == "Gene"){
          
          DF <- gen_dat()
          
          
        }else if (gvar_y() == "Meta") {
          
          DF <- meta_dat()
          
        }
        
        colnames(DF)
        
      })
      
      
      
      
      df <- reactive({Xproj$a()%>% 
          filter( meta.definition %in% input$genecor_samp) %>% 
          filter(!is.na(input$Gene1)) %>% 
          filter(!is.na(input$Gene2))
      })
      
      ## to make reactive facet_grid function 
      
      facet_cat <- reactive({paste("~", input$Facet)})
      
      down_width <- reactive({paste(".", input$download_width)})
      
      
      
      
      ## ggiraph package is used to form interactive graph from ggplot 2 
      
      
      gene_graph <- reactive({
        
        dataset = df()
        
        p = ggplot(data = dataset) 
        
        {if(input$notification == T) p <- p + geom_point_interactive(aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]],size= dataset[[input$Gene3]], alpha= dataset[[input$Gene4]],color = dataset[[input$Gene5]],
                                                                     tooltip = paste0( "Gender: ",dataset[["meta.gender"]] , "<br/>",
                                                                                       "Patient ID: ", dataset[["meta.patient"]], "<br/>",
                                                                                       "Race: ", dataset[["meta.race"]])))}
        
        
        
        {if(input$notification == F) p <- p + geom_point_interactive(aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]], alpha= dataset[[input$Gene3]]), size= dataset[[input$Gene4]], color = dataset[[input$Gene5]])}       
        
        
        {if(input$formula & input$regline_type == 'linear') p <- p + stat_cor(mapping = aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]]), family = "Arial", size = 7, geom = "label")}        
        
        # {if(input$genecor_regline) p <- p + stat_smooth(mapping = aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]]), method = "loess", color = input$col2)}
        
        {if(input$regline_type == "linear" & input$genecor_regline) p <- p + stat_smooth(mapping = aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]]), method = "lm", color = input$col2, level = input$conf_lev)}   
        
        {if(input$regline_type == "loess" & input$genecor_regline) p <- p + stat_smooth(mapping = aes(x = .data[[input$Gene1]], y = .data[[input$Gene2]]), method = "loess", color = input$col2, span = input$span, level = input$conf_lev)}   
        
        {if(input$facet) p <- p + facet_wrap(facet_cat(), labeller = as_labeller(dataset, default=label_wrap_gen(16)), scales = "free_y")+ 
            
            theme(aspect.ratio = 1)}   
        
        p <- p + theme_bw(base_size = 16)
        
        legend_title1 <- input$Gene3
        legend_title2 <- input$Gene4
        legend_title3 <- input$Gene5
        
        p <- p + scale_size(legend_title1, range = c(1, 10))
        
        p <- p + scale_alpha_continuous(legend_title2, range = c(0.1, 1))
        
        p <- p + scale_colour_continuous(legend_title3, type = getOption("ggplot2.continuous.colour"))
        
        p <- p + theme(text=element_text(family="Arial", face="bold", size=input$text_size))
        
        girafe(ggobj = p, 
               width_svg = input$gene_width , height_svg = input$gene_height,
               options = list(
                 opts_sizing(rescale = TRUE, width = 0.15),
                 opts_hover_inv(css = "opacity:0.1;"),
                 opts_zoom(max = 5)
               )
               
        )
        
        
      })
      
      
      observeEvent(input$do, {
        
        
        
        
        output$plot <- renderGirafe({
          
          
          req(input$do)
          input$do
          
          isolate({
            
            req(iv$is_valid())
            
            validate(
              need(input$genecor_samp, "Choose at least one sample type"),
              need(input$Gene1, "Don't forget to choose a variable for x-axis"),
              need(input$Gene2, "Don't forget to choose a variable for y-axis"),
              need(input$gene_width, "Don't forget to choose the width of the plot"),
              need(input$gene_height, "Don't forget to choose the height of the plot"),
              need(input$text_size, "Don't forget to choose the font size")
            )
            
            if (input$Facet < 0  && input$facet == T ) {
              validate("Choose a faceting parameter")
            }
            
            gene_graph()
            
          })
          
        })
        

        
      })
      
    }
  )
}
