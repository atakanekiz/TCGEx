library(dplyr)
library(shiny)
library(ggpubr)
library(data.table)

gene_vs_cat_ui <- function(id, label, choices){
  
  ns <- NS(id)
  
  tagList(
    
    ui <- fluidPage(
      
      sidebarPanel(
        
        
        selectizeInput(ns("exprs_samptyp"), multiple=T,
                       "Select sample types",
                       choices=NULL, # will be updated dynamically
                       options=list(placeholder = "eg. Primary solid tumor")),
        
        hr(),
        
        
        
        selectizeInput(ns("cat_plotvar"), multiple=F,
                       "Select categorical variable (x-axis)",
                       choices=NULL, # will be updated dynamically
                       options = list(placeholder = "eg. meta.tumor_stage")),
        
        
        selectizeInput(ns("cat_plotvarsubset"), multiple=T,
                       "Select categories to include in the plot",
                       choices=NULL, # will be updated dynamically
                       options = list(placeholder = "eg. meta.tumor_stage")),
        
        
        selectizeInput(ns("num_plotvar"), multiple=F,
                       "Select numerical gene/feature (y-axis)",
                       choices=NULL, # will be updated dynamically
                       options= list(placeholder = "eg. CD8A",
                                     plugins = list('restore_on_backspace'))),
        
        selectizeInput(ns("facet_plotvar"), multiple=F,
                       "Select faceting variable (optional)",
                       choices=NULL, # will be updated dynamically
                       options = list(placeholder = "eg. meta.gender")),
        
        checkboxInput(inputId = ns("exprs_stats"), "Show statistics?", F),
        
        actionBttn(inputId = ns("cat_gene_run"), 
                   label = "Analyze",
                   style = "unite",
                   block = TRUE,
                   color = "primary"),
        br(),
        downloadButton(ns("downloadPlot4"), "Download plot", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
        
        #help section UI
        
        introjsUI(),
        actionButton(ns("genecat_help"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
        
        textOutput(ns("filewarning_box")) ,     
        
        width = 3
        
      ),
      
      mainPanel(
        
        
        plotOutput(outputId = ns("exprs_plot")),
        
        conditionalPanel(           ## Before changing this section, panel was including every conditional panel. It can be reverted if needed
          
          h3("Graphing options"),
          
          condition = "output.exprs_plot",
          
          ns=ns,
          
          
          column(4, 
                 
                 selectInput(inputId = ns("exprs_plotpal"), "Select plotting colors", multiple = F, 
                             selectize = T, 
                             selected = "jco",
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
                                 "Brown-Green" = "BrBG")) 
                             
                 ),
                 
                 checkboxInput(inputId = ns("exprs_rotatex"), "Rotate x-axis labels?", F)
                 
                 
          ), 
          
          
          column(4,
                 
                 selectInput(inputId = ns("exprs_plotadd"), "Add features to plot", 
                             multiple = F, choices = c("Jitter" = "jitter", 
                                                       "Mean" = "mean",
                                                       "Mean+SE" = "mean_se",
                                                       "None" = "none"),
                             selected = "none")),
          
          conditionalPanel(
            
            condition = "input.exprs_plotadd !== 'none'",
            ns=ns,
            
            numericInput(inputId = ns("exprs_plotaddsize"), "Size of the points", 
                         value = 1, min = 0.5, max=3, step=0.5),
            
            selectInput(inputId = ns("exprs_plotaddcolor"), "Color of the points",
                        choices = c("black", "blue", "red", "gold2"),
                        selected = "black")
          )),
        
        column(4, 
               
               conditionalPanel(
                 
                 
                 condition="input.exprs_stats ",
                 
                 ns=ns,
                 
                 selectInput(inputId = ns("exprs_statmethod"), "Select statistical test", 
                             choices = c("T-test"="t.test",
                                         "Wilcoxon"="wilcox.test"), selected = "t.test"),
                 
                 
                 selectInput(inputId = ns("exprs_statlabel"), "Select how to show p-values",
                             choices = c("p.adj symbol" = "p.signif.adj", 
                                         "p.adj numeric" = "p.adj",
                                         "p symbol"="p.signif",
                                         "p numeric"="p.format"),
                             selected = "p.signif.adj"),
                 
                 radioButtons(inputId = ns("select_comp"), "Select comparisons",
                              choiceNames = c("Single reference",
                                              "Pairwise"),
                              choiceValues = c("ref", "pairs"), selected = "ref"),
                 
                 conditionalPanel(
                   
                   condition = "input.select_comp == 'ref' ",
                   ns=ns,
                   
                   selectizeInput(inputId = ns("exprs_statref"), "Select reference", multiple=F,
                                  choices=NULL, # will be updated dynamically
                                  options = list(placeholder = "eg. meta.tumor_stage")),
                   
                   
                   
                 ),
                 
                 conditionalPanel(
                   
                   condition = "input.select_comp == 'pairs'",
                   ns=ns,
                   
                   selectizeInput(inputId = ns("exprs_statpairs"), "Select group pairs", 
                                  multiple=T,  
                                  choices=NULL, #update dynamically
                                  options = list(maxItems=2)
                   ),
                   
                   
                   actionButton(inputId = ns("exprs_addcomp"), "Add to comparisons"),
                   
                   actionButton(inputId = ns("exprs_resetcomp"), "Reset comparisons"),
                   
                   p("Selected pairwise comparisons"),
                   
                   verbatimTextOutput(outputId = ns("exprs_statpairlist")
                                      
                   ),
                   
                   
                   
                   
                 ))
        )  
      )
    )
  )
}





gene_vs_cat_server <- function(id,Xproj){
  
  moduleServer(id,function(input, output, session){
    
    
    output$filewarning_box <- renderText({
      
      if (!is.null(Xproj$fileInfost())) {
        shinyalert("Warning!", "Please ensure that the column names containing clinic information in your data start with 'meta.' . Otherwise, you will not be able to use this analysis module.") }
    })
    
    
    genecat_steps <- reactive({
      
      
      return(
        
        data.frame(
          
          element = paste0("#", session$ns(c(NA, 
                                             "exprs_samptyp + .selectize-control", 
                                             "cat_plotvar + .selectize-control",
                                             "cat_plotvarsubset + .selectize-control",
                                             "num_plotvar + .selectize-control",
                                             "facet_plotvar + .selectize-control",
                                             "exprs_stats",
                                             "downloadPlot4"))),
          
          intro = c(
            "This is feature-to-metadata visualization module. You can select various data subsets and see how numeric features such as gene expression differs among these. Continue the tutorial to see how the module works.",
            "Select sample type(s) you would like to analyze here.",
            "Define which data subsets you want to include in the analysis.",
            "Here, you can select the categorical variable to place on the x-axis",
            "Select a numeric feature to plot in the y-axis",
            "You can also select faceting variables to examine how feature-vs-category relationship might be different in various data subsets (eg. males and females).",
            "Click here if you would like to show statistics on the graph. You can also define reference group or specific pair-wise comparisons in the main panel.",
            "You can download the plot by clicking here."
          )
          
        )
        
      )
    })
    
    
    observeEvent(input$genecat_help, {
      
      introjs(session, options = list(steps = genecat_steps()) )
      
    })
    
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session,
                                  "exprs_samptyp",
                                  choices = Xproj$a()$meta.definition,
                                  server = T)})
    
    # Remove cols which is including more than 10 levels and numeric
    
    hidden_cols_cat_plotvar<-reactive({c("meta.definition","meta.treatments")})
    
    meta_cols_cat_plotvar <- reactive({colnames(Xproj$a())[grep("^meta\\.", colnames(Xproj$a()))]})
    
    sel_cols_cat_plotvar <- reactive({meta_cols_cat_plotvar()[unlist(lapply(Xproj$a()[, meta_cols_cat_plotvar(), with = FALSE], function(x) length(levels(x)))) < 10]})
    
    numeric_cols_cat_plotvar <- reactive({unlist(lapply(Xproj$a()[, sel_cols_cat_plotvar(), with = FALSE], function(x) is.numeric(x)))})
    
    available_cols_cat_plotvar <- reactive({setdiff(sel_cols_cat_plotvar(), c(hidden_cols_cat_plotvar(), sel_cols_cat_plotvar()[numeric_cols_cat_plotvar()]))})
    
    
    
    # Remove cols which is including more than 10 levels and numeric for facet_plotvar
    
    hidden_cols_facet_plotvar<-reactive({c("meta.definition", "meta.treatments")})
    
    meta_cols_facet_plotvar <- reactive({colnames(Xproj$a())[grep("^meta\\.", colnames(Xproj$a()))]})
    
    sel_cols_facet_plotvar <- reactive({meta_cols_facet_plotvar()[unlist(lapply(Xproj$a()[, meta_cols_facet_plotvar(), with = FALSE], function(x) length(levels(x)))) < 10]})
    
    numeric_cols_facet_plotvar <- reactive({unlist(lapply(Xproj$a()[, sel_cols_facet_plotvar(), with = FALSE], function(x) is.numeric(x)))})
    
    available_cols_facet_plotvar <- reactive({setdiff(sel_cols_facet_plotvar(), c(hidden_cols_facet_plotvar(), sel_cols_facet_plotvar()[numeric_cols_facet_plotvar()]))})
   
     #   # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session,
                                  "cat_plotvar", selected="",
                                  choices = available_cols_cat_plotvar(),
                                  server = T)})
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session,
                                  "cat_plotvarsubset", selected=levels(factor(Xproj$a()[[input$cat_plotvar]])),
                                  choices = levels(factor(Xproj$a()[[input$cat_plotvar]])),
                                  server = T)})
    
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session,
                                  "num_plotvar", selected="",
                                  choices = colnames(Xproj$a()%>% select(-starts_with("meta."))),
                                  server = T)})
    
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session,
                                  "facet_plotvar", selected="",
                                  choices = available_cols_facet_plotvar(),
                                  server = T)})
    
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session,
                                  "exprs_statref", selected="",
                                  choices = levels(factor(Xproj$a()[[input$cat_plotvar]])),
                                  server = T)})
    
    
    # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
    observe({updateSelectizeInput(session,
                                  "exprs_statpairs", selected="",
                                  choices = input$cat_plotvarsubset,
                                  server = T)})
    
    
    exprs_plotdat <- reactive({
      
      validate(need(input$exprs_samptyp, "Select sample type"),
               need(input$cat_plotvar, "Select categorical variable for x-axis"),
               need(input$num_plotvar, "Select numerical variable for y-axis"),
               need(input$cat_gene_run, "Click 'Analyze' button to generate the graph"),
               if(input$exprs_stats){                              ## Validation is active only when statistics button is chosen
                 need(input$exprs_statref, "Specify the reference group or choose pairwise combinations below to show statistical comparisons (or, simply unclick 'Show statistics?' box on the left panel to plot graphs without statistics)")}
      )
      
      if(input$facet_plotvar == "") facetvar <- NULL else facetvar <- input$facet_plotvar
      
      sel_cols <- c("meta.definition", input$num_plotvar, input$cat_plotvar, facetvar)
      
      df <- Xproj$a()[, ..sel_cols]
      
      df <- df[meta.definition %in% input$exprs_samptyp, ]
      
      df <- df[get(input$cat_plotvar) %in% input$cat_plotvarsubset, ]
      
      df
      
    })
    
    complist <- reactiveVal(NULL)
    
    
    observeEvent(input$exprs_addcomp, {
      
      updatedcomps <- append(complist(), list(input$exprs_statpairs))
      
      complist(updatedcomps)
      
      
      # WARNING Input to asJSON(keep_vec_names=TRUE) is a named vector.
      observe({updateSelectizeInput(session, 
                                    "exprs_statpairs", selected="",
                                    choices = input$cat_plotvarsubset,
                                    server = T)})
      
    })
    
    
    observeEvent(input$exprs_resetcomp, complist(NULL))
    
    
    output$exprs_statpairlist <- renderPrint({
      
      complist()
      
    })
    
    
    
    stat_layer <-reactive({ 
      validate(need(input$exprs_samptyp, "Select sample type"),
               need(input$cat_plotvar, "Select categorical variable for x-axis"),
               need(input$num_plotvar, "Select numerical variable for y-axis"),
               need(input$cat_gene_run, " Click generate button to generate the graph"),
               if(input$exprs_stats){                              ## Validation is active only when statistics button is chosen
                 need(input$exprs_statref, "Select the reference group below")}
      )
      
      
      form <- as.formula(paste0(input$num_plotvar, "~", input$cat_plotvar))
      
      if(input$facet_plotvar == "") facetvar <- NULL else facetvar <- input$facet_plotvar
      
      
      if(input$select_comp == "pairs") refgrp = NULL else refgrp =input$exprs_statref
      
      
      if(is.null(facetvar)) df = exprs_plotdat() else {
        
        df = subset(exprs_plotdat(), ave(get(input$num_plotvar),
                                         get(facetvar),
                                         get(input$cat_plotvar), FUN = length)>1)
        
      }
      
      
      # Data validation to check if there are enough observations for the t-test
      validate(need(nrow(df) >= 3, "Not enough observations for t-test"))
      
      # to ensure there are enough unique groups in the categorical variable for t-test to work
      validate(need(length(unique(df[[input$cat_plotvar]])) >= 2, "Not enough unique groups in the categorical variable for t-test"))
      
      
      res <-compare_means(form, df,
                          method=input$exprs_statmethod,
                          p.adjust.method = input$exprs_padjmethod,
                          group.by = facetvar, 
                          ref.group = refgrp)
      
      
      res$p.signif.adj <- ifelse(res$p.adj <= 0.0001, "****",
                                 ifelse(res$p.adj <= 0.001, "***",
                                        ifelse(res$p.adj <= 0.01, "**",
                                               ifelse(res$p.adj <= 0.05, "*", "ns"))))
      
      res$y.position <- max(exprs_plotdat()[[input$num_plotvar]], na.rm = T)
      
      
      
      if(input$select_comp == "pairs"){
        
        
        reslev <- level_extractor(res)
        
        selpos <- unlist(lapply(complist(), function(x) which(reslev == paste(levels(factor(x)), collapse = "_"))))
        
        
        res2 <- res[selpos, ]
        
        return(res2)
        
      } else  return(res)   # return stats data frame
      
      
      
    })
    
    
    val3 <- reactiveValues()
    
    output$exprs_plot <- renderPlot({
      if(input$facet_plotvar == "") facetvar <- NULL else facetvar <- input$facet_plotvar
      
      pp <- ggboxplot(exprs_plotdat(), input$cat_plotvar, input$num_plotvar, 
                      fill = input$cat_plotvar,
                      facet.by = facetvar,
                      palette =  input$exprs_plotpal,
                      add = input$exprs_plotadd,
                      add.params = list(size=input$exprs_plotaddsize,
                                        color=input$exprs_plotaddcolor),
                      outlier.shape=NA,
                      font.x=14, font.y=14, font.tickslab = 12, # Adjust axis label and tick font size
                      panel.labs.font = list(size=16)) # Adjust panel label font size
      
      if(input$exprs_rotatex) pp <- pp + rotate_x_text(angle=45)
      
      if(input$exprs_stats) {
        pp <- pp + stat_pvalue_manual(stat_layer(),
                                      label=input$exprs_statlabel, 
                                      step.increase = 0.1, 
                                      step.group.by = facetvar)
        pp <- pp + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
        
        # Move "ns" labels to the top of the plot
        pp <- pp + theme(plot.caption = element_text(hjust = 0.5))
      }
      
      # Determine the number of facets
      num_facets <- ifelse(is.null(facetvar), 0, length(unique(exprs_plotdat()[[facetvar]])))
      
      # Dynamically adjust legend settings based on the number of facets
      if (num_facets > 0) {
        if (num_facets > 1) {
          font_size_legend <- 12
          legend_title_size <- 14
          legend_position <- "bottom"
        } else {
          font_size_legend <- 14
          legend_title_size <- 16
          legend_position <- "bottom"
        }
        
        # Adjust font size based on the number of facets
        font_size_x <- switch(num_facets,
                              1, 10,
                              2, 9,
                              3, 8,
                              4, 7,
                              5, 6,
                              6, 5,
                              7, 5,
                              8, 4,
                              9, 4,
                              10, 4,
                              12)
        
        font_size_y <- font_size_x
        
        # Adjust the size of facet panels based on the number of facets
        panel_spacing <- switch(num_facets,
                                1, 0.2,
                                2, 0.25,
                                3, 0.3,
                                4, 0.35,
                                5, 0.4,
                                6, 0.45,
                                7, 0.5,
                                8, 0.55,
                                9, 0.6,
                                10, 0.65,
                                12, 0.7)
        
        pp <- pp + theme(legend.text = element_text(size = font_size_legend),       # Font size for legend text
                         legend.title = element_text(size = legend_title_size),     # Font size for legend title
                         axis.text.x = element_text(size = font_size_x),            # Font size for x-axis labels
                         axis.text.y = element_text(size = font_size_y),            # Font size for y-axis labels
                         legend.position = legend_position,                         # Move the legend to the bottom of the plot
                         plot.title = element_text(size = 18),                       # Increase the title font size
                         plot.margin = unit(c(1, 1, 2, 1), "lines"),                 # Adjust the bottom margin
                         legend.margin = margin(0, 0, 0, 0),                         # Reduce the legend margin
                         legend.box.margin = margin(0, 0, 0, 0),                     # Reduce the legend box margin
                         legend.spacing = unit(0.1, "lines"),                        # Reduce the legend spacing
                         strip.text = element_text(size = font_size_x),              # Adjust facet text size
                         panel.spacing = unit(panel_spacing, "lines"))               # Adjust the size of facet panels
      } else {
        pp <- pp + theme(legend.text = element_text(size = 14),                    
                         legend.title = element_text(size = 16),
                         axis.text.x = element_text(size = 12),                      
                         axis.text.y = element_text(size = 12),
                         legend.position = "bottom",                                
                         plot.title = element_text(size = 18),                      
                         plot.margin = unit(c(1, 1, 2, 1), "lines"),                 
                         legend.margin = margin(0, 0, 0, 0),                        
                         legend.box.margin = margin(0, 0, 0, 0),                     
                         legend.spacing = unit(0.1, "lines"))                        
      }
      
      val3$pp <- pp
      
      print(pp)
    })
    
    
    
    output$downloadPlot4 <- downloadHandler(
      filename = function() {
        paste("Meta_Data_plot.png")
      },
      content = function(file) {
        
        png(file)
        print(val3$pp)
        dev.off()
        
      })
    
  }
  )
}
