#Source
source("select_data.R")
source("heatmap.R")
source("ml_glm.R")
source("roc.R")
source("gene_vs_gene.R")
source("cox_ph.R")
source("pca.R")
source("gene_vs_category.R")
source("kaplan_meier.R")
source("gsea.R")
source("gene_correlations.R")

#Library and Theme

library(shinyjs)
library(shinydashboardPlus)
library(fresh)

library(shiny)
library(shinydashboardPlus)
library(fresh)

mytheme <- create_theme(
  adminlte_color(
    light_blue = "#434C5E"
  ),
  adminlte_sidebar(
    width = "300px",
    dark_bg = "#D8DEE9",
    dark_hover_bg = "#81A1C1",
    dark_color = "#2E3440" 
  ),
  adminlte_global(
    content_bg = "#FFF",
    box_bg = "#D8DEE9", 
    info_box_bg = "#D8DEE9"
  )
)


ui <- navbarPage(  
  
  id = "TCGEx",
    
    title = div(
  tags$a(href = "https://tcgex.iyte.edu.tr/"

         # ,tags$img(src = 'images/left_logo.png', align = 'middle', height = "50px", width = "130px")
         )
),
  
  theme = mytheme,

  tags$head(
    tags$script(src = "https://www.googletagmanager.com/gtag/js?id=G-HWCPP52NZ2"),
    tags$script(
      HTML(
        "window.dataLayer = window.dataLayer || [];
        function gtag(){dataLayer.push(arguments);}
        gtag('js', new Date());
        gtag('config', 'G-HWCPP52NZ2');"
      )
    )
  ),

## gtag('config', 'G-HWCPP52NZ2') # initial tracker for single machine
## gtag('config', 'G-DVM1837XKH') # tracker for srv-1
## gtag('config', 'G-D7TFB6PP0Y') # tracker for srv-2
## gtag('config', 'G-BNNRL1KBD8') # tracker for srv-3
## gtag('config', 'G-0D1PDQNTF9') # tracker for srv-4


# tabsetPanel(id="app_tabset",
#             type = "tabs", ############################################################
  
  tabPanel(
    # "TCGEx",
     tags$img(src = "images/left_logo.png", align = 'top', height = "30px", width = "80px"),
    includeHTML("home.html"),
    tags$script(src = "plugins/scripts.js"),
    tags$head(
      tags$link(rel = "stylesheet", 
                type = "text/css", 
                href = "plugins/font-awesome-4.7.0/css/font-awesome.min.css"),
      tags$link(rel = "icon", 
                type = "image/png", 
                href = "images/logo_tab_icon.png"))
  ),
  tabPanel(
    "DATA SELECTION",
    fluidPage(h1("Please select cancer project(s)"), select_data_ui("seldata")
  )),
  tabPanel(
    "KAPLAN-MEIER",
    fluidPage(h1("Kaplan-Meier Survival Analysis"), km_ui("km"))
  ),
  tabPanel(
    "COX-PH",
    fluidPage(h1("Cox Proportional Hazards Survival Analysis"), cox_ui("cox"))
  ),
  tabPanel(
    "BOXPLOT",
    fluidPage(h1("Boxplot Metadata Analysis"), gene_vs_cat_ui("boxplot"))
  ),
  tabPanel(
    "SCATTERPLOT",
    fluidPage(h1("Scatterplot Correlation Analysis"), gene_vs_gene_ui("scatterplot"))
  ),
  tabPanel(
    "CORRELATED GENES",
    fluidPage(h1("Correlated Gene Table Analysis"), gene_cor_UI("genecor"))
  ),
  tabPanel(
    "HEATMAP",
    fluidPage(h1("Heatmap Analysis"), heatmap_ui("heatmap"))
  ),
  tabPanel(
    "GSEA",
    fluidPage(h1("Gene Sets Enrichment Analysis (GSEA)"), gsea_ui("gsea"))
  ),
  tabPanel(
    "ROC",
    fluidPage(h1("Receiver Operating Characteristic (ROC) Analysis"), roc_ui("roc"))
  ),
  tabPanel(
    "PCA",
    fluidPage(h1("Principal Compenent Analysis (PCA)"), pca_ui("pca"))
  ),
  tabPanel(
    "MACHINE LEARNING",
    fluidPage(ml_ui("ml"))
  ),
  tabPanel(
    "ABOUT",
    fluidPage(includeHTML("about.html"), shinyjs::useShinyjs())
  ),
  
# ), ################################################################### tabsetpanel 

  footer = includeHTML("footer.html")
)


#Main Server

server <- function(input, output, session) {
  
  # # User React
  # lastInteractionTime <- reactiveVal(NULL)
  # 
  # # Update last interaction
  # observeEvent(input, {
  #   lastInteractionTime(Sys.time())
  # })
  # 
  # # Close app after 5 min without interaction
  # autoCloseApp <- function() {
  #   if (!is.null(lastInteractionTime())) {
  #     if (as.numeric(difftime(Sys.time(), lastInteractionTime(), units = "secs")) >= 300) {
  #       session$close()
  #     }
  #   }
  # }
  # 
  # # Start the timer and constantly check the interaction time  
  # observe({
  #   invalidateLater(1000) # Start the timer every 1 second
  #   autoCloseApp() # Check Auto Close
  # })
  
  Xproj<-reactiveValues()
  
  select_data_server("seldata",Xproj=Xproj)

  
  observeEvent(input$TCGEx, {
    if(input$TCGEx == "MACHINE LEARNING") {
      data_prep_ml_server("ml",Xproj=Xproj)
      df <- data_prep_ml_server("ml",Xproj=Xproj)
      ml_main_server("ml",regress_data  = df,Xproj=Xproj)
    }
  })

  observeEvent(input$TCGEx, {
    if(input$TCGEx == "PCA") {pca_server("pca",Xproj=Xproj)}
  })

  observeEvent(input$TCGEx, {
    if(input$TCGEx == "ROC") {roc_server("roc",Xproj=Xproj)}
  })
  
  
  observeEvent(input$TCGEx, {
    if(input$TCGEx == "SCATTERPLOT") {gene_vs_gene_server("scatterplot",Xproj=Xproj)}
  })
  
  observeEvent(input$TCGEx, {
    if(input$TCGEx == "COX-PH") {cox_server("cox",Xproj=Xproj)}
  })
 
  
  observeEvent(input$TCGEx, {
    if(input$TCGEx == "BOXPLOT") {gene_vs_cat_server("boxplot",Xproj=Xproj)}
  })
  
  
  observeEvent(input$TCGEx, {
    if(input$TCGEx == "KAPLAN-MEIER") {km_server("km",Xproj=Xproj)}
    })
  
  
  observeEvent(input$TCGEx, {
    if(input$TCGEx == "HEATMAP") {heatmap_server("heatmap",Xproj=Xproj)}
  })
  
  
  observeEvent(input$TCGEx, {
    if(input$TCGEx == "GSEA") {gsea_server("gsea",Xproj=Xproj)}
  })
  
  
  observeEvent(input$TCGEx, {
    if(input$TCGEx == "CORRELATED GENES") {
      gene_cor_tb_server("genecor",Xproj=Xproj)
      gene_cor_pl_server("genecor",Xproj=Xproj)
    }
  })
  
 
  
}

#run the app

shinyApp(ui = ui, server = server)
