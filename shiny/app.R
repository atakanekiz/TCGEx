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
  
  title = div(
  tags$a(href = "https://tcgex.iyte.edu.tr/"
         
         ,tags$img(src = 'images/left_logo.png', align = 'middle', height = "50px", width = "130px"))
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
  
  tabPanel(
    "HOMEPAGE",
    # tags$img(src = "images/Untitled.png", height = "27px", width = "90px"),
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
    fluidPage(h1("Please select the cancer data"), select_data_ui("module"))
  ),
  tabPanel(
    "KAPLAN-MEIER",
    fluidPage(h1("Kaplan-Meier Survival Analysis "), km_ui("module"))
  ),
  tabPanel(
    "COX-PH",
    fluidPage(h1("Cox Proportional Hazards Survival Analysis"), cox_ui("module"))
  ),
  tabPanel(
    "METADATA ANALYSIS",
    fluidPage(h1("MetaData Analysis"), gene_vs_cat_ui("module"))
  ),
  tabPanel(
    "CORRELATION ANALYSIS",
    fluidPage(h1("Correlation Analysis"), gene_vs_gene_ui("module"))
  ),
  tabPanel(
    "CORRELATED GENE TABLE",
    fluidPage(h1("Correlated Gene Table Analysis "), gene_cor_UI("module"))
  ),
  tabPanel(
    "HEATMAP",
    fluidPage(h1("Heatmap Analysis "), heatmap_ui("module"))
  ),
  tabPanel(
    "GSEA",
    fluidPage(h1("Gene Sets Enrichment Analysis (GSEA) "), gsea_ui("module"))
  ),
  tabPanel(
    "ROC",
    fluidPage(h1("Receiver Operating Characteristic (ROC) Analysis "), roc_ui("module"))
  ),
  tabPanel(
    "PCA",
    fluidPage(h1("Principal Compenent Analysis (PCA)"), pca_ui("module_pca"))
  ),
  tabPanel(
    "MACHINE LEARNING",
    fluidPage(ml_ui("ml"))
  ),
  tabPanel(
    "ABOUT",
    fluidPage(includeHTML("about.html"), shinyjs::useShinyjs())
  ),
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
  
  select_data_server("module",Xproj=Xproj)
  
  heatmap_server("module",Xproj=Xproj)
  
  data_prep_ml_server("ml",Xproj=Xproj)
  df <- data_prep_ml_server("ml",Xproj=Xproj)
  ml_main_server("ml",regress_data  = df,Xproj=Xproj)
  
  pca_server("module_pca",Xproj=Xproj)
  
  roc_server("module",Xproj=Xproj)
  
  gene_vs_gene_server("module",Xproj=Xproj)
  
  cox_server("module",Xproj=Xproj)
  
  gene_vs_cat_server("module",Xproj=Xproj)
  
  km_server("module",Xproj=Xproj)
  
  heatmap_server("module",Xproj=Xproj)
  
  gsea_server("module",Xproj=Xproj)
  
  gene_cor_tb_server("module",Xproj=Xproj)
  gene_cor_pl_server("module",Xproj=Xproj)
  
}

#run the app

shinyApp(ui = ui, server = server)
