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


#User Interface

ui <- dashboardPage(skin = "black", 
                    title = "TCGExplorer",
                    dashboardHeader(
                      title = tags$img(src='images/tcgex_logo_2.jpg', align = 'middle', height = '50px')
                      ),
                    dashboardSidebar(
                      minified = F,
                      sidebarMenu(
                        menuItem("HOMEPAGE",icon = icon("home"),tabName="home"),
                        menuItem("DATA SELECTION",icon = icon("dna"),tabName="select_data"),
                        menuItem("KAPLAN-MEIER",icon = icon("disease"),tabName="km"),
                        menuItem("COXPH",icon = icon("biohazard"),tabName="cox"),
                        menuItem("GENE-TO-GENE",icon = icon("chart-line"),tabName="gene_vs_gene"),
                        menuItem("GENE-TO-CATEGORY",icon = icon("chart-bar"),tabName="gene_vs_cat"),
                        menuItem("GENE CORRELATIONS",icon = icon("magnifying-glass"),tabName="gene_cor"),
                        menuItem("HEATMAP",icon = icon("microscope"),tabName="heatmap"),
                        menuItem("GSEA",icon = icon("chart-line"),tabName="gsea"),
                        menuItem("PCA",icon = icon("circle-nodes"),tabName="pca"),
                        menuItem("ROC",icon = icon("capsules"),tabName="roc"),
                        menuItem("MACHINE LEARNING",icon = icon("atom"),tabName="ml"),
                        menuItem("ABOUT",icon = icon("address-card"),tabName="about")
                        
                      )
                    ),
                    dashboardBody(
                                  tags$head(
                                            includeHTML(("google-analytics.html"))
                                  ),
                                  
                                  tags$div(
                                    id = "banner",
                                    style = "position: fixed; bottom: 0; left: 0; right: 0; height: 45px; background-color: #AAA4A3; color: white; padding: 10px; width: 100%; z-index: 998;",
                                    tags$p(
                                      style = "font-size: 15px; margin: 0; display: inline-block;",
                                      "We use Google Analytics to give you the best experience on our site and analyze traffic.",
                                      tags$a(
                                        href = "https://policies.google.com/technologies/cookies",
                                        style = "color: white; font-weight: bold; margin-left: 10px; text-decoration: underline;",
                                        "Learn more"
                                      )
                                    ),
                                    tags$button(
                                      id = "close-button",
                                      "I got it!",
                                      style = "float: right; background-color: #ffffff; color: #AAA4A3; border: none; border-radius: 5px; padding: 5px; font-size: 14px; margin-left: 10px; font-weight: bold;",
                                      onclick = "this.parentNode.style.display = 'none';")
                                  ),
                                  
                                  use_theme(mytheme),
                                  tabItems(
                                    tabItem(tabName = "about",
                                            includeHTML("about.html"),
                                            shinyjs::useShinyjs()),
                                    tabItem(tabName="select_data",
                                            fluidPage(h1("Please Select The Cancer Data"),
                                                      select_data_ui("module"))),
                                    tabItem(tabName="heatmap",
                                            fluidPage(h1("Heatmap Analysis "),
                                                      heatmap_ui("module"))),
                                    tabItem(tabName="ml",
                                            fluidPage(h1("Machine Learning Algorithms Analysis"),
                                                      ml_ui("ml"))),
                                    tabItem(tabName="pca",
                                            fluidPage(h1("Principal Compenent Analysis (PCA)"),
                                                      pca_ui("module_pca"))),
                                    tabItem(tabName="roc",
                                            fluidPage(h1("Receiver Operating Characteristic (ROC) Curve Analysis "),
                                                      roc_ui("module"))),
                                    tabItem(tabName="gene_vs_gene",
                                            fluidPage(h1("Feature Correlation Analysis"),
                                                      gene_vs_gene_ui("module"))),
                                    tabItem(tabName="cox",
                                            fluidPage(h1("Cox Proportional Hazards Regression"),
                                                      cox_ui("module"))),
                                    tabItem(tabName="gene_vs_cat",
                                            fluidPage(h1("Metadata Analysis"),
                                                      gene_vs_cat_ui("module"))),
                                    tabItem(tabName="km",
                                            fluidPage(h1("Kaplan Meier Survival Analysis "),
                                                      km_ui("module"))),
                                    tabItem(tabName="gsea",
                                            fluidPage(h1("Gene Sets Enrinchment Analysis (GSEA) "),
                                                      gsea_ui("module"))),
                                    tabItem(tabName="gene_cor",
                                            fluidPage(h1("Correlated Genes Analysis "),
                                                      gene_cor_UI("module"))),
                                    tabItem(tabName="home",
                                             includeHTML("home.html"),
                                             tags$script(src = "plugins/scripts.js"),
                                             tags$head(
                                               tags$link(rel = "stylesheet", 
                                                         type = "text/css", 
                                                         href = "plugins/font-awesome-4.7.0/css/font-awesome.min.css"),
                                               tags$link(rel = "icon", 
                                                         type = "image/png", 
                                                         href = "images/logo_icon.png")))
                                   ),
                                  dashboardFooter(
                                      includeHTML("footer.html")
                                      ) 
                                  
                                    
                                  ))

#Main Server

server <- function(input, output, session) {
  
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