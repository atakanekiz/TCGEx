library(dplyr)
library(Hmisc)
library(janitor)
library(plotly)
library(DT)
library(shinybusy)
library(rintrojs)
library(shinyvalidate)
library(shinyWidgets)
library(stats)
library(openxlsx)
library(corrplot)
library(RColorBrewer)
gene_cor_UI <- function(id) {
  ns <- NS(id)
  tagList(
  
    navbarPage(
      "TCGExplorer",
      tabPanel(
        "Corr Table",
        sidebarPanel(
          
          add_busy_spinner(
            spin = "cube-grid",
            position = "top-right",
            color = "#01303f",
            margins = c(300, 500),
            height = "60px",
            width = "60px"),
          
          selectizeInput(ns("genecor_samp2"), multiple=T,
                         "*Please select sample types",
                         choices=NULL,
                         options=list(placeholder = "eg. Primary solid tumor")),
          
          
          
          selectizeInput(
            inputId = ns("p_gene"),
            label = "*Please select primary gene",
            choices = NULL,
            options=list(placeholder = "eg. TSPAN6")
            
          ),
          numericInput(
            inputId = ns("top_high"),
            label = "How many genes with the highest positive correlation should be shown?",
            value = 10 ,
            min = 2
          ),
          numericInput(
            inputId = ns("top_low"),
            label = "How many genes with the highest negative correlation should be shown?",
            value = 10 ,
            min = 2
          ),
          radioButtons(ns("corr"), "Select correlation calculation method", choices = c("pearson", "spearman")),
          
          
          
          
          actionBttn(inputId = ns("act"),
                     label = "Generate",
                     style = "unite",
                     block = TRUE,
                     color = "primary"),
          br(),
          downloadButton(ns("downloadDat"), "Download results", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
          br(),
          br(),
          #help section UI
          introjsUI(),
          actionButton(ns("intro4"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7")
          
          
        ),
        
        mainPanel(
          textOutput(ns("note")),
          verbatimTextOutput(ns("p_text")),
          DTOutput(ns("g_table"))
        )
        
        
      ),
      tabPanel(
        "Corr Plot",
        sidebarPanel(
          selectizeInput(ns("genecor_samp3"), multiple=T,
                         "*Please select sample types",
                         choices=NULL,
                         options=list(placeholder = "eg. Primary solid tumor")),
          radioButtons(ns("corr2"), "Select correlation calculation method", choices = c("pearson", "spearman")),
          checkboxInput(ns("coef"), "Show correlation coefficient on the plot", value = TRUE),
          radioButtons(ns("gen_sel"), "Select gene selection method", choices = c("Choose manually", "Upload File")),
          
          conditionalPanel(
            condition = "input.gen_sel == 'Choose manually' ", ns = ns,
            selectizeInput(
              inputId = ns("p_gene2"),
              label = "*Select at least 2 genes",
              choices = NULL,
              multiple = T,
              options=list(placeholder = "eg. TSPAN6")
              
            )
            
          ),
          
          conditionalPanel(
            
            condition = "input.gen_sel == 'Upload File' ", ns = ns,
            
            fileInput(inputId = ns("corr_up"),
                      label = "Upload your gene set as .csv file",
                      accept = c("text/csv", "text/comma-separated-values,text/plain",
                                 ".csv")),  
          ),
          
          actionBttn(inputId = ns("act2"),
                     label = "Generate",
                     style = "unite",
                     block = TRUE,
                     color = "primary"),
          br(),
          #help section UI
          introjsUI(),
          actionButton(ns("intro5"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7")
          
          
        
        ),
        mainPanel(
          plotOutput(ns("corre_plot"),
                     width = "125%",
                     height = "600px")
        )
      )
    )
    
    
    
    
    
    
    
  )
}

data_prep <- function(df){
  
  df_prep <- df %>% 
    select(!contains("meta.")) %>% 
    select(where(is.numeric))
  
  df_prep <- na.omit(df_prep)
  
  df_prep <- remove_constant(df_prep)
  
  df_prep
  
}



gene_cor_tb_server <- function(id,Xproj) {
  moduleServer(
    id,
    function(input, output, session) {
      
      
      
      ## help section server
      
      ns <- session$ns
      
      inf <- reactive({
          
            data.frame(
                
                element = paste0("#", session$ns(c(NA, "genecor_samp2 + .selectize-control","p_gene + .selectize-control ", "top_high", "top_low", "corr"))),
                
                intro = paste(c(
                  "This is the Gene Correlation Analysis module. In this first tab, you can select a gene and tabulate its top positively and negatively correlated genes. You can also visualize correlations in the second tab above. Continue the tutorial to learn the features of this module. <b>NOTE:</b> Since pair-wise correlation is calculated genome-wide, this analysis can take a some time.",
                  "Select sample type (eg. primary and/or metastatic) to focus the analysis on specific patient subsets.",
                  "Select your gene of interest here.",
                  "Here you can change the number of top positively correlated genes shown.",
                  "You can also change the number of top negatively correlated genes shown.",
                  "Specify how the correlation should be calculated."

                ))
              )
              
           
            
         
        })
        
       
        observeEvent(input$intro4, {
        
        introjs(session, options = list(steps = inf() ) )
        
      })
      
      #shinyvalidate
      
        num_iv <- InputValidator$new()
        num_iv$add_rule("top_high", ~ if (input$top_high < 2  & !anyNA(input$top_high)) "Choice must be equal to or greater than 2")
        num_iv$add_rule("top_low", ~ if (input$top_low < 2 & !anyNA(input$top_low)) "Choice must be equal to or greater than 2")  
        
        iv <- InputValidator$new()
        iv$add_validator(num_iv)
        
        iv$enable()
        
        
      ##inputs
      
      observe({updateSelectizeInput(session, "genecor_samp2",choices = as.data.frame(Xproj$a())[["meta.definition"]], server = T)})
      observe({updateSelectizeInput(session, 'p_gene', choices = colnames(gene_choices()), server = TRUE, selected = "")})
      
      

      
      
      ## input choices for gene selection
      
      gene_choices <- reactive({
        
        Xproj$a() %>%
          select(!starts_with("meta.")) %>% 
          select(where(is.numeric))
        
      })
      
      
      ## data preparation 
      
      
      pre_data <- reactive({
        pre <- as.data.frame(Xproj$a()) %>% 
        filter(meta.definition %in% input$genecor_samp2) %>%
        select(!starts_with("meta.")) %>% 
        select(where(is.numeric))
        
        pre
        
      })
      
      
      f_data_1 <- reactive({data_prep(pre_data())})
      
      f_data_2 <- reactive({
        
        f_data_1() %>% 
           
           select(-input$p_gene)
        
        
      })
      
     
      
      cor_dat <-  reactive({
       
       data_cor_a <- cor(as.matrix(f_data_1()[,input$p_gene]),as.matrix(f_data_2()), method = input$corr)
       
       data_cor_b <- sapply(f_data_2(), FUN=function(x, y) cor.test(x, y, method = input$corr)$p.value, y=f_data_1()[[input$p_gene]])
        
        data_cor_a = t(data_cor_a)
        
        data_cor_a = as.data.frame(data_cor_a)
        
        data_cor_b = as.data.frame(data_cor_b)
        
        
        
        
        colnames(data_cor_a)[1]  <- "correlation_coefficient"
        colnames(data_cor_b)[1]  <- "p_value"
        
        
        data_cor = cbind(data_cor_b, data_cor_a)
        
        data_cor <- data_cor %>% mutate(across(is.numeric, signif, digits=3))
        
        data_cor$Genes <- rownames(data_cor)
        
        data_cor
      })
      
     
     ##high gene selection 
     
     high_genes <- reactive ({ cor_dat() %>%                                      
       arrange(desc(correlation_coefficient)) %>% 
       head(n= input$top_high)
     })
     
     ##low gene selection 
     
      low_genes <- reactive({cor_dat() %>%                                      
       arrange(desc(correlation_coefficient)) %>% 
       tail(n= input$top_low)
      })
        
      ## high-low table prep
      
      hl_table <- reactive(rbind(high_genes(), low_genes()))
     
     
      p_table <- reactive ({ pre_data() %>% 
        select(hl_table()$Genes)
      
      })
      
      p_table_2 <- reactive ({data.frame(zero_patient = colSums(p_table()==0, na.rm = T)) })
      
        
        gene_table <- reactive({
        
          gene_tab <- cbind(hl_table(), p_table_2())
        
          gene_tab <- gene_tab %>% 
            relocate(Genes)%>% 
            mutate_if(is.numeric, signif, digits=3)
          
        })
     
    
        
      
       observeEvent(input$act, {
         

         
         ##table formation
         
         output$g_table <- renderDataTable({
           
           req(input$act)
           input$act
           
           isolate({
             
             req(iv$is_valid()) 
             
             validate(
               need(input$genecor_samp2, 'Please select at least one sample type'),
               need(input$p_gene, 'Please select a gene'),
               need(input$top_high, 'Please choose how many top positive correlators should be shown'),
               need(input$top_low, 'Please choose how many top negative correlators should be shown')
             )
             
             gene_table()}) 
           
         })
 
 
         
         ##message 
        
         output$p_text <- renderText({
           
           req(input$act)
           input$act
           
           isolate({

             req(iv$is_valid())
             
             validate(
               need(input$genecor_samp2, ''),
               need(input$p_gene, ''),
               need(input$top_high, ''),
               need(input$top_low, '')
             )
             
             sum(pre_data()[[input$p_gene]] == 0, na.rm = T)})
           
         })
        
         
         
         
         
         
         ##note
          
         
         output$note <- renderText({
           
           req(input$act)
           input$act
           
           isolate({
             req(iv$is_valid())
             validate(
               need(input$genecor_samp2, ''),
               need(input$p_gene, ''),
               need(input$top_high, ''),
               need(input$top_low, '')
             )
             
             "Number of samples where the chosen gene is not expressed"
           })
           
         })
         
         
         #download 
         
         output$downloadDat <- downloadHandler(
           filename = "gene_correlations.xlsx",
           content = function(file) {
             write.xlsx(gene_table(), file,colnames = TRUE,
                        rownames = F, append = FALSE, showNA = TRUE)
           }
         )
         
         
      })
      
    }
  )
}



gene_cor_pl_server <- function(id,Xproj) {
  moduleServer(
    id,
    function(input, output, session) {
    
      ns <- session$ns
      
      cor_help <- reactive({
        
        data.frame(
          
          element = paste0("#", session$ns(c(NA, "genecor_samp3 + .selectize-control","corr2", "gen_sel" ))),
          
          intro = paste(c(
            "Here, you can visualize gene-to-gene correlations. Continue the tutorial to learn the features of this module.",
            "Select sample type (eg. primary and/or metastatic) to focus the analysis on specific patient subsets.",
            "Specify how the correlation should be calculated.",
            "You can choose genes manually or you can upload the gene list as a csv file. You can create csv file in spreadsheet software by providing your genes in the first column. Please only use gene symbols and ensure that you are using the correct capitalization."
          ))
        )
        

      })
      
      
      observeEvent(input$intro5, {
        
        introjs(session, options = list(steps = cor_help() ) )
        
      })
      
      
      #shinyvalidate
      
     iv <- InputValidator$new()
     iv$add_rule("p_gene2", ~ if (length(input$p_gene2) == 1  & !anyNA(length(input$p_gene2)) & input$gen_sel == 'Choose manually') "You must choose at least 2 genes")
     
      iv$enable()
      
      observe({updateSelectizeInput(session, "genecor_samp3",choices = as.data.frame(Xproj$a())[["meta.definition"]], server = T)})
      observe({updateSelectizeInput(session, 'p_gene2', choices = colnames(gene_choices_2()), server = TRUE, selected = "")})
      
      ## input choices for gene selection
      
      gene_choices_2 <- reactive({
        
        Xproj$a() %>%
          select(!starts_with("meta.")) %>% 
          select(where(is.numeric))
        
      })
      
      
      #gene list upload

      g_list <- reactive({
        
        req(input$corr_up)

        glist <- read.table(input$corr_up$datapath, header = FALSE, sep = ",", quote = "\"'",
                            dec = ".")

        colnames(glist)[1] <- "genes"

        glist

      })
      
      #correlation plot data 
      
      pre_cor_dat <- reactive({
        if(input$gen_sel == 'Choose manually'){
          
          cor_df_a <- Xproj$a() %>% 
            filter(meta.definition %in% input$genecor_samp3)%>% 
            select(input$p_gene2) %>% 
            drop_na() 
          
          
          cor_df_a
          
        }else if(input$gen_sel == 'Upload File'){
        
          cor_df_b  <- Xproj$a() %>% 
            filter(meta.definition %in% input$genecor_samp3)   
          
          cor_df_b <- cor_df_b [,intersect(g_list()[["genes"]], colnames(cor_df_b)),with=F]
          
          cor_df_b <- cor_df_b %>% drop_na() 
          
          cor_df_b
        }
        
      }) 
      
      #correlation calculation
      
      corre_dat <- reactive({
        
        cor_pd <- cor(pre_cor_dat(), method = input$corr2)
        
      })
      
      #plot 
      
      corr_pl <- reactive({
        
        if(input$coef){
          
          corrplot(corre_dat(), 
                   method= "circle",
                   order="hclust", 
                   type = 'upper', 
                   addCoef.col = "black",
                   tl.col="black", 
                   tl.srt=45, 
                   number.digits = 2,
                   col=brewer.pal(n=8, name="BrBG"))
        }else {
          
          corrplot(corre_dat(), 
                   method= "circle",
                   order="hclust", 
                   type = 'upper', 
                   addCoef.col = NULL ,
                   tl.col="black", 
                   tl.srt=45, 
                   number.digits = 2,
                   col=brewer.pal(n=8, name="BrBG"))
        }
        

        
      })
      
      
      observeEvent(input$act2, {
        
          output$corre_plot <- renderPlot({
            
            
            req(input$act2)
            input$act2
            
            isolate({
              req(iv$is_valid())
              validate(
                need(input$genecor_samp3, "Choose at least 1 sample type")
                
              )
              
              if(input$gen_sel == 'Choose manually'){
                
                validate(need(input$p_gene2, "Choose at least 2 genes"))
              }
              
              if(input$gen_sel == 'Upload File'){
                
                validate(need(input$corr_up, "Please upload your file"))
              }
              
              corr_pl()
              
            })
          
        })
          
      })
      
    }
  )
}