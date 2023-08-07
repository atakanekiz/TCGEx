#Load Libraries

library(dplyr)
library(shiny)
library(shinydashboard)
library(ggplot2)
library(broom)
library(DT)
library(msigdbr)
library(data.table)
library(janitor)
library(shinybusy)
library(fgsea)
library(gridExtra)
library(plyr)
library(openxlsx)
library(grid)
library(readxl)

#Load Sources

source("functions/gene_ranker.R")
source("functions/gsea_replot.R")
source("functions/plotGseaTable2.R")
source("functions/top_plotter.R")


##GSEA ANALYSIS

gsea_ui <- function(id, label, choices) {
  
  ns <- NS(id)
  
  tagList(
    
    ui <- fluidPage(
      add_busy_spinner(
        spin = "cube-grid",
        position = "top-right",
        color = "#01303f",
        margins = c(300, 500),
        height = "60px",
        width = "60px"),
      
      sidebarPanel(
        
        
        selectizeInput(inputId = ns ("gsea_samptyp"), 
                       multiple=T,
                       label = "1. Select sample types",
                       choices=NULL, # will be updated dynamically
                       options=list(placeholder = "eg. Primary solid tumor")),
        hr(),
        
        selectizeInput(inputId = ns("gsea_feat"),
                       label = "2. Select feature",
                       choices = NULL, # will be updated dynamically
                       options=list(placeholder = "eg. IFNG or meta.tissue_type",
                                    plugins = list('restore_on_backspace'))),
        
        conditionalPanel(
          
          condition = "output.gsea_var_status",
          ns=ns, 
          
          numericInput(inputId = ns("gsea_high_cutoff"), "High cutoff percent (Sample Group Left)", 50, 1, 100),
          
          numericInput(inputId = ns("gsea_low_cutoff"), "Low cutoff percent (Reference Group Right)", 50, 1, 100)
        ),
        
        conditionalPanel(
          
          condition = "output.gsea_var_status2",
          ns=ns, 
          selectizeInput(inputId = ns ("gsea_sel_feat_meta_groups"), "Select samples", multiple=T,
                         choices = NULL, # will be updated dynamically
                         options = list(placeholder = "eg. female,male", maxItems = 2)),
          
          selectizeInput(inputId = ns ("gsea_sel_feat_meta_groups2"), "Select Sample ID  ", multiple=F,
                         choices = NULL, # will be updated dynamically
                         options = list(placeholder = "eg. male"))
          
        ),
        
        hr(),
        
        numericInput(inputId = ns("nperm"), 
                     label = tags$span(
                       "3. nPerm Value",
                       tags$i(
                         class = "glyphicon glyphicon-info-sign",
                         style = "color:#0072B2;",
                         title = "If you increase the nPerm Value, the analysis takes more time"
                       )),
                     min = 1000, 
                     step = 100, 
                     value = 1000),
        
        radioButtons(ns("gsea_gene_sets"), "Choose gene set collection", choices = c("MSigDB", "Custom Gene Set")),
        
        conditionalPanel(
          
          condition = "input.gsea_gene_sets == 'MSigDB'",
          ns=ns,
          
          
          selectizeInput(ns("gsea_cat"), "Please select an MSigDB Collection", choices = c("Hallmark gene sets (H)" = "H",
                                                                                          "Positional gene sets (C1)" = "C1",
                                                                                          "Curated gene sets (C2)" = "C2",
                                                                                          "Regulatory target gene sets (C3)" = "C3",
                                                                                          "Computational gene sets (C4)" = "C4",
                                                                                          "Ontology gene sets (C5)" = "C5" ,
                                                                                          "Oncogenic gene sets (C6)" = "C6",
                                                                                          "Immunologic gene sets (C7)" = "C7",
                                                                                          "Cell type signature gene sets (C8)" = "C8")),
          conditionalPanel(condition = "input.gsea_cat == 'C2'|input.msigdb_setnames_response =='C3'|
                                      input.gsea_cat =='C4'|
                                      input.gsea_cat =='C5'|
                                      input.gsea_cat =='C7' ", ns = ns, 
                           selectizeInput(ns("gsea_subcat"),"Please select a subcategory" ,choices = c(""))
          ),
          
          
          radioButtons(ns("individual"), "Show", choices = c("Top Pathways", "Specific Pathway"), selected = "Specific Pathway"),
          
          conditionalPanel(
            
            condition = 'input.individual == "Specific Pathway"',
            ns=ns,
            
            
            selectizeInput(ns("gsea_pathway"), 
                           "*Please select a pathway", 
                           options=list(placeholder = "eg. HALLMARK_DNA_REPAIR"),
                           choices = NULL)
            
          )
          
          
          
        ),
        
        conditionalPanel(
          
          condition = "input.gsea_gene_sets == 'Custom Gene Set'",
          ns=ns,
          
          fileInput(inputId = ns("gset_up"),
                    label = tags$span(
                      "Please upload your xlsx/xls file.",
                      tags$i(
                      class = "glyphicon glyphicon-info-sign",
                      style = "color:#0072B2;",
                      title = "The xlsx/xls file should contain two unnamed columns: the first column should contain the gene set name, and the second column should contain human gene names. Each gene should be associated with a gene set (ie. no missing data), and multiple gene sets can be provided in one file."
                    ),tags$br(),
                    a(href="sample_gsea_input.xlsx", "Sample Input File", download=NA, target="_blank")),
                    accept = c(".xls", ".xlsx" # "text/csv", "text/comma-separated-values,text/plain", ".csv"
                               )),  
          radioButtons(ns("individual_2"), "Show", choices = c("Top Pathways", "Specific Pathway"), selected = "Top Pathways"),
          
          conditionalPanel(
            
            condition = 'input.individual_2 == "Specific Pathway"',
            ns=ns,
            
            
            selectizeInput(ns("cre_path"), 
                           "*Please select a pathway", 
                           options=list(placeholder = "eg. EXAMPLE_PATHWAY"),
                           choices = NULL)
            
          )
          
          
        ),
        
        
        
        actionBttn(inputId = ns("gsea_run"), 
                   label = "Analyze",
                   style = "unite",
                   block = TRUE,
                   color = "primary"),
        br(),
        
        numericInput(inputId = ns("gsea_width"), "Choose the width of the plot to download", min = 100, step = 10, value = 480),
        numericInput(inputId = ns("gsea_height"), "Choose the height of the plot to download", min = 100, step = 10, value = 480),
        
        downloadButton(ns("g_downloadPlot"), "Download GSEA Plot", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
        br(),
        br(),
        downloadButton(ns("download_lead"), "Download Leading Edge Genes ", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
        br(),
        br(),
        downloadButton(ns("downloadData"), "Download ranked data", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
        br(),
        br(),
        downloadButton(ns("download_fgsea"), "Download GSEA result", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
        br(),
        br(),
        #help section UI
        
        introjsUI(),
        actionButton(ns("intro3"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
        
        width = 3
        
      ),
      
      mainPanel(
        
        plotOutput(outputId = ns("gsea_plot")),
        
        DTOutput(ns("gene_text"))
        
      )
    )
  )
} 

gsea_server <- function(id,Xproj) {
  moduleServer(id,function(input, output, session){
    
    ns <- session$ns
    
    ## msigdb_database reading
    
    msigdb_gene_sets =  reactive({readRDS(paste0("genesets/", "msigdb_long", ".rds"))})
    
    ## help section server
    
    help_sec <- reactive({
      
      
      
      if(input$gsea_gene_sets == 'Custom Gene Set') {
        
        data.frame(
          
          element = paste0("#", session$ns(c(NA, "gsea_samptyp + .selectize-control", "gsea_feat + .selectize-control ", "nperm", "gsea_gene_sets", "individual_2"))),
          
          intro = paste(c(
            "This is the gene set enrichment analysis (GSEA) module. Here, you can categorize samples based on a custom criteria and examine whether previously defined or user-provided gene sets are enriched in either of the data subset. Continue with the tutorial to learn more about this module. <b>NOTE:</b> This module can take some time depending on the data set and user selections.",
            "You can select which sample types should be included in the analysis (eg. primary and/or metastatic).",
            "GSEA is performed between two groups of data. If you would like to perform GSEA for a categorical clinical feature, you are expected to select two data subsets and define one of them as the 'sample' for the analysis (the other one will become reference). If you would like to perform GSEA for a numerical feature such as gene expression, then you can categorize samples based on gene expression values as 'high' and 'low' through user-defined quantiles. Setting high and low cutoff to 50, will categorize gene expression at the median value. You can set these numbers to 25 to compare the top 25% expressors to bottom 25% expressors.", 
            "You can increase the number of permutations for preliminary estimation of P-values. If you see NA's in enrichment scores and p-values during your analyses, you can try increasing this number. Note that the calculation time will be longer accordingly.",
            "You can perform GSEA using previously defined gene sets from the <a href='https://www.gsea-msigdb.org/gsea/msigdb/'>Molecular Signatures Database (MSigDB)</a> or provide your own gene sets",
            "If you select the 'Top Pathways' option, pathways with the highest and lowest enrichments are shown </i>(they may not always be significant!)</i>. If you choose 'Specific Pathway', the enrichment plot is prepared only for the pathway you select."
          ))
          
        )
        
      }else if (input$gsea_gene_sets == 'MSigDB'){
        
        data.frame(
          
          element = paste0("#", session$ns(c(NA, "gsea_samptyp + .selectize-control", "gsea_feat + .selectize-control ", "nperm","gsea_gene_sets" ,"gsea_cat+ .selectize-control ", "individual"))),
          
          intro = paste(c(
            "This is the gene set enrichment analysis (GSEA) module. Here, you can categorize samples based on a custom criteria and examine whether previously defined or user-provided gene sets are enriched in either of the data subset. <b>NOTE:</b> This module can take some time depending on the data set. Continue with the tutorial to learn more about this module.",
            "You can select which sample types should be included in the analysis (eg. primary and/or metastatic).",
            "GSEA is performed between two groups of data. If you would like to perform GSEA for a categorical clinical feature, you are expected to select two data subsets and define one of them as the 'sample' for the analysis (the other one will become reference). If you would like to perform GSEA for a numerical feature such as gene expression, then you can categorize samples based on gene expression values as 'high' and 'low' through user-defined quantiles. Setting high and low cutoff to 50, will categorize gene expression at the median value. You can set this numbers to 25 to compare the top 25% expressors to bottom 25% expressors.",
            "You can increase the number of permutations for preliminary estimation of P-values.",
            "You can perform GSEA using previously defined gene sets from the <a href='https://www.gsea-msigdb.org/gsea/msigdb/'>Molecular Signatures Database (MSigDB)</a> or provide your own gene sets",
            "Select MSigDB gene set collection to use in the analyses <b>NOTE:</b> Since some of these collections (such as curated and GO gene sets) contain a large number of entries, analysis can take some time.",
            "If you select the 'Top Pathways' option, pathways with the highest and lowest enrichments are shown </i>(please note that they may not always be significant)</i>. If you choose 'Specific Pathway', the enrichment plot is prepared only for the pathway you select."
          ))
          
        )
        
      }
      
      
      
      
    })
    
    observeEvent(input$intro3, {
      
      introjs(session, options = list(steps = help_sec() ) )
      
    })
    
    #shinyvalidate
    observe(req(input$gsea_feat))
    
    panel_iv <- InputValidator$new()
    panel_iv$condition(~ input$gsea_feat %in% colnames(con_dat()))  
    panel_iv$add_rule("gsea_sel_feat_meta_groups", ~ if (length(input$gsea_sel_feat_meta_groups) < 2 & !anyNA(input$gsea_sel_feat_meta_groups)) "Please select one more data subsets")
    
    numeric_iv <- InputValidator$new()
    numeric_iv$condition(~ input$gsea_feat %in% colnames(con_dat2()))
    numeric_iv$add_rule("gsea_high_cutoff", ~ if (input$gsea_high_cutoff + input$gsea_low_cutoff > 100 & !anyNA(input$gsea_low_cutoff) & !anyNA(input$gsea_high_cutoff)) "For proper categorization, high and low cutoffs can't exceed 100 when added together")
    numeric_iv$add_rule("gsea_low_cutoff", ~ if (input$gsea_high_cutoff + input$gsea_low_cutoff > 100 & !anyNA(input$gsea_low_cutoff) & !anyNA(input$gsea_high_cutoff)) "For proper categorization, high and low cutoffs can't exceed 100 when added together")
    
    nperm_iv <- InputValidator$new()
    nperm_iv$add_rule("nperm", ~ if (input$nperm < 1000 & !anyNA(input$nperm)) "nPerm Value must be greater than or equal to 1000")
    
    
    iv <- InputValidator$new()
    iv$add_validator(panel_iv)
    iv$add_validator(numeric_iv)
    iv$add_validator(nperm_iv)
    iv$enable()
    
    ##choice reactive for reference ID   
    
    reference <- reactive({
      
      input$gsea_sel_feat_meta_groups[input$gsea_sel_feat_meta_groups != input$gsea_sel_feat_meta_groups2]
      
    })  
    
    
    
    ## uploaded data 
    
    gdata <- reactive({
      
      req(input$gset_up)
      
      gene_s <- as.data.frame(read_excel(input$gset_up$datapath, sheet = 1, col_names = F))
      

      validate(
        need(
          {if(length(colnames(gene_s))== 2) TRUE else FALSE},
          "This file does not contain the appropriate number of columns.Please ensure that the first column contains gene set name and the second column contains gene names"))
      
      colnames(gene_s) <- c("geneset_name", "genes")
      
      validate(
        need(
          {if(class(gene_s$geneset_name) == "character" & class(gene_s$genes) == "character")TRUE else FALSE}, 
          "The columns must contain only character inputs.Please ensure that the first column contains gene set names and the second column contains gene names"))
      
      custom_sets <- unique(gene_s$geneset_name)
      
      custom_list <- list()
      
      
      for(i in custom_sets) {
        
        custom_list[[i]] <- gene_s %>%
          filter(geneset_name == i) %>% 
          pull(genes)
        
        
      }
      
      custom_list
      
    })
    
    
    
    ## data frames for condition
    
    con_dat <- reactive({
      
      req(Xproj$a())
      
      Xproj$a() %>% select(-where(is.numeric))
      
    })
    
    con_dat2 <- reactive({
      
      Xproj$a() %>% select(where(is.numeric))
      
    })
    
    observe({updateSelectizeInput(session, 
                                  "gsea_samptyp",
                                  choices = Xproj$a()$meta.definition, 
                                  server = T)})
    
    observe({updateSelectizeInput(session, 
                                  "gsea_feat", selected="",
                                  choices = colnames(Xproj$a()), 
                                  server = T)})
    
    observe({updateSelectizeInput(session,
                                  "gsea_sel_feat_meta_groups", selected = "", 
                                  choices = levels(as.factor(Xproj$a()[[input$gsea_feat]])),
                                  server = T)})
    
    observe({updateSelectizeInput(session,
                                  "gsea_sel_feat_meta_groups2", selected = "", 
                                  choices = input$gsea_sel_feat_meta_groups,
                                  server = T)})
    
    
    observe({
      req(input$gsea_cat)
      
      pc_subcat = names(msigdb_gene_sets()[[input$gsea_cat]])
      if (length(pc_subcat) > 1)  {
        updateSelectizeInput(session,'gsea_subcat', choices = pc_subcat , server = TRUE)
      } 
    })
    
    
    
    
    observe({
      if(input$gsea_cat %in% c("C2","C3","C4","C5","C7")) {
        
        gsea_path = names(msigdb_gene_sets()[[input$gsea_cat]][[input$gsea_subcat]])
        updateSelectizeInput(session,'gsea_pathway', choices = gsea_path , server = TRUE)
        
        
      } else {
        
        gsea_path = names(msigdb_gene_sets()[[input$gsea_cat]][[]])
        updateSelectizeInput(session,'gsea_pathway', choices = gsea_path , server = TRUE)
      }
      
    })
    

    
    observe({updateSelectizeInput(session, 
                                  "cre_path",selected="",
                                  choices = names(gdata()) , 
                                  server = T)})
    

    # Prompt conditional panel if variable is a gene
    output$gsea_var_status <- reactive({
      
      is.numeric(Xproj$a()[[input$gsea_feat]])
      
    })
    
    outputOptions(output, "gsea_var_status", suspendWhenHidden = FALSE)
    
    # Prompt conditional panel if variable is a categorical variable
    output$gsea_var_status2 <- reactive({
      
      is.character(Xproj$a()[[input$gsea_feat]]) | is.factor(Xproj$a()[[input$gsea_feat]])
      
    })
    
    outputOptions(output, "gsea_var_status2", suspendWhenHidden = FALSE)
    
    ## Data preparation for gsea_plotter
    
    gsea_dat <- reactive({ 
      
      dat <- Xproj$a()
      
      # De-duplicate patients. Expand on this to drop certain sample groups when duplicates are present
      dat <- dat[!duplicated(meta.patient), ]
      
      dat <- dat[meta.definition %in% input$gsea_samptyp, ]
      
      
      if(is.numeric(dat[[input$gsea_feat]])){
        
        
        dat <- dat[ , .SD, .SDcols = is.numeric]
        
        gsea_feat_zero_count <<- sum(dat[[input$gsea_feat]] == 0, na.rm=T)
        gsea_feat_na_count <<- sum(is.na(dat[[input$gsea_feat]]))
        
        dat[, sample := ifelse(dat[[input$gsea_feat]] >= quantile(dat[[input$gsea_feat]], (100-input$gsea_high_cutoff)/100, na.rm = T), "high",
                               ifelse(dat[[input$gsea_feat]] <= quantile(dat[[input$gsea_feat]], input$gsea_low_cutoff/100, na.rm = T), "low","empty"))]
        
        dat <- dat[, (input$gsea_feat):=NULL]
        
        
        dat$sample <- as.character(dat$sample)
        
        dat <- select(dat, -starts_with("meta."))
        
        
        dat = as.data.frame(dat)
        
        return(dat)
        
      }else {
        
        
        dat <- dat[get(input$gsea_feat) %in% input$gsea_sel_feat_meta_groups, ]
        
        
        dat <- setnames(dat, (input$gsea_feat) , "sample")
        
        dat$sample <- as.character(dat$sample)
        
        dat <- dat[, c(.(sample = sample),.SD), .SDcols = is.numeric]
        
        dat <- select(dat, -starts_with("meta."))
        
        dat = as.data.frame(dat)
        
        return(dat)
        
      }
      
    })
    
    
    
    
    plot_ind <- reactive({
      
      if(input$individual == "Specific Pathway" && input$gsea_gene_sets == 'MSigDB'){
        
        input$gsea_pathway
        
      }else if(input$individual == "Top Pathways" && input$gsea_gene_sets == 'MSigDB'){
        
        NULL
        
      }else if(input$individual_2 == "Specific Pathway" && input$gsea_gene_sets == 'Custom Gene Set'){
        
        input$cre_path
        
      }else if(input$individual_2 == "Top Pathways" && input$gsea_gene_sets == 'Custom Gene Set'){
        
        NULL
        
      }
      
    })  
    
    #Preranked genes
    
    preranked_genes <- reactive({
      
      if(is.numeric(Xproj$a()[[input$gsea_feat]])) {
        
        p <- gene_ranker(exprs = gsea_dat(),
                         sample_id = "high",
                         reference_id = "low",
                         save_ranked_object = T,
                         method = 's2n',
                         verbose = T)
        
        p
        
      }else {
        
        p <-  gene_ranker(
          gsea_dat(),
          sample_id = input$gsea_sel_feat_meta_groups2 , 
          reference_id = reference(),
          save_ranked_object = T,
          method = 's2n',
          verbose = T
        )
        
        p
        
      }
    })
    
    # ranked data as data table 
    
    
    ranked_data <- reactive({
      
      
      df <- ldply (preranked_genes(), data.table)
      
      colnames(df)<- c("Genes", "Ranks")
      
      df <- df %>% 
        mutate_if(is.numeric, signif, digits=3) %>% 
          mutate(Gene_Order = seq_along(df$Ranks)) %>% 
               relocate(Gene_Order) %>% 
                  select(c("Gene_Order", "Genes"))
      
      
      df
      
      
    })
    
    
    
    
    
    #gene set 
    
    gene_set <- reactive({
      
      if (input$gsea_gene_sets == 'MSigDB') {
        
        if(input$individual == "Specific Pathway"){
          
         
          if(input$gsea_cat %in% c("C2","C3","C4","C5","C7")) {
            
            gene_set <- list(msigdb_gene_sets()[[input$gsea_cat]][[input$gsea_subcat]][[input$gsea_pathway]])
            
            names(gene_set) <- input$gsea_pathway
            
            gene_set
            
            
          } else {
            
            gene_set <- list(msigdb_gene_sets()[[input$gsea_cat]][[]][[input$gsea_pathway]])
            
            names(gene_set) <- input$gsea_pathway
            
            gene_set
          }
          
        }else if(input$individual == "Top Pathways") {
          
          if(input$gsea_cat %in% c("C2","C3","C4","C5","C7")) {
            
            gene_set <- msigdb_gene_sets()[[input$gsea_cat]][[input$gsea_subcat]]
            
            gene_set
            
            
          } else {
            
            gene_set <- msigdb_gene_sets()[[input$gsea_cat]][[]]
            
            gene_set
          }
          
        }
        
        
      } else if (input$gsea_gene_sets == "Custom Gene Set"){
        
        gene_set <- gdata()
        
      }
      
    })
    
    
    #sample_ id
    
    sample_id <- reactive({
      
      if(is.numeric(Xproj$a()[[input$gsea_feat]])) {
        
        # print("high")
        
        return("high")
        
      }else {
        
        input$gsea_sel_feat_meta_groups2
        
      }
      
    })
    
    #reference_id 
    
    reference_id <- reactive({
      
      if(is.numeric(Xproj$a()[[input$gsea_feat]])) {
        
        # print("low")
        
        return("low")
        
      }else {
        
        reference()
        
      }
      
    })
    
    
    #fgsea reactive
    
    res <- reactive({
    
    fgsea(pathways = gene_set(), stats = preranked_genes(), minSize = 1, maxSize = ncol(gsea_dat())-3,  nPermSimple = input$nperm)
    
    }) 
    
    gsea_res <- reactive({
      
      g_res = data.frame(res())
      
      g_res = select(g_res, c("pathway", "pval", "padj", "log2err", "ES", "NES"))
      
    })
  
    
    
    
    # plot 
    
    gsea_plot <- reactive({
      
      plot_individual = plot_ind()
      
      verbose = T 
      
      set.seed(123)
      
      if (is.null(plot_individual)) {
        
        
        plot_title <- paste(sample_id(), "vs", reference_id())
        
        
        
        plot_grob <- top_plotter(gsea_results = res(), ranked_genes = preranked_genes(), gene_set = gene_set(), top_n = 10, gseaParam = 1,
                                 plot_title = plot_title, do.plot = F)
        grid.arrange(plot_grob)
        
        
      } else {
        
        
        if (verbose)
          message(paste("Plotting", plot_ind()))
        
        annot_pval <- signif(as.numeric(res()[res()$pathway == plot_ind(), "pval"]), digits = 2)
        annot_NES <- signif(as.numeric(res()[res()$pathway == plot_ind(), "NES"]), digits = 2)
        annot_ES <- signif(as.numeric(res()[res()$pathway == plot_ind(), "ES"]), digits = 2)
        x_pos <- length(preranked_genes())/4
        
        
        annot_text <- paste("p: ", annot_pval, "\nNES:", annot_NES)
        
        
        
        plot_grob <- plotEnrichment(pathway = gene_set()[[plot_ind()]], stats = preranked_genes()) + 
          labs(title = plot_ind(), subtitle = paste(sample_id(),"vs", reference_id() 
                                                
          )) + 
          annotate("text", x = x_pos, y = annot_ES/2, label = annot_text, colour = "black",
                   size = 4, fontface = 2) + 
          theme(plot.title = element_text(size = 20, hjust = 0.5), 
                plot.subtitle = element_text(size = 12, hjust = 0.5))
        
        print(plot_grob)
        
        
        
        
      }
      
      
    })
    
    
    
    ind <- reactive(input$individual)
    ind_2 <- reactive(input$individual_2)
    ges <- reactive(input$gsea_gene_sets)
    
    leading_genes <- reactive({

      if(ind() == "Specific Pathway"  && ges() == 'MSigDB' ){

        dt <- res() %>%
          filter(pathway == input$gsea_pathway) %>%
          select(c("pathway", "leadingEdge"))


        dt <- unnest(dt, cols = c(leadingEdge))

        dt

      }else if(ind_2() == "Specific Pathway" && ges() == "Custom Gene Set" ){

        dt <- res() %>%
          filter(pathway == input$cre_path) %>%
          select(c("pathway", "leadingEdge"))


        dt <- unnest(dt, cols = c(leadingEdge))

        dt

      }else if(ind() == "Top Pathways"  && ges() == 'MSigDB') {

        high_enrich <- res() %>%
            arrange(desc(NES)) %>%
            head(n= 5)

        low_enrich <- res() %>%
          arrange(desc(NES)) %>%
          tail(n= 5)

        dt <- rbind(high_enrich, low_enrich)

        dt <- select(dt, c("pathway", "leadingEdge"))


        dt <- unnest(dt, cols = c(leadingEdge))

        dt

      }else if(ind_2() == "Top Pathways" && ges() == "Custom Gene Set") {



        dt <- res() %>%
          select(c("pathway", "leadingEdge"))


        dt <- unnest(dt, cols = c(leadingEdge))

        dt

      }



    })

   

    
    # GSEA analysis plot and gene list
    
    observeEvent(input$gsea_run, {

     #gene list

      output$gene_text <- renderDT({

        req(input$gsea_run)
        input$gsea_run

        isolate({

          req(iv$is_valid())
          req(input$gsea_samptyp)
          req(input$gsea_feat)
          req(input$nperm)
          req(input$gsea_width)
          req(input$gsea_height)

          if(input$gsea_feat %in% colnames(con_dat())){
              req(input$gsea_sel_feat_meta_groups)
              req(input$gsea_sel_feat_meta_groups2)
          }

          if(input$gsea_feat %in% colnames(con_dat2())){
              req(input$gsea_high_cutoff)
              req(input$gsea_low_cutoff)
          }

          if(input$gsea_gene_sets == "Custom Gene Set"){
            
            req(input$gset_up)
            
            if (input$individual_2 == "Specific Pathway") {
          
              req(input$cre_path)
            }
          }
            
            if(input$gsea_gene_sets == 'MSigDB'){
      
            req(input$gsea_cat)
              
            if(input$gsea_cat %in% c("C2","C3","C4","C5","C7")) {
                
                req(input$gsea_subcat)
                
              }
            
            if (input$individual == "Specific Pathway") {

              req(input$gsea_pathway)
              
            }
            
          }
        

          leading_genes()
          



        })

      })

      
      
      
      #download 
      
      
      
      output$downloadData <- downloadHandler(
        filename = "preranked_genes.xlsx",
        content = function(file) {
          write.xlsx(ranked_data(), file, colnames = TRUE,
                     rownames = F, append = FALSE, showNA = TRUE)
          
          
        } 
      )
      
      
      output$download_fgsea <- downloadHandler(
        filename = "GSEA_result.xlsx",
        content = function(file) {
          write.xlsx(gsea_res(), file, colnames = TRUE,
                     rownames = F, append = FALSE, showNA = TRUE)
          
          
        } 
      )
      
      output$download_lead <- downloadHandler(
        filename = "leading_edge_genes.xlsx",
        content = function(file) {
          write.xlsx(leading_genes(), file, colnames = TRUE,
                     rownames = F, append = FALSE, showNA = TRUE)


        }
      )
      
      #plot 
      output$gsea_plot <- renderPlot({
        
        req(input$gsea_run)
        input$gsea_run
        
        isolate({
          
          req(iv$is_valid()) 
          
          validate(
            
            need(input$gsea_samptyp, 'Please select at least one sample type'),
            need(input$gsea_feat, 'Feature choice is required'),
            need(input$nperm, 'nPerm Value is required'), 
            need(input$gsea_width, 'Please pick a number although you will not download the plot'),
            need(input$gsea_height, 'Please pick a number although you will not download the plot')
          )
          
          
          
          if(input$gsea_feat %in% colnames(con_dat())){
            
            validate(
              
              need(input$gsea_sel_feat_meta_groups, 'Samples are required'),
              need(input$gsea_sel_feat_meta_groups2, 'Sample ID must be assigned')
              
            )
            
          }
          
          if(input$gsea_feat %in% colnames(con_dat2())){
            
            validate(
              need(input$gsea_high_cutoff, "High cutoff percent is required"),
              need(input$gsea_low_cutoff, "Low cutoff percent is required")
            )
          }
          
          
          if(input$gsea_gene_sets == "Custom Gene Set"){
            
            validate(need(input$gset_up, "Don't forget to upload your xlsx/xls file"))
  
           
             if (input$individual_2 == "Specific Pathway") {
              validate(
                need(input$cre_path, "Pathway choice is required")
                )
            }
          }
          
          if(input$gsea_gene_sets == 'MSigDB'){
            
            validate(need(input$gsea_cat, 'Please select gene set'))
            
            if(input$gsea_cat %in% c("C2","C3","C4","C5","C7")) {
              
              validate(
                need(input$gsea_subcat, "Subcategory choice is needed"))
              
            }
            
            if (input$individual == "Specific Pathway") {
              validate(
                need(input$gsea_pathway, "Pathway choice is needed"))
              
            }
            
          }
          
          gsea_plot()
          
          
        })
      })
      
      
      
      
      output$g_downloadPlot <- downloadHandler(
        filename = function() {
          paste("GSEA_plot.png")
        },
        content = function(file) {
          
          png(file, width = input$gsea_width, height = input$gsea_height)
          print(gsea_plot())
          dev.off()
          
        })
      
      
    })
  })
}