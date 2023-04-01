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


#Load Sources

source("functions/gene_ranker.R")
# source("functions/df_extractor.R")
# source("functions/gene_grapher.R")
source("functions/gsea_replot.R")
source("functions/plotGseaTable2.R")
source("functions/top_plotter.R")
# source("functions/data_importer.R")
# source("functions/data_transformer.R")

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
          
          numericInput(inputId = ns("gsea_high_cutoff"), "High cutoff percent", 50, 1, 100),
          
          numericInput(inputId = ns("gsea_low_cutoff"), "Low cutoff percent", 50, 1, 100)
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
        
        radioButtons(ns("gsea_gene_sets"), "Choose gene set collection", choices = c("MSigDB", "Custom Gene Set")),
      
        conditionalPanel(
          
          condition = "input.gsea_gene_sets == 'MSigDB'",
          ns=ns,
          
          selectizeInput(ns("gsea_cat"), "Please select one of the following MSigDB collections", 
                         choices = NULL, 
                         options=list(placeholder = "eg. Hallmark")
          ),
          
          radioButtons(ns("individual"), "Show", choices = c("Top Pathways", "Specific Pathway")),
          
          conditionalPanel(
            
            condition = 'input.individual == "Specific Pathway"',
            ns=ns,
            
            
            selectizeInput(ns("path"), 
                           "*Please select a pathway", 
                           options=list(placeholder = "eg. HALLMARK_DNA_REPAIR"),
                           choices = NULL)
            
          )
          
          
          
        ),
        
        conditionalPanel(
         
          condition = "input.gsea_gene_sets == 'Custom Gene Set'",
          ns=ns,
          
          fileInput(inputId = ns("gset_up"),
                    label = "Upload your gene set as .csv file",
                    accept = c("text/csv", "text/comma-separated-values,text/plain",
                               ".csv")),
          
          #'[#########################################################################################################################]
          #'[Please add hover info for file formatting. Sample code is below]
          #'[fileInput(..., label = tags$span(
          #'["3. Please upload your csv file.", 
          #'[tags$i(
          #'[   class = "glyphicon glyphicon-info-sign", 
          #'[   style = "color:#0072B2;",
          #'[   title = "The csv file should contain two unnamed columns: the first column should contain the gene set name, and the second column should contain gene names. Each gene should be associated with a gene set (ie. no missing data), and multiple gene sets can be provided in one file."
          #'[))
          #'[#########################################################################################################################]    
          
          
          
          radioButtons(ns("individual_2"), "Show", choices = c("Top Pathways", "Specific Pathway")),
          
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
                   label = "Perform GSEA",
                   style = "unite",
                   block = TRUE,
                   color = "primary"),
        br(),
        downloadButton(ns("g_downloadPlot"), "Download GSEA Plot", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
        br(),
        br(),
        downloadButton(ns("downloadData"), "Download ranked data", style="color: #eeeeee; background-color: #01303f; border-color: #01303f"),
        br(),
        br(),
        #help section UI
        
        introjsUI(),
        actionButton(ns("intro3"), "App Tutorial", style="color: #FFFFFF; background-color: #81A1C1; border-color: #02a9f7"),
        
        
      ),
      
      mainPanel(
        
        plotOutput(outputId = ns("gsea_plot")),
        
        DTOutput(ns("gene_text"))
        
        # tableOutput(ns("test"))
        
    
      )
    )
  )
} 

gsea_server <- function(id,Xproj) {
  moduleServer(id,function(input, output, session){
    
    ns <- session$ns
    
    
    
    ## help section server
    
    help_sec <- reactive({
      
      
      
      if(input$gsea_gene_sets == 'Custom Gene Set') {
        
        data.frame(
          
          element = paste0("#", session$ns(c(NA, "gsea_samptyp + .selectize-control", "gsea_feat + .selectize-control ", "gsea_gene_sets", "individual_2"))),
          
          intro = paste(c(
            "This is the gene set enrichment analysis (GSEA) module. Here, you can categorize samples based on a custom criteria and examine whether previously defined or user-provided gene sets are enriched in either of the data subset. <b>NOTE:</b> This module can take some time depending on the data set. Continue with the tutorial to learn more about this module.",
            "You can select which sample types should be included in the analysis (eg. primary and/or metastatic).",
            "GSEA is performed between two groups of data. If you would like to perform GSEA for a categorical clinical feature, you are expected to select two data subsets and define one of them as the 'sample' for the analysis (the other one will become reference). If you would like to perform GSEA for a numerical feature such as gene expression, then you can categorize samples based on gene expression values as 'high' and 'low' through user-defined quantiles. Setting high and low cutoff to 50, will categorize gene expression at the median value. You can set this numbers to 25 to compare the top 25% expressors to bottom 25% expressors.", 
            "You can perform GSEA using previously defined gene sets from the <a href='https://www.gsea-msigdb.org/gsea/msigdb/'>Molecular Signatures Database (MSigDB)</a> or provide your own gene sets",
            "If you select the 'Top Pathways' option, pathways with the highest and lowest enrichments are shown </i>(they may not always be significant!)</i>. If you choose 'Specific Pathway', the enrichment plot is prepared only for the pathway you select."
          ))
          
        )
        
      }else if (input$gsea_gene_sets == 'MSigDB'){
        
        data.frame(
          
          element = paste0("#", session$ns(c(NA, "gsea_samptyp + .selectize-control", "gsea_feat + .selectize-control ", "gsea_gene_sets" ,"gsea_cat+ .selectize-control ", "individual"))),
          
          intro = paste(c(
            "This is the gene set enrichment analysis (GSEA) module. Here, you can categorize samples based on a custom criteria and examine whether previously defined or user-provided gene sets are enriched in either of the data subset. <b>NOTE:</b> This module can take some time depending on the data set. Continue with the tutorial to learn more about this module.",
            "You can select which sample types should be included in the analysis (eg. primary and/or metastatic).",
            "GSEA is performed between two groups of data. If you would like to perform GSEA for a categorical clinical feature, you are expected to select two data subsets and define one of them as the 'sample' for the analysis (the other one will become reference). If you would like to perform GSEA for a numerical feature such as gene expression, then you can categorize samples based on gene expression values as 'high' and 'low' through user-defined quantiles. Setting high and low cutoff to 50, will categorize gene expression at the median value. You can set this numbers to 25 to compare the top 25% expressors to bottom 25% expressors.",
            "You can perform GSEA using previously defined gene sets from the <a href='https://www.gsea-msigdb.org/gsea/msigdb/'>Molecular Signatures Database (MSigDB)</a> or provide your own gene sets",
            "Select MSigDB gene set collection to use in the analyses <b>NOTE:</b> Since some of these collections contain a large number of gene sets, analysis can take some time.",
            "If you select the 'Top Pathways' option, pathways with the highest and lowest enrichments are shown </i>(they may not always be significant!)</i>. If you choose 'Specific Pathway', the enrichment plot is prepared only for the pathway you select."
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
    
    iv <- InputValidator$new()
    iv$add_validator(panel_iv)
    iv$add_validator(numeric_iv)
    iv$enable()
    
##choice reactive for reference ID   
    
  reference <- reactive({
    
    input$gsea_sel_feat_meta_groups[input$gsea_sel_feat_meta_groups != input$gsea_sel_feat_meta_groups2]
    
  })  
         
## function for pathway choice input
    
  pathway_choice <- reactive({
  
  if(input$gsea_cat == "hallmark") {
    
   load(file='genesets/msigdb_hallmark.rda')
    
    names(msigdb_hallmark)
    
  } else if(input$gsea_cat == "go") {
    
     load(file='genesets/msigdb_go.rda')
    
    names(msigdb_go)
    
  }else if(input$gsea_cat == "curated") {
    
     load(file='genesets/msigdb_curated.rda')
    
    names(msigdb_curated)
    
  }else if(input$gsea_cat == "immune") {
    
   load(file='genesets/msigdb_immune.rda')
    
    names(msigdb_immune)
    
  }else if(input$gsea_cat == "motif") {
    
    load(file='genesets/msigdb_motif.rda')
    
    names(msigdb_motif)
    
  }else if(input$gsea_cat == "all") {
    
    load(file='genesets/msigdb_all.rda')
    
    names(msigdb_all)
    
  }
    
})
    
  
  ## uploaded data 
  
  gdata <- reactive({
    
    req(input$gset_up)
    
    gene_s <- read.table(input$gset_up$datapath, header = FALSE, sep = ";", quote = "\"'",
                         dec = ".")  
    
    colnames(gene_s) <- c("geneset_name", "genes")
    
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
    
    
    observe({updateSelectizeInput(session, 
                                  "path",selected="",
                                  choices = pathway_choice() , 
                                  server = T)})
    
    observe({updateSelectizeInput(session, 
                                  "cre_path",selected="",
                                  choices = names(gdata()) , 
                                  server = T)})
    
    observe({updateSelectizeInput(session, 
                                  "gsea_cat", selected="",
                                  choices = list(`Gene Sets` = list("Hallmark" = "hallmark",
                                                                    "GO" = "go",
                                                                    "Curated" = "curated", 
                                                                    "Immune" = "immune",
                                                                    "Motif" = "motif",
                                                                    "All Gene Sets" = "all")
                                  ), 
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
      
      input$path
      
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
    
    df <- df %>% mutate_if(is.numeric, signif, digits=3)
    
  })
  
  
  
  
  # gene set reactive
  
  gene_s <- reactive({
    
    if(input$gsea_gene_sets == 'Custom Gene Set') {
      
      "Custom Gene Set"
      
    }else if (input$gsea_gene_sets == 'MSigDB'){
      
      input$gsea_cat
      
    }
    
  })
  
  
  #gene set 
  
  gene_set <- reactive({
    
    if (input$gsea_cat == "hallmark" && input$gsea_gene_sets == 'MSigDB') {
      
      load(file='genesets/msigdb_hallmark.rda')
      
      gene_set <-msigdb_hallmark
      
    } else if (input$gsea_cat == "go" && input$gsea_gene_sets == 'MSigDB') {
      
      load(file='genesets/msigdb_go.rda')
      
      gene_set <- msigdb_go
      
    } else if (input$gsea_cat == "curated" && input$gsea_gene_sets == 'MSigDB') {
      
      load(file='genesets/msigdb_curated.rda')
      
      gene_set <- msigdb_curated
      
    } else if (input$gsea_cat == "immune" && input$gsea_gene_sets == 'MSigDB') {
      
      load(file='genesets/msigdb_immune.rda')
      
      gene_set <- msigdb_immune
      
    } else if (input$gsea_cat == "motif" && input$gsea_gene_sets == 'MSigDB') {
      
      load(file='genesets/msigdb_motif.rda')
      
      gene_set <- msigdb_motif
      
    } else if (input$gsea_cat == "all" && input$gsea_gene_sets == 'MSigDB') {
      
      load(file='genesets/msigdb_all.rda')
      
      gene_set <- msigdb_all
      
    }else if (input$gsea_gene_sets == "Custom Gene Set"){
      
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
    
    exprs = gsea_dat()
    
    minSize = 1 
    
    maxSize = ncol(exprs)-3
    
    fgsea(pathways = gene_set(), stats = preranked_genes(), minSize = minSize, maxSize = maxSize)
    
    # assign("gsea_res", res, .GlobalEnv)
    
  }) 
  
  
  # plot 
  
  gsea_plot <- reactive({
    
    exprs = gsea_dat()
    
    ranked_genes = preranked_genes()
    
    # pos_marker = NULL 
    
    # neg_marker = NULL
    
    sample_id = sample_id()
    
    # sample_cluster = NULL
    
    reference_id = reference_id()
    
    # reference_cluster = NULL
    
    gene_set = gene_set()
    
    minSize = 1 
    
    maxSize = ncol(exprs)-3
    
    top_n = 10
    
    gseaParam = 1
    
    plot_individual = plot_ind()
    
    append_title = T
    
    top_plots_title = T
    
    verbose = T 
    
    annot_text_color = "black"
      
    annot_text_size = 4
    
    annot_text_fontface = 2
    
    res = res()
    
    set.seed(123)
    
    if (is.null(plot_individual)) {
      
      if (top_plots_title == T) {
        
        # arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = paste(pos_marker, collapse = "."), neg = paste(neg_marker,
        #                                                                                                                               collapse = "."))
        # select_non_null <- !sapply(arg_list, function(x) {
        #   identical(x, "")
        # })
        
        
        
        main_title <- paste(sample_id, "vs", reference_id)
        # plot_subtitle <- paste(names(arg_list[select_non_null]), arg_list[select_non_null], sep = ": ", collapse = "__")
        
        # plot_title <- paste0(main_title, "\n", plot_subtitle)
        
        plot_title <- main_title
        
      } else {
        plot_title = ""
      }
      
      plot_grob <- top_plotter(gsea_results = res, ranked_genes = ranked_genes, gene_set = gene_set, top_n = top_n, gseaParam = gseaParam,
                               plot_title = plot_title, do.plot = F)
      grid.arrange(plot_grob)
      
      
    } else {
      
      
      hits <- c(grep(plot_individual, res$pathway, ignore.case = T, value = T))
       
      assign('hits', hits, .GlobalEnv)
      # 
      # 
      # if (length(hits) > 1) {
      #   multiple_hits <- t(t(hits))
      #   colnames(multiple_hits) <- "Multiple pathway matches"
      #   rownames(multiple_hits) <- c(1:length(hits))
      #   print(multiple_hits)
      #   num <- as.numeric(readline(prompt = "Multiple pathways are found. Select a number from the list above "))
      #   
      #   
      #   while (!num %in% 1:length(hits)) {
      #     num <- as.numeric(readline(prompt = paste0("Please pick a number between 1 and ", length(hits), ":    ")))
      #   }
      #   assign("num", num, .GlobalEnv)
      #   
      #   annot_padj <- signif(as.numeric(res[res$pathway == hits[num], "padj"]), digits = 2)
      #   annot_NES <- signif(as.numeric(res[res$pathway == hits[num], "NES"]), digits = 2)
      #   annot_ES <- signif(as.numeric(res[res$pathway == hits[num], "ES"]), digits = 2)
      #   
      #   # grob<- grobTree(textGrob(paste('adj.p: ', annot_padj, '\nNES:', annot_NES), x= 0.1, y=annot_ES, hjust = 0, gp = gpar(col='red',
      #   # fontsize=3, fontface='italic')))
      #   
      #   annot_text <- paste("adj.p: ", annot_padj, "\nNES:", annot_NES)
      #   
      #   if (append_title == F) {
      #     
      #     # plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) + labs(title = hits[num]) + annotation_custom(grob)
      #     
      #     plot_grob <- plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) + labs(title = hits[num]) + annotate("text",
      #                                                                                                                             x = x_pos, y = annot_ES/2, label = annot_text, colour = annot_text_color, size = annot_text_size, fontface = annot_text_fontface) +
      #       theme(plot.title = element_text(size = 5, hjust = 0.5))
      #     print(plot_grob)
      #     
      #   } else {
      #     
      #     # arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = pos_marker, neg = neg_marker)
      #     # select_non_null <- !sapply(arg_list, function(x) {
      #     #   identical(x, "")
      #     # })
      #     # select_non_null2 <- !sapply(arg_list, is.null)
      #     # select_non_null <- as.logical(select_non_null * select_non_null2)
      #     # 
      #     # plot_subtitle <- paste(names(arg_list[select_non_null]), arg_list[select_non_null], sep = ": ", collapse = "__")
      #     
      #     
      #     
      #     # plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) + labs(title = paste0(hits[num], sample_id, ' vs ', reference_id),
      #     # subtitle = plot_subtitle)+ annotation_custom(grob)
      #     
      #     plot_grob <- plotEnrichment(pathway = gene_set[[hits[num]]], stats = ranked_genes) + 
      #       labs(title = hits[num], subtitle = paste(sample_id, "vs", reference_id 
      #                                                # plot_subtitle
      #                                                )) + 
      #       annotate("text", x = x_pos, y = annot_ES/2, label = annot_text, colour = annot_text_color,
      #                                                                                                                                                                                  size = annot_text_size, fontface = annot_text_fontface) + theme(plot.title = element_text(size = 10, hjust = 0.5), plot.subtitle = element_text(size = 6,
      #                                                                                                                                                                                                                                                                                                                                  hjust = 0.5))
      #     print(plot_grob)
      #     
      #    
      #     
      #   }
      #   
      #   
      #   
      #   ##########################
      #   
      # } else {
        
        if (verbose)
          message(paste("Plotting", hits))
        
        annot_padj <- signif(as.numeric(res[res$pathway == hits, "padj"]), digits = 2)
        annot_NES <- signif(as.numeric(res[res$pathway == hits, "NES"]), digits = 2)
        annot_ES <- signif(as.numeric(res[res$pathway == hits, "ES"]), digits = 2)
        x_pos <- length(ranked_genes)/4
        
        # grob<- grobTree(textGrob(paste('adj.p: ', annot_padj, '\nNES:', annot_NES), x= 0.1, y=annot_ES, hjust = 0, gp = gpar(col='red',
        # fontsize=3, fontface='italic')))
        
        annot_text <- paste("adj.p: ", annot_padj, "\nNES:", annot_NES)
        
        
        if (append_title == F) {
          
          # plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) + labs(title = hits) + annotation_custom(grob)
          
          plot_grob <- plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) + labs(title = hits) + annotate("text", x = x_pos,
                                                                                                                        y = annot_ES/2, label = annot_text, colour = annot_text_color, size = annot_text_size, fontface = annot_text_fontface)
          print(plot_grob)
          
        } else {
          
          # arg_list <- list(samp_clu = sample_cluster, ref_clu = reference_cluster, pos = paste(pos_marker, collapse = "."), neg = paste(neg_marker,
          #                                                                                                                               collapse = "."))
          # 
          # select_non_null <- !sapply(arg_list, function(x) {
          #   identical(x, "")
          # })
          # select_non_null2 <- !sapply(arg_list, is.null)
          # select_non_null <- as.logical(select_non_null * select_non_null2)
          # 
          # plot_subtitle <- paste(names(arg_list[select_non_null]), arg_list[select_non_null], sep = ": ", collapse = " ")
          
          # plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) + labs(title = paste0(hits, sample_id, ' vs ', reference_id), subtitle =
          # plot_subtitle)+ annotation_custom(grob)
          
          plot_grob <- plotEnrichment(pathway = gene_set[[hits]], stats = ranked_genes) + 
            labs(title = hits, subtitle = paste(sample_id,"vs", reference_id 
                                                # plot_subtitle
                                                )) + 
            annotate("text", x = x_pos, y = annot_ES/2, label = annot_text, colour = annot_text_color,
                     size = annot_text_size, fontface = annot_text_fontface) + 
            theme(plot.title = element_text(size = 10, hjust = 0.5), 
                  plot.subtitle = element_text(size = 6, hjust = 0.5))
          
          print(plot_grob)
        }
        
      # }
      
    }
    
    
  })
  
  
  

    # GSEA analysis plot and gene list
    
    observeEvent(input$gsea_run, {
      
      #gene list
      
      output$gene_text <- renderDataTable({
        
        req(input$gsea_run)
        input$gsea_run
        
        isolate({
          
          req(iv$is_valid()) 
          
          validate(

            need(input$gsea_samptyp, ''),
            need(input$gsea_feat, '')
          )

          
          if(input$gsea_feat %in% colnames(con_dat())){
            
            validate(
              
              need(input$gsea_sel_feat_meta_groups, ''),
              need(input$gsea_sel_feat_meta_groups2, '')
              
            )
            
          }
          
          if(input$gsea_feat %in% colnames(con_dat2())){
            
            validate(
              need(input$gsea_high_cutoff, ""),
              need(input$gsea_low_cutoff, '')
            )
          }
          
          
          if(input$gsea_gene_sets == 'MSigDB'){
            
            validate(
              need(input$gsea_cat, '')
            )
            if (input$individual == "Specific Pathway") {
              validate(
                need(input$path, ""))
            }
            
          }
          
          
          if(input$gsea_gene_sets == "Custom Gene Set"){
            
            validate(need(input$gset_up, ""))
            
            if (input$individual_2 == "Specific Pathway") {
              validate(
                need(input$cre_path, ""))
            }
            
          }
     

          ranked_data()
          
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
      
      
      #plot 
        output$gsea_plot <- renderPlot({
          
          req(input$gsea_run)
          input$gsea_run
          
          isolate({
            
            req(iv$is_valid()) 
            
            validate(

              need(input$gsea_samptyp, 'Please select at least one sample type'),
              need(input$gsea_feat, 'Feature choice is required')
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
              
              validate(need(input$gset_up, "Don't forget to upload your csv file"))
              
              if (input$individual_2 == "Specific Pathway") {
                validate(
                  need(input$cre_path, "Pathway choice is required"))
              }
            }
            
              if(input$gsea_gene_sets == 'MSigDB'){
              
             validate(need(input$gsea_cat, 'Please select gene set'))
              
              if (input$individual == "Specific Pathway") {
                validate(
                  need(input$path, "Pathway choice is needed"))
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
              
              png(file)
              print(gsea_plot())
              dev.off()
              
            })
          
      
          
          
          
        
        

 
        
        
      })
  })
}
