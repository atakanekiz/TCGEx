# Dockerfile

# R tabanlı resmi imajı kullan

FROM rocker/shiny:latest

WORKDIR /app

COPY . /app

# Gerekli R paketlerini yükleyin

RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libv8-dev \
    && R -e "install.packages(c('heatmaply','grid','plyr','colourpicker','gridExtra','corrplot','stats','shinyjs','fresh','openxlsx','readxl','shinydashboardPlus','survminer', 'survival', 'data.table', 'dplyr', 'Hmisc', 'janitor', 'plotly', 'shinybusy', 'rintrojs', 'shinyvalidate', 'ggpubr', 'shiny', 'ggplot2', 'shinydashboard', 'ggiraph', 'extrafont', 'broom', 'msigdbr', 'pheatmap', 'RColorBrewer','shinyWidgets', 'tidyverse', 'glmnet', 'shinythemes', 'factoextra', 'kableExtra', 'plotROC'))"

RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install(c('edgeR', 'fgsea','ComplexHeatmap'))"

# Uygulamanın çalışmasını sağla

EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/app/app.R', host = '0.0.0.0', port = 3838)"]
