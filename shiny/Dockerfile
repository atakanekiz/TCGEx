FROM rocker/shiny:4.2.2

RUN apt-get update -qq && \ 
  apt-get install -y --no-install-recommends \
    cmake \
    gsfonts \
    imagemagick \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libicu-dev \
    libjpeg-dev \
    libmagick++-dev \
    libpng-dev \
    libssl-dev \
    libtiff-dev \
    libxml2-dev \
    make \
    pandoc \
    zlib1g-dev \
    && apt-get clean \ 
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('heatmaply','shinyBS','grid','plyr','colourpicker','gridExtra','corrplot','stats','shinyjs','fresh','openxlsx','readxl','shinydashboardPlus','survminer', 'survival','DT','data.table', 'dplyr', 'Hmisc', 'janitor', 'plotly', 'shinybusy', 'rintrojs', 'shinyvalidate', 'ggpubr', 'shiny', 'ggplot2', 'shinydashboard', 'ggiraph', 'extrafont', 'broom', 'msigdbr', 'pheatmap', 'RColorBrewer','shinyWidgets', 'tidyverse', 'glmnet', 'shinythemes', 'factoextra', 'kableExtra', 'plotROC'))"
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install(c('edgeR', 'fgsea','ComplexHeatmap'))"

WORKDIR /app   # project folder should be empty

COPY . /app

#wget ..... googledrv.gz  > projects/

#gzip googledrv.gz

#rm googledrv 


EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/app/app.R', host = '0.0.0.0', port = 3838)"]