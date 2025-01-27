
library(dplyr)


path <- "C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/Masters Thesis/R/projects/TCGEx/shiny/projects"

rds_files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)


for (file in rds_files) {
  
  data <- readRDS(file)
  
  if ("meta.gender" %in% colnames(data)) {
    colnames(data)[colnames(data) == "meta.gender"] <- "meta.sex"
  }
  
  saveRDS(data, file)
}

#
colnames(`SKCM-cpm_TCGA`)[colnames(`SKCM-cpm_TCGA`) == "meta.gender"] <- "meta.sex"

saveRDS(`SKCM-cpm_TCGA`, "SKCM-cpm_TCGA.rds")
#





library(dplyr)

directory_path <- "C:/Users/memre/OneDrive - Izmir Yuksek Teknoloji Enstitusu/IYTE/Masters Thesis/R/projects/TCGEx/shiny/projects"

rds_files <- list.files(directory_path, pattern = "\\.rds$", full.names = TRUE)

results <- data.frame(
  File = character(),
  MetaSexExists = logical(),
  MetaGenderExists = logical(),
  stringsAsFactors = FALSE
)

for (file in rds_files) {
  data <- tryCatch(readRDS(file), error = function(e) NULL)
  
  if (is.null(data)) {
    warning(paste("Cannot read file:", file))
    next
  }
  
  meta_sex_exists <- "meta.sex" %in% colnames(data)
  meta_gender_exists <- "meta.gender" %in% colnames(data)
  
  results <- results %>% 
    add_row(File = basename(file),
            MetaSexExists = meta_sex_exists,
            MetaGenderExists = meta_gender_exists)
}
print(results)
