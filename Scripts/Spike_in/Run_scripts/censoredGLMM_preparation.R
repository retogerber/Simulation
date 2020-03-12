# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Spike_in/Run_scripts")


################################################################################
## load packages
devtools::load_all("/home/reto/polybox/ETH/Master_Thesis/Code/censcyt")

library(flowCore)
library(dplyr)
library(reshape2)
library(purrr)
library(censcyt)
library(diffcyt)
library(CATALYST)
library(SummarizedExperiment)


main_dir <- "../../.."
data_type <- c("Preprocessed","PositiveSignal")
filename_pattern <- c("_pp","_pp_ps")
for (iii in 1:2) {
  data_dir <- paste0(main_dir,"/Data/Real/", data_type[iii], "/")
  data_dir_meta <- paste0(main_dir,"/Data/Real/")
  MetaDataFile <- paste0(data_dir_meta, "MetaDataTrain.csv")
  ChannelsFile <- paste0(data_dir_meta, "FlowCAPchannels.csv")
  
  MetaData <- as_tibble(read.csv(MetaDataFile, sep = ",", header = TRUE)) %>% 
    dplyr::filter(!is.na(Status)) %>% 
    select(-Prediction) %>% 
    mutate(patient_id = as.factor(seq(1:length(Status))), Status = as.factor(Status))
  
  cat("Preprocessing:\n")
  MetaData_melt <- melt(MetaData, id.vars = c("Status","Survival.Time","patient_id"), value.name = "sample_name", variable.name = "condition")
  # sort according to time
  MetaData_melt <- MetaData_melt %>% 
    # mutate(TrueSurvivalTime = ifelse(Status == 0, Survival.Time*2, Survival.Time)) %>% 
    arrange(Survival.Time)
  experiment_info <- MetaData_melt %>% 
    mutate(sample_id = as.factor(paste0("Sa_", gsub(".fcs","",sample_name)))) %>% 
    select(-sample_name) #%>% 
  
  Channels <- as_tibble(read.csv(ChannelsFile, sep = ",", header = TRUE, stringsAsFactors = FALSE))
  marker_info <- Channels %>% 
    mutate(marker_class = if_else((X == "Not used") | (Channel.Name %in% c("Time","FSC-A","FSC-H","SSC-A")), "none","type")) %>% 
    transmute(channel_name = Channel.Name, marker_name = if_else(Reagent != "", Reagent, channel_name), marker_class) #%>% 
  
  
  # subset of processed files
  files <- experiment_info %>%
    mutate(rowid = row_number()) %>%
    select(sample_id) %>%
    transmute(filename = paste0(data_dir,as.character(gsub("Sa_","",sample_id)), filename_pattern[iii], ".fcs")) %>%
    as_vector()
  
  stopifnot(all(basename(files) %in% list.files(data_dir)))
  
  # subset of experiment_info
  experiment_info <- experiment_info %>% 
    mutate(rowid = row_number()) %>% 
    select(-rowid) %>% 
    mutate(patient_id = droplevels(as.factor(patient_id)),
          sample_id = droplevels(as.factor(sample_id)),
          condition = droplevels(as.factor(condition)))
  
  # experiment_info <- experiment_info %>% select(-patient_id)
  
  design <- createDesignMatrix(
    experiment_info, cols_design = c("Survival.Time","condition")
  )
  # design
  contrast <- createContrast(c(0, 1, 0))
  cat("dimension of contrast matrix == design matrix: ")
  nrow(contrast) == ncol(design)
  
  cols_to_include <- marker_info %>% transmute(marker_class_b = marker_class == "type") %>% as_vector()
  
  d_flowSet <- read.flowSet(
    files, transformation = FALSE, truncate_max_range = FALSE
  )
  
  marker_info <- marker_info %>% dplyr::rename(channel = channel_name, antigen = marker_name, class = marker_class)
  daf <- prepData(x = d_flowSet, 
                  panel = marker_info, 
                  md = experiment_info %>% mutate(file = files), 
                  panel_cols = list(channel = "channel",
                                    antigen = "antigen",
                                    class = "class"),
                  md_cols = list(file = "file",
                                 id = "sample_id",
                                 factors = c("condition","Survival.Time")),
                  cofactor = NULL)
  daf <- CATALYST::cluster(daf, xdim = 10, ydim = 10, maxK = 20, seed = 123)
  
  saveRDS(daf, paste0(main_dir,"/Data/Real/RDS_files/", data_type[iii] ,"/daf", filename_pattern[iii] ,"_clustered.rds"))

    # transform to diffcyt format
  for (clustering_to_use in c("meta10","meta20","som100")) {
  # clustering_to_use <- "meta20"
  code_id <- colData(daf)$cluster_id
  cluster_id <- metadata(daf)$cluster_codes[, clustering_to_use][code_id]
  stopifnot(length(cluster_id) == nrow(colData(daf)), 
            length(code_id) == nrow(colData(daf)))
  colData(daf)$code_id <- code_id
  colData(daf)$cluster_id <- cluster_id
  clustering_name <- clustering_to_use
  cs_by_s <- split(seq_len(ncol(daf)), colData(daf)$sample_id)
  cs <- unlist(cs_by_s[metadata(daf)$experiment_info$sample_id])
  es <- t(assays(daf)[["exprs"]])[cs, , drop = FALSE]
  d_se <- SummarizedExperiment(assays = list(exprs = es), 
                               rowData = colData(daf)[cs, ], colData = rowData(daf), 
                               metadata = metadata(daf))
  
  d_counts <- calcCounts(d_se)
  d_medians <- calcMedians(d_se)
  
  da_formula1 <- createFormula(experiment_info, 
                               cols_fixed = c("Survival.Time","condition"),
                               cols_random = c("patient_id" ,"sample_id"))
  experiment_info_re_time <- experiment_info
  experiment_info_re_time["Survival.Time"] <- experiment_info_re_time["Survival.Time"] *(1/max(experiment_info_re_time["Survival.Time"]))
  da_formula1$data <- experiment_info_re_time
  da_formula1$formula <- formula(y ~ Surv(Survival.Time,Status) + condition + (1|patient_id) + (1|sample_id))
  da_formula1$data$Status <- as.integer(levels(da_formula1$data$Status))[as.integer(da_formula1$data$Status)]
  da_formula1$data$condition <- ifelse(da_formula1$data$condition == "Stim", 
                                       rep(1,length(da_formula1$data$condition)),
                                       rep(0,length(da_formula1$data$condition)))
  levels(da_formula1$data$sample_id) <- gsub("Sa_","",levels(da_formula1$data$sample_id))
  da_formula1$data$sample_id <- as.integer(levels(da_formula1$data$sample_id)[as.integer(da_formula1$data$sample_id)])
  formula <- da_formula1
  
  tmplist <- list(d_counts,formula,contrast)
  saveRDS(tmplist, paste0(main_dir,"/Data/Real/RDS_files/", data_type[iii] ,"/res", 
                          filename_pattern[iii] ,
                          "_",
                          clustering_to_use,
                          "_counts_formula_contrast.rds"))
  
  }
}
