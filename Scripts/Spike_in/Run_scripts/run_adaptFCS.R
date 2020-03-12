### run on cluster with number of cores:
nc_wanted <- 50

# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Spike_in/Run_scripts")

# Raw data is only present in compressed form, to rerun please first extract the 
# files


################################################################################
## load packages
library(parallel)
library(reshape2)
library(tidyverse)
library(censcyt)
library(data.table)

data_dir <- "../../../Data/Real"

nc <- max(min(nc_wanted,detectCores()-2),1)


data_dir_raw <- paste(getwd(),data_dir,"Raw",sep="/")
data_dir_pp <- paste(dirname(data_dir_raw),"Preprocessed",sep = "/")
data_dir_ps <- paste(dirname(data_dir_raw),"PositiveSignal",sep = "/")
data_dir_ps_2 <- paste(dirname(data_dir_raw),"PositiveSignalControl",sep = "/")
  
if(!dir.exists(data_dir_raw)) abort("Raw data not found, check directory")
if(!dir.exists(data_dir_pp)) {dir.create(data_dir_pp)}
if(!dir.exists(data_dir_ps)) {dir.create(data_dir_ps)}

## preprocessing
filenames <- list.files(data_dir_raw, pattern = ".fcs", full.names=TRUE)
par_walk(filenames,
        .f = fcs_preprocessing,
        outdir = data_dir_pp,
        compensate = TRUE,
        cofactor = 150,
        filter_extremes = FALSE,
        live_cd3_gate = FALSE,
        mc.cores = nc)

## adapt FCS (introduce true positive signal)
filenames_pp <- list.files(data_dir_pp, pattern = ".fcs", full.names=TRUE)

MetaData <- as_tibble(read.csv(paste(dirname(data_dir_raw),"MetaDataTrain.csv",sep="/"), sep = ",", header = TRUE)) %>% 
  dplyr::filter(!is.na(Status)) %>% 
  select(-c(Prediction)) 

MetaDataMelt <- melt(MetaData,id.vars = c("Status", "Survival.Time"))

sfit <- survival::survfit(survival::Surv(MetaDataMelt$Survival.Time, as.numeric(MetaDataMelt$Status))~1)
sumsfit <- summary(sfit)
dt = data.table(Time = sumsfit$time, id = seq_along(sumsfit$time))
setkey(dt, Time)
id_match <- dt[J(MetaDataMelt$Survival.Time), roll = TRUE, which = TRUE]
MetaDataMelt <- MetaDataMelt  %>%
  mutate(survival_time_rel =  (summary(sfit)$surv[id_match]), 
         filebasename = str_replace(value,"\\.","_pp."), 
         filebasename_ps = str_replace(value,"\\.","_pp_ps."))

marker_names <-  c("IFNg","TNFa","CD4","CD27","CD107a","CD154","CD3","CCR7","IL2","CD8","CD57","CD45RO","Vivid/CD14")

# make random sampling in 'adapt_fcs' reproducible
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
s <- .Random.seed
mc.reset.stream()
# run parallel
par_pwalk(MetaDataMelt[ ,c("survival_time_rel","filebasename")], 
          ~ adapt_fcs(filename = paste(data_dir_pp,..2,sep="/"),
                      phenotype = c("CD3+CD4+"),
                      outdir = data_dir_ps,
                      marker_names = marker_names,
                      channel_marker_ref = filenames_pp[2],
                      type = "frequency", 
                      fn = function(x,time_relative){x*time_relative},
                      time_relative = ..1),
          mc.cores = nc,
          mc.set.seed = TRUE)


# run again to get new thresholds and population sizes
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
s <- .Random.seed
mc.reset.stream()
par_pwalk(MetaDataMelt[ ,c("survival_time_rel","filebasename_ps")], 
          ~ adapt_fcs(filename = paste(data_dir_ps,..2,sep="/"),
                      phenotype = c("CD3+CD4+"),
                      outdir = data_dir_ps_2,
                      marker_names = marker_names,
                      channel_marker_ref = filenames_pp[2],
                      type = "frequency", 
                      fn = function(x,time_relative){x*time_relative},
                      time_relative = ..1,
                      dry_run = TRUE),
          mc.cores = nc,
          mc.set.seed = TRUE)

