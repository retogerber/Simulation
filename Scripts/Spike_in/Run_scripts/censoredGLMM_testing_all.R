# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Spike_in/Run_scripts")


################################################################################
## load packages
devtools::load_all("/home/reto/polybox/ETH/Master_Thesis/Code/censcyt")



main_dir <- "../../.."
data_type <- c("Preprocessed","PositiveSignal")
filename_pattern <- c("_pp","_pp_ps")
survival_time_scale <- c("original","relative","log")
clusterings <- c("meta10","meta20","som100")
clus_pats <- expand.grid(filename_pattern,clusterings,stringsAsFactors = FALSE)


cat("prepare data ...")
data_comb <- purrr::map2(rep(seq_along(filename_pattern),length(clusterings)),rep(clusterings,each=length(filename_pattern)), function(iii,ii){
  tmplist <- readRDS(paste0(main_dir,"/Data/Real/RDS_files/", data_type[iii] ,"/res", 
         filename_pattern[iii] ,
         "_",
         ii,
         "_counts_formula_contrast.rds"))

  formula_original <- tmplist[[2]]
  formula_original$data$Survival.Time <- tmplist[[2]]$data$Survival.Time*5855  
  formula_relative <- tmplist[[2]]
  formula_relative$data$Survival.Time <- tmplist[[2]]$data$Survival.Time
  formula_log <- tmplist[[2]]
  formula_log$data$Survival.Time <- log(tmplist[[2]]$data$Survival.Time*5855+abs(min(tmplist[[2]]$data$Survival.Time*5855))+1)
  
  return(list(conditions = c(data_type[iii],ii),
              d_counts = tmplist[[1]],
              formula = list(original = formula_original,
                             relative = formula_relative,
                             log = formula_log),
              contrast = tmplist[[3]]))
})
cat("done\nstart imputations ...")
conditions <- expand.grid(
  scale = survival_time_scale[1], 
  x = seq_len(dim(clus_pats)[1]), 
  method_est = c("cc","rs","km","pmm","mrl","ppd"),
  stringsAsFactors= FALSE)
nc_cores <- 72
res_list <- par_pmap(conditions, function(scale,x, method_est){
  tmplist <- data_comb[[x]]
  d_counts <- tmplist[["d_counts"]]
  contrast <- tmplist[["contrast"]]
  formula <- tmplist[["formula"]][[scale]]
  formula$data$Survival.Time <-  formula$data$Survival.Time + abs(min( formula$data$Survival.Time))+1
  rm(tmplist)
  list(result = testDA_censoredGLMM(d_counts = d_counts, formula = formula,
                                            contrast = contrast,
                                            method_est = method_est,
                                            verbose = FALSE, m = 50,
                                            print_diagnostics = FALSE, 
                                            nr_cores = 2),
       conditions = c(scale, data_type = clus_pats[x,1], n_clusters = clus_pats[x,2],method_est))
}, mc.cores = nc_cores)
cat("done\nsave results ...")
saveRDS(
  res_list,
  paste0(main_dir,"/Data/Real/RDS_files/", data_type[iii], "/res_list_complete_testDA_censoredGLMM.rds"))

cat("done\n")
