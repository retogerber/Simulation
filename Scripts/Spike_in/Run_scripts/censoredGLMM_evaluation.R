# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Spike_in/Run_scripts")


################################################################################
## load packages
devtools::load_all("/home/reto/polybox/ETH/Master_Thesis/Code/censcyt")
suppressPackageStartupMessages(library(ggplot2))
library(tibble)
library(magrittr)
suppressPackageStartupMessages(library(diffcyt))
suppressPackageStartupMessages(library(CATALYST))
suppressPackageStartupMessages(library(SummarizedExperiment))


main_dir <- "../../.."
data_type <- c("Preprocessed","PositiveSignal")
iii <- 1
filename_pattern <- c("_pp","_pp_ps")
survival_time_scale <- c("original", "cluster_dependent")
# survival_time_scale <- c("original","relative","log")
adapt_counts <- FALSE
plot_dir <- "/home/reto/polybox/ETH/Master_Thesis/plots/Framework/PositiveSignal/diffcyt_Framework"

p_melt_ls <- purrr::map(c(FALSE), function(adapt_counts){
  # res_list <-
  #   purrr::map2(rep(seq_along(data_type), each = length(survival_time_scale)),
  #               rep(seq_along(survival_time_scale), length(data_type)), 
  #               function(iii, ii) {
  #     readRDS(
  #       paste0(
  #         main_dir,
  #         "/FlowCap/results/",
  #         data_type[iii] ,
  #         "/censdiffcyt_Framework/res_list_scale_",
  #         survival_time_scale[ii],
  #         "_",
  #         filename_pattern[iii],
  #         "_testDA_censoredGLM",
  #         ifelse(adapt_counts,"_adapted_survivaltime",""),
  #         ".rds"
  #       )
  #     )
  # })
  res_list <- readRDS(paste0(main_dir,"/Data/Real/RDS_files/", data_type[iii], "/res_list_complete_testDA_censoredGLMM.rds"))
  
  if (adapt_counts){
    true_diff_cluster_nr_ls <- c(9,9,3,3)
  } else {
    true_diff_cluster_nr_ls <- list(c(NA),c(1,4),c(NA),c(1,4))
  }
  p_melt_ls <- purrr::map(seq_along(res_list), function(x){
    conditions_res <- res_list[[x]]
    # print(true_diff_cluster_nr_ls[[x]])
    # print(x)
    # true_diff_cluster_nr <- true_diff_cluster_nr_ls[[x]]
    if (class(conditions_res)=="character"){
      return(NULL)
    }
    p_names <-  unlist(purrr::map(seq_along(conditions_res), function(x){
      conditions_res$conditions[3]
    }))
    # p_ls <- purrr::map(seq_along(conditions_res), function(x){
    p_ls <- matrix(tryCatch(rowData(conditions_res$result)$p_adj,error=function(e){NA}),ncol=1)
    # colnames(p_ls) <- conditions_res$conditions[3]
    tibble(pvals = p_ls[,1], 
           scale = conditions_res$conditions[1], 
           clustering = conditions_res$conditions[2], 
           imputation = conditions_res$conditions[3],
           method_est = conditions_res$conditions[4])
    # list(pvals = p_ls,conditions = tibble(conditions_res$conditions))
      # tmp
      # positive_signal_result[[tmpind]]$conditions[3]
    })
    # tmp_cond <- conditions_res[[1]]$conditions[c(1,2)]
    # p_comb <-  dplyr::bind_cols(p_ls[!is.na(p_ls)]) %>% 
    #   dplyr::mutate(diff = FALSE, 
    #                 adapt_counts = factor(adapt_counts,levels = c(FALSE,TRUE),labels = c("spike-in before clustering","spike-in after clustering")),
    #                 data_type = factor(data_type[as.integer(tmp_cond[2])], levels = c("Preprocessed","PositiveSignal"),labels = c("Preprocessed","PositiveSignal")),
    #                 survival_type = tmp_cond[1]
    #                 )#,
    #                 # condition = paste("adapt_counts:",adapt_counts," -- ",data_type[as.integer(tmp_cond[2])]))
    # 
    # if (!is.na(true_diff_cluster_nr)) {
    #   p_comb[true_diff_cluster_nr, "diff"] <- TRUE
    # }
    # colnames(p_comb) <- c(p_names[!is.na(p_ls)],"diff","adapt_counts","data_type","survival_type")
    # # print(p_comb)
    # p_melt <- reshape2::melt(p_comb, id.vars = c("diff","adapt_counts","data_type","survival_type"))
  #   return(p_melt)
  # })
  return(dplyr::bind_rows(p_melt_ls))
  
})
p_melt_comb <- dplyr::bind_rows(p_melt_ls) #%>% 
  # dplyr::filter((survival_type != "cluster_dependent") | (adapt_counts != "spike-in after clustering"))

plt1 <- ggplot(p_melt_comb,aes(method_est,pvals)) + 
  geom_jitter(width = 0.2) + 
  labs(x = "", y = "p-value", color = "Real difference") +
  ylim(c(0, 1)) +
  facet_wrap(clustering ~ imputation)
plt1
ggsave(paste0(plot_dir,"spike_in_test_results.png") ,plt1, width = 20, height = 20, units = "cm")






























# mice.impute.tester <- function(y, ry, x = NULL, wy = NULL, ...){
#   print("y")
#   print(str(y))
#   print("ry")
#   print(str(ry))
#   print("y[ry]")
#   print(str(y[ry]))
#   print("x")
#   print(head(x))
#   print("wy")
#   print(str(wy))
#   return(runif(sum(wy)))
# }
# tmptmp <- mice(micedata, maxit = 1, m = 1, seed = 1, method = "tester")
# 
# tmpdata <- data_comb[[2]]$formula[[1]]$data
# library(mice)
# survival_time_na <- tmpdata$Survival.Time
# survival_time_na[tmpdata$Status==0] <- NA
# survival_time_compl <- tmpdata$Survival.Time
# survival_time_compl[tmpdata$Status==1] <- 1e5
# 
# micedata <- data.frame(Survival.Time_na = survival_time_na/max(tmpdata$Survival.Time), 
#                        Survival.Time_compl = survival_time_compl/max(tmpdata$Survival.Time),
#                        Status = as.factor(tmpdata$Status),
#                        sample_id = as.factor(tmpdata$sample_id),
#                        patient_id = as.factor(tmpdata$patient_id),
#                       counts = assay(data_comb[[2]]$d_counts)[3,]/colSums(assay(data_comb[[2]]$d_counts)),
#                       weights = colSums(assay(data_comb[[2]]$d_counts)))
# micedata$counts_rel <- micedata$counts/micedata$weights
# # md.pattern(micedata)
# imp <- mice(micedata, m = 1, seed = 1)
# stripplot(imp, Survival.Time, pch = 19, xlab = "Imputation number")
# fit <- with(imp, lme4::glmer(counts ~ Survival.Time + (1|patient_id) + (1|sample_id), family = "binomial", weights = micedata$weights))
# summary(pool(fit))

# i <- 1
# for (i in 1:10){
#   coxdata <- data.frame(Survival.Time = tmpdata$Survival.Time, 
#                         Status = tmpdata$Status, 
#                         counts = assay(data_comb[[2]]$d_counts)[i,])
#   
#   
#   coxfit <- survival::coxph(survival::Surv(Survival.Time,Status) ~ counts, coxdata)
#   cat("p-value for row ", i, " :\t",summary(coxfit)$logtest[3],"\n")
# }



# preprocessed_signal_result <- res_list[[1]]
# positive_signal_result <- res_list[[4]]
# 
# tmpind <- 10
# print(preprocessed_signal_result[[tmpind]]$conditions)
# topTable(preprocessed_signal_result[[tmpind]]$result)
# 
# tmpind <- 10
# print(positive_signal_result[[tmpind]]$conditions)
# topTable(positive_signal_result[[tmpind]]$result)

