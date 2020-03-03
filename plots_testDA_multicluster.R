library(tidyverse)
library(yardstick)
library(iCOBRA)
library(cowplot)
library(ggpubr)
# dir_save <- "/home/reto/polybox/ETH/Master_Thesis/Portmacquarie/tmp/"
dir_save <- "/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/"
plot_dir <- "/home/reto/polybox/ETH/Master_Thesis/plots/Simulations/testDA/"
test_results <- readRDS(paste0(dir_save,"test_results_testDA_multicluster.rds"))


# grouping_quo <- dplyr::quos(cluster_id,simType, formula, impType, n_dat, censoring_rate,
#                             b0_True, b1_True,b2_True,mi_rep, cov_dep_cens,
#                             error_variance, variance_fixeff, variance_raneff)
params <- purrr::map(test_results[[2]], ~ .x[c("mrl","rs","km","pmm","cc","ppd")]) 
plts <- purrr::map(seq(3,36,by=3),function(x){
  params_one <- params[[x]]
  params_one <- purrr::map(params_one, function(x){
    # print(class(x))
    if(class(x)=="tbl_df"){
      return(x)
    }else{
      return(NULL)
    }
  })
  params_one_comb <- purrr::map(params_one, ~bind_rows(.x)) %>% bind_rows()
    # purrr::reduce(rbind)
  # params <- bind_rows(bind_rows(test_results[[2]][[1]][c("mrl","rs","km","pmm")]),bind_rows(test_results[[2]][[2]][c("mrl","rs","km","pmm")]))
  
  tib2plt <- params_one_comb %>% 
    mutate(proportion_censored = n_cens/n_dat,
           truestru = as.numeric(cluster_id) %% 20 == 1,
           class = as.factor(truestru),
           impType = factor(impType, levels = c("cc","pmm","km", "rs","mrl","ppd")),
           simType = factor(simType, levels = c("lm","glm","glmer")),
           b1_True = as.factor(b1_True), 
           n_dat= as.factor(n_dat),
           cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR")))
  
  
  cobratib <- tib2plt %>%
    mutate(id = seq_len(dim(tib2plt)[1])) #%>%
    # filter(censoring_rate==0.3, cov_dep_cens=="MAR") #%>%
  #   select(!!!grouping_quo,n_cens,p_adj,impType)
  # cobratib_melt <- spread(cobratib,impType,p_adj)
  
  pval <- data.frame(cc = cobratib$p_val[cobratib$impType=="cc"])
  padj <- data.frame(cc = cobratib$p_adj[cobratib$impType=="cc"])
  for (impType in unique(cobratib$impType)) {
    pval[[impType]] <- cobratib$p_val[cobratib$impType==impType]
    padj[[impType]] <- cobratib$p_adj[cobratib$impType==impType]
  }
  rownames(pval) <- paste0("I_",cobratib$id[cobratib$impType=="cc"])
  rownames(padj) <- paste0("I_",cobratib$id[cobratib$impType=="cc"])

  truth <- data.frame(status = as.integer(cobratib$class[cobratib$impType=="cc"])-1,
                      cov_dep_cens = cobratib$cov_dep_cens[cobratib$impType=="cc"],
                      sim_id = cobratib$sim_id[cobratib$impType=="cc"],
                      row.names = paste0("I_",cobratib$id[cobratib$impType=="cc"])
                                    )
  cobradata <- COBRAData(pval = pval, truth = truth,padj = padj)
  
  
  # data(cobradata_example)
  cobraperf <- calculate_performance(cobradata,
                                     binary_truth = "status",
                                     aspects = c("fdrtpr", "fdrtprcurve",
                                                 "tpr", "roc"),
                                     thrs = c(0.01, 0.05, 0.1))#, splv = "sim_id")
  
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", 
                                     facetted = TRUE)
  # plot_tpr(cobraplot)
  plt1 <- plot_fdrtprcurve(cobraplot,pointsize = 1,linewidth = 0.5) + 
    labs(title=paste0("n",unique(params_one_comb$n_dat),
                      "_c",100*unique(params_one_comb$censoring_rate),
                      "_",ifelse(unique(params_one_comb$cov_dep_cens)==1,"MAR","MCAR"),
                      "_",unique(params_one_comb$transform_fn)))
  return(plt1)
})
plt_comb <- ggarrange(plotlist = plts)#,nrow=1,ncol=2)
plt_comb
ggsave(paste0(plot_dir,"tpr_fdr.png"), plt_comb, width = 100, height = 100, units = "cm")

