# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Simulated/Plot_scripts")


################################################################################
## load packages and data
library(tidyverse)
library(yardstick)
library(iCOBRA)
# library(cowplot)
library(ggpubr)
library(scales)
dir_save <- "../../../Data/Simulated/"
plot_dir <- "../../../Plots/Simulated/"


test_results <- readRDS(paste0(dir_save,"test_results_testDA_multicluster_10reps.rds"))

################################################################################
## processing and plotting

params <- purrr::map(test_results[[2]], ~ .x[c("mrl","rs","km","pmm","cc","ppd")]) 
conditions <- purrr::map(params, function(dfs){
  ind <- which(unlist(purrr::map(dfs,~is.data.frame(.x))))[1]
  if (is.na(ind)) return(tibble::tibble(n_dat=NA,censoring_rate=NA,cov_dep_cens=NA,transform_fn=NA))
  dfs[[ind]] %>% dplyr::select(n_dat,censoring_rate,cov_dep_cens,transform_fn) %>% unique()
}) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(ID=seq_along(n_dat)) %>% 
  dplyr::arrange(n_dat,censoring_rate,transform_fn,cov_dep_cens)

conditions %>% 
  select(-ID) %>% 
  group_by(n_dat,censoring_rate,cov_dep_cens,transform_fn) %>% 
  tally()

all_impTypes <- c("cc","km","mrl","pmm","ppd","rs")
transform_fn <- c("identity","log_positive","boxcox_positive")
n_dat <- c(10,20,50)
censrate <- c(0.3,0.5,0.7)
cov_dep_cens <- c(0)
arg_list <- as.list(expand.grid(n_sam=n_dat,
                                censrate=censrate,
                                cens_mech=cov_dep_cens,
                                trans_fn=transform_fn))
attributes(arg_list) <- NULL
arg_df <- data.frame(id = seq_along(arg_list[[1]]),
                     n_sam=arg_list[[1]],
                     censrate=arg_list[[2]],
                     cens_mech=arg_list[[3]],
                     trans_fn=arg_list[[4]])
cond_filt <- purrr::pmap(arg_list,function(n_sam,censrate,cens_mech,trans_fn){
  conditions %>% 
    dplyr::filter(
      censoring_rate==censrate,
      cov_dep_cens==cens_mech,
      transform_fn==trans_fn,
      n_dat==n_sam,
    ) %>% 
    dplyr::select(ID) %>% 
    purrr::as_vector()
})

# number of clusters per simulation condition 180 clusters per simulation times 10 repetitions
length_theo <- 1800
plts <-   purrr::map(cond_filt,function(cond_n_dat){
  params_one_comb <- purrr::map(cond_n_dat,function(x){
  params_one <- params[[x]]
  params_one <- purrr::map(seq_along(params_one), function(x){
    if(class(params_one[[x]])[1]=="tbl_df"){
      return(params_one[[x]])
    }else{
      NA_tib <- params_one$rs %>% 
        mutate(p_val=NA,p_adj=NA,impType=names(params_one)[x])
      return(NA_tib)
    }
  })
   purrr::map(params_one, ~bind_rows(.x)) %>% bind_rows()
  }) %>% 
    dplyr::bind_rows()
  
  tib2plt <- params_one_comb %>% 
    mutate(proportion_censored = n_cens/n_dat,
           truestru = stringr::str_detect(cluster_id,"T"),
           class = as.factor(truestru),
           impType = factor(impType, levels = c("cc","pmm","km", "rs","mrl","ppd")),
           simType = factor(simType, levels = c("lm","glm","glmer")),
           b1_True = as.factor(b1_True), 
           n_dat= as.factor(n_dat),
           cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR")))
  
  
  cobratib <- tib2plt %>%
    mutate(id = seq_len(dim(tib2plt)[1])) 
  cobratib$impType <- droplevels(cobratib$impType)
  tmp_impTypes <- c()
  pval <- data.frame(rs = cobratib$p_val[cobratib$impType=="rs"])
  padj <- data.frame(rs = cobratib$p_adj[cobratib$impType=="rs"])
  for (impType in unique(cobratib$impType)) {
    na_ind <- is.na(cobratib$p_adj[cobratib$impType==impType])
    print(paste(impType,":",sum(na_ind),"NA's"))
    if (sum(na_ind)<length_theo){
      pval[[impType]] <- cobratib$p_val[cobratib$impType==impType]
      padj[[impType]] <- cobratib$p_adj[cobratib$impType==impType]
      tmp_impTypes <- c(tmp_impTypes,impType)
    }
  }
  tmp_impTypes <- sort(tmp_impTypes)
  rownames(pval) <- paste0("I_",cobratib$id[cobratib$impType=="rs"])
  rownames(padj) <- paste0("I_",cobratib$id[cobratib$impType=="rs"])
  
  truth <- data.frame(status = as.integer(cobratib$class[cobratib$impType=="rs"])-1,
                      cov_dep_cens = cobratib$cov_dep_cens[cobratib$impType=="rs"],
                      sim_id = cobratib$sim_id[cobratib$impType=="rs"],
                      row.names = paste0("I_",cobratib$id[cobratib$impType=="rs"])
  )
  cobradata <- COBRAData(pval = pval, truth = truth,padj = padj)
  cobraperf <- calculate_performance(cobradata,
                                     binary_truth = "status",
                                     aspects = c("fdrtpr", "fdrtprcurve",
                                                 "tpr", "roc"),
                                     thrs = c(0.01, 0.05, 0.1))
  print(tmp_impTypes)
  color_ind <- unlist(purrr::map(tmp_impTypes, ~which(.x==all_impTypes)))
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = c(RColorBrewer::brewer.pal(length(all_impTypes),"Dark2")[color_ind]), 
                                     facetted = TRUE)
  return(cobraplot)
})


################################################################################
## censoring
tmp_plt_ind <- arg_df %>% 
  filter(n_sam==20, trans_fn=="log_positive") 
tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
    plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)  + 
      scale_x_continuous(breaks = c(0.005,0.01,0.05,0.1,1),limits = c(0.005,1),trans = "log10") +
      labs(title=paste0("censoring rate: ",100*tmp_plt_ind$censrate[i],"%")) + 
      theme(aspect.ratio=1) + 
      guides(colour = guide_legend(nrow = 1))
})
plt_comb <- ggarrange(plotlist = tmp_plts,nrow=1,ncol=3, common.legend = TRUE)#,ncol=2)
plt_comb <- annotate_figure(plt_comb, 
                fig.lab = "Censoring rate",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_multicluster_effect_cens.png"), plt_comb, width = 30, height = 13, units = "cm")

################################################################################
## size
tmp_plt_ind <- arg_df %>% 
  filter(censrate==0.7, trans_fn=="log_positive") %>% 
  arrange(desc(n_sam))
tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
  plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)  + 
    scale_x_continuous(breaks = c(0.005,0.01,0.05,0.1,1),limits = c(0.005,1),trans = "log10") +
    labs(title=paste0("Size: ",tmp_plt_ind$n_sam[i])) + 
    theme(aspect.ratio=1) + 
    guides(colour = guide_legend(nrow = 1))
})
plt_comb <- ggarrange(plotlist = tmp_plts,nrow=1,ncol=3, common.legend = TRUE)#,ncol=2)
plt_comb <- annotate_figure(plt_comb, 
                            fig.lab = "Sample size",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_multicluster_effect_size.png"), plt_comb, width = 30, height = 13, units = "cm")

################################################################################
## transformation
tmp_plt_ind <- arg_df %>% 
  filter(censrate==0.7, n_sam==20) %>% 
  arrange(desc(n_sam))
tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
  plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)  + 
    scale_x_continuous(breaks = c(0.005,0.01,0.05,0.1,1),limits = c(0.005,1),trans = "log10") +
    labs(title=paste0("Transformation: ",tmp_plt_ind$trans_fn[i])) + 
    theme(aspect.ratio=1) + 
    guides(colour = guide_legend(nrow = 1))
})
plt_comb <- ggarrange(plotlist = list(tmp_plts[[2]],tmp_plts[[3]],tmp_plts[[1]]),nrow=1,ncol=3, common.legend = TRUE)#,ncol=2)
plt_comb <- annotate_figure(plt_comb, 
                            fig.lab = "Data transformation",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_multicluster_effect_transform.png"), plt_comb, width = 30, height = 13, units = "cm")

