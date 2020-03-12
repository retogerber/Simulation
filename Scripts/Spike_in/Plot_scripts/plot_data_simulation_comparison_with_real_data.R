library(dplyr)
devtools::load_all("/home/reto/polybox/ETH/Master_Thesis/Code/censcyt")

main_dir <- "/home/reto/polybox/ETH/Master_Thesis/Portmacquarie"
data_type <- c("Preprocessed","PositiveSignal")
filename_pattern <- c("_pp","_pp_ps")
survival_time_scale <- c("original","relative","log")
clusterings <- c("meta10","meta20","som100")
clus_pats <- expand.grid(filename_pattern,clusterings,stringsAsFactors = FALSE)


cat("prepare data ...")
data_comb <- purrr::map2(rep(seq_along(filename_pattern),length(clusterings)),rep(clusterings,each=length(filename_pattern)), function(iii,ii){
  
  tmplist <- readRDS(paste0(main_dir,"/FlowCap/results/", data_type[iii] ,"/censdiffcyt_Framework/res", filename_pattern[iii] ,"_",
                            ii,"_counts_formula_contrast.rds"))
  
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
tmptib <- tibble::tibble(
  y = SummarizedExperiment::assay(data_comb[[2]]$d_counts)[1,]/ SummarizedExperiment::rowData(data_comb[[2]]$d_counts)$n_cells[1],
  x = data_comb[[1]]$formula$original$data$Survival.Time,
  I = data_comb[[1]]$formula$original$data$Status
)



set.seed(123)
cens_scales <- scales_for_censoring(
  censoring_val = c(0.3,0.4,0.5,0.6,0.7,0.8),
  log_ratio_val = c(0),
  seq_to_test = c(seq(0.01,0.99,length.out = 50), seq(1,10,by=0.5)),
  n = 10000,
  shapes = list(shape_x = 0.5,
                shape_c1 = 1,
                shape_c2 = 0.5),
  scale_x = 0.25)


formula_cens_glmm <- formula(y~Surv(X,I)+z+(1|r))
set.seed(1234)
data_sim <- simulate_data(n = 382, formula = formula_cens_glmm, type = "glmer",
                          error_variance = 0,
                          b = c(-6,0,0),
                          weibull_params = list(X = list(shape = 1, scale = 0.25),
                                                C = list(shape = 1, scale = cens_scales$C1[5])),
                          variance_fixeff = 0.5,
                          variance_raneff = 1,
                          n_levels_raneff = NULL,
                          n_levels_fixeff = 2,
                          censoring_dependent_on_covariate = FALSE,
                          weibull_params_covariate_dependent_censoring = list(shape = 0.5, scale = 0.09),
                          transform_fn = "identity",
                          verbose = TRUE)
b0 <- -5
b1 <- -0.5
data_sim$z <- factor(data_sim$z,labels = c(0,1))
maxx <- max(max(data_sim$y),max(tmptib$y))
data_sim$TrVal <- data_sim$TrVal/max(data_sim$TrVal)*max(tmptib$x)
plt_sim <- ggplot(data_sim) +
  geom_point(aes(TrVal,y,color=factor(I,levels = c(1,0),labels = c("observed","censored"))),size=0.5)+
  labs(x="x",color="",title = "Simulated")+
  ylim(c(0,maxx)) +
  geom_label(aes(x=0.8*max(TrVal),y=maxx,label=paste(round((1-sum(data_sim$I)/dim(data_sim)[1])*100,1),"% Censored"))) + 
  guides(color = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values=scales::hue_pal()(3)[2:1])
plt_sim_hist <- ggplot(data_sim) +
  geom_histogram(aes(TrVal,fill="Total"),alpha=0.3,binwidth = 100)+
  geom_histogram(aes(TrVal,fill=factor(I,levels = c(1,0),labels = c("observed","censored"))),alpha=0.5,position = "identity",binwidth = 70)+
  labs(x="x",fill="")  +
  ylim(c(0,52))
plt_example <- ggplot(tmptib) +
  geom_point(aes(x,y,color=factor(I,levels=c(1,0),labels = c("observed","censored"))),size=0.5)+
  labs(color="",title = "Real")+
  ylim(c(0,maxx)) +
  geom_label(aes(x=0.8*max(x),y=maxx,label=paste(round((1-sum(tmptib$I)/dim(tmptib)[1])*100,1),"% Censored"))) + 
  guides(color = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values=scales::hue_pal()(3)[2:1])
plt_example_hist <- ggplot(tmptib) +
  geom_histogram(aes(x,fill="Total"),alpha=0.3,binwidth = 100)+
  geom_histogram(aes(x,fill=factor(I,levels = c(1,0),labels = c("observed","censored"))),alpha=0.5,position = "identity",binwidth = 70)+
  labs(x="x",fill="") +
  ylim(c(0,52))
  
plts <- ggpubr::ggarrange(plt_sim,plt_example,common.legend = TRUE,align = "v")
hists <- ggpubr::ggarrange(plt_sim_hist,plt_example_hist,common.legend = TRUE,align = "v")
aligned_plts <- ggpubr::ggarrange(plts, hists,align = "v",nrow=2)
ggsave("/home/reto/polybox/ETH/Master_Thesis/plots/Simulations/comparison_real_simulated_data.png",aligned_plts,width = 20,height = 20,units = "cm")
