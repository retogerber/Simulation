# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Simulated/Plot_scripts")


################################################################################
## load packages and data
library(tidyverse)
library(ggpubr)
library(iCOBRA)

dir_save <- "../../../Data/Simulated/"
plot_dir <- "../../../Plots/Simulated/"
methods <- c("cc","pmm","rs","km","mrl","ppd")
methods_no_cc <- c("pmm","rs","km","mrl","ppd")
test_results_null <- readRDS(paste0(dir_save,"test_results_normal_null.rds"))

test_results_null_proc <- test_results_null[[2]] %>% 
  mutate(impType = factor(impType, levels = c("cc","pmm","km", "rs","mrl","ppd")),
         simType = factor(simType, levels = c("lm","glm","glmer")),
         cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR")))

test_results <- readRDS(paste0(dir_save,"test_results_normal.rds"))

test_results_proc <- test_results[[2]] %>% 
  mutate(impType = factor(impType, levels = c("cc","pmm","km", "rs","mrl","ppd")),
         simType = factor(simType, levels = c("lm","glm","glmer")),
         cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR")))

# fdr_satis_mult %>% 
#   group_by(method) %>% 
#   summarise(mean = mean(satis)) %>% 
#   arrange(method)


### kulback leibler divergence
kulback_leibler_divergence <- function(x){
  tmp_cut <- cut_width(x,width = 0.1,center = 0.05)
  tmp_probs <- table(tmp_cut)/length(tmp_cut)
  theo_probs <- rep(0.1,10)
  map(seq_along(tmp_probs),function(i){
    tmp_probs[i]*log(tmp_probs[i]/theo_probs[i])
  }) %>% 
    reduce(sum)
}
### MSE
mse_fn <- function(x,b1_True_val=-1.5){
  (1/length(x))*sum((x-b1_True_val)^2)
}


cond_ls <- 
  list(
################################################################################
## kulbach leibler divergence
    g7_kls_n10 = list(
      data = test_results_null_proc %>%  filter(n_dat == 10),
      func = kulback_leibler_divergence,
      x = "pval",
      lvls = c(-Inf,0.05,0.1,Inf),
      name = "g7_n10",
      methods = methods_no_cc
    ),
    g7_kls_n20 = list(
      data = test_results_null_proc %>%  filter(n_dat == 20),
      func = kulback_leibler_divergence,
      x = "pval",
      lvls = c(-Inf,0.05,0.1,Inf),
      name = "g7_n20",
      methods = methods
    ),
    g7_kls_n50 = list(
      data = test_results_null_proc %>%  filter(n_dat == 50),
      func = kulback_leibler_divergence,
      x = "pval",
      lvls = c(-Inf,0.05,0.1,Inf),
      name = "g7_n50",
      methods = methods
    )
)

################################################################################


template <- list(
  mse_b1_template_name = list(
    data = data.frame(),
    func = function(x){mse_fn(x,-0.5)},
    x = "b1",
    lvls = c(-Inf,0.05,0.1,Inf),
    name = "mse_b1_template_name",
    methods = methods
  ),
  mean_pval_template_name = list(
    data = data.frame(),
    func = mean,
    x = "pval",
    lvls = c(-Inf,0.05,0.1,Inf),
    name = "mean_pval_template_name",
    methods = methods
  ),
  max_pval_template_name = list(
    data = data.frame(),
    func = max,
    x = "pval",
    lvls = c(-Inf,0.05,0.1,Inf),
    name = "max_pval_template_name",
    methods = methods
  )
)


create_list_entry <- function(name, data, methods = methods, templ = template){
  templ_out <- map(1:3, function(i){
    templ[[i]]$data <- data
    templ[[i]]$name <- 
      stringr::str_replace(templ[[i]]$name,"template_name",name)
    templ[[i]]
  })
  names(templ_out) <- stringr::str_replace_all(names(templ),"template_name",name)
  return(templ_out)
}

add_list_entry <- function(name, data, methods = methods, templ = template, tmp_cond_ls = cond_ls){
  tmp_entry <- create_list_entry(name = name, data = data, methods = methods, templ = templ)
  tmp_entry_names <- names(tmp_entry)
  for (i in 1:3) {
    tmp_cond_ls[[tmp_entry_names[i]]] <- tmp_entry[[i]]
  }
  return(tmp_cond_ls)
}

cond_ls <- add_list_entry("g1_identity", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="identity"
                            ))
cond_ls <- add_list_entry("g1_boxcox", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="boxcox_positive"
                            ))
cond_ls <- add_list_entry("g1_log", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ))

test_results_all_cens_all_n_dat <- test_results_proc %>%   
  filter(error_variance==0,
         variance_fixeff == 0.5, 
         variance_raneff == 0.1,
         simType == "glmer",
         b1_True == -0.5,
         cov_dep_cens=="MCAR", 
         transform_fn =="log_positive"
  )
cond_ls <- add_list_entry("g2_n10", 
                          test_results_all_cens_all_n_dat %>%   
                            filter(n_dat==10, censoring_rate==0.7),
                          methods = methods_no_cc)
cond_ls <- add_list_entry("g2_n20", 
                          test_results_all_cens_all_n_dat %>%   
                            filter(n_dat==20, censoring_rate==0.7))
cond_ls <- add_list_entry("g2_n50", 
                          test_results_all_cens_all_n_dat %>%   
                            filter(n_dat==50, censoring_rate==0.7))

cond_ls <- add_list_entry("g3_c30", 
                          test_results_all_cens_all_n_dat %>%   
                            filter(n_dat==20, censoring_rate==0.3))
cond_ls <- add_list_entry("g3_c50", 
                          test_results_all_cens_all_n_dat %>%   
                            filter(n_dat==20, censoring_rate==0.5))
cond_ls <- add_list_entry("g3_c70", 
                          test_results_all_cens_all_n_dat %>%   
                            filter(n_dat==20, censoring_rate==0.7))

cond_ls <- add_list_entry("g4_MCAR", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ))
cond_ls <- add_list_entry("g4_MAR", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ))

cond_ls <- add_list_entry("g5_low", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ))
cond_ls <- add_list_entry("g5_high", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 1,
                                   simType == "glmer",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ))

cond_ls <- add_list_entry("g6_GLMM", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ))
cond_ls <- add_list_entry("g6_GLM", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glm",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ))
cond_ls <- add_list_entry("g6_LM", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "lm",
                                   b1_True == -0.5,
                                   censoring_rate==0.7,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ))

cond_ls <- add_list_entry("g8_-0.5", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -0.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ))
template_b1_1 <- template
template_b1_1$mse_b1_template_name$func <- function(x){mse_fn(x,-1)}
cond_ls <- add_list_entry("g8_-1", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -1,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ),
                          templ = template_b1_1)
template_b1_1.5 <- template
template_b1_1.5$mse_b1_template_name$func <- function(x){mse_fn(x,-1.5)}
cond_ls <- add_list_entry("g8_-1.5", 
                          test_results_proc %>%   
                            filter(error_variance==0,
                                   variance_fixeff == 0.5, 
                                   variance_raneff == 0.1,
                                   simType == "glmer",
                                   b1_True == -1.5,
                                   cov_dep_cens=="MCAR", 
                                   n_dat == 20,
                                   censoring_rate == 0.7,
                                   transform_fn =="log_positive"
                            ),
                          templ = template_b1_1.5)
# cond <- cond_ls$fdr_satis_mult
sum_data <- map(cond_ls, function(cond){
  # method <- "cc"
  tmp_tib <- map(cond[["methods"]],function(method){
    data <- cond[["data"]][cond[["data"]]$impType==method,][[cond[["x"]]]]
    tibble(method_names = method, 
           !!cond[["name"]] := do.call(cond[["func"]],args=list(x=data)))
  }) %>% 
    bind_rows() %>% 
    rename(!!paste0(cond[["name"]],"_val") := !!cond[["name"]]) %>% 
    mutate(!!paste0(cond[["name"]]) := cut(!!sym(paste0(cond[["name"]],"_val")), 
                                                    breaks = cond[["lvls"]], 
                                                    labels = c(3,2,1)),
           !!paste0(cond[["name"]]) := as.numeric(as.character(!!sym(paste0(cond[["name"]])))))
}) %>% 
  reduce(full_join) %>% 
  mutate(method_names = factor(method_names,levels = c("cc","pmm","rs","km","mrl","ppd")))

keywords <- c("g1_identity", "g1_boxcox", "g1_log",
              "g2_n10","g2_n20","g2_n50",
              "g3_c30", "g3_c50", "g3_c70", 
              "g4_MAR", "g4_MCAR",
              "g5_low", "g5_high",
              "g6_GLMM", "g6_GLM", "g6_LM",
              "g8_-0.5", "g8_-1", "g8_-1.5")
sum_data_prep <- map(keywords, function(keyword){
  sum_data %>% select(method_names,contains(keyword)) %>% 
    mutate(!!paste0("mse_b1_",keyword) := as.numeric(as.character(!!sym(paste0("mse_b1_",keyword)))),
           !!paste0("mean_pval_",keyword) := as.numeric(as.character(!!sym(paste0("mean_pval_",keyword)))),
           !!paste0("max_pval_",keyword) := as.numeric(as.character(!!sym(paste0("max_pval_",keyword))))) %>% 
    mutate(!!paste0(keyword) := 
             round(0.5 * !!sym(paste0("mse_b1_",keyword)) + 
                     0.25 * !!sym(paste0("mean_pval_",keyword)) + 
                     0.25 * !!sym(paste0("max_pval_",keyword)))) %>% 
    select(method_names,paste0(keyword))
}) %>% 
  reduce(full_join)  %>% 
  full_join(sum_data %>% 
  select(-contains(keywords)))


sum_data_gath <- sum_data_prep %>% select(-contains("_val")) %>% 
  mutate(total_score = rowMeans(.[-1],na.rm = TRUE)) %>% 
  arrange(total_score) %>% 
  mutate(method_names = factor(method_names,levels = .data$method_names)) %>% 
  select(-total_score) %>% 
  gather( "test", "score",-method_names) %>% 
  mutate(score=factor(score,levels = c(3,2,1), labels = c("\nGood\n","\nIntermediate\n","\nPoor\n"))) 

sum_data_gath_gr <- sum_data_gath %>% 
  mutate(group = stringr::str_extract(test,"^g[:digit:]{1}"),
         test = stringr::str_sub(test,start = 4))
labels_to_use <- c("Transformation", 
                   "Sample \n Size (n)",
                   "Censoring \n rate (c)", 
                   "Censoring \n Type",
                   "Variance \n Random \n Effect",
                   "Regression \n Type",
                   "P-values \n Null", 
                   "Regression \n Coefficient")
for (i in seq_along(labels_to_use)) {
  sum_data_gath_gr$group <- stringr::str_replace_all(sum_data_gath_gr$group,paste0("g",i),labels_to_use[i])
}
sum_data_gath_gr$group <- factor(sum_data_gath_gr$group, 
                                 levels = labels_to_use)

plt <- ggplot(sum_data_gath_gr ) + 
  geom_tile(aes(method_names,test,fill=score), width=0.9, height=0.9,color = "black") + 
  scale_fill_manual(values=RColorBrewer::brewer.pal(11,"RdYlBu")[c(9,5,3)], 
                      na.translate=FALSE) +
  # theme_bw() + 
  facet_grid(group~., scales = "free", space = "free") +
  theme(line = element_blank(),
        legend.title = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  # scale_fill_discrete() + 
  labs(x=element_blank(),y=element_blank())
plt
ggsave(paste0(plot_dir,"summary.png"), plt, width = 20, height = 25, units = "cm")

