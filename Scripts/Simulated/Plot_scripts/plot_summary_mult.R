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

cobra_list <- readRDS(paste0(dir_save,"multicluster_cobra_data.rds"))
cobra_data <- cobra_list$cobradata
cobra_meta <- cobra_list$arg_df

template <- 
  list(
################################################################################
## Error corrections   
  fdr_satis_template_name = list(
    data = data.frame(),
    func = mean,
    x = "satis",
    lvls = c(-Inf,1/3+0.001,2/3+0.001,Inf),
    name = "Error correction_template_name",
    methods = methods
  ),
  tpr_template_name = list(
    data = data.frame(),
    func = function(x){1-median(x)},
    x = "TPR",
    lvls = c(-Inf,0.2,0.4,Inf), # limits of median of 0.8 and 0.6
    name = "TPR_template_name",
    methods = methods
  ),
  fdr_template_name = list(
    data = data.frame(),
    func = function(x){median(x)},
    x = "FDR",
    lvls = c(-Inf,0.2,0.4,Inf),
    name = "FDR_template_name",
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

add_list_entry_simplified <- function(name, filtering_list, data=cobra_data, meta=cobra_meta, methods = methods, templ = template, tmp_cond_ls = cond_ls){
  tmp_id <- meta %>% 
    filter(trans_fn == filtering_list[["trans_fn"]], 
           n_sam==filtering_list[["n_sam"]],
           censrate==filtering_list[["censrate"]],
           cens_mech==filtering_list[["cens_mech"]]) %>% 
    select(id) %>% 
    as_vector()
  fdr_satis <- map(tmp_id, ~ fdrtpr(data[[.x]])[ ,c("method","satis","FDR","TPR")]) %>% 
    reduce(bind_rows) %>% 
    mutate(satis = ifelse(satis == "yes",0,1)) %>% 
    rename(impType = method)
  return(add_list_entry(name = name, data = fdr_satis, methods = methods, templ = templ, tmp_cond_ls = tmp_cond_ls))
}

cond_ls <- list()
cond_ls <- add_list_entry_simplified("g1_identity",list(trans_fn="identity",
                                    n_sam=20,
                                    censrate=0.7,
                                    cens_mech=0))
cond_ls <- add_list_entry_simplified("g1_log",list(trans_fn="log_positive",
                                                         n_sam=20,
                                                         censrate=0.7,
                                                         cens_mech=0))
cond_ls <- add_list_entry_simplified("g1_boxcox",list(trans_fn="boxcox_positive",
                                                    n_sam=20,
                                                    censrate=0.7,
                                                    cens_mech=0))

cond_ls <- add_list_entry_simplified("g2_10",list(trans_fn="log_positive",
                                                         n_sam=10,
                                                         censrate=0.7,
                                                         cens_mech=0))
cond_ls <- add_list_entry_simplified("g2_20",list(trans_fn="log_positive",
                                                    n_sam=20,
                                                    censrate=0.7,
                                                    cens_mech=0))
cond_ls <- add_list_entry_simplified("g2_50",list(trans_fn="log_positive",
                                                       n_sam=50,
                                                       censrate=0.7,
                                                       cens_mech=0))

cond_ls <- add_list_entry_simplified("g3_0.3",list(trans_fn="log_positive",
                                                    n_sam=20,
                                                    censrate=0.3,
                                                    cens_mech=0))
cond_ls <- add_list_entry_simplified("g3_0.5",list(trans_fn="log_positive",
                                                    n_sam=20,
                                                    censrate=0.5,
                                                    cens_mech=0))
cond_ls <- add_list_entry_simplified("g3_0.7",list(trans_fn="log_positive",
                                                    n_sam=20,
                                                    censrate=0.7,
                                                    cens_mech=0))

# cond <- cond_ls[[3]]
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
              "g2_10","g2_20","g2_50",
              "g3_0.3", "g3_0.5", "g3_0.7")
sum_data_prep <- map(keywords, function(keyword){
  sum_data %>% select(method_names,contains(keyword)) %>% 
    mutate(!!paste0("Error correction_",keyword) := as.numeric(as.character(!!sym(paste0("Error correction_",keyword)))),
           !!paste0("TPR_",keyword) := as.numeric(as.character(!!sym(paste0("TPR_",keyword)))),
           !!paste0("FDR_",keyword) := as.numeric(as.character(!!sym(paste0("FDR_",keyword))))) %>% 
    mutate(!!paste0(keyword) := 
             round(0.5 * !!sym(paste0("Error correction_",keyword)) + 
                     0.25 * !!sym(paste0("TPR_",keyword)) + 
                     0.25 * !!sym(paste0("FDR_",keyword)))) %>% 
    select(method_names,paste0(keyword))
}) %>% 
  reduce(full_join)  %>% 
  full_join(sum_data %>% 
              select(-contains(keywords)))



sum_data_gath <- sum_data %>% select(-contains("_val")) %>% 
  mutate(total_score = rowMeans(.[-1],na.rm = TRUE)) %>% 
  arrange(total_score) %>% 
  mutate(method_names = factor(method_names,levels = .data$method_names)) %>% 
  select(-total_score) %>% 
  gather( "test", "score",-method_names) %>% 
  mutate(score=factor(score,levels = c(3,2,1), labels = c("\nGood\n","\nIntermediate\n","\nPoor\n"))) 

sum_data_gath_gr <- sum_data_gath %>% 
  mutate(group = stringr::str_extract(test,"g[:digit:]{1}"),
         group_inner = stringr::str_extract(test,"g[:digit:]{1}_[:graph:]+") %>% stringr::str_sub(start=4),
         test = stringr::str_extract(test,"[[:graph:]]+_") %>% stringr::str_sub(end=-5))
labels_to_use <- c("Transf.",
                   "Sample \n Size ",
                   "Censoring \n rate ")

labels_to_use_inner <- list(Transformation = c("identity","log","boxcox"))
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
  facet_grid(group+group_inner~., scales = "free", space = "free") +
  theme(line = element_blank(),
        legend.title = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill="white"),
        text = element_text(size=15)) +
  # scale_fill_discrete() + 
  labs(x=element_blank(),y=element_blank())
plt
ggsave(paste0(plot_dir,"summary_multicluster.png"), plt, width = 20, height = 25, units = "cm")

