# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts")


################################################################################
## load packages and data
library(tidyverse)
library(ggpubr)
dir_save <- "../../../Data/Simulated/"
plot_dir <- "../../../Plots/Simulated/"


res <- readRDS(paste0(dir_save,"/res_plot_mi_reps.rds"))
mi_reps <- c(2,3,4,5,7,10,20,30,40,50,75,100,125,150,200,300,400,500)

################################################################################
## plots

res_df <- bind_rows(res) %>% 
  rename("p-value":=pval)
res_long <- pivot_longer(data=res_df,cols=c("estimate","std_err","df","p-value","user_cpu","system_cpu","elapsed_time"))  %>% 
  arrange(method_est,mi_rep) %>% 
  mutate(sim_cond = factor(sim_id,levels = 1:6,labels=c("c0.3-MCAR","c0.5-MCAR","c0.7-MCAR","c0.3-MAR","c0.5-MAR","c0.7-MAR")))


plt0 <- ggplot(expand.grid(censoring=c(0.3,0.5,0.7),mdm = c("MCAR","MAR"))) +
  geom_tile(aes(x = factor(censoring), y = factor(mdm),fill = factor(1:6))) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(4,"Greens")[-1],RColorBrewer::brewer.pal(4,"Blues")[-1])) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(line = element_blank(),rect = element_blank(),legend.position = "none",axis.text = element_text(size=10)) +
  labs(x=element_blank(),y=element_blank())

plt1 <- ggplot(res_long %>% filter(name %in% c("p-value"))) +
  geom_point(aes(mi_rep,value,color=sim_cond))+
  geom_line(aes(mi_rep,value,color=sim_cond)) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(4,"Greens")[-1],RColorBrewer::brewer.pal(4,"Blues")[-1])) +
  scale_x_continuous(breaks = c(0,mi_reps[-c(1,2,3,4,5,6,7,8,9,11,13)])) +
  scale_y_continuous(breaks = c(0,0.05,0.25,0.5,0.75,1)) +
  facet_wrap(name~method_est,ncol = 5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size=0.2),
        panel.grid.major.x = element_line(size=0.3),
        legend.position = "none") +
  guides(color=guide_legend(nrow=2,byrow = TRUE)) +
  labs(x="",y="",color="Simulation\nCondition")

plt2 <- ggplot(res_long %>% filter(name %in% c("estimate"))) +
  geom_point(aes(mi_rep,value,color=sim_cond))+
  geom_line(aes(mi_rep,value,color=sim_cond)) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(4,"Greens")[-1],RColorBrewer::brewer.pal(4,"Blues")[-1])) +
  scale_x_continuous(breaks = c(0,mi_reps[-c(1,2,3,4,5,6,7,8,9,11,13)])) +
  scale_y_continuous(breaks = c(0,2,4,6)) +
  facet_wrap(name~method_est,ncol = 5) +
  geom_hline(yintercept = 1) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size=0.2),
        panel.grid.major.x = element_line(size=0.3),
        legend.position = "none") +
  labs(x="Multiple Imputation repetitions",y="",color="Simulation\nCondition")
plt_empty <- ggplot(res_long) + geom_blank() + theme_bw() + theme(line = element_blank(),rect = element_blank())
plt_comb <- ggarrange(ggarrange(plt_empty,plt0,plt_empty,ncol = 3),plt1,plt2,nrow=3, heights = c(1,5,5))
ggsave(paste0(plot_dir,"pval_estimate_mi_reps.png"),plt_comb, width = 12,height = 8)

