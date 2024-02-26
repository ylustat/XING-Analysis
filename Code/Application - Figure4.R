library(ggplot2)
library(tidyverse)
library(data.table)
library(gg.gap)
library(ggnewscale)

load("transmedcisfull.RData")
fold <- 2
calculate_total_eff_of_tissues <- function(total_eff) {
  rd <- total_eff %>% apply(1,function(x) {
    if(sum(!is.na(x))==0) return(NA)
    pos <- x[x>0]; neg <- x[x<0]
    if(sign(x[which.max(abs(x))])>0){
      grp = pos
    } else {
      grp = neg
    }
    grp <- abs(grp)
    n <- sum(grp > (1/fold) * max(grp,na.rm = T),na.rm = T)
    return(n)
  }) 
  return(rd)
}

# trans-eQTL, total effect
df_tab[[1]]$cis_eQTL_num[is.na(df_tab[[1]]$cis_eQTL_num)] <- 0
ind_keep_noindirect <- df_tab[[1]] %>% 
  filter(trans_gene == 1 & cis_gene == 0) %>% 
  select(probe,variant_id,ind) %>%
  distinct_at(vars(probe, variant_id),.keep_all = T)

indirect_exp_cross <- lapply(df_tab,function(x) x$indirect_exp) %>% 
  do.call("cbind",.) %>% apply(1,function(x) max(abs(x)))
ind_keep_indirect <- df_tab[[1]] %>% 
  mutate(indirect_exp_cross = indirect_exp_cross) %>% 
  filter(trans_gene == 1 &  cis_gene == 1) %>% 
  select(probe,variant_id,indirect_exp_cross,ind) %>%
  group_by(probe,variant_id) %>% 
  slice_max(abs(indirect_exp_cross),n = 1, with_ties = F) %>% 
  select(-indirect_exp_cross)

ind_keep <- rbind(ind_keep_noindirect,ind_keep_indirect)
df_tab_sub <- df_tab %>% lapply(function(x) {
  x <- x %>% filter(ind %in% ind_keep$ind)
  return(x)
})
total_eff <- lapply(df_tab_sub, function(x) x$base_full) %>% do.call("cbind",.)
indirect_eff <- lapply(df_tab_sub, function(x) x$indirect_exp) %>% do.call("cbind",.)
trans_num <- df_tab_sub[[1]]$trans_num
trans_gene_indicator <- df_tab_sub[[1]]$trans_gene
cis_gene_indicator <- df_tab_sub[[1]]$cis_gene
cis_eQTL_num_indicator <- df_tab_sub[[1]]$cis_eQTL_num
total_eff <- total_eff[trans_gene_indicator==1,]
trans_num <- trans_num[trans_gene_indicator==1]

cis_gene_indicator <- cis_gene_indicator[trans_gene_indicator==1]
cis_eQTL_num_indicator <- cis_eQTL_num_indicator[trans_gene_indicator==1]
number_tissue_total_eff <- calculate_total_eff_of_tissues(total_eff)
number_tissue_indirect_eff <- calculate_total_eff_of_tissues(indirect_eff)


source('~/Downloads/gg.gap-master/R/gg.gap.R')
binwid <- .33
trans_num_category <- cut(trans_num,breaks = c(-Inf,1,4,Inf),labels = c(1,2,3)) %>% as.numeric()
fig <- data.frame(total = number_tissue_total_eff, 
                  indirect = number_tissue_indirect_eff,
                  sign.num = trans_num_category) %>%
  pivot_longer(cols = total:indirect,values_to = "eff",names_to = "type") %>% 
  group_by(type,eff,sign.num) %>% summarise(n=n()) %>% na.omit() %>%
  group_by(type,sign.num) %>% mutate(prop = n/sum(n))
fig <- ggplot()+
  geom_col(data = fig %>% filter(type == "total"), aes(x = sign.num-binwid/2, y = n, fill = eff),
           width = binwid,position = "stack",size=0.1)+
  scale_fill_gradient(low="white",high = "mediumpurple",
                      name="# tissues with shared\ntotal trans-eQTL effect",breaks=1:4*2)+
  new_scale_fill()+
  geom_col(data = fig %>% filter(type == "indirect"), aes(x = sign.num+binwid/2, y = n, fill = eff),
           width = binwid,position = "stack",size=0.1)+
  scale_fill_gradient(low="white",high = "dodgerblue1",
                      name="# tissues with shared\nindirect trans-eQTL effect",breaks=1:4*2)+
  new_scale_fill()+
  geom_col(data = fig %>% group_by(type,sign.num) %>% summarise(sums=sum(n)) %>% filter(type == "total"),
           aes(x = sign.num-binwid/2, y = sums), fill = "white",
           width = binwid,position = "stack",size=.15,alpha=0,col="black")+
  scale_fill_gradient(low="white",high = "gray")+
  new_scale_fill()+
  geom_col(data = fig %>% group_by(type,sign.num) %>% summarise(sums=sum(n)) %>% filter(type == "indirect"),
           aes(x = sign.num+binwid/2, y = sums), fill = "white",
           width = binwid,position = "stack",size=0.15,alpha=0,col="black")+
  scale_fill_gradient(low="white",high = "blue")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white"))+
  scale_x_continuous(breaks = 1:3,labels=c("1","2-4","5 or more")) +
  labs(x = "# tissues with non-zero trans-eQTL effect",
       y = "# trans-eQTLs")+
  theme(legend.key = element_rect(colour='grey32'))+
  geom_line(aes(x = c(1-binwid,1+binwid),y = c(150,150)), colour = "grey15",linetype = "dashed")+
  geom_line(aes(x = c(1-binwid,1+binwid),y = c(1000,1000)), colour = "grey15",linetype = "dashed")
legend <- cowplot::get_legend(fig)
ggfig <- gg.gap1(plot=fig,
                 segments=c(150,1000),
                 tick_width = c(50,9000),
                 rel_heights=c(0.5,3.5,0,2),
                 ylim=c(0,20000))+theme(axis.text=element_text(size=18),
                                        axis.title=element_text(size=18),
                                        axis.title.y = element_text(size = 24))+
  theme(plot.margin = unit(c(0,0,.0,.5), "cm"))
ggdraw(legend)
ggfig



  
  
  # trans mQTL, total effect
  ind_keep_noindirect <- df_tab[[1]] %>% 
    filter(trans_gene == 0 & cis_cpg == 0) %>% 
    select(probe,variant_id,ind) %>%
    distinct_at(vars(probe, variant_id),.keep_all = T)
  
  indirect_meth_cross <- lapply(df_tab,function(x) x$indirect_meth) %>% 
    do.call("cbind",.) %>% apply(1,function(x) max(abs(x)))
  ind_keep_indirect <- df_tab[[1]] %>% 
    mutate(indirect_meth_cross = indirect_meth_cross) %>% 
    filter(trans_gene == 0 &  cis_cpg == 1) %>% 
    select(probe,variant_id,indirect_meth_cross,ind) %>%
    group_by(probe,variant_id) %>% 
    slice_max(abs(indirect_meth_cross),n = 1, with_ties = F) %>% 
    select(-indirect_meth_cross)
  
  ind_keep <- rbind(ind_keep_noindirect,ind_keep_indirect)
  
  df_tab_sub <- df_tab %>% lapply(function(x) {
    x <- x %>% filter(ind %in% ind_keep$ind)
    return(x)
  })
  total_eff <- lapply(df_tab_sub, function(x) x$base_full) %>% do.call("cbind",.)
  indirect_eff <- lapply(df_tab_sub, function(x) x$indirect_meth) %>% do.call("cbind",.)
  trans_num <- df_tab_sub[[1]]$trans_num
  trans_gene_indicator <- df_tab_sub[[1]]$trans_gene
  cis_gene_indicator <- df_tab_sub[[1]]$cis_cpg
  total_eff <- total_eff[trans_gene_indicator==0,]
  trans_num <- trans_num[trans_gene_indicator==0]
  indirect_eff <- indirect_eff[trans_gene_indicator==0,]
  cis_gene_indicator <- cis_gene_indicator[trans_gene_indicator==0]
  number_tissue_total_eff <- calculate_total_eff_of_tissues(total_eff)
  number_tissue_indirect_eff <- calculate_total_eff_of_tissues(indirect_eff)
  
  binwid = 0.35
  trans_num_category <- cut(trans_num,breaks = c(-Inf,1,4,Inf),labels = c(1,2,3)) %>% as.numeric()
  fig <- data.frame(total = number_tissue_total_eff, 
                    indirect = number_tissue_indirect_eff,
                    sign.num = trans_num_category) %>%
    pivot_longer(cols = total:indirect,values_to = "eff",names_to = "type") %>% 
    group_by(type,eff,sign.num) %>% summarise(n=n()) %>% na.omit() %>%
    group_by(type,sign.num) %>% mutate(prop = n/sum(n))
  fig <- ggplot()+
    geom_col(data = fig %>% filter(type == "total"), aes(x = sign.num-binwid/2, y = n, fill = eff),
             width = binwid,position = "stack",size=0.1)+
    scale_fill_gradient(low="white",high = "mediumpurple",name="# tissues with shared\ntotal trans-mQTL effect",breaks=1:4*2)+
    new_scale_fill()+
    geom_col(data = fig %>% filter(type == "indirect"), aes(x = sign.num+binwid/2, y = n, fill = eff),
             width = binwid,position = "stack",size=0.1)+
    scale_fill_gradient(low="white",high = "dodgerblue1",name="# tissues with shared\nindirect trans-mQTL effect",1:4*2)+
    new_scale_fill()+
    geom_col(data = fig %>% group_by(type,sign.num) %>% summarise(sums=sum(n)) %>% filter(type == "total"),
             aes(x = sign.num-binwid/2, y = sums), fill = "white",
             width = binwid,position = "stack",size=0.15,alpha=0,col="black")+
    scale_fill_gradient(low="white",high = "gray")+
    new_scale_fill()+
    geom_col(data = fig %>% group_by(type,sign.num) %>% summarise(sums=sum(n)) %>% filter(type == "indirect"),
             aes(x = sign.num+binwid/2, y = sums), fill = "white",
             width = binwid,position = "stack",size=0.15,alpha=0,col="black")+
    scale_fill_gradient(low="white",high = "blue")+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill="white"))+
    scale_x_continuous(breaks = 1:3,labels=c("1","2-4","5 or more")) +
    labs(x = "# tissues with non-zero trans-mQTL effect",
         y = "# trans-mQTLs")+
    theme(legend.key = element_rect(colour='grey32'))+
    geom_line(aes(x = c(1-binwid,1+binwid),y = c(80,80)), colour = "grey15",linetype = "dashed")+
    geom_line(aes(x = c(1-binwid,1+binwid),y = c(1000,1000)), colour = "grey15",linetype = "dashed")
  
  ggfig <- gg.gap1(plot=fig,
                   segments=c(80,1000),
                   tick_width = c(40,7000),
                   rel_heights=c(0.5,3.5,0,2),
                   ylim=c(0,15000))+theme(axis.text=element_text(size=18),
                                          axis.title=element_text(size=18),
                                          axis.title.y = element_text(size = 24))+
    theme(plot.margin = unit(c(0,0,.0,.5), "cm"))
  legend <- cowplot::get_legend(fig)
  ggfig
