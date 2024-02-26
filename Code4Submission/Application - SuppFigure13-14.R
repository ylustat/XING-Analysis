library(ggplot2)
library(tidyverse)
library(data.table)

df <- fread("/Users/luyihao/Downloads/transexpmed.txt")
fdr_level <- max(df$exp_mediator_pval[which(qvalue::qvalue(df$exp_mediator_pval)$qvalue<.05)])

df %>% 
  filter(trans_num == 1) %>% 
  mutate(significance = ifelse(exp_mediator_pval<fdr_level,1,0)) %>% 
  mutate(cis_exp_col = paste0(trans_num, significance)) %>% 
  ggplot() +
  geom_point(aes(x = trans_mediated_by_e, 
                 y = -log10(exp_mediator_pval), 
                 fill = cis_exp_col),pch=21,size = 2.5)+
  theme_classic() +
  labs(y = expression('-log'[10]*'(P-value of mediation)'),
       x = "% reduction in trans-eQTL effects") +
  scale_fill_manual(values = c("grey80","royalblue1","grey45","red1"),name="Type of cis-mediator",
                    labels=c("mQTL-specific","mQTL-specific","mQTL with eQTL effect","mQTL-specific")) +
  guides(color="none",fill="none")+theme(strip.background = element_blank(),
                                         strip.text.x = element_blank(),
                                         panel.spacing = unit(3, "lines"))+
  geom_vline(xintercept = 0, size = 1, col = "grey", linetype = "dashed", alpha = 0.5)


df %>% 
  filter(trans_num >= 2) %>% 
  mutate(significance = ifelse(exp_mediator_pval<fdr_level,1,0)) %>% 
  mutate(cis_exp_col = paste0(trans_num, significance)) %>% 
  ggplot() +
  geom_point(aes(x = trans_mediated_by_e, 
                 y = -log10(exp_mediator_pval), 
                 fill = cis_exp_col),pch=21,size = 2.5)+
  theme_classic() +
  labs(y = expression('-log'[10]*'(P-value of mediation)'),
       x = "% reduction in trans-eQTL effects") +
  scale_fill_manual(values = c("grey80","darkorange","grey45","red1"),name="Type of cis-mediator",
                    labels=c("mQTL-specific","mQTL-specific","mQTL with eQTL effect","mQTL-specific")) +
  guides(color="none",fill="none")+theme(strip.background = element_blank(),
                                         strip.text.x = element_blank(),
                                         panel.spacing = unit(3, "lines"))+
  xlim(c(-0.5,1.2))+
  geom_vline(xintercept = 0, size = 1, col = "grey", linetype = "dashed", alpha = 0.5)




df <- fread("/Users/luyihao/Downloads/transmethmed.txt")

fdr_level <- max(a$meth_mediator_pval[which(qvalue::qvalue(a$meth_mediator_pval)$qvalue<.05)])
df %>% 
  filter(trans_num == 1) %>% 
  mutate(significance = ifelse(meth_mediator_pval<fdr_level,1,0)) %>% 
  mutate(cis_cpg_col = paste0(trans_num, significance)) %>% 
  ggplot() +
  geom_point(aes(x = trans_mediated_by_m, 
                 y = -log10(meth_mediator_pval), 
                 fill = cis_cpg_col),pch=21,size = 2.5)+
  theme_classic() +
  labs(y = expression('-log'[10]*'(P-value of mediation)'),
       x = "% reduction in trans-mQTL effects") +
  scale_fill_manual(values = c("grey80","royalblue1","grey45","red1"),name="Type of cis-mediator",
                    labels=c("mQTL-specific","mQTL-specific","mQTL with eQTL effect","mQTL-specific")) +
  guides(color="none",fill="none")+theme(strip.background = element_blank(),
                                         strip.text.x = element_blank(),
                                         panel.spacing = unit(3, "lines"))+
  xlim(c(-1,1))+
  geom_vline(xintercept = 0, size = 1, col = "grey", linetype = "dashed", alpha = 0.5)

df %>% 
  filter(trans_num >= 2) %>% 
  mutate(significance = ifelse(meth_mediator_pval<fdr_level,1,0)) %>% 
  mutate(cis_cpg_col = paste0(trans_num, significance)) %>% 
  ggplot() +
  geom_point(aes(x = trans_mediated_by_m, 
                 y = -log10(meth_mediator_pval), 
                 fill = cis_cpg_col),pch=21,size = 2.5)+
  theme_classic() +
  labs(y = expression('-log'[10]*'(P-value of mediation)'),
       x = "% reduction in trans-mQTL effects") +
  scale_fill_manual(values = c("grey80","darkorange","grey45","red1"),name="Type of cis-mediator",
                    labels=c("mQTL-specific","mQTL-specific","mQTL with eQTL effect","mQTL-specific")) +
  guides(color="none",fill="none")+theme(strip.background = element_blank(),
                                         strip.text.x = element_blank(),
                                         panel.spacing = unit(3, "lines"))+
  xlim(c(-1,1))+
  geom_vline(xintercept = 0, size = 1, col = "grey", linetype = "dashed", alpha = 0.5)




