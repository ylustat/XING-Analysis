#==================================================================#
#==================== codes to reproduce figures ==================#
#==================================================================#
library(tidyverse)
library(ggplot2)
library(data.table)
# Figure 2a. varied sample size
df <- fread("ResultFigure2a.txt")
df %>% 
  pivot_longer(ALG1:hteqtl) %>% 
  mutate(name = factor(name, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  ggplot(aes(x = n1, y = value, col = name))+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x = "Sample size",y = "AUC")+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+ 
  theme(legend.text.align = 0)+
  scale_y_continuous(limits = c(0.74,0.95))

# Figure 2b. cross-context correlation
df <- fread("ResultFigure2b.txt")
df %>% 
  pivot_longer(ALG1:hteqtl) %>% 
  mutate(name = factor(name, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  ggplot(aes(x = r12, y = value, col = name))+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x = "Correlation across contexts",y = "AUC")+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  theme(legend.text.align = 0)+
  scale_y_continuous(limits = c(0.75,0.95))

# Figure 2c. cross-data correlation
df <- fread("ResultFigure2c.txt")
df %>%
  pivot_longer(ALG1:hteqtl) %>% 
  mutate(name = factor(name, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  ggplot(aes(x = r21, y = value, col = name))+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x = "Correlation across data types",y = "AUC")+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  theme(legend.text.align = 0)+
  scale_y_continuous(limits = c(0.72,0.85))

# Figure 2d. proportion of uncorrelated trios
df <- fread("ResultFigure2d.txt")
df %>% 
  pivot_longer(ALG1:hteqtl) %>% 
  mutate(name = factor(name, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  ggplot(aes(x = prop_uncorrelated, y = value, col = name))+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x = "Proportion of tested units without\ncontext-shared information",y = "AUC")+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  theme(legend.text.align = 0)+
  scale_y_continuous(limits = c(0.58,0.9))


#==================================================================#
#=============== data generation and simulations ==================#
#==================================================================#
# See Simulation - SuppFigure1-7.R