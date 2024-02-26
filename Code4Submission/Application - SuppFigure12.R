library(tidyverse)
library(ggplot2)
library(data.table)
library(ggrepel)
library(scales)
library(gg.gap)
library(ggbreak)

criterion <- "top_SNP"
h2 <- read.csv(paste0("/Users/luyihao/Downloads/remove.",criterion,"_h2.result.qtl.csv"))
h2 <- na.omit(h2)
colnames(h2)[1:5] <- c("index","file","h2.x","h2.eQTL","h2.mQTL")
h2$Type[h2$Type == "d"] <- "Disease\n(N=31)"
h2$Type[h2$Type == "t"] <- "Trait\n(N=19)"

# violin plot
h2 <- h2 %>% 
  mutate(`Trans-eQTL hotspot` = (100 * (1 - h2.eQTL/h2.x)), 
         `Trans-mQTL hotspot` = (100 * (1 - h2.mQTL/h2.x)),
         `eQTL.mean.change` = (100 * (1 - h2.eQTL/h2.x)/eQTL_num_hotspot), 
         `mQTL.mean.change` = (100 * (1 - h2.mQTL/h2.x)/mQTL_num_hotspot),
         e.abs.diff = h2.x - h2.eQTL,
         m.abs.diff = h2.x - h2.mQTL,
         e.abs.diff.perSNP = (h2.x - h2.eQTL)/eQTL_num_hotspot,
         m.abs.diff.perSNP = (h2.x - h2.mQTL)/mQTL_num_hotspot)

# eQTL hotspot
gg1 <- h2 %>%
  pivot_longer(cols = c(`eQTL.mean.change`,`mQTL.mean.change`),values_to = "h2_diff",names_to = "type") %>%
  filter(type == "eQTL.mean.change") %>%
  ggplot(aes(x = Type,y = h2_diff)) +
  geom_violin(aes(fill=Type),outlier.shape = NA,width=0.95,alpha=0.5) +
  geom_boxplot(outlier.shape = NA,width=0.1,fill="gold") +
  facet_wrap(~type,scales = "free")+
  labs(x = "", y = "") + # Change of heritability per hotspot region (%)
  theme_classic()+scale_fill_manual(values=c("#f0650e", "#0091ff"),label = c("Disease","Trait")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none")
trans <- function(x) ifelse(x > 0.8,0.8+(x-0.8)/6, x)
inv <- function(x) ifelse(x > 0.8, 6*(x-0.8) + 0.8, x)
my_trans <- trans_new("my_trans", trans, inv)
gg1 <- gg1 + 
  scale_y_break(c(0.8,1.4), scales = 1.5)+
  scale_y_continuous(trans = my_trans,
                     breaks = c(seq(0, 0.8, by = 0.4),seq(0.8,8,1.6)))

print(gg1)

# mQTL hotspot
gg1 <- h2 %>%
  pivot_longer(cols = c(`eQTL.mean.change`,`mQTL.mean.change`),values_to = "h2_diff",names_to = "type") %>%
  filter(type == "mQTL.mean.change") %>%
  ggplot(aes(x = Type,y = h2_diff)) +
  geom_violin(aes(fill=Type),outlier.shape = NA,width=0.95,alpha=0.5) +
  geom_boxplot(outlier.shape = NA,width=0.1,fill="gold") +
  facet_wrap(~type,scales = "free")+
  labs(x = "", y = "") + # Change of heritability per hotspot region (%)
  theme_classic()+scale_fill_manual(values=c("#f0650e", "#0091ff"),label = c("Disease","Trait")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none")
trans <- function(x) ifelse(x > 0.4,0.4+(x-0.4)/6, x)
inv <- function(x) ifelse(x > 0.4, 6*(x-0.4) + 0.4, x)
my_trans <- trans_new("my_trans", trans, inv)
gg1 <- gg1 + scale_y_continuous(trans = my_trans,
                                breaks = c(seq(0, 0.4, by = 0.2),seq(0.4,3.6,0.8)))+
  scale_y_break(c(0.8,1.2), scales = 1.5)
print(gg1)




