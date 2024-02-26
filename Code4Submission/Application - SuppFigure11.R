library(tidyverse)
library(data.table)
gwas <- fread("/Users/luyihao/Downloads/gwas_catalog_v1.0-associations_e104_r2021-11-22.tsv")

gwas$`DISEASE/TRAIT` <- gwas$`DISEASE/TRAIT` %>% 
  str_replace_all(" ","_") %>% str_replace_all("/","_") %>% str_replace_all("[(]|[)]","_") %>% 
  str_replace_all("[']","")

trait.used <- read.csv("/Users/luyihao/Downloads/trait_disease_used.csv")
colnames(trait.used) <- c("","","","DISEASE","num.SNP")
gwas <- gwas %>% filter(`DISEASE/TRAIT` %in% trait.used$DISEASE)

SNP.removed <- matrix(NA,nrow = 80,ncol = 4)
for (i in 1:80) {
  print(i)
  gwas.sub <- gwas %>% filter(`DISEASE/TRAIT` == trait.used$`DISEASE`[i]) %>% 
    dplyr::select(SNPS,REGION,CHR_ID,CHR_POS,
                  `RISK ALLELE FREQUENCY`,`P-VALUE`,
                  `OR or BETA`,`95% CI (TEXT)`) %>% 
    mutate(CHR_POS = as.numeric(CHR_POS)) %>% 
    arrange(CHR_ID,CHR_POS)
  
  gwas.sub <- gwas.sub[!duplicated(gwas.sub$SNPS),]
  
  total <- clump(dat = gwas.sub, 
                 SNP_col = "SNPS", 
                 pval_col = "P-VALUE",
                 clump_kb = 1, 
                 clump_r2 = 0.9999, 
                 pop="EUR")
  
  gwas.sub.ind.0.1 <- clump(dat = gwas.sub, SNP_col = "SNPS", pval_col = "P-VALUE",clump_kb = 25,
                            clump_r2 = 0.1, pop="EUR")
  gwas.sub.ind.0.3 <- clump(dat = gwas.sub, SNP_col = "SNPS", pval_col = "P-VALUE",clump_kb = 25,
                            clump_r2 = 0.3, pop="EUR")
  gwas.sub.ind.0.5 <- clump(dat = gwas.sub, SNP_col = "SNPS", pval_col = "P-VALUE",clump_kb = 25,
                            clump_r2 = 0.5, pop="EUR")
  gwas.sub.ind.0.7 <- clump(dat = gwas.sub, SNP_col = "SNPS", pval_col = "P-VALUE",clump_kb = 25, 
                            clump_r2 = 0.7, pop="EUR")
  gwas.sub.ind.0.1 <-gwas.sub.ind.0.3 <- gwas.sub.ind.0.5 <- gwas.sub.ind.0.7
  SNP.removed[i,] <- c(nrow(total) - nrow(gwas.sub.ind.0.1),
                       nrow(total) - nrow(gwas.sub.ind.0.3),
                       nrow(total) - nrow(gwas.sub.ind.0.5),
                       nrow(total) - nrow(gwas.sub.ind.0.7))
}

colnames(SNP.removed) <- c(0.1,0.3,0.5,0.7)
SNP.removed <- SNP.removed %>% apply(2,function(x) x/trait.used$num.SNP)
SNP.removed %>% as.data.frame() %>% 
  pivot_longer(cols = `0.1`:`0.7`,
               names_to = "threshold",
               values_to = "Proportion of SNPs removed") %>% 
  mutate(threshold = factor(threshold)) %>% 
  ggplot(aes(y = `Proportion of SNPs removed`,x = threshold,fill = threshold)) + 
  geom_violin(width=0.7) +
  geom_boxplot(outlier.shape = NA,width=0.25,fill = "gold") +
  theme_clean() + labs(x = expression(r^2 ~ threshold ~ of ~ clumping)) +
  guides(fill="none")
