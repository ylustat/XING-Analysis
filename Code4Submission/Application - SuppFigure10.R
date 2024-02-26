library(tidyverse)
library(ggplot2)
library(CCA)
library(ggtext)
library(glue)
library(ggpattern)
setwd("/Users/yihaolu/Downloads/Research/XING_revision/")
############################################################
############# Expression Deconvolution #####################
############################################################
load("Data/data_for_cell_type.RData")
Z1 <- data_for_cell_type$Z1; Z2 <- data_for_cell_type$Z2
X_full_dat1 <- data_for_cell_type$X1; X_full_dat2 <- data_for_cell_type$X2
tissues <- str_remove_all(colnames(Z1),"eQTL_")
colnames(X_full_dat1) <- tissues
load("Data/merge_mean_cell_type.RData")
res$res_exp <- res$res_exp %>% data.frame() %>% dplyr::select(-c("P.value","Correlation","RMSE"))
colnames(res$res_exp) <- c("Naive B Cells",
                           "Memory B Cells",
                           "Plasma Cells",
                           "CD8 T Cells",
                           "Naive CD4 T Cells",
                           "Resting Memory CD4 T Cells",
                           "Activated Memory CD4 T Cells",
                           "Follicular Helper T Cells",
                           "Regulatory T Cells (Tregs)",
                           "Gamma Delta T Cells",
                           "Resting NK Cells",
                           "Activated NK Cells",
                           "Monocytes",
                           "M0 Macrophages",
                           "M1 Macrophages",
                           "M2 Macrophages",
                           "Resting Dendritic Cells",
                           "Activated Dendritic Cells",
                           "Resting Mast Cells",
                           "Activated Mast Cells",
                           "Eosinophils",
                           "Neutrophils")

CC <- 4
sd1 <- apply(X_full_dat1, 2, sd)
sd2 <- apply(X_full_dat2, 2, sd)
mean1 <- colMeans(X_full_dat1)
mean2 <- colMeans(X_full_dat2)
ccas <- cc(X_full_dat1, X_full_dat2)
Ahat <- ccas$xcoef
Bhat <- ccas$ycoef
U <- X_full_dat1 %*% Ahat
V <- X_full_dat2 %*% Bhat
X_est1 <- U[,1:CC] %*% ginv(Ahat[,1:CC])
X_est2 <- V[,1:CC] %*% ginv(Bhat[,1:CC])

X_est1 <- X_est1 %*% diag(sd1) + t(mean1 %*% matrix(1,nrow = 1, ncol = nrow(X_full_dat1)))
X_est2 <- X_est2 %*% diag(sd2) + t(mean2 %*% matrix(1,nrow = 1, ncol = nrow(X_full_dat2)))

X_res1 <- X_full_dat1 - X_est1   # Residual for data 1 after the CCA
X_res2 <- X_full_dat2 - X_est2   # Residual for data 2 after the CCA

pca_obj <- prcomp(scale(X_res1))
a=cor(pca_obj$rotation,res$res_exp[tissues,])[1,]
ind <- which(colnames(res$res_exp) %in% c("B.cells.naive","T.cells.CD8","T.cells.CD4.naive","NK.cells.resting","Monocytes","Dendritic.cells.resting","Eosinophils","Neutrophils","Plasma.cells"))
y = pca_obj$rotation[,1:3]
x = data.matrix(res$res_exp[tissues,])
# Calculate p-value
pval <- apply(y, 2, function(y1) apply(x,2,function(x1) cor.test(y1,x1)$p.value)) %>% apply(1,min)

# Create a function to annotate symbols based on p-value
get_pval_annotation <- function(pval){
  if (pval < 0.001){
    return('***')
  } else if (pval < 0.01){
    return('**')
  } else if (pval < 0.05){
    return('*')
  } else if (pval < 0.1){
    return('.')
  } else {
    return('')
  }
}

# Calculate correlation
data <- apply(y, 2, function(y1) apply(x,2,function(x1) cor(y1,x1))) %>% 
  apply(1,function(x) c(max(abs(x)),which.max(abs(x)))) %>% 
  t() %>% 
  data.frame() %>% 
  mutate(pval_annotation = sapply(pval, get_pval_annotation)) %>%  # Get p-value annotation
  rename(adjR2 = "X1",whichrank = "X2") %>% 
  mutate(cells = rownames(.)) %>% 
  arrange(whichrank, desc(adjR2)) %>% 
  mutate(cells = factor(cells, levels = unique(cells))) %>% 
  mutate(adjR2 = ifelse(adjR2 < 0, 0, adjR2)) %>% 
  filter(pval_annotation!="")
data_exp <- data
# Generate bar plot
highlight = function(x, pat, color="black", family="") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

ggplot(data, aes(x = cells, y = adjR2, fill = factor(whichrank)))+
  geom_col(width = 0.5, alpha = 0.8, col = "black")+
  geom_text(aes(label=pval_annotation), vjust=-0.25,size=5) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Cell", y = "|Correlation|", fill="")+
  scale_fill_manual(values = c("red2","gold","steelblue"),labels = c("PC"["exp"]~1,"PC"["exp"]~2,"PC"["exp"]~3))+
  expand_limits(y = max(data$adjR2) * 1.1)+
  scale_x_discrete(labels= function(x) highlight(x, "Macrophages|Neutrophils|Memory CD4 T Cells|Dendritic Cells|Monocytes", "purple4")) +
  theme(axis.text.x=element_markdown(),
        text = element_text(size = 14),
        plot.margin = margin(l = 20))


############################################################
############# DNA methylation Deconvolution ################
############################################################
tissues <-  c("Ovary","Breast_Mammary_Tissue","Prostate","Testis","Kidney_Cortex","Lung","Whole_Blood","Colon_Transverse","Muscle_Skeletal")
colnames(X_res2) <- tissues
load("Data/merge_mean_cell_type.RData")
colnames(res$res_dnam) <- c("Epithelial Cells",
                            "Fibroblasts",
                            "Immature Cells",
                            "Monocytes",
                            "Dendritic Cells",
                            "Macrophages",
                            "Neutrophils",
                            "Eosinophils",
                            "Regulatory T Cells (Tregs)",
                            "Naive T Cells",
                            "Memory T Cells",
                            "CD8 T Cells",
                            "Natural Killer Cells",
                            "B Cells")


res$res_exp <- res$res_exp %>% data.frame() %>% dplyr::select(-c("P.value","Correlation","RMSE"))

pca_obj <- prcomp(scale(X_res2))
ind <- which(colnames(res$res_exp) %in% c("B.cells.naive","T.cells.CD8","T.cells.CD4.naive","NK.cells.resting","Monocytes","Dendritic.cells.resting","Eosinophils","Neutrophils","Plasma.cells"))

y = pca_obj$rotation[,1:3]
x = data.matrix(res$res_dnam[tissues,])
# Calculate p-value
pval <- apply(y, 2, function(y1) apply(x,2,function(x1) cor.test(y1,x1)$p.value)) %>% apply(1,min)

# Create a function to annotate symbols based on p-value
get_pval_annotation <- function(pval){
  if (pval < 0.001){
    return('***')
  } else if (pval < 0.01){
    return('**')
  } else if (pval < 0.05){
    return('*')
  } else if (pval < 0.1){
    return('.')
  } else {
    return('')
  }
}

significance_labels <- data.frame(
  Significance = c("***", "**", "*", "."),
  Description = c("p < 0.001", "p < 0.01", "p < 0.05", "p < 0.1")
)

# Calculate correlation
data <- apply(y, 2, function(y1) apply(x,2,function(x1) cor(y1,x1))) %>% 
  apply(1,function(x) c(max(abs(x)),which.max(abs(x)))) %>% 
  t() %>% 
  data.frame() %>% 
  mutate(pval_annotation = sapply(pval, get_pval_annotation)) %>%  # Get p-value annotation
  rename(adjR2 = "X1",whichrank = "X2") %>% 
  mutate(cells = rownames(.)) %>% 
  arrange(whichrank, desc(adjR2)) %>% 
  mutate(cells = factor(cells, levels = unique(cells))) %>% 
  mutate(adjR2 = ifelse(adjR2 < 0, 0, adjR2)) %>% 
  filter(pval_annotation!="")
# Generate bar plot
highlight = function(x, pat, color="black", family="") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

ggplot(data, aes(x = cells, y = adjR2, fill = factor(whichrank)))+
  geom_col(width = 0.5, alpha = 0.8, col = "black")+
  geom_col_pattern(width = 0.5, alpha = 0.8, col = "black",
                   pattern_type = "stripe",
                   pattern_colour = "black",
                   pattern_size = 0.1,
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.02,
                   pattern_key_scale_factor = 1) +
  geom_text(aes(label=pval_annotation), vjust=-0.25,size=5) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Cell", y = "|Correlation|", fill="")+
  scale_fill_manual(values = c("red2","gold","steelblue"),labels = c("PC"["DNAm"]~1,"PC"["DNAm"]~2,"PC"["DNAm"]~3))+
  expand_limits(y = max(data$adjR2) * 1.1)+
  scale_y_continuous(breaks = c(0:5)*0.2)+
  scale_x_discrete(labels= function(x) highlight(x, "Macrophages|Neutrophils|Memory T Cells|Dendritic Cells|Monocytes", "purple4")) +
  theme(axis.text.x=element_markdown(),
        text = element_text(size = 13),
        plot.margin = margin(l = 20))




