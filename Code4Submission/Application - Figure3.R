library(ggplot2)
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dendextend)
tissue.names <- read.csv("tissue_used.csv")
########################################################################
############################## trans-association #########################
########################################################################
load("enrichmenttransasso.RData")
result_order <- apply(result,2,function(x) as.numeric(cut(length(x) - rank(x) + 1,breaks=c(0,1,28))))
result_order2 <- apply(result,2,function(x) as.numeric(cut(length(x) - rank(x) + 1,breaks=c(0,1,2,28))))
result_order3 <- apply(result,2,function(x) as.numeric(cut(length(x) - rank(x) + 1,breaks=c(0,2,3,28))))

rownames(result_order) <- rownames(result_order2) <- 
  rownames(result_order3) <- rownames(result)
colnames(result_order) <- colnames(result_order2) <- 
  colnames(result_order3) <- colnames(result)

result_order <- t(result_order)
result_order2 <- t(result_order2)
result_order3 <- t(result_order3)

rownames(result_order) <-
  paste0(1:length(rownames(result_order)),": ",rownames(result_order))

row_dend = hclust(dist(result_order))
col_dend = hclust(dist(t(result_order)))

mycols <- colorRamp2(breaks = c(1, 2), 
                     colors = c("red", "white"))

h1 <- Heatmap(result_order,
              name = "order",
              row_names_gp = gpar(fontsize = 7),
              rect_gp = gpar(col = "grey"),
              cluster_rows = row_dend,
              col = mycols)
result_order <- result_order[row_order(h1),column_order(h1)]
colOrder <- order(apply(result_order,2,function(x) ifelse(length(which(x==1))==0,100,which(x==1))))

h2 <- Heatmap(result_order,
              name = "order",
              row_names_gp = gpar(fontsize = 7),
              rect_gp = gpar(col = "grey"),
              col = mycols,cluster_rows = F,cluster_columns = F,
              column_order = colOrder)

result <- result[,row_order(h1)]
result <- result[column_order(h1),][column_order(h2),]
result_order <- result_order[,rownames(result)]

################################################
result <- result[order(rownames(result)),]
result <- t(result)
result_order <- result_order[,order(colnames(result_order))]
result_order <- result_order %>% as.data.frame() %>% 
  arrange_all()

rownames(result_order) <- rownames(result_order) %>% 
  str_split(": ") %>% 
  sapply(function(x) x[2])

result <- result[rownames(result_order),]

result <- as.matrix(result)
result_order <- as.matrix(result_order)

result_tmp <- result
colnames(result_tmp) <- tissue.names$Tissue

result_order <- result_order[rownames(result),colnames(result)]
result_order2 <- result_order2[rownames(result),colnames(result)]
result_order3 <- result_order3[rownames(result),colnames(result)]

result_order_melt <- melt(result_order) %>% filter(value == 1)
result_order_melt2 <- melt(result_order2) %>% filter(value == 2)
result_order_melt3 <- melt(result_order3) %>% filter(value == 2)

result_tmp_melt <- melt(t(result_tmp))
colnames(result_tmp_melt)[1:2] <- c("Tissue","Trait/Disease")

gg1 <- result_tmp_melt %>% 
  ggplot() +
  geom_tile(aes(y=`Trait/Disease`,x=Tissue, 
                fill=value),
            col="grey")+
  scale_fill_gradient2(low="#4575B4", high="#D73027",
                       limit = c(-6,6),breaks=c(-4,0,4),
                       name="Scaled proportions")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),text = element_text(size=10)) +
  labs(y = "Trait/Disease") +
  geom_rect(data = result_order_melt2, fill = NA, col = "black",linetype="dotted",
            size = 0.4,
            aes(xmin = as.numeric(Var2) - 0.5, ymin = as.numeric(Var1) - 0.5,
                xmax = as.numeric(Var2) + 0.5, ymax = as.numeric(Var1) + 0.5)) +
  geom_rect(data = result_order_melt, fill = NA, col = "black",
            aes(xmin = as.numeric(Var2) - 0.5, ymin = as.numeric(Var1) - 0.5, 
                xmax = as.numeric(Var2) + 0.5, ymax = as.numeric(Var1) + 0.5)) +
  geom_rect(data = result_order_melt3, fill = NA, col = "black", linetype="dotted",
            size = 0.4,
            aes(xmin = as.numeric(Var2) - 0.5, ymin = as.numeric(Var1) - 0.5, 
                xmax = as.numeric(Var2) + 0.5, ymax = as.numeric(Var1) + 0.5))

gg1

trans_result_order <- result_order
trans_result_tmp <- result_tmp


########################################################################
############################## cis-association #########################
########################################################################
load("/Users/luyihao/Downloads/enrichmentcisasso.RData")
result_order <- apply(result,2,function(x) as.numeric(cut(length(x) - rank(x) + 1,breaks=c(0,1,28))))

rownames(result_order) <- rownames(result)
colnames(result_order) <- colnames(result)

result_order <- t(result_order)

rownames(result_order) <- paste0(1:length(rownames(result_order)),": ",rownames(result_order))

result <- t(result)

result_tmp <- result
rownames(result_order) <- rownames(result_order) %>% 
  str_split(": ") %>% 
  sapply(function(x) x[2])

colnames(result_order)[colnames(result_order) == "Brain Spinal cord cervical c-1"] <-
  colnames(result_tmp)[colnames(result_tmp) == "Brain Spinal cord cervical c-1"] <- 
  "Brain Spinal cord cervical c.1"
result_order <- result_order[rownames(trans_result_order),colnames(colnames(trans_result_tmp))]
colnames(result_tmp) <- tissue.names$Tissue
result_tmp <- result_tmp[rownames(trans_result_tmp),colnames(trans_result_tmp)]
result_order_melt <- melt(result_order) %>% filter(value == 1)
result_tmp_melt <- melt(t(result_tmp))
trans_result_order_melt <- melt(trans_result_order) %>% filter(value == 1)
colnames(result_tmp_melt)[1:2] <- c("Tissue","Trait/Disease")


gg1 <- result_tmp_melt  %>% ggplot() +
  geom_tile(aes(y=`Trait/Disease`,x=Tissue, fill=value),col="grey")+
  scale_fill_gradient2(low="#4575B4", high="#D73027", 
                       limit = c(-6,6),breaks=c(-4,0,4),
                       name="Scaled proportions")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size=10)) +
  labs(y = "Trait/Disease") +
  geom_rect(data = result_order_melt2, fill = NA, col = "black",linetype="dotted",
            size = 0.4,
            aes(xmin = as.numeric(Var2) - 0.5, ymin = as.numeric(Var1) - 0.5, 
                xmax = as.numeric(Var2) + 0.5, ymax = as.numeric(Var1) + 0.5)) +
  geom_rect(data = result_order_melt3, fill = NA, col = "black",linetype="dotted",
            size = 0.4,
            aes(xmin = as.numeric(Var2) - 0.5, ymin = as.numeric(Var1) - 0.5, 
                xmax = as.numeric(Var2) + 0.5, ymax = as.numeric(Var1) + 0.5)) +
  geom_rect(data = trans_result_order_melt, fill = NA, col = "black",
            size = 0.4,
            aes(xmin = as.numeric(Var2) - 0.5, ymin = as.numeric(Var1) - 0.5, 
                xmax = as.numeric(Var2) + 0.5, ymax = as.numeric(Var1) + 0.5),
            size = 0.3)

gg1
