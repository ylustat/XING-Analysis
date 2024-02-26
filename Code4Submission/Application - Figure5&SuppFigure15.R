library(ggplot2)
library(tidyverse)
library(data.table)
load("countsspatialdegs.RData")
for (trait in c("autism","scz")) {
  diff_mat <- (result[[trait]] - result$overall)
  rownames(diff_mat) <- 
    paste0("Sample ",1:12)
  
  # quantile clip for plot
  upper_max <- tail(quantile(diff_mat,na.rm = T,1:50*0.02),2)[1]
  lower_max <- head(quantile(diff_mat,na.rm = T,1:50*0.02),2)[1]
  diff_mat[diff_mat>upper_max] <- upper_max
  diff_mat[diff_mat<lower_max] <- lower_max
  
  if(trait == "autism") {
    diff <- diff_mat[,c(2,5,6,7,1,3,4)]
    diff_mat_order <- apply(diff,2,function(x) {
      y <- x
      y[x>0] <- 1
      y[x<0] <- 0
      return(y)
    })
    diff_mat_order <- diff_mat_order %>% data.frame() %>% arrange_all(desc)
    colnames(diff_mat_order) <- str_replace(colnames(diff_mat_order),"[.]"," ")
    
    diff <- diff[rev(rownames(diff_mat_order)),] %>% melt()
  } else if (trait == "scz") {
    diff <- (diff_mat)[,c(2,5,7,1,3,4,6)]
    diff_mat_order <- apply(diff,2,function(x) {
      y <- x
      y[x>0] <- 1
      y[x<0] <- 0
      return(y)
    })
    diff_mat_order <- diff_mat_order %>% data.frame() %>% arrange_all(desc)
    colnames(diff_mat_order) <- str_replace(colnames(diff_mat_order),"[.]"," ")
    
    diff <- diff[rev(rownames(diff_mat_order)),] %>% melt() 
  }
  gg1 <- diff %>% 
    ggplot() +
    geom_tile(aes(y=Var1,x=Var2, fill=value),col="grey",size=0.5,alpha=0.8)+
    theme_classic()+
    scale_fill_gradient2(high = "red",low = "white",
                         name="Difference") +
    labs(x = "Layer", y = "Sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),text = element_text(size=10)) +
    theme(legend.position = "none",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  print(gg1)
  
}