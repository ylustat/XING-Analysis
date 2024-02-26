#==================================================================#
#==================== codes to reproduce figures ==================#
#==================================================================#
library(ggplot2)
library(tidyverse)
library(data.table)
df <- fread("ResultFigureS9a.txt")
df <- rbind(df, 
            df %>% filter(cc == 20) %>% mutate(cc = 9.5) %>% mutate(X1 = X1),
            df %>% filter(cc == 6) %>% mutate(cc = 9.5) %>% mutate(X1 = X1 - 0.03))
df %>% 
  mutate(h2 = factor(h2, levels = c(0.2,0.15,0.1))) %>% 
  ggplot(aes(x = cc, y = X1, col = h2, fill = h2))+
  geom_line(size=1.5,alpha=0.5)+geom_point()+
  geom_point(size=2,pch=21,col="black",stroke=0.5, alpha=0.9)+
  scale_fill_manual(values = c("red","pink","orange"), labels = c(0.2,0.15,0.1))+
  scale_color_manual(values = c("red","pink","orange"), labels = c(0.2,0.15,0.1))+
  theme_classic() +
  labs(x = "# CC", y = "AUC", col = expression(theta['\u2113']), fill = expression(theta['\u2113']))+
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(color = "black", linewidth = 0.5),
    legend.key = element_blank(),
    legend.key.size = unit(0.4, 'cm')
  )+
  scale_y_continuous(limits = c(0.6,0.92))+
  theme(axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank())+
  ggbreak::scale_x_break(breaks=c(6.2,14),
                         scales = 0.7,
                         ticklabels=seq(20,40,by=10),
                         space=1)+
  scale_x_continuous(breaks = 0:10)+
  theme(axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_blank())



df <- fread("ResultFigureS9b.txt")
df <- rbind(df, 
            df %>% filter(cc == 20) %>% mutate(cc = 9.5) %>% mutate(X1 = X1) %>% mutate(X2 = X2),
            df %>% filter(cc == 6) %>% mutate(cc = 9.51) %>% mutate(X1 = X1) %>% mutate(X2 = X2))

df %>% 
  ggplot(aes(x = cc, y = X1, col = factor(h2), fill = factor(h2)))+
  geom_line(size=1.5,alpha=0.5)+geom_point()+
  geom_point(size=2,pch=21,col="black",stroke=0.5, alpha=0.9)+
  scale_fill_manual(values = c("orange","pink","red"), labels = c(0.1,0.15,0.2))+
  scale_color_manual(values = c("orange","pink","red"), labels = c(0.1,0.15,0.2))+
  theme_classic() +
  labs(x = "# PC", y = "AUC", col = expression(theta['\u2113']), fill = expression(theta['\u2113']))+
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(color = "black", linewidth = 0.5),
    legend.key = element_blank(),
    legend.key.size = unit(0.4, 'cm')
  )+
  scale_y_continuous(limits = c(0.6,0.92))+
  theme(axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank())+
  ggbreak::scale_x_break(breaks=c(6.2,16),
                         scales = 0.7,
                         ticklabels=seq(20,40,by=10),
                         space=1)+
  scale_x_continuous(breaks = 0:10)+
  theme(axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_blank())

#==================================================================#
#=============== data generation and simulations ==================#
#==================================================================#
library(ggplot2)
library(tidyverse)
library(plotROC)
library(RColorBrewer)
library(patchwork)
library(CCA)
library(MASS)
library(pROC)
library(pbapply)
library(ggbreak)
conflicted::conflicts_prefer(mvtnorm::rmvnorm)
setwd("/Users/yihaolu/Downloads/Research/XING_revision/Scripts/")
source('NG_simulation_formal_function_cri.R')
source('Simulation_XING.R')
source('generate_correlated_binary.R')
load("params.RData")
args <- c(400,400,10,0.5,0.5,0.5,0.5,
          0.5,0.96,0,
          20,20,20,20,0,
          0.5,0,1,2,2,
          2,0.1,0,0)
n1 = as.numeric(args[1]) # sample size of data1
n2 = as.numeric(args[2]) # sample size of data2
g = as.numeric(args[3]) # number of genes (M = g * 50)
r11 = as.numeric(args[4]) # correlation within data 1
r12 = as.numeric(args[5]) # correlation within data 2
r21 = as.numeric(args[6]) # correlation 1 across data 1 and data 2
r22 = as.numeric(args[7]) # correlation 2 across data 1 and data 2
h = sqrt(as.numeric(args[8])) # input heritability (h2)
threshold = as.numeric(args[9]) # logit threshold
ld = as.numeric(args[10]) # LD correlation among SNPs
k11 = as.numeric(args[11])
k12 = as.numeric(args[12])
k21 = as.numeric(args[13])
k22 = as.numeric(args[14])
sigmaerror = as.numeric(args[15])
h_second = sqrt(as.numeric(args[16]))
prop_uncorrelated = as.numeric(args[17])
structured = as.numeric(args[18]) # 1: true, 0: false
cc0 = as.numeric(args[19])
pc1 = as.numeric(args[20])
pc2 = as.numeric(args[21])
pi_init = as.numeric(args[22])
proportion_ld = as.numeric(args[23])
ld_second = as.numeric(args[24])
corr_pattern = "exch"
#==================================================================#
# parameter setting
## dimension parameter
M <- numeric() # number of marginal tests
N <- 2 # number of studies
K <- c(k11+k12,k21+k22) # number of conditions for each study


#==================================================================#
# data generation
repind <- 1234
df <- list()
for (repind in 1:3) {
  print(repind)
  res <- list()
  cnt <- 1
  for (h in sqrt(seq(0.2,0.1,-0.05))) {
    res[[cnt]] <- pblapply(seq(0.3,0.5,0.1)[3],function(r11) {
      r12 <- r11
      r21 <- r11
      ## sample size parameter
      n <- c(n1, n2) # sample size of genotype
      l <- c(10, 10) # number of blocks for the cis-SNPs of each gene
      k <- c(50, 50) # number of SNPs in each block
      g <- c(g, g) # number of genes
      M <- unique(g * l * k) # number of marginal tests in each condition
      
      ## correlation parameter
      r <- ld # indepedent: 0, moderate LD: 0.4, high LD: 0.8
      r_in_data <- c(r11, r12) # co-regulation correlation across tissues
      r_cross_data <- r21 # co-regulation correlation across omics
      R <- cor_mat_gen(r1 = rep(r_in_data,2), r2 = r_cross_data, cut = c(k11,k21),K=K)
      ## heritability
      # h = list(data1 = rep(c(h,h_second),times = c(k11,k12)),
      #          data2 = rep(c(h,h_second),times = c(k21,k22)))
      
      ## generate latent status
      set.seed(repind*round(sum(as.numeric(args))))
      gamma <- rmvbin(M, margprob = rep(1-threshold,nrow(R)), bincorr = R,simulvals = params) 
      gamma <- list(gamma[,1:K[1]],gamma[,(K[1]+1):sum(K)])
      
      ## generate the phenotype
      phenotype <- genotype <- beta <- lapply(g, function(x) rep(list(NA),x))
      for (ni in 1:N) {
        cis_size <- l[ni] * k[ni]
        for (i in 1:g[ni]) {
          # set.seed(i)
          # genotype
          genotype[[ni]][[i]] <- genotype_gen(n[ni], block = l[ni], block_size = k[ni], r)
          genotype[[ni]][[i]] <- scale(genotype[[ni]][[i]])
          # coefficient
          num_signal <- sqrt(mean(colSums(gamma[[ni]][((i-1) * cis_size + 1):(i * cis_size),])))
          if(structured) {
            dimindex <- c(0,cumsum(K))
            R_sub = R[(dimindex[ni]+1):(dimindex[ni+1]),(dimindex[ni]+1):(dimindex[ni+1])]
            beta[[ni]][[i]] <- structured_beta_gen(cis_size,K[ni],rep(h,K[ni])/num_signal,R_sub) # self-defined heritability
          } else {
            beta[[ni]][[i]] <- beta_gen(cis_size,K[ni],h/num_signal) # self-defined heritability
          }
          
          # error term
          epsilon <- matrix(rnorm(n[ni] * K[ni],0,sqrt(1-h^2)),nrow = n[ni])
          # generate phenotype
          phenotype[[ni]][[i]] <- genotype[[ni]][[i]] %*%
            (beta[[ni]][[i]] * gamma[[ni]][((i-1) * cis_size + 1):(i * cis_size),])+
            epsilon
        }
      }
      
      ## generate the summary statistics
      Z <- list(matrix(NA,nrow=M,ncol=K[1]),
                matrix(NA,nrow=M,ncol=K[2]))
      beta_est <- se_est <- pval_est <-
        list(matrix(NA, nrow = M, ncol = K[1]),
             matrix(NA, nrow = M, ncol = K[2]))
      
      # Loop through each sample (ni) in the dataset (1 to N)
      for (ni in 1:N) {
        cis_size <- l[ni] * k[ni]
        
        # Apply the function in parallel for each gene (gi)
        result <- pbapply::pblapply(1:g[ni], function(gi) {
          res_gi <- list(Z = matrix(NA, nrow = cis_size, ncol = K[ni]),
                         beta_est = matrix(NA, nrow = cis_size, ncol = K[ni]),
                         se_est = matrix(NA, nrow = cis_size, ncol = K[ni]),
                         pval_est = matrix(NA, nrow = cis_size, ncol = K[ni]))
          
          for (Mi in 1:cis_size) { # snp
            for (ti in 1:K[ni]) { # tissue
              x = genotype[[ni]][[gi]][, Mi]
              y = phenotype[[ni]][[gi]][, ti]
              model <- simplelm(x = x, y = y)
              res_gi$Z[Mi, ti] <- model[3]
              res_gi$beta_est[Mi, ti] <- model[1]
              res_gi$se_est[Mi, ti] <- model[2]
              res_gi$pval_est[Mi, ti] <- model[4]
            }
          }
          return(res_gi)
        })  # Use 20 cores for parallel computation
        
        # Combine the results from the parallel computation
        for (gi in 1:g[ni]) {
          idx <- (gi - 1) * cis_size
          Z[[ni]][idx + 1:(cis_size), ] <- result[[gi]]$Z
          beta_est[[ni]][idx + 1:(cis_size), ] <- result[[gi]]$beta_est
          se_est[[ni]][idx + 1:(cis_size), ] <- result[[gi]]$se_est
          pval_est[[ni]][idx + 1:(cis_size), ] <- result[[gi]]$pval_est
        }
      }
      
      
      ## get true beta
      beta <- lapply(beta, function(x) do.call("rbind",x))
      true_beta <- lapply(1:N, function(x) {
        beta[[x]] * gamma[[x]]
      })
      
      # shuffle
      if(prop_uncorrelated > 0){
        shuffle <- round(M * prop_uncorrelated)
        shuffle_index <- lapply(1:ncol(Z[[1]]),function(x) sample(1:shuffle,shuffle))
        for (i in 1:ncol(Z[[1]])) {
          Z[[1]][1:shuffle,i] <- Z[[1]][shuffle_index[[i]],i]
          gamma[[1]][1:shuffle,i] <- gamma[[1]][shuffle_index[[i]],i]
        }
      }
      
      cov_data1 <- solve(cov(Z[[1]]))
      cov_data2 <- solve(cov(Z[[2]]))
      Alg1_data1 <- Alg1(betahat = Z[[1]],
                         Lambda = cov_data1,
                         iterT = 50,
                         bound=1e-8, pi_init=pi_init,vk_init=0.1)
      Alg1_data2 <- Alg1(betahat = Z[[2]],
                         Lambda = cov_data2,
                         iterT = 50,
                         bound=1e-8, pi_init=pi_init,vk_init=0.1)
      Alg2_data1 <- Alg2_ind1(betahat = Z[[1]],
                              Lambda = cov_data1,
                              PC = pc1,
                              results_alg1 = Alg1_data1,
                              eps_thresh=1e-2,
                              iterT=50)
      Alg2_data2 <- Alg2_ind1(betahat = Z[[2]],
                              Lambda = cov_data2,
                              PC = pc2,
                              results_alg1 = Alg1_data2,
                              eps_thresh=1e-2,
                              iterT=50)
      performance_PC <- pbapply::pblapply(c(0,1,2,3,4,5,6,20,30,40), function(pc1) {
        pc2 <- pc1
        Alg4_mQTL_eQTL <- Alg4(betahat1 = Z[[1]],
                               betahat2 = Z[[2]],
                               Lambda1 = cov_data1,
                               Lambda2 = cov_data2,
                               CC = cc0,
                               PC1 = pc1,
                               PC2 = pc2,
                               results_alg1_dat1=Alg1_data1,
                               results_alg1_dat2=Alg1_data2,
                               results_alg2_dat1 = Alg2_data1,
                               results_alg2_dat2 = Alg2_data2,
                               iterT = 30,
                               bound=1e-4)
        return(c(auc(c(gamma[[1]]), c(Alg4_mQTL_eQTL$alphajkm2_dat1),quiet=T),
                 auc(c(gamma[[2]]), c(Alg4_mQTL_eQTL$alphajkm2_dat2),quiet=T)))
      },cl=5)

      return(performance_PC)
    })
    cnt <- cnt+1
  }
  df[[repind]] <- res %>% 
    lapply(function(x) lapply(x,function(y) do.call("rbind",y) %>% data.frame())) %>% 
    lapply(function(y) {
      mapply(function(x,ind) data.frame(x,r = ind), y, seq(0.3,0.5,0.1)[3],SIMPLIFY = F)
    }) %>% 
    lapply(function(x) do.call("rbind",x)) %>% 
    mapply(function(x,y) data.frame(x,h2=y), ., seq(0.2,0.1,-0.05),SIMPLIFY = F) %>% 
    do.call("rbind",.)
}
df <- lapply(df, function(x) {
  x$cc <- c(0,1,2,3,4,5,6,20,30,40)
  return(x)
})
df <- Reduce("+",df)/length(df)





