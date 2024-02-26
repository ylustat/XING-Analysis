#==================================================================#
#==================== codes to reproduce figures ==================#
#==================================================================#
library(ggplot2)
library(tidyverse)
library(data.table)
res <- fread("/Users/luyihao/Downloads/XING_revision_r3/Data4Submission/ResultFigureS8.txt")
res %>%
  ggplot(aes(x = ld, y = ., col = factor(theta), fill = factor(theta)))+
  geom_line(size = 2)+geom_point()+
  geom_point(size=6,pch=21,col="black",stroke=0.5)+
  scale_fill_manual(values = brewer.pal(4,"Set2"), labels = c(0.01,0.05,0.15,0.25))+
  scale_color_manual(values = brewer.pal(4,"Set2"), labels = c(0.01,0.05,0.15,0.25))+
  theme_classic() +
  labs(x = "Pairwise correlation among SNPs within the same block", y = "AUC", col = expression(theta["\u2113"]), fill = expression(theta["\u2113"]))+
  theme(
    text = element_text(size = 25),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(color = "black", size = 0.5),
    legend.key = element_blank()
  )+
  scale_y_continuous(limits = c(0.5,0.85))

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
setwd("/Users/yihaolu/Downloads/Research/XING_revision/Scripts/")
source('NG_simulation_formal_function_cri.R')
source('Simulation_XING.R')
args <- c(300,300,60,0.3,0.3,
          0.3,0.3,0.5,0.9,0,
          7,3,7,3,0,
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
res <- list()
cnt <- 1
for (h in sqrt(c(seq(0.25,0.05,-0.1),0.01))) {
  res[[cnt]] <- pblapply(seq(0,0.4,0.05),function(ld) {
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
    h = list(data1 = rep(c(h,h_second),times = c(k11,k12)),
             data2 = rep(c(h,h_second),times = c(k21,k22)))
    ## generate latent status
    set.seed(repind*round(sum(as.numeric(args))))
    X <- latent_X_gen(M,K,R)
    gamma <- latent_gamma_gen(X,threshold = threshold)
    phenotype <- rep(list(rep(list(NA),g[1])),N)
    genotype <- rep(list(rep(list(NA),g[1])),N)
    beta <- rep(list(rep(list(NA),g[1])),N)
    for (ni in 1:N) {
      cis_size <- l[ni] * k[ni]
      for (i in 1:g[ni]) {
        set.seed(i)
        # genotype
        if(proportion_ld > 0) {
          l1 <- round(l[ni] * (1 - proportion_ld))
          l2 <- l[ni] - l1
          geno1 <- genotype_gen(n[ni], block = l1, block_size = k[ni], r, corr_pattern = corr_pattern)
          geno2 <- genotype_gen(n[ni], block = l2, block_size = k[ni], ld_second, corr_pattern = corr_pattern)
          genotype[[ni]][[i]] <- cbind(geno1,geno2)
        } else {
          genotype[[ni]][[i]] <- genotype_gen(n[ni], block = l[ni], block_size = k[ni], r, corr_pattern = corr_pattern)
        }
        genotype[[ni]][[i]] <- scale(genotype[[ni]][[i]])
        # coefficient
        num_signal <- sqrt(mean(colSums(gamma[[ni]][((i-1) * cis_size + 1):(i * cis_size),])))
        if(structured) {
          dimindex <- c(0,cumsum(K))
          beta[[ni]][[i]] <- structured_beta_gen(cis_size,K[ni],h[[ni]]/num_signal,
                                                 R[(dimindex[ni]+1):(dimindex[ni+1]),(dimindex[ni]+1):(dimindex[ni+1])]) # self-defined heritability
        } else {
          beta[[ni]][[i]] <- beta_gen(cis_size,K[ni],h[[ni]]/num_signal) # self-defined heritability
        }
        # error term
        epsilon <- mvtnorm::rmvnorm(n[ni],
                                    mean = rep(0,K[ni]),
                                    sigma = diag(1-h[[ni]]^2))
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
                            iterT=20)
    Alg2_data2 <- Alg2_ind1(betahat = Z[[2]],
                            Lambda = cov_data2,
                            PC = pc2,
                            results_alg1 = Alg1_data2,
                            eps_thresh=1e-2,
                            iterT=20)
    performance_PC <- pbapply::pblapply(2, function(pc1) {
      print("Start Alg4 ...")
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
                             iterT = 20,
                             bound=1e-4,
                             sparse = F)
      return(c(auc(c(gamma[[1]]), c(Alg4_mQTL_eQTL$alphajkm2_dat1),quiet=T)))
    })
    return(performance_PC)
  },cl=5)
  cnt <- cnt+1
}

res <- mapply(function(y,ind) {
  do.call("rbind",lapply(y,function(x) x[[1]])) %>%
    data.frame(ld = seq(0,0.4,0.05),theta=ind)
},res,c(seq(0.25,0.05,-0.1),0.01),SIMPLIFY = F) %>%
  do.call("rbind",.) 