#==================================================================#
#==================== codes to reproduce figures ==================#
#==================================================================#
library(ggplot2)
library(tidyverse)
library(data.table)
setwd("/Users/luyihao/Downloads/XING_revision_r3/Data4Submission")
# Figure S1
df <- fread("ResultFigureS1.txt")
df %>%
  pivot_longer(cols = ALG4:metasoft,names_to = "method",values_to = "probability") %>%
  mutate(method = factor(method, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  mutate(method = factor(method),
         line_type=factor(as.integer(method)%%2)) %>%
  ggplot() +
  geom_roc(aes(d = as.integer(truevalue), m = probability,
               color = method),size = 0.4, labels = F,n.cuts = 10) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") + theme_classic()+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  xlim(0,0.2)+
  theme(legend.text.align = 0)+
  scale_y_continuous(breaks = seq(0,1, by = 0.2))

# Figure S2
df <- fread("ResultFigureS2.txt")
df %>%
  pivot_longer(cols = ALG4:metasoft,names_to = "method",values_to = "probability") %>%
  mutate(method = factor(method, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  mutate(method = factor(method),
         line_type=factor(as.integer(method)%%2)) %>%
  ggplot() +
  geom_roc(aes(d = as.integer(truevalue), m = probability,
               color = method),size = 0.4, labels = F,n.cuts = 10) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") + theme_classic()+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  xlim(0,0.2)+
  theme(legend.text.align = 0)+
  scale_y_continuous(breaks = seq(0,1, by = 0.2), limits = c(0,0.9))

# Figure S3
load("ResultFigureS3.txt")
df <- fread("ResultFigureS3.txt")
df %>%
  filter(Data == "Data1") %>%
  mutate(k2=k21+k22) %>%
  dplyr::select(k2, `X-ING`:mash) %>%
  pivot_longer(cols = `X-ING`:mash,names_to = "Method") %>%
  mutate(Method = factor(Method, levels = c("X-ING","mash")),
         line_type=factor(as.integer(Method)%%2)) %>%
  ggplot(aes(x = k2, y = value, group = Method,col=Method,shape=Method,linetype=Method))+
  geom_point()+geom_line()+
  theme_classic()+
  labs(x = "Dimension of data A",y = "AUC of data B")+
  scale_color_manual(values=c("red","purple"),
                     name="Method")

df %>%
  filter(Data == "Data2") %>%
  mutate(k2=k21+k22) %>%
  dplyr::select(k2, `X-ING`:mash) %>%
  pivot_longer(cols = `X-ING`:mash,names_to = "Method") %>%
  mutate(Method = factor(Method, levels = c("X-ING","mash")),
         line_type=factor(as.integer(Method)%%2)) %>%
  ggplot(aes(x = k2, y = value, group = Method,col=Method,shape=Method,linetype=Method))+
  geom_point()+geom_line()+
  theme_classic()+
  labs(x = "Dimension of data A",y = "AUC of data A")+
  scale_color_manual(values=c("red","purple"),
                     name="Method")

# Figure S4a
df <- fread("ResultFigureS4a.txt")
df %>%
  pivot_longer(ALG1:hteqtl) %>% 
  mutate(name = factor(name, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  ggplot(aes(x = h2, y = value, col = name))+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x = "Percentage of variance explained by\nsimulated predictors (unstructured effects)",y = "AUC")+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  theme(legend.text.align = 0)+
  scale_y_continuous(limits = c(0.67,0.93))


# Figure S4b
df <- fread("ResultFigureS4b.txt")
df %>%
  pivot_longer(ALG1:hteqtl) %>% 
  mutate(name = factor(name, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  ggplot(aes(x = h2, y = value, col = name))+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x = "Percentage of variance explained by\nsimulated predictors (structured effects)",y = "AUC")+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  theme(legend.text.align = 0)+
  scale_y_continuous(limits = c(0.67,0.93))


# Figure S4c
df <- fread("ResultFigureS4c.txt")
df %>% 
  pivot_longer(ALG1:hteqtl) %>% 
  mutate(name = factor(name, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  ggplot(aes(x = h_second, y = value, col = name))+
  geom_point()+
  geom_line()+
  theme_classic()+
  labs(x = "Percentage of variance explained by simulated\npredictors for the second group of contexts",y = "AUC")+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  theme(legend.text.align = 0)+
  scale_y_continuous(limits = c(0.7,0.9))+
  scale_x_continuous(breaks = seq(0.05,0.3,0.05))

# Figure S5
df <- fread("ResultFigureS5.txt")
df %>% 
  pivot_longer(cols = `X-ING`:mash, names_to = "method") %>% 
  mutate(method=factor(method,levels = c("X-ING","mash"))) %>% 
  ggplot() + 
  geom_col(aes(x = factor(R2), y = value, group = method, fill = method),
           position = "dodge",alpha=0.5,col="black",width=0.6) + # ,width=0.06
  theme_classic() +
  scale_fill_manual(values=c("red","purple"),name="Method") + 
  labs(x = "Sample size",y = "RMSE of truly non-null effects")



# Figure S6
df <- fread("ResultFigureS6.txt")
df %>%
  pivot_longer(ALG1:hteqtl,names_to = "method",values_to = "time") %>% 
  mutate(method = factor(method, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>% 
  ggplot(aes(x = g, y = time, col = method))+
  geom_line(size=1.5,alpha=0.7)+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  theme_classic()+
  theme(legend.text.align = 0)+
  labs(x = "Number of tested units", y = "Computation Time (second)")

# Figure S7
df <- fread("ResultFigureS7.txt")
df %>% 
  pivot_longer(ALG1:hteqtl,names_to = "method",values_to = "time") %>% 
  mutate(method = factor(method, levels=c("ALG4","ALG2","mash","hteqtl","ALG1","metasoft"))) %>%
  ggplot(aes(x = k11+3, y = time*10, col = method))+
  geom_line(size=1.5,alpha=0.7)+
  scale_color_manual(values=c("red","orange","purple","skyblue","pink","black"),
                     labels = c(expression('X-ING'),expression('X-ING'[single-omics]),
                                'mash','MT-eQTL',expression('X-ING'[starting]),'Metasoft'),
                     name="Method")+
  theme_classic()+
  theme(legend.text.align = 0)+
  coord_cartesian(ylim = c(0,5000))+
  labs(x = "Number of contexts", y = "Computation Time (second)")
#==================================================================#
#=============== data generation and simulations ==================#
#==================================================================#
args <- commandArgs(trailingOnly = T)
n1 = as.numeric(args[1]) # sample size of data1
n2 = as.numeric(args[2]) # sample size of data2
g = as.numeric(args[3]) # number of genes (M = g * 50)
r11 = as.numeric(args[4]) # correlation within data 1
r12 = as.numeric(args[5]) # correlation within data 2
r21 = as.numeric(args[6]) # correlation 1 across data 1 and data 2
prop_uncorrelated = as.numeric(args[7]) # correlation 2 across data 1 and data 2
h = sqrt(as.numeric(args[8])) # input heritability (h2)
threshold = as.numeric(args[9]) # logit threshold
ld = as.numeric(args[10]) # LD correlation among SNPs
k11 = as.numeric(args[11])
k12 = as.numeric(args[12])
k21 = as.numeric(args[13])
k22 = as.numeric(args[14])
repind = as.numeric(args[15])

h_second = h
structured = 0
cc0 = 2
pc1 = 2
pc2 = 2
pi_init = 0.1
proportion_ld = 0
ld_second = 0

print(Sys.time())
library(ggplot2)
library(tidyverse)
library(plotROC)
library(RColorBrewer)
library(patchwork)
setwd("/scratch/t.phs.yihaolu/XING_NG/simulation")
source('Revision_NG_simulation_formal_function_cri.R')
source('generate_correlated_binary.R')
load("params.RData")
#==================================================================#
# parameter setting
## dimension parameter
M <- numeric() # number of marginal tests
N <- 2 # number of studies
K <- c(k11+k12,k21+k22) # number of conditions for each study

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
#==================================================================#
# data generation
#==================================================================#
auc_list <- list()
cat("This is the ", repind, " time\n")
token <- paste0(c(args,repind),collapse = "_")
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
  list(matrix(NA,nrow=M,ncol=K[1]),
       matrix(NA,nrow=M,ncol=K[2]))
for (ni in 1:N) {
  cis_size <- l[ni] * k[ni]
  for (gi in 1:g[ni]) { # gene
    print(gi)
    for (Mi in 1:cis_size) { # snp
      for (ti in 1:K[ni]) { # tissue
        x = genotype[[ni]][[gi]][,Mi]
        y = phenotype[[ni]][[gi]][,ti]
        model <- simplelm(x = x, y = y)
        Z[[ni]][(gi-1)*cis_size+Mi,ti] <- model[3]
        beta_est[[ni]][(gi-1)*cis_size+Mi,ti] <- model[1]
        se_est[[ni]][(gi-1)*cis_size+Mi,ti] <- model[2]
        pval_est[[ni]][(gi-1)*cis_size+Mi,ti] <- model[4]
      }
    }
  }
}


# get true beta
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

## run our algorithm
library(MASS)
library(pROC)
library(ashr)
library(mashr)
source("Revision_CCA.R")
source('Revision_AlgorithmsHeir1.R')
source('Revision_HierarchicalLLRAlg3-4_1.R')
print("Start Alg1 ... data 1")
cov_data1 <- solve(cov(Z[[1]]))
cov_data2 <- solve(cov(Z[[2]]))

Alg1_data1 <- Alg1(betahat = Z[[1]],
                   Lambda = cov_data1,
                   iterT = 10,
                   bound=1e-8, pi_init=pi_init,vk_init=0.1)

Alg1_data2 <- Alg1(betahat = Z[[2]],
                   Lambda = cov_data2,
                   iterT = 10,
                   bound=1e-8, pi_init=pi_init,vk_init=0.1)

Alg2_data1 <- Alg2_ind1(betahat = Z[[1]],
                        Lambda = cov_data1,
                        PC = cc0+pc1,
                        results_alg1 = Alg1_data1,
                        eps_thresh=1e-2,
                        iterT=5)
Alg2_data2 <- Alg2_ind1(betahat = Z[[2]],
                        Lambda = cov_data2,
                        PC = cc0+pc2,
                        results_alg1 = Alg1_data2,
                        eps_thresh=1e-2,
                        iterT=5)

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


data_mashr1 <- mash_set_data(beta_est[[1]], se_est[[1]], V = cor(beta_est[[1]]))
U.1 = cov_canonical(data_mashr1)
m.1 = mash(data_mashr1, U.1, output_lfdr=T)


meta_input_data1 <- matrix(NA, nrow = M, ncol = K[1] * 2)
meta_input_data1[,2*(1:K[1])-1] <- beta_est[[1]]
meta_input_data1[,2*(1:K[1])] <- se_est[[1]]
meta_input_data1 <- data.frame(SNP = paste0("rs",1:M),meta_input_data1)
fwrite(meta_input_data1, file = paste0("/scratch/t.phs.yihaolu/LLR/simulation/Metasoft/input/",token,"input_data1.txt"), sep = "\t", col.names = F)

system(paste0("cd /scratch/t.phs.yihaolu/LLR/simulation/Metasoft/;java -jar Metasoft.jar -input ../Metasoft/input/",
              token,"input_data1.txt -mvalue true -mvalue_p_thres 0.05 -mvalue_prior_sigma 0.001 -output ../Metasoft/input/",token,"_data1"))
result_data1 <- fread(paste0("/scratch/t.phs.yihaolu/LLR/simulation/Metasoft/input/",token,"_data1"))
pp1 <- data.matrix(result_data1[,(17+K[1]):(17+2*K[1]-1)]);pp1[is.na(pp1)] <- 0

### HT-eQTL
setwd("/scratch/t.phs.yihaolu/XING_NG/simulation")
source('Revision_HT_eQTL.R')
# clip for robustness
Z <- lapply(Z, function(x) {
  x[x>25] <- 25
  x[x<(-25)] <- (-25)
  return(x)
})
# to speed up computation.
ht_pp1_first <- HT_eQTL(Z[[1]][,1:k11])
ht_pp1_second <- HT_eQTL(Z[[1]][,(k11+1):(k11+k12)])
ht_pp1 <- cbind(ht_pp1_first,ht_pp1_second)

#==================================================================#
## result evaluation
#==================================================================#
dynamic.index1 <- which((rowSums(gamma[[1]])>=2 & rowSums(gamma[[1]])<=5)|(rowSums(gamma[[1]])==0))
auc_df <- data.frame(data1 = c(auc(c(gamma[[1]]), c(Alg1_data1$alphajkm1),quiet=T),
                               auc(c(gamma[[1]]), c(Alg2_data1$alphajkm1),quiet=T),
                               auc(c(gamma[[1]]), c(Alg4_mQTL_eQTL$alphajkm2_dat1),quiet=T),
                               auc(c(gamma[[1]]), c(1-m.1$result$lfdr),quiet=T),
                               auc(c(gamma[[1]]), c(pp1),quiet=T),
                               auc(c(gamma[[1]]), c(ht_pp1),quiet=T)),
                     data1.dynamic=c(auc(c(gamma[[1]][dynamic.index1,]),
                                         c(Alg1_data1$alphajkm1[dynamic.index1,]),quiet=T),
                                     auc(c(gamma[[1]][dynamic.index1,]),
                                         c(Alg2_data1$alphajkm1[dynamic.index1,]),quiet=T),
                                     auc(c(gamma[[1]][dynamic.index1,]),
                                         c(Alg4_mQTL_eQTL$alphajkm2_dat1[dynamic.index1,]),quiet=T),
                                     auc(c(gamma[[1]][dynamic.index1,]),
                                         c(1-m.1$result$lfdr[dynamic.index1,]),quiet=T),
                                     auc(c(gamma[[1]][dynamic.index1,]),
                                         c(pp1[dynamic.index1,]),quiet=T),
                                     auc(c(gamma[[1]][dynamic.index1,]),
                                         c(ht_pp1[dynamic.index1,]),quiet=T)))

rmse_df <- data.frame(Method = c("LLR","mash","Marginal model"),
                      data1 = c(rmse(true_beta[[1]],(Alg4_mQTL_eQTL$mujkm2_dat1 * se_est[[1]])),
                                rmse(true_beta[[1]],m.1$result$PosteriorMean),
                                rmse(true_beta[[1]],beta_est[[1]])),
                      data1.true = c(rmse(true_beta[[1]][gamma[[1]]!=0],(Alg4_mQTL_eQTL$mujkm2_dat1 * se_est[[1]])[gamma[[1]]!=0]),
                                     rmse(true_beta[[1]][gamma[[1]]!=0],m.1$result$PosteriorMean[gamma[[1]]!=0]),
                                     rmse(true_beta[[1]][gamma[[1]]!=0],beta_est[[1]][gamma[[1]]!=0])),
                      data1.dynamic = c(rmse(true_beta[[1]][dynamic.index1,],(Alg4_mQTL_eQTL$mujkm2_dat1 * se_est[[1]])[dynamic.index1,]),
                                        rmse(true_beta[[1]][dynamic.index1,],m.1$result$PosteriorMean[dynamic.index1,]),
                                        rmse(true_beta[[1]][dynamic.index1,],beta_est[[1]][dynamic.index1,])))
result <- list(parameter = args, auc = auc_df, rmse = rmse_df)








