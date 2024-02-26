args <- commandArgs(trailingOnly = T)
chunk <- as.numeric(args[1])
size <- as.numeric(args[2])
bound <- as.character(args[3])
type <- as.character(args[4])

library(tidyverse)
library(dplyr)
library(data.table)
setwd("/scratch/t.phs.yihaolu/LLR/result")
# load cis-results
pc0 = 2; result <- list()
for (chr in 1:22) {
  print(chr)
  load(paste("pi_e0.1pi_m0.001vk0.01",
             "chr",chr,"pc",pc0,"RData",sep = "."))
  result[[chr]] <- result_list
  load(paste0("/scratch/t.phs.yihaolu/GTEx_data/merge_data/merge_eQTL_mQTL_chr",chr,"_zscore.RData"))
  result[[chr]][[3]] <- merge_eQTL_mQTL %>% na.omit() %>% dplyr::select(Gene, cpg, variant_id)
  result[[chr]][[4]] <- merge_eQTL_mQTL %>% na.omit() %>% dplyr::select(contains(".x"))
  result[[chr]][[5]] <- merge_eQTL_mQTL %>% na.omit() %>% dplyr::select(contains(".y"))
}
mQTL_result <- lapply(result, function(x) data.frame(x[[1]])) %>% rbindlist()
eQTL_result <- lapply(result, function(x) data.frame(x[[2]])) %>% rbindlist()
variant_result <- lapply(result, function(x) x[[3]]) %>% rbindlist()

variant_result <- variant_result %>% mutate(ind=1:n())
cis_eQTL_sign <- which(rowSums(eQTL_result>0.8)>=1)
cis_mQTL_sign <- which(rowSums(mQTL_result>0.8)>=1)
variant_result <- variant_result[union(cis_eQTL_sign,cis_mQTL_sign),]
eQTL_result <- eQTL_result[union(cis_eQTL_sign,cis_mQTL_sign),]
mQTL_result <- mQTL_result[union(cis_eQTL_sign,cis_mQTL_sign),]
colnames(mQTL_result) <- colnames(eQTL_result) <-
  c("Breast_Mammary_Tissue","Kidney_Cortex","Muscle_Skeletal",
    "Prostate","Testis","Colon_Transverse","Lung","Ovary","Whole_Blood")
cis_eQTL_tissue_index <- apply(eQTL_result,1,function(x) which(x>0.8))
cis_mQTL_tissue_index <- apply(mQTL_result,1,function(x) which(x>0.8))
cis_eQTL_tissue_num <- lengths(cis_eQTL_tissue_index)
cis_mQTL_tissue_num <- lengths(cis_mQTL_tissue_index)
cis_eQTL_tissue = sapply(cis_eQTL_tissue_index,function(x) paste(x,collapse = "+"))
cis_mQTL_tissue = sapply(cis_mQTL_tissue_index,function(x) paste(x,collapse = "+"))

cis_info <- variant_result %>%
  mutate(eQTL_tissue=cis_eQTL_tissue,
         mQTL_tissue=cis_mQTL_tissue,
         cis_eQTL_num = cis_eQTL_tissue_num,
         cis_mQTL_num = cis_mQTL_tissue_num)

cis_eQTL_info <- cis_info %>%
  select(Gene,variant_id,eQTL_tissue,cis_eQTL_num) %>%
  filter(cis_eQTL_num >= 1) %>%
  distinct_all()
colnames(cis_eQTL_info) <- c("cis_gene","variant_id","cis_eQTL_tissue","cis_eQTL_num")
cis_mQTL_info <- cis_info %>%
  select(cpg,variant_id,mQTL_tissue,cis_mQTL_num) %>%
  filter(cis_mQTL_num >= 1) %>%
  distinct_all()
colnames(cis_mQTL_info) <- c("cis_cpg","variant_id","cis_mQTL_tissue","cis_mQTL_num")
# load trans-results
disease_trait_list <- read.csv("/scratch/t.phs.yihaolu/LLR/trait_disease_used.csv")
list_to_use = disease_trait_list[,4]
setwd("/scratch/t.phs.yihaolu/external_data/mediation/result/zscore_from_pval")
res <- list()
for (i in 1:length(list_to_use)) {
  filename <- paste0(list_to_use[i],".sign.5e-06.RData")
  if(file.exists(filename)) {
    load(filename)
    res[[i]] <- trans_res
    res[[i]][[1]] <- res[[i]][[1]] %>% mutate(trait = list_to_use[i])
  }
}

info <- res %>% lapply(function(x) x[[1]]) %>% do.call("rbind",.)
eQTL <- res %>% lapply(function(x) x[[2]]) %>% do.call("rbind",.)
mQTL <- res %>% lapply(function(x) x[[3]]) %>% do.call("rbind",.)

trans_eQTL_tissue_index <- apply(eQTL,1,function(x) which(x>0.8))
trans_mQTL_tissue_index <- apply(mQTL,1,function(x) which(x>0.8))
trans_eQTL_tissue = sapply(trans_eQTL_tissue_index,function(x) paste(x,collapse = "+"))
trans_mQTL_tissue = sapply(trans_mQTL_tissue_index,function(x) paste(x,collapse = "+"))
trans_eQTL_tissue_num <- lengths(trans_eQTL_tissue_index)
trans_mQTL_tissue_num <- lengths(trans_mQTL_tissue_index)
trans_info <- info %>% mutate(eQTL_tissue=trans_eQTL_tissue,
                              mQTL_tissue=trans_mQTL_tissue,
                              eQTL_num = trans_eQTL_tissue_num,
                              mQTL_num = trans_mQTL_tissue_num)

merge_info <- trans_info %>%
  left_join(cis_eQTL_info, by = "variant_id") %>%
  left_join(cis_mQTL_info, by = "variant_id") %>%
  mutate(trans_gene = ifelse(eQTL_num>0,1,0),
         trans_cpg = ifelse(mQTL_num>0,1,0))
trans_eQTL_info <- merge_info %>%
  filter(eQTL_num > 0) %>%
  select(gene,variant_id,eQTL_tissue,
         cis_gene,cis_cpg,cis_eQTL_tissue,cis_mQTL_tissue,
         eQTL_num,cis_eQTL_num,cis_mQTL_num) %>%
  distinct_all()
cpg_gene <- fread("/scratch/t.phs.yihaolu/LLR/trans_mQTL/cpg_gene_cor/gene_cpg_df.txt")
cpg_gene <- cpg_gene %>% group_by(gene) %>% slice_max(abs(rho),n=1,with_ties = F)
trans_mQTL_info <- merge_info %>%
  filter(mQTL_num > 0) %>%
  select(cpg,variant_id,mQTL_tissue,
         cis_gene,cis_cpg,cis_eQTL_tissue,cis_mQTL_tissue,
         mQTL_num,cis_eQTL_num,cis_mQTL_num) %>%
  distinct_all() %>%
  filter(!(is.na(cis_gene) & is.na(cis_cpg))) %>%
  filter(!(cpg %in% cpg_gene$cpg))

colnames(trans_mQTL_info)[8] <- colnames(trans_eQTL_info)[8] <- "trans_num"
colnames(trans_eQTL_info)[1:5] <- colnames(trans_mQTL_info)[1:5] <-
  c("probe","variant_id","tissue","cisgene","ciscpg")

full_mat <- rbind(trans_eQTL_info,trans_mQTL_info)

full_mat <- full_mat %>% ungroup() %>%
  mutate(trans_gene = c(rep(1,nrow(trans_eQTL_info)),
                        rep(0,nrow(trans_mQTL_info)))) %>%
  mutate(cis_gene = ifelse(is.na(cisgene),0,1), cis_cpg = ifelse(is.na(ciscpg),0,1))

full_mat <- full_mat %>%
  separate(col = "variant_id", into = c("SNP.chr","SNP.pos"),remove=F)
fwrite(full_mat, file = paste0("/scratch/t.phs.yihaolu/external_data/mediation/",bound,".full_mat_wide_format_z_p_at_least_1.txt"))

full_mat <- trans_mQTL_info
full_mat <- full_mat %>% ungroup() %>%
  mutate(trans_gene = c(rep(0,nrow(trans_mQTL_info)))) %>%
  mutate(cis_gene = ifelse(is.na(cisgene),0,1), cis_cpg = ifelse(is.na(ciscpg),0,1)) %>%
  separate(col = "variant_id", into = c("SNP.chr","SNP.pos"),remove=F)
fwrite(full_mat, file = paste0("/scratch/t.phs.yihaolu/external_data/mediation/",bound,".full_mat_wide_format_z_p_at_least_1_extra.txt"))

full_mat <- fread(paste0("/scratch/t.phs.yihaolu/external_data/mediation/",bound,".full_mat_wide_format_",type,"_at_least_1_extra.txt"))
if((1+(chunk-1)*size)>nrow(full_mat)) {
  stop("EXCEED MAX.")
}

setwd("/scratch/t.phs.yihaolu/external_data/mediation/")
tissue_used <- c("Breast_Mammary_Tissue","Kidney_Cortex","Muscle_Skeletal",
                 "Prostate","Testis","Colon_Transverse","Lung","Ovary","Whole_Blood")
expmat <- covmat <- methmat <- methcovmat <- list()
for (i in 1:length(tissue_used)) {
  print(i)
  expmat[[i]] <- fread(paste0("GTEx_Analysis_v8_eQTL_expression_matrices/",tissue_used[i],".v8.normalized_expression.bed.gz"))
  covmat[[i]] <- fread(paste0("GTEx_Analysis_v8_eQTL_covariates/",tissue_used[i],".v8.covariates.txt"))
  methmat[[i]] <- fread(paste0("/gpfs/data/linchen-lab/Yihao/mQTL_mapping/Phenotypes/",tissue_used[i],".bed.gz"))
  methmat[[i]] <- methmat[[i]] %>% filter(gene_id %in% union(unique(full_mat$ciscpg),unique(full_mat$probe)))
  methcovmat[[i]] <- fread(paste0("/gpfs/data/linchen-lab/Yihao/mQTL_mapping/Covariates/",tissue_used[i],".covariates.txt"))
}
names(expmat) <- names(covmat) <- names(methmat) <- names(methcovmat) <- tissue_used
covmat <- lapply(covmat,function(x) x[,-1] %>% data.matrix() %>% t())
methcovmat <- lapply(methcovmat,function(x) x[,-1] %>% data.matrix() %>% t())

result_mat_sub <- full_mat[(1+(chunk-1)*size):min(chunk*size,nrow(full_mat)),]

if(nrow(result_mat_sub) == 0) next 
genotype = fread("raw.use.xing.correct")

SNP_list <- unique(result_mat_sub$variant_id)
geno_id <- genotype$FID
genotype_sub <- SNP_list %>% sapply(function(x) which(str_starts(colnames(genotype),x)))
genotype_sub <- data.matrix(genotype[,..genotype_sub])
rownames(genotype_sub) <- geno_id
final_result <- rep(list(matrix(NA,nrow = nrow(result_mat_sub),ncol = 11)),9)


for (i in 1:nrow(result_mat_sub)) {
  if (i%%100 == 0) print(i)
  genoind <- which(str_starts(colnames(genotype_sub),result_mat_sub[i,]$variant_id))
  x0 <- genotype_sub[,genoind]
  tryCatch({
    for (tissue_index in 1:9) {
      if(result_mat_sub$cis_gene[i]==1) {
        cis_gene <- expmat[[tissue_index]] %>% 
          filter(gene_id == result_mat_sub$cisgene[i])
        cis_gene <- cis_gene[,-c(1:4)] %>% data.matrix() %>% t()
      }
      if(result_mat_sub$cis_cpg[i]==1) {
        cis_cpg <- methmat[[tissue_index]] %>% 
          filter(gene_id == result_mat_sub$ciscpg[i])
        cis_cpg <- cis_cpg[,-c(1:4)] %>% data.matrix() %>% t()
      }
      if(result_mat_sub$trans_gene[i] == 0){
        trans_probe <- methmat[[tissue_index]] %>% 
          filter(gene_id == result_mat_sub$probe[i])
        trans_probe <- trans_probe[,-c(1:4)] %>% data.matrix() %>% t()
      } else {
        trans_probe <- expmat[[tissue_index]] %>%
          filter(gene_id == result_mat_sub$probe[i])
        trans_probe <- trans_probe[,-c(1:4)] %>% data.matrix() %>% t()
      }
      
      
      base_model <- model_meth <- model_exp <- model_both <-
        exp_mediator_pval <- meth_mediator_pval <- 
        trans_pval <- cis_exp_pval <- cis_meth_pval <- 
        both_mediator_pval <- NA
      if(result_mat_sub$trans_gene[i] == 1){
        sample_overlap <- Reduce(intersect,list(names(x0),
                                                rownames(covmat[[tissue_index]]),
                                                rownames(trans_probe)))
        x <- x0[sample_overlap] %>% data.matrix()
        trans_probe_reg <- trans_probe[sample_overlap,] %>% data.matrix()
        cov_used <- covmat[[tissue_index]][sample_overlap,] %>% data.matrix()
        model <- summary(lm(trans_probe_reg ~ x + cov_used))[[4]]
        trans_pval <- model[2,4]
        base_full <- model[2,1]
      }
      if(result_mat_sub$trans_gene[i] == 0){
        sample_overlap <- Reduce(intersect,list(names(x0),
                                                rownames(methcovmat[[tissue_index]]),
                                                rownames(trans_probe)))
        x <- x0[sample_overlap] %>% data.matrix()
        trans_probe_reg <- trans_probe[sample_overlap,] %>% data.matrix()
        cov_used <- methcovmat[[tissue_index]][sample_overlap,] %>% data.matrix()
        model <- summary(lm(trans_probe_reg ~ x + cov_used))[[4]]
        trans_pval <- model[2,4]
        base_full <- model[2,1]
      }
      if(result_mat_sub$cis_gene[i] == 1){
        sample_overlap <- Reduce(intersect,list(names(x0),
                                                rownames(covmat[[tissue_index]]),
                                                rownames(cis_gene)))
        x <- x0[sample_overlap] %>% data.matrix()
        cis_gene_reg <- cis_gene[sample_overlap,] %>% data.matrix()
        cov_used <- covmat[[tissue_index]][sample_overlap,] %>% data.matrix()
        cis_exp_pval <- summary(lm(cis_gene_reg ~ x + cov_used))[[4]][2,4]
      }
      if(result_mat_sub$cis_cpg[i] == 1){
        sample_overlap <- Reduce(intersect,list(names(x0),
                                                rownames(methcovmat[[tissue_index]]),
                                                rownames(cis_cpg)))
        x <- x0[sample_overlap] %>% data.matrix()
        cis_cpg_reg <- cis_cpg[sample_overlap,] %>% data.matrix()
        cov_used <- methcovmat[[tissue_index]][sample_overlap,] %>% data.matrix()
        cis_meth_pval <- summary(lm(cis_cpg_reg ~ x + cov_used))[[4]][2,4]
      }
      if(result_mat_sub$cis_gene[i] == 1& result_mat_sub$trans_gene[i] == 1) {
        sample_overlap <- Reduce(intersect,list(names(x0),rownames(cis_gene),
                                                rownames(covmat[[tissue_index]]),
                                                rownames(trans_probe)))
        if(length(sample_overlap) >= 10) {
          x <- x0[sample_overlap] %>% data.matrix()
          cis_gene_reg <- cis_gene[sample_overlap,] %>% data.matrix()
          trans_probe_reg <- trans_probe[sample_overlap,] %>% data.matrix()
          cov_used <- covmat[[tissue_index]][sample_overlap,] %>% data.matrix()
          
          base_model <- summary(lm(trans_probe_reg ~ x + cov_used))[[4]][2,1]
          model_exp <- summary(lm(trans_probe_reg ~ x + cis_gene_reg + cov_used))[[4]][2,1]
          exp_mediator_pval <- summary(lm(trans_probe_reg ~ cis_gene_reg + x + cov_used))[[4]][2,4]
        }
        
        
      }
      if(result_mat_sub$cis_cpg[i] == 1 & result_mat_sub$trans_gene[i] == 0) {
        sample_overlap <- Reduce(intersect,list(names(x0),rownames(cis_cpg),
                                                rownames(methcovmat[[tissue_index]]),
                                                rownames(trans_probe)))
        if(length(sample_overlap) >= 10) {
          x <- x0[sample_overlap] %>% data.matrix()
          cis_cpg_reg <- cis_cpg[sample_overlap,] %>% data.matrix()
          trans_probe_reg <- trans_probe[sample_overlap,] %>% data.matrix()
          cov_used <- methcovmat[[tissue_index]][sample_overlap,] %>% data.matrix()
          
          base_model <- summary(lm(trans_probe_reg ~ x + cov_used))[[4]][2,1]
          model_meth <- summary(lm(trans_probe_reg ~ x + cis_cpg_reg + cov_used))[[4]][2,1]
          meth_mediator_pval <- summary(lm(trans_probe_reg ~ cis_cpg_reg + x + cov_used))[[4]][2,4]
        }
        
        
      }
      if(result_mat_sub$cis_gene[i] == 1 & result_mat_sub$trans_gene[i] == 0) {
        sample_overlap <- Reduce(intersect,list(names(x0),
                                                rownames(cis_gene),
                                                rownames(methcovmat[[tissue_index]]),
                                                rownames(trans_probe)))
        if(length(sample_overlap) >= 10) {
          x <- x0[sample_overlap] %>% data.matrix()
          cis_gene_reg <- cis_gene[sample_overlap,] %>% data.matrix()
          trans_probe_reg <- trans_probe[sample_overlap,] %>% data.matrix()
          cov_used <- methcovmat[[tissue_index]][sample_overlap,] %>% data.matrix()
          
          base_model <- summary(lm(trans_probe_reg ~ x + cov_used))[[4]][2,1]
          model_exp <- summary(lm(trans_probe_reg ~ x + cis_gene_reg + cov_used))[[4]][2,1]
          exp_mediator_pval <- summary(lm(trans_probe_reg ~ cis_gene_reg + x + cov_used))[[4]][2,4]
          
        }
        
      }
      if(result_mat_sub$cis_cpg[i] == 1 & result_mat_sub$trans_gene[i] == 1) {
        sample_overlap <- Reduce(intersect,list(names(x0),
                                                rownames(cis_cpg),
                                                rownames(methcovmat[[tissue_index]]),
                                                rownames(trans_probe)))
        if(length(sample_overlap) >= 10) {
          x <- x0[sample_overlap] %>% data.matrix()
          cis_cpg_reg <- cis_cpg[sample_overlap,] %>% data.matrix()
          trans_probe_reg <- trans_probe[sample_overlap,] %>% data.matrix()
          cov_used <- methcovmat[[tissue_index]][sample_overlap,] %>% data.matrix()
          
          base_model <- summary(lm(trans_probe_reg ~ x + cov_used))[[4]][2,1]
          model_meth <- summary(lm(trans_probe_reg ~ x + cis_cpg_reg + cov_used))[[4]][2,1]
          meth_mediator_pval <- summary(lm(trans_probe_reg ~ cis_cpg_reg + x + cov_used))[[4]][2,4]
        }
      }
      if(result_mat_sub$cis_gene[i] == 1 & result_mat_sub$cis_cpg[i] == 1) {
        sample_overlap <- Reduce(intersect,list(names(x0),
                                                rownames(cis_gene),rownames(trans_probe),
                                                rownames(cis_cpg),rownames(trans_probe)))
        if(length(sample_overlap) >= 10) {
          x <- x0[sample_overlap] %>% data.matrix()
          cis_gene_reg <- cis_gene[sample_overlap,] %>% data.matrix()
          cis_cpg_reg <- cis_cpg[sample_overlap,] %>% data.matrix()
          trans_probe_reg <- trans_probe[sample_overlap,] %>% data.matrix()
          cov_used <- methcovmat[[tissue_index]][sample_overlap,] %>% data.matrix()
          cov_used <- cov_used[,-c(3:5)]
          
          # base_model <- summary(lm(trans_probe ~ x + cov_used))[[4]][2,1]
          # model_exp <- summary(lm(trans_probe ~ x + cis_gene + cov_used))[[4]][2,1]
          # model_meth <- summary(lm(trans_probe ~ x + cis_cpg + cov_used))[[4]][2,1]
          model_both <- summary(lm(trans_probe_reg ~ x + cis_gene_reg + cis_cpg_reg + cov_used))[[4]][2,1]
          
          x[is.na(x)] <- mean(x,na.rm=T)
          both_mediator_model1 <- (lm(trans_probe_reg ~ cov_used))
          both_mediator_model2 <- (lm(trans_probe_reg ~ cis_gene_reg + cis_cpg_reg + x + cov_used))
          
          both_mediator_pval <- anova(both_mediator_model2,both_mediator_model1)$`Pr(>F)`[2]
        }
        
      }
      
      final_result[[tissue_index]][i,] <- c(base_model, model_exp, model_meth, model_both, 
                                            exp_mediator_pval, meth_mediator_pval, both_mediator_pval,
                                            trans_pval,cis_exp_pval,cis_meth_pval,base_full)
    }
  },silent = TRUE, error = function(x) return(paste0("error message ",i)))
  
}
final_result <- lapply(final_result, function(x) {
  colnames(x) <- c("base_model","model_exp","model_meth","model_both",
                   "exp_mediator_pval","meth_mediator_pval","both_mediator_pval",
                   "trans_pval","cis_exp_pval","cis_meth_pval","base_full")
  return(x)
})

final_result <- lapply(final_result, function(x) data.frame(result_mat_sub, x))
save(final_result, file = paste0("calculate_eff/eff_result_chr",chunk,
                                 "_",size,"_",bound,"_",type,"_at_least_1.RData"))

# script
cd /scratch/t.phs.yihaolu/external_data/mediation/script
bound=5e-06
type=z_p
for chunk in {1..700}; do
for size in 1000; do echo "
#PBS -l mem=10gb
cd /scratch/t.phs.yihaolu/external_data/mediation
Rscript calculate_eff.R $chunk $size $bound $type">calculate_eff$chunk.$size.$bound.$type.sh
qsub calculate_eff$chunk.$size.$bound.$type.sh; done; done



# collect
library(tidyverse)
library(dplyr)
library(data.table)
setwd("/scratch/t.phs.yihaolu/external_data/mediation/calculate_eff/")
df <- list();count=1
bound <- "5e-06";size=1000
type = "b_se"
type = "z_p"
for (chunk in c(1:1000)) {
  filename <-  paste0("eff_result_chr",chunk,
                      # "_",size,"_",bound,"_b_se.RData")
                      "_",size,"_",bound,"_",type,"_at_least_1.RData")
  if(file.exists(filename)) {
    load(filename)
    df[[count]] <- final_result
  } else {
    print(chunk)
    # system(paste0("qsub ../script/calculate_eff",chunk,".",size,".",bound,".sh"))
  }
  count <- count+1
}
df <- df[lengths(df)>0]
df_list <- list()
for (i in 1:9) {
  print(i)
  df_list[[i]] <- lapply(df, function(x) x[[i]]) %>% do.call("rbind",.)
}

df_tab <- lapply(df_list, function(x) {
  x <- x %>% mutate(trans_mediated_by_e = 1 - model_exp/base_full,
                    trans_mediated_by_m = 1 - model_meth/base_full,
                    trans_mediated_by_em = 1 - model_both/base_full)
  return(x)})

save(df_tab, file = paste0("df_tab",bound,"_iter_",type,".RData"))
