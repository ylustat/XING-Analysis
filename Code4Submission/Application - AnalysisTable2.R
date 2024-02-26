# may also need the 'XING_revision.R' file. This is to find top/secondary mQTLs and run XING.
# count lead mQTL
library(tidyverse)
library(data.table)
setwd("/Users/yihaolu/Downloads/Research/XING_revision/Data/")
filenames <- dir()
filenames <- filenames[str_detect(filenames, ".regular.perm.fdr.txt")]
tissues <- str_replace_all(filenames, ".regular.perm.fdr.txt", "")
cpg_list <- lapply(filenames,fread) %>% 
  lapply(function(x) x %>% filter(qval < 0.05)) %>% 
  lapply(function(x) x$cpg_id)
cpg_list %>% unlist() %>% unique() %>% length() # 286152
names(cpg_list) <- tissues
# count secondary mQTL
library(tidyverse)
library(data.table)
args <- commandArgs(trailingOnly = T)
i <- as.numeric(args[1])
tissue <- c("BreastMammaryTissue", "ColonTransverse", "KidneyCortex", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")[i]
setwd("/scratch/t.phs.yihaolu/map_mQTL")
df <- fread(paste0(tissue,".mQTLs.conditional.txt.gz"))
df <- df %>% filter(str_detect(V1,":1")) %>% 
  group_by(V1) %>% slice_min(V7,n=1)
fwrite(df, file = paste0(tissue,".mQTLs.secondary.txt"))

cd /scratch/t.phs.yihaolu/map_mQTL/script
for i in {1..9}; do echo "
#PBS -l mem=50gb
cd /scratch/t.phs.yihaolu/map_mQTL
Rscript slice_secondary.R $i" > slice_$i.sh
qsub slice_$i.sh; done


cd /scratch/t.phs.yihaolu/map_mQTL/script
for tissue in BreastMammaryTissue ColonTransverse KidneyCortex Lung MuscleSkeletal Ovary Prostate Testis WholeBlood; do echo "
#PBS -l mem=5gb
cd /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/results
cp $tissue.cor.stats.txt /scratch/t.phs.yihaolu/map_mQTL/
cd /scratch/t.phs.yihaolu/map_mQTL/
cut -d $'\t' -f1,2,5 $tissue.cor.stats.txt > $tissue.cor.stats.simple.txt
awk 'BEGIN{FS=OFS=\"\t\"} NR==1 {print; next} !(\$1 in max) || abs(\$3) > abs(max[\$1]) { max[\$1] = \$3; line[\$1] = \$0 } END {for (i in line) print line[i]} function abs(v) {return v < 0 ? -v : v}' $tissue.cor.stats.simple.txt > $tissue.cor.stats.max.txt
" > slice_$tissue.sh
qsub slice_$tissue.sh; done

# Write the first file, including its header, to the output file
cd /scratch/t.phs.yihaolu/map_mQTL/
  cat BreastMammaryTissue.cor.stats.max.txt > combined.txt

# For the remaining files, skip the first line (the header) before appending to the output file
tail -n +2 ColonTransverse.cor.stats.max.txt >> combined.txt
tail -n +2 Lung.cor.stats.max.txt >> combined.txt
tail -n +2 MuscleSkeletal.cor.stats.max.txt >> combined.txt
tail -n +2 Ovary.cor.stats.max.txt >> combined.txt
tail -n +2 Prostate.cor.stats.max.txt >> combined.txt
tail -n +2 Testis.cor.stats.max.txt >> combined.txt
tail -n +2 WholeBlood.cor.stats.max.txt >> combined.txt

cd /scratch/t.phs.yihaolu/map_mQTL/script
echo "
#PBS -l mem=5gb
cd /scratch/t.phs.yihaolu/map_mQTL/
awk 'BEGIN{FS=OFS=\"\t\"} NR==1 {print; next} !(\$1 in max) || abs(\$3) > abs(max[\$1]) { max[\$1] = abs(\$3); line[\$1] = \$0 } END {for (i in line) print line[i]} function abs(v) {return v < 0 ? -v : v}' combined.txt | sort -k1 > combined.uniq.txt" > uniq.file.sh
qsub uniq.file.sh
# merge data
library(tidyverse)
library(data.table)
setwd("/Users/yihaolu/Downloads/Research/XING_revision/Data/")
filenames <- dir()
filenames <- filenames[str_detect(filenames, ".mQTLs.secondary.txt")]
tissues <- str_replace_all(filenames, ".mQTLs.secondary.txt", "")
secondary_cpg_list <- lapply(filenames,fread) %>%
  lapply(function(x) x %>% group_by(V1) %>% slice_min(V7,n=1,with_ties = T)) %>%
  do.call("rbind",.) %>%
  mutate(V1 = str_replace_all(V1,":1","")) %>%
  select(V1,V2) %>%
  distinct_all() %>%
  mutate(type = 1)
colnames(secondary_cpg_list) <- c("cpg_id","variant_id","type")

filenames <- dir()
filenames <- filenames[str_detect(filenames, ".regular.perm.fdr.txt")]
tissues <- str_replace_all(filenames, ".regular.perm.fdr.txt", "")
top_cpg_list <- lapply(filenames,fread) %>%
  lapply(function(x) x %>% filter(qval < 0.05)) %>%
  do.call("rbind",.) %>%
  select(cpg_id, variant_id) %>%
  distinct_all() %>%
  mutate(type = 2)

cpg_list <- rbind(top_cpg_list,secondary_cpg_list)

cpg_list <- cpg_list %>%
  group_by(cpg_id,variant_id) %>%
  slice_min(type, n = 1)
fwrite(cpg_list, file = "cpg_list.txt")

gene_cpg_cor <- fread("combined.uniq.txt")
colnames(gene_cpg_cor) <- unlist(tail(gene_cpg_cor,1)); gene_cpg_cor <- gene_cpg_cor[1:(nrow(gene_cpg_cor) - 1),]
gene_cpg_df <- gene_cpg_cor %>%
  inner_join(cpg_list, by = c("cpg" = "cpg_id")) %>%
  select(cpg,gene,variant_id,type) %>%
  relocate(gene, cpg, variant_id) %>%
  mutate(Chr = str_split(variant_id,"_") %>% sapply(function(x) x[1]))
fwrite(gene_cpg_df, file = "gene_cpg_info.txt")



library(data.table)
library(tidyverse)
setwd("/scratch/t.phs.yihaolu/XING_revision/map_eQTL/output")
# step 1: get gene-CpG list
gene_cpg_df <- fread("gene_cpg_info/gene_cpg_info.txt")
gene_cpg_df <- gene_cpg_df %>%
  split(list(gene_cpg_df$Chr))
for (i in 1:22) {
  print(i)
  fwrite(gene_cpg_df[[i]], file = paste0("gene_cpg_info/all_cpg_",names(gene_cpg_df)[i],"_gene_cpg_info.txt"),sep = "\t")
}
# step 2: subset mQTL data only keep cpg in gene-CpG list
cd /scratch/t.phs.yihaolu/XING_revision/map_mQTL/script
for chr in {1..22}; do
for tissue in Ovary BreastMammaryTissue Prostate Testis KidneyCortex Lung WholeBlood ColonTransverse MuscleSkeletal; do echo "
#PBS -l mem=5gb
cd /scratch/t.phs.yihaolu/XING_revision/map_mQTL/
awk -F\"\t\" 'NR==FNR{values[\$2]; next} \$1 in values' /scratch/t.phs.yihaolu/XING_revision/map_eQTL/output/gene_cpg_info/all_cpg_chr${chr}_gene_cpg_info.txt output/full_mQTL/${tissue}_chr${chr}.txt > filtered_mQTL/all_cpg_filtered_${tissue}_chr${chr}.txt" > fastqtl_${tissue}_${chr}.sh
qsub fastqtl_${tissue}_${chr}.sh; done; done
# step 3: merge multi-tissue mQTL
args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
stat = as.character(args[2])
library(data.table)
library(tidyverse)
setwd("/scratch/t.phs.yihaolu/XING_revision/map_mQTL")
tissue_list <- c("Ovary","BreastMammaryTissue","Prostate","Testis","KidneyCortex","Lung","WholeBlood","ColonTransverse","MuscleSkeletal")
prefix <- "all_cpg_"
df <- list()
for (i in 1:9) {
  print(i)
  df[[i]] <- fread(paste0("filtered_mQTL/",prefix,"filtered_",tissue_list[i],"_chr",chr,".txt")) 
}

if(stat == "pval") {
  df <- lapply(df,function(x) x %>% select(V1:V2,V7) %>% distinct_all())) %>%
  Reduce(function(x,y) inner_join(x,y,by="token"), .) %>% 
  distinct_all()
}

if(stat == "beta") {
  df <- lapply(df,function(x) x %>% select(V1:V2,V8) %>% distinct_all()) %>% 
    Reduce(function(x,y) inner_join(x,y,by=c("V1","V2")), .) %>% 
    distinct_all()
}

if(stat == "se") {
  df <- lapply(df,function(x) x %>% select(V1:V2,V9) %>% distinct_all()) %>% 
    Reduce(function(x,y) inner_join(x,y,by=c("V1","V2")), .) %>% 
    distinct_all()
}

colnames(df) <- c("cpg","SNP",tissue_list)
fwrite(df, file = paste0("merge_tissue/",prefix,"merge_",stat,"_chr",chr,".txt"))

cd /scratch/t.phs.yihaolu/XING_revision/map_mQTL/script
for chr in {1..22}; do 
for stat in beta se pval; do echo "
#PBS -l mem=100gb
cd /scratch/t.phs.yihaolu/XING_revision/map_mQTL/
Rscript merge_mQTL.R $chr $stat" > merge_${chr}_${stat}.sh
qsub merge_${chr}_${stat}.sh; done; done


library(data.table)
library(tidyverse)
setwd("/scratch/t.phs.yihaolu/map_mQTL")
gene_cpg_cor <- fread("combined.uniq.txt")
tissue_list1 <- c("Ovary","Breast_Mammary_Tissue","Prostate","Testis","Kidney_Cortex","Lung","Whole_Blood","Colon_Transverse","Muscle_Skeletal",
                  "Adipose_Subcutaneous","Brain_Spinal_cord_cervical_c-1","Nerve_Tibial","Adipose_Visceral_Omentum","Brain_Substantia_nigra","Adrenal_Gland","Pancreas","Artery_Aorta","Cells_Cultured_fibroblasts","Pituitary","Artery_Coronary","Cells_EBV-transformed_lymphocytes","Artery_Tibial","Colon_Sigmoid","Skin_Not_Sun_Exposed_Suprapubic","Brain_Amygdala","Skin_Sun_Exposed_Lower_leg","Brain_Anterior_cingulate_cortex_BA24","Esophagus_Gastroesophageal_Junction","Small_Intestine_Terminal_Ileum","Brain_Caudate_basal_ganglia","Esophagus_Mucosa","Spleen","Brain_Cerebellar_Hemisphere","Esophagus_Muscularis","Stomach","Brain_Cerebellum","Heart_Atrial_Appendage","Brain_Cortex","Heart_Left_Ventricle","Thyroid","Brain_Frontal_Cortex_BA9","Uterus","Brain_Hippocampus","Liver","Vagina","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Brain_Putamen_basal_ganglia")
df <- lapply(tissue_list1[1:9], function(x) fread(paste0("GTEx_Analysis_v8_eQTL/",x,".v8.signif_variant_gene_pairs.txt.gz"))) %>%
  do.call("rbind",.)
eQTL_list <- df %>% select(variant_id, gene_id) %>% distinct_all()
fwrite(eQTL_list, file = "eQTL_list.txt")
eQTL_list <- fread("eQTL_list.txt")

cpg_list <- fread("cpg_list.txt")
closest <- lapply(1:22, function(chr) fread(paste0("closest_gene_cpg_",chr,".txt"))) %>% do.call("rbind",.)
cpg_list <- cpg_list %>%
  inner_join(closest, by = c("cpg_id")) 
fwrite(cpg_list,"most_near_gene_cpg_info.txt")

# split lookup table
cd /scratch/t.phs.yihaolu/XING_revision/map_eQTL/
  wget https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz
zcat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz|awk -F"\t" '{print > ("lookup/lookup_" $2".txt")}'

# merge eQTL and mQTL, need snp, gene, cpg, ref allele, alt allele, eQTL eff/se/pval, mQTL eff/se/pval
args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
stat = as.character(args[2])
library(data.table)
library(tidyverse)
setwd("/scratch/t.phs.yihaolu/XING_revision/map_mQTL")
prefix = "all_cpg_"
gene_info <- fread("/home/yihaolu/XING/gene_info_GTEx_v8.txt") %>% select(gene_id,gene_symbol,gene_type)
eQTL <- fread(paste0("../map_eQTL/output/merge_tissue/merged_",stat,"_",chr,".txt"))
colnames(eQTL) <- paste0("eQTL_",colnames(eQTL))

mQTL <- fread(paste0("merge_tissue/",prefix,"merge_",stat,"_chr",chr,".txt"))
colnames(mQTL) <- paste0("mQTL_",colnames(mQTL))

gene_cpg_df <- fread(paste0("/scratch/t.phs.yihaolu/map_mQTL/most_near_gene_cpg_info.txt"))

eQTL <- eQTL %>% 
  inner_join(gene_cpg_df, by = c("eQTL_gene_id" = "gene_id","eQTL_variant_id" = "variant_id")) 
eQTL <- eQTL %>% 
  inner_join(mQTL, by = c("cpg_id" = "mQTL_cpg", "eQTL_variant_id" = "mQTL_SNP"))
eQTL <- eQTL %>%
  relocate(cpg_chr, cpg_id)

colnames(eQTL)[1:4] <- c("Chr","cpg","gene_id","variant_id")


lookup <- fread(paste0("../map_eQTL/lookup/lookup_chr",chr,".txt"),h=F)
cn <- fread("../map_eQTL/lookup/lookup_chr.txt")
colnames(lookup) <- colnames(cn)
eQTL <- eQTL %>% 
  inner_join(gene_info, by = "gene_id") %>% 
  inner_join(lookup, by = c("variant_id")) %>% 
  relocate(chr,ref,alt,rs_id_dbSNP151_GRCh38p7,gene_type) %>% 
  select(-Chr)
eQTL <- eQTL %>% select(-c("variant_id_b37","num_alt_per_site","variant_pos"))

# eQTL <- eQTL %>% filter(gene_type == "protein_coding")
fwrite(eQTL, paste0("processed_mQTL/nearby_",prefix,"processed_eQTL_mQTL_",stat,"_chr",chr,".txt"))


cd /scratch/t.phs.yihaolu/map_mQTL/script
for stat in beta se pval; do
for chr in {1..22}; do echo "
#PBS -l mem=45gb
cd /scratch/t.phs.yihaolu/XING_revision/map_mQTL/
Rscript merge_eQTL_mQTL.R $chr $stat" > merge_${chr}_${stat}.sh
qsub merge_${chr}_${stat}.sh; done; done




# perform XING
library(MASS)
library(CCA)
library(RGCCA)
library(pbapply)
library(data.table)
library(tidyverse)
setwd("/Users/yihaolu/Downloads/Research/XING_revision/Data")
source("../Scripts/AlgorithmsHeir1.R")
source("../Scripts/HierarchicalLLRAlg3-4_1.R")
tissue_list1 <- c("Ovary","Breast_Mammary_Tissue","Prostate","Testis","Kidney_Cortex","Lung","Whole_Blood","Colon_Transverse","Muscle_Skeletal",
                  "Adipose_Subcutaneous","Brain_Spinal_cord_cervical_c-1","Nerve_Tibial","Adipose_Visceral_Omentum","Brain_Substantia_nigra","Adrenal_Gland","Pancreas","Artery_Aorta","Cells_Cultured_fibroblasts","Pituitary","Artery_Coronary","Cells_EBV-transformed_lymphocytes","Artery_Tibial","Colon_Sigmoid","Skin_Not_Sun_Exposed_Suprapubic","Brain_Amygdala","Skin_Sun_Exposed_Lower_leg","Brain_Anterior_cingulate_cortex_BA24","Esophagus_Gastroesophageal_Junction","Small_Intestine_Terminal_Ileum","Brain_Caudate_basal_ganglia","Esophagus_Mucosa","Spleen","Brain_Cerebellar_Hemisphere","Esophagus_Muscularis","Stomach","Brain_Cerebellum","Heart_Atrial_Appendage","Brain_Cortex","Heart_Left_Ventricle","Thyroid","Brain_Frontal_Cortex_BA9","Uterus","Brain_Hippocampus","Liver","Vagina","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Brain_Putamen_basal_ganglia")
tissue_list1 <- c("Ovary","Breast_Mammary_Tissue","Prostate","Testis","Kidney_Cortex","Lung","Whole_Blood","Colon_Transverse","Muscle_Skeletal")
tissue_list1 <- c("Ovary","Breast_Mammary_Tissue","Prostate","Testis","Kidney_Cortex","Lung","Whole_Blood","Colon_Transverse","Muscle_Skeletal",
                  "Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Artery_Aorta","Artery_Coronary","Artery_Tibial","Colon_Sigmoid","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Heart_Atrial_Appendage","Brain_Cortex","Heart_Left_Ventricle","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia")
tissue_list2 <- c("Ovary","BreastMammaryTissue","Prostate","Testis","KidneyCortex","Lung","WholeBlood","ColonTransverse","MuscleSkeletal")

conflicted::conflict_prefer("select","dplyr")
conflicted::conflict_prefer("filter","dplyr")
beta_list <- se_list <- pval_list <- list()
for (chr in 1:22) {
  print(chr)
  beta_list[[chr]] <- fread(paste0("nearby_all_cpg_processed_eQTL_mQTL_beta_chr",chr,".txt"))
  se_list[[chr]] <- fread(paste0("nearby_all_cpg_processed_eQTL_mQTL_se_chr",chr,".txt"))
  pval_list[[chr]] <- fread(paste0("nearby_all_cpg_processed_eQTL_mQTL_pval_chr",chr,".txt"))
}

beta_list <- beta_list %>% do.call("rbind",.) %>% 
  select(rs_id_dbSNP151_GRCh38p7,cpg,gene_symbol,type,variant_id,
         paste0(
           "eQTL_",tissue_list1
         ),
         paste0(
           "mQTL_",tissue_list2
         ))
se_list <- se_list %>% do.call("rbind",.) %>% 
  select(rs_id_dbSNP151_GRCh38p7,cpg,gene_symbol,type,variant_id,
         paste0(
           "eQTL_",tissue_list1
         ),
         paste0(
           "mQTL_",tissue_list2
         ))
pval_list <- pval_list %>% do.call("rbind",.) %>% 
  select(rs_id_dbSNP151_GRCh38p7,cpg,gene_symbol,type,variant_id,
         paste0(
           "eQTL_",tissue_list1
         ),
         paste0(
           "mQTL_",tissue_list2
         ))

info_full <- beta_list %>% 
  select(rs_id_dbSNP151_GRCh38p7,gene_symbol,cpg,type,variant_id)

drop_index <- unique(which(is.na(beta_list)|is.na(se_list)|se_list==0,arr.ind=T)[,1])
if(length(drop_index)>0){
  beta_list <- beta_list[-drop_index,]
  se_list <- se_list[-drop_index,]
  pval_list <- pval_list[-drop_index,]
}

info <- beta_list %>% 
  select(rs_id_dbSNP151_GRCh38p7,gene_symbol,cpg,type,variant_id) %>% 
  mutate(n = 1:n())
##################################################################################################################
beta_list <- lapply(c("eQTL","mQTL"),function(x) beta_list[info$n,] %>% select(contains(x)) %>% data.matrix())
se_list <- lapply(c("eQTL","mQTL"),function(x) se_list[info$n,] %>% select(contains(x)) %>% data.matrix())
pval_list <- lapply(c("eQTL","mQTL"),function(x) pval_list[info$n,] %>% select(contains(x)) %>% data.matrix())
zscore_list <- mapply(function(x,y) x/y, beta_list, se_list,SIMPLIFY = F)
zscore_list <- append(zscore_list,list(info))

cc0 <- 4; pc1 <- pc2 <- 4
t1 <- 4; t2 <- 4; v1 <- 0.1; v2 <- 0.1

cov_data1 <- solve(cov(zscore_list[[1]][rowSums(abs(zscore_list[[1]])>t1)==0,]))
cov_data2 <- solve(cov(zscore_list[[2]][rowSums(abs(zscore_list[[2]])>t2)==0,]))
Alg1_data1 <- Alg1(betahat = zscore_list[[1]],
                   Lambda = cov_data1,
                   iterT = 10,
                   bound=1e-8, pi_init=colMeans(pval_list[[1]]<1e-5),vk_init=2*v1*colMeans(pval_list[[1]]<1e-5))
Alg1_data2 <- Alg1(betahat = zscore_list[[2]],
                   Lambda = cov_data2,
                   iterT = 10,
                   bound=1e-8, pi_init=colMeans(pval_list[[2]]<1e-5),vk_init=2*v2*colMeans(pval_list[[2]]<1e-5))
Alg2_data1 <- Alg2_ind1(betahat = zscore_list[[1]],
                        Lambda = cov_data1,
                        PC = pc1,
                        results_alg1 = Alg1_data1,
                        eps_thresh=1e-2,
                        iterT=5)
Alg2_data2 <- Alg2_ind1(betahat = zscore_list[[2]],
                        Lambda = cov_data2,
                        PC = pc2,
                        results_alg1 = Alg1_data2,
                        eps_thresh=1e-2,
                        iterT=5)
Alg4_mQTL_eQTL <- Alg4(betahat1 = zscore_list[[1]],
                       betahat2 = zscore_list[[2]],
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

pp1 <- Alg4_mQTL_eQTL$alphajkm2_dat1; pp2 <- Alg4_mQTL_eQTL$alphajkm2_dat2
# total number of SNP-CpG pairs
tot_pair <- info %>% nrow()
info$cpg %>% unique() %>% length()
# number of SNP-CpG pairs with non-zero effect in >=2 tissues
sign_pair <- sum(rowSums(Alg4_mQTL_eQTL$alphajkm2_dat2>0.8)>1)
# number of SNP-CpG pairs with non-zero effect in >=2 tissues 
coocc_pair <- sum(rowSums(Alg4_mQTL_eQTL$alphajkm2_dat2>0.8)>1 & rowSums(Alg4_mQTL_eQTL$alphajkm2_dat1>0.8)>0)
tot_pair; sign_pair; coocc_pair; sign_pair/tot_pair; coocc_pair/sign_pair

# unique SNPs
tot_SNP <- info$variant_id %>% unique() %>% length()
# number of SNP-CpG pairs with non-zero effect in >=2 tissues
sign_SNP <- info$variant_id[rowSums(Alg4_mQTL_eQTL$alphajkm2_dat2>0.8)>1] %>% unique() %>% length()
# number of SNP-CpG pairs with non-zero effect in >=2 tissues 
coocc_SNP <- info$variant_id[rowSums(Alg4_mQTL_eQTL$alphajkm2_dat2>0.8)>1 & rowSums(Alg4_mQTL_eQTL$alphajkm2_dat1>0.8)>0] %>% 
  unique() %>% length()
tot_SNP; sign_SNP; coocc_SNP; sign_SNP/tot_SNP; coocc_SNP/sign_SNP

Alg4_mQTL_eQTL <- Alg4_mQTL_eQTL[c("alphajkm2_dat1","alphajkm2_dat2","Lq_iter_dat1","Lq_iter_dat2","mujkm2_dat1","mujkm2_dat2")]
Alg4_mQTL_eQTL[3:6] <- NA
Alg4_mQTL_eQTL[3:4] <- pval_list
result <- append(Alg4_mQTL_eQTL, list(info))
i = 2
save(result, file = paste0("lead_mQTL_eQTL_",i,".RData"))



library(data.table)
library(tidyverse)
setwd("/scratch/t.phs.yihaolu/map_mQTL")
args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
x <- paste0("chr",chr)

gene_info <- fread("/home/yihaolu/XING/gene_info_GTEx_v8.txt") %>% select(Chr,gene_id,TSS)
cpg_info <- fread("cpg_info.txt") %>% select(`#chr`,gene_id,start)
colnames(gene_info) <- c("gene_chr","gene_id","gene_pos")
colnames(cpg_info) <- c("cpg_chr","cpg_id","cpg_pos")
cpg_list <- fread("cpg_list.txt")

closest <- cpg_info %>%
  filter(cpg_chr == x) %>%
  inner_join(gene_info %>%
               filter(gene_chr == x),
             by = c("cpg_chr" = "gene_chr")) %>%
  mutate(distance = abs(cpg_pos - gene_pos)) %>%
  filter(distance<1.5e6) %>%
  group_by(cpg_id) %>%
  slice_min(distance,n=1)
fwrite(closest, file = paste0("closest_gene_cpg_",chr,".txt"))

cd /scratch/t.phs.yihaolu/XING_revision/map_mQTL/script
for chr in {1..22}; do echo "
#PBS -l mem=100gb
cd /scratch/t.phs.yihaolu/XING_revision/map_mQTL/
Rscript find_closest.R $chr" > merge_${chr}.sh
qsub merge_${chr}.sh; done