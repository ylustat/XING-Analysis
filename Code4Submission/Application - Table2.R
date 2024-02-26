# download mQTL data
cd /scratch/t.phs.yihaolu/
  ./gdrive download 1FIWN9lTPIeiHiYCmx0ufoQz1TT1_qcpU --path /scratch/t.phs.yihaolu/external_data/external_mQTLs

zcat FUSION.all.tsv.gz|awk '{print > "FUSION_chr"$2".txt"}' -
  
  for i in {1..22}; do
echo "
  #PBS -l mem=15gb
  cd /scratch/t.phs.yihaolu/external_data/external_mQTLs
  cat FUSION_chr$i.txt | sed 's/|/ /' | awk '{print \$1\"\t\"\$9\"\t\"\$19}' > FUSION_chr${i}_simple.txt
  " > cut_data_${i}.sh
qsub cut_data_${i}.sh
done

rm cut_data_*.sh*
  
# download eQTL data
#PBS -l mem=15gb
cd /scratch/t.phs.yihaolu/external_data/external_eQTLs
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
zcat 2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | awk -F'\t' '{ print $1, $2, $3, $9}' > eQTLGen.txt
awk '{print > "eQTLGen_chr"$3".txt"}' eQTLGen.txt

cd /scratch/t.phs.yihaolu/external_data/external_eQTLs
vi process_eQTL.sh
qsub process_eQTL.sh

# download lookup data
cd /scratch/t.phs.yihaolu/external_data/
  wget https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz

zcat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz | awk -F'\t' '{ print $1, $2, $7 }' > lookup.txt
grep -v -E '(\.|chrX|chrY)' lookup.txt > lookup_clean.txt
mv lookup_clean.txt lookup.txt
awk '{print > "lookup_"$2".txt"}' lookup.txt


awk 'FNR==NR {a[$3]; next} $3 in a' ../lookup_chr1.txt FUSION_chr1_simple.txt > FUSION_chr1_simple_intersect.txt
awk 'NR==FNR{values[$3]++; next} $3 in values' ../lookup_chr1.txt FUSION_chr1_simple.txt > FUSION_chr1_simple_intersect.txt

###########################################################################################################################################################
## mQTL REPLICATION
###########################################################################################################################################################
args <- commandArgs(TRUE)
file = as.numeric(args[1])
chr = as.numeric(args[2])

library(data.table)
library(tidyverse)
library(pbapply)
setwd("/scratch/t.phs.yihaolu/XING_revision/")
home_path <- "/home/yihaolu/XING/"
external_path <- "/scratch/t.phs.yihaolu/external_data/"
gene_info <- fread(paste0(home_path,"gene_info_GTEx_v8.txt")) %>% select(gene_id, Chr, gene_symbol)
gene_info$gene_id <- str_split(gene_info$gene_id,"[.]") %>% lapply(function(x) x[1]) %>% unlist()

load(paste0("result/Upload/lead_mQTL_eQTL",file,".RData"))
intersect_data <- list()

print(chr)
info <- result[[7]] %>% ungroup() %>%
  mutate(ind = 1:n()) %>%
  inner_join(gene_info, by = "gene_symbol") %>%
  filter(Chr == paste0("chr",chr))

FUSION <- fread(paste0(external_path, "external_mQTLs/FUSION_chr",chr,"_simple.txt"))
colnames(FUSION) <- c("CPG","P","SNP")

FUSION <- FUSION[which(FUSION$CPG %in% info$cpg),]
FUSION <- FUSION[which(FUSION$SNP %in% info$rs_id_dbSNP151_GRCh38p7),]

intersect_data <- info %>%
  inner_join(FUSION, by = c("rs_id_dbSNP151_GRCh38p7" = "SNP", "cpg" = "CPG"))

fwrite(intersect_data, file = paste0("result/intersect_data_chr",chr,"_file",file,".txt"))

cd /scratch/t.phs.yihaolu/XING_revision/map_eQTL/script
for file in 2; do
for chr in {1..22}; do echo "
#PBS -l mem=15gb
cd /scratch/t.phs.yihaolu/XING_revision/
Rscript replicate_mQTL.R ${file} $chr" > simplify_${file}_${chr}.sh
qsub simplify_${file}_${chr}.sh; done; done



library(data.table)
library(tidyverse)
library(pbapply)
for (file in c(2)) {
  setwd("/scratch/t.phs.yihaolu/XING_revision/")
  home_path <- "/home/yihaolu/XING/"
  external_path <- "/scratch/t.phs.yihaolu/external_data/"
  gene_info <- fread(paste0(home_path,"gene_info_GTEx_v8.txt")) %>% select(gene_id, Chr, gene_symbol)
  gene_info$gene_id <- str_split(gene_info$gene_id,"[.]") %>% lapply(function(x) x[1]) %>% unlist()
  
  load(paste0("result/Upload/lead_mQTL_eQTL_",file,".RData"))
  
  intersect_data <- list()
  for (chr in 1:22) {
    if(file.exists(paste0("result/intersect_data_chr",chr,"_file",file,".txt"))) {
      intersect_data[[chr]] <- fread(paste0("result/intersect_data_chr",chr,"_file",file,".txt"))
    }
  }
  intersect_data <- do.call("rbind",intersect_data)
  
  # overall replication
  print(1 - qvalue::qvalue_truncp(intersect_data$P[which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) == 0)])$pi0)
  print(1 - qvalue::qvalue_truncp(intersect_data$P[which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) >= 1)])$pi0)
  
  mean(intersect_data$P[which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) >= 0)]<6e-7)
  mean(intersect_data$P[which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) >= 1)]<6e-7)
  mean(intersect_data$P[which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) >= 2)]<6e-7)
  
  # cross-omics replication comparison
  # intersect_data <- intersect_data %>% filter(type == 2)
  subPP_e <- result$alphajkm2_dat1[intersect_data$ind,]
  subPP_m <- result$alphajkm2_dat2[intersect_data$ind,]
  index_mQTL <- apply(subPP_m, 2, function(x) which(x > 0.8))
  index_eQTL <- which(rowSums(subPP_e>0.8)>=1)
  
  index_multi <- which(rowSums(subPP_m > 0.8)>1)
  index_mQTL <- lapply(index_mQTL, function(x) intersect(x,index_multi))
  
  rbind(lapply(index_mQTL, function(x) setdiff(x, index_eQTL)) %>%
          sapply(function(x) 1 - qvalue::qvalue_truncp(intersect_data$P[x])$pi0),
        
        lapply(index_mQTL, function(x) intersect(x, index_eQTL)) %>%
          sapply(function(x) 1 - qvalue::qvalue_truncp(intersect_data$P[x])$pi0),
        
        lapply(index_mQTL, function(x) setdiff(x, index_eQTL)) %>%
          sapply(function(x) mean(intersect_data$P[x] < 6e-7)),
        
        lapply(index_mQTL, function(x) intersect(x, index_eQTL)) %>%
          sapply(function(x) mean(intersect_data$P[x] < 6e-7)))
  
  mean(intersect_data$P[setdiff(index_mQTL %>% unlist() %>% unique(),
                                index_eQTL %>% unlist() %>% unique())] < 6e-7)
  mean(intersect_data$P[intersect(index_mQTL %>% unlist() %>% unique(),
                                  index_eQTL %>% unlist() %>% unique())] < 6e-7)
  # first column: all
  index_input_mQTL <- apply(result[[4]][intersect_data$ind,], 2, function(x) which(x < 1e-5))
  index_input_multi <- which(rowSums(result[[4]][intersect_data$ind,]<1e-5)>1)
  index_input_mQTL <- lapply(index_input_mQTL, function(x) setdiff(x,index_input_multi))
  
  index_XING_mQTL <- apply(subPP_m, 2, function(x) which(x > 0.8))
  index_multi <- which(rowSums(subPP_m > 0.8)>1)
  index_XING_specific_mQTL <- lapply(index_XING_mQTL, function(x) setdiff(x,index_multi)) %>%
    lapply(function(x) intersect(x,index_eQTL))
  index_XING_multi_mQTL <- lapply(index_XING_mQTL, function(x) intersect(x,index_multi))
  index_mash <- apply(result[[5]][intersect_data$ind,], 2, function(x) which(x < 0.05))
  
  single_tissue_single_omic <- lapply(index_XING_mQTL, function(x) setdiff(x,index_multi)) %>%
    lapply(function(x) setdiff(x,index_eQTL))
  single_tissue_multi_omic <- lapply(index_XING_mQTL, function(x) setdiff(x,index_multi)) %>%
    lapply(function(x) intersect(x,index_eQTL))
  multi_tissue_single_omic <- lapply(index_XING_mQTL, function(x) intersect(x,index_multi)) %>%
    lapply(function(x) setdiff(x,index_eQTL))
  multi_tissue_multi_omic <- lapply(index_XING_mQTL, function(x) intersect(x,index_multi)) %>%
    lapply(function(x) intersect(x,index_eQTL))
  
  rbind(
    sapply(single_tissue_single_omic, function(x) mean(intersect_data$P[x] < 6e-7)),
    sapply(single_tissue_multi_omic, function(x) mean(intersect_data$P[x] < 6e-7)),
    sapply(multi_tissue_single_omic, function(x) mean(intersect_data$P[x] < 6e-7)),
    sapply(multi_tissue_multi_omic, function(x) mean(intersect_data$P[x] < 6e-7))
  ) %>%
    apply(1,paste0,collapse=",")
  
  paste0("(",
         rbind(
           sapply(single_tissue_single_omic, function(x) (sum(intersect_data$P[x] < 6e-7))),
           sapply(single_tissue_multi_omic, function(x) (sum(intersect_data$P[x] < 6e-7))),
           sapply(multi_tissue_single_omic, function(x) (sum(intersect_data$P[x] < 6e-7))),
           sapply(multi_tissue_multi_omic, function(x) (sum(intersect_data$P[x] < 6e-7)))
         ),"/",
         rbind(
           sapply(single_tissue_single_omic, function(x) length(x)),
           sapply(single_tissue_multi_omic, function(x) length(x)),
           sapply(multi_tissue_single_omic, function(x) length(x)),
           sapply(multi_tissue_multi_omic, function(x) length(x))
         ),")"
  ) %>% matrix(nrow = 4) %>% apply(1,paste0,collapse=",")
  
  length(which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) >= 2))
  
  sum(intersect_data$P[which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) >= 0)]<6e-7)
  mean(intersect_data$P[which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) >= 0)]<6e-7)
  mean(intersect_data$P[which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) >= 1)]<6e-7)
  mean(intersect_data$P[which(rowSums(result$alphajkm2_dat2[intersect_data$ind,]>0.8) >= 2)]<6e-7)
  mean(intersect_data$P[which(rowSums(result[[5]][intersect_data$ind,]<0.05) >= 2)]<6e-7)
  
  setdiff(index_input_mQTL[[1]],index_XING_mQTL[[1]]) %>% length()
  setdiff(index_XING_mQTL[[1]],index_input_mQTL[[1]]) %>% length()
  
  tissue <- "MuscleSkeletal"
  top <- fread(paste0("/scratch/t.phs.yihaolu/map_mQTL/",tissue,".regular.perm.fdr.txt"))
  sec <- fread(paste0("/scratch/t.phs.yihaolu/map_mQTL/",tissue,".mQTLs.secondary.txt"))
  sec$V1 <- str_replace(sec$V1,":1","")
  intersect_data[rowSums(subPP_m>0.8)>=1,] %>%
    inner_join(top, by = c("variant_id"="variant_id","cpg"="cpg_id")) %>%
    summarise(mean(P<6e-7),sum(P<6e-7))
  intersect_data[rowSums(subPP_m>0.8)>=1,] %>%
    inner_join(sec, by = c("variant_id"="V2","cpg"="V1")) %>%
    summarise(mean(P<6e-7),sum(P<6e-7))
  
  info %>%
    inner_join(top, by = c("variant_id"="variant_id","cpg"="cpg_id"))
  info %>%
    inner_join(sec, by = c("variant_id"="V2","cpg"="V1"))
}

cd /scratch/t.phs.yihaolu/XING_revision/map_eQTL/script
echo "
#PBS -l mem=15gb
cd /scratch/t.phs.yihaolu/XING_revision/
Rscript run_cis.R" > run_cis.sh
qsub run_cis.sh

###########################################################################################################################################################
## mQTL REPLICATION - GoDMC
###########################################################################################################################################################
args <- commandArgs(TRUE)
file = as.numeric(args[1])

library(data.table)
library(tidyverse)
library(pbapply)
setwd("/scratch/t.phs.yihaolu/XING_revision/")
home_path <- "/home/yihaolu/XING/"
external_path <- "/scratch/t.phs.yihaolu/external_data/"
gene_info <- fread(paste0(home_path,"gene_info_GTEx_v8.txt")) %>% select(gene_id, Chr, gene_symbol)
gene_info$gene_id <- str_split(gene_info$gene_id,"[.]") %>% lapply(function(x) x[1]) %>% unlist()

load(paste0("result/Upload/lead_mQTL_eQTL_",file,".RData"))
intersect_data <- list()
info <- result[[7]] %>% ungroup() %>%
  mutate(ind = 1:n())
# for (chr in 1:22) {
#   print(chr)
#   GoDMC <- fread(paste0(external_path, "external_mQTLs/GoDMC_chr",chr,".txt"))
#   GoDMC <- GoDMC %>% select(cpg,pval,rsid)
#   colnames(GoDMC) <- c("CPG","P","SNP")
#   fwrite(GoDMC, file = paste0(external_path,"external_mQTLs/GoDMC_simple_chr",chr,".txt"))
# }
GoDMC <- list()
for (chr in 1:22) {
  print(chr)
  GoDMC[[chr]] <- fread(paste0(external_path, "external_mQTLs/GoDMC_simple_chr",chr,".txt"))
}
GoDMC <- do.call("rbind",GoDMC)

GoDMC <- GoDMC[which(GoDMC$CPG %in% info$cpg),]
GoDMC <- GoDMC[which(GoDMC$SNP %in% info$rs_id_dbSNP151_GRCh38p7),]

GoDMC <- GoDMC %>% filter(P<6e-7)

subPP_e <- result$alphajkm2_dat1[info$ind,]
subPP_m <- result$alphajkm2_dat2[info$ind,]
index_XING_mQTL <- apply(subPP_m, 2, function(x) which(x > 0.8))
index_multi <- which(rowSums(subPP_m > 0.8)>1)
index_eQTL <- which(rowSums(subPP_e>0.8)>=1)

single_tissue_single_omic <- lapply(index_XING_mQTL, function(x) setdiff(x,index_multi)) %>%
  lapply(function(x) setdiff(x,index_eQTL))
single_tissue_multi_omic <- lapply(index_XING_mQTL, function(x) setdiff(x,index_multi)) %>%
  lapply(function(x) intersect(x,index_eQTL))
multi_tissue_single_omic <- lapply(index_XING_mQTL, function(x) intersect(x,index_multi)) %>%
  lapply(function(x) setdiff(x,index_eQTL))
multi_tissue_multi_omic <- lapply(index_XING_mQTL, function(x) intersect(x,index_multi)) %>%
  lapply(function(x) intersect(x,index_eQTL))

intersect_data <- info %>%
  inner_join(GoDMC, by = c("rs_id_dbSNP151_GRCh38p7" = "SNP", "cpg" = "CPG"))

nrow(intersect_data)/nrow(info)
nrow(intersect_data[rowSums(subPP_m[intersect_data$n,]>0.8)>=2,])/nrow(info[rowSums(subPP_m>0.8)>=2,])

(sum(intersect_data$ind %in% (single_tissue_single_omic %>% unlist() %>% unique()))/
    length(single_tissue_single_omic %>% unlist() %>% unique()))
(sum(intersect_data$ind %in% (single_tissue_multi_omic %>% unlist() %>% unique()))/
    length(single_tissue_multi_omic %>% unlist() %>% unique()))
(sum(intersect_data$ind %in% (multi_tissue_single_omic %>% unlist() %>% unique()))/
    length(multi_tissue_single_omic %>% unlist() %>% unique()))
(sum(intersect_data$ind %in% (multi_tissue_multi_omic %>% unlist() %>% unique()))/
    length(multi_tissue_multi_omic %>% unlist() %>% unique()))

rbind(
  sapply(single_tissue_single_omic, function(x) (sum(intersect_data$ind %in% x)/length(x))),
  sapply(single_tissue_multi_omic, function(x) (sum(intersect_data$ind %in% x)/length(x))),
  sapply(multi_tissue_single_omic, function(x) (sum(intersect_data$ind %in% x)/length(x))),
  sapply(multi_tissue_multi_omic, function(x) (sum(intersect_data$ind %in% x)/length(x)))
) %>% apply(1,paste0,collapse=",")

paste0("(",
       rbind(
         sapply(single_tissue_single_omic, function(x) (sum(intersect_data$ind %in% x))),
         sapply(single_tissue_multi_omic, function(x) (sum(intersect_data$ind %in% x))),
         sapply(multi_tissue_single_omic, function(x) (sum(intersect_data$ind %in% x))),
         sapply(multi_tissue_multi_omic, function(x) (sum(intersect_data$ind %in% x)))
       ),"/",
       rbind(
         sapply(single_tissue_single_omic, function(x) length(x)),
         sapply(single_tissue_multi_omic, function(x) length(x)),
         sapply(multi_tissue_single_omic, function(x) length(x)),
         sapply(multi_tissue_multi_omic, function(x) length(x))
       ),")"
) %>% matrix(nrow = 4) %>% apply(1,paste0,collapse=",")

fwrite(intersect_data, file = paste0("result/intersect_data_chr",chr,"_file",file,".txt"))

cd /scratch/t.phs.yihaolu/XING_revision/map_eQTL/script
for file in 2; do
for chr in {1..22}; do echo "
#PBS -l mem=15gb
cd /scratch/t.phs.yihaolu/XING_revision/
Rscript replicate_mQTL.R ${file} $chr" > simplify_${file}_${chr}.sh
qsub simplify_${file}_${chr}.sh; done; done
