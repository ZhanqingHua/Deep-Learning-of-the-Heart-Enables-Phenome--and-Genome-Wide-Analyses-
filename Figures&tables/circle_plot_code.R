
library(data.table)

vcf <- fread("/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas/categorical-20002-1075_hg38.vcf", skip = "#CHROM")

setnames(vcf, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))

vcf[, PVAL := sub(".*PVAL=([^;]+).*", "\\1", INFO)]

vcf_out <- vcf[, .(
  ID = paste0(CHROM, ":", POS, ":", REF, ":", ALT),
  chr = gsub("chr", "", CHROM),
  pos = POS,
  ref = REF,
  alt = ALT,
  PVAL = as.numeric(PVAL)
)]

fwrite(vcf_out, "/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas/categorical-20002-1075_hg38_withID.tsv", sep="\t")




library(data.table)
library(CMplot)

wdir <- "/data/circ_ukb/AHua/project/GWAS/regenie/height/metaanalysis"
cad_file <- "/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas/categorical-20002-1075_hg38_withID.tsv"

features <- c("X47_MajorAxisLength", "X47_Sphericity", "X47_Maximum2DDiameterColumn","X47_SurfaceArea", "X47_SurfaceVolumeRatio","X47_Flatness")

overall_df <- NULL
for (feat in features) {
  f <- paste0(wdir, "/height_", feat, "_all_metal1.txt")
  df <- fread(f, header=TRUE)
  
  
  df[, c("chr","pos","ref","alt") := tstrsplit(MarkerName, ":", fixed=TRUE)]
  df[, chr := as.integer(gsub("chr","", chr))]
  df[, pos := as.integer(pos)]
  
  df_sub <- df[, .(SNP = MarkerName, Chromosome = chr, Position = pos, P = `P-value`)]
  setnames(df_sub, "P", feat)
  
  if (is.null(overall_df)) {
    overall_df <- df_sub
  } else {
    overall_df <- merge(overall_df, df_sub, by=c("SNP","Chromosome","Position"))
  }
}

cat("Merged features:", ncol(overall_df)-3, "\n")



cad <- fread(cad_file, header=TRUE)
cad_sub <- cad[, .(
  SNP = ID,
  Chromosome = as.integer(chr),
  Position = pos,
  CAD = as.numeric(PVAL)
)]


cad_sub[, Chromosome := as.integer(Chromosome)]
overall_df[, Chromosome := as.integer(Chromosome)]

df_cad <- merge(overall_df, cad_sub, by=c("SNP","Chromosome","Position"))

cat("Merged columns:", colnames(df_cad), "\n")


setwd("/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas")
CMplot(df_cad, type="p", plot.type="c", col=c("#ADD8E6"), chr.labels=paste("Chr", 1:22, sep=""), r=0.4, outward=TRUE, cir.chr.h=1.3, chr.den.col="black", threshold=5e-7, ylab=c("X47_MajorAxisLength", "X47_Sphericity", "X47_Maximum2DDiameterColumn","X47_SurfaceArea", "X47_SurfaceVolumeRatio","X47_Flatness","CAD"), file="pdf", dpi=300, file.output=TRUE, verbose=TRUE, width=10, height=10, amplify=TRUE, signal.col="red", signal.cex=1, file.name="height_features_cad_circle_x47")














library(data.table)
library(CMplot)

# 路径
wdir <- "/data/circ_ukb/AHua/project/GWAS/regenie/height/metaanalysis"
cad_file <- "/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas/categorical-20002-1075_hg38_withID.tsv"

# features 列表
features <- c("X45_LeastAxisLength", "X45_Sphericity", "X45_Elongation","X45_Flatness", "X45_Maximum2DDiameterColumn")

# === 1. 读入所有 height feature GWAS 文件并 merge ===
overall_df <- NULL
for (feat in features) {
  f <- paste0(wdir, "/height_", feat, "_all_metal1.txt")
  df <- fread(f, header=TRUE)
  
  # 假设文件里有列名: MarkerName, P-value
  df[, c("chr","pos","ref","alt") := tstrsplit(MarkerName, ":", fixed=TRUE)]
  df[, chr := as.integer(gsub("chr","", chr))]
  df[, pos := as.integer(pos)]
  
  df_sub <- df[, .(SNP = MarkerName, Chromosome = chr, Position = pos, P = `P-value`)]
  setnames(df_sub, "P", feat)
  
  if (is.null(overall_df)) {
    overall_df <- df_sub
  } else {
    overall_df <- merge(overall_df, df_sub, by=c("SNP","Chromosome","Position"))
  }
}

cat("Merged features:", ncol(overall_df)-3, "\n")



# CAD GWAS
cad <- fread(cad_file, header=TRUE)

# 构造 CAD 子表
cad_sub <- cad[, .(
  SNP = ID,
  Chromosome = as.integer(chr),
  Position = pos,
  CAD = as.numeric(PVAL)
)]


cad_sub[, Chromosome := as.integer(Chromosome)]
overall_df[, Chromosome := as.integer(Chromosome)]


df_cad <- merge(overall_df, cad_sub, by=c("SNP","Chromosome","Position"))

cat("Merged columns:", colnames(df_cad), "\n")


setwd("/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas")
CMplot(df_cad, type="p", plot.type="c", col=c("#ADD8E6"), chr.labels=paste("Chr", 1:22, sep=""), r=0.4, outward=TRUE, cir.chr.h=1.3, chr.den.col="black", threshold=5e-7, ylab=c("X45_LeastAxisLength", "X45_Sphericity", "X45_Elongation","X45_Flatness", "X45_Maximum2DDiameterColumn","CAD"), file="pdf", dpi=300, file.output=TRUE, verbose=TRUE, width=10, height=10, amplify=TRUE, signal.col="red", signal.cex=1, file.name="height_features_cad_circle")






library(data.table)

# 读 VCF，忽略掉 header 之前的 "##" 行
vcf <- fread("/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas/categorical-20002-1075_hg38.vcf", skip = "#CHROM")

# 给列起名字
setnames(vcf, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))

# 从 INFO 里提取 PVAL
vcf[, PVAL := sub(".*PVAL=([^;]+).*", "\\1", INFO)]

# 构造 ID 和整理成 TSV 样式
vcf_out <- vcf[, .(
  ID = paste0(CHROM, ":", POS, ":", REF, ":", ALT),
  chr = gsub("chr", "", CHROM),
  pos = POS,
  ref = REF,
  alt = ALT,
  PVAL = as.numeric(PVAL)
)]

# 保存成 TSV
fwrite(vcf_out, "/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas/categorical-20002-1075_hg38_withID.tsv", sep="\t")



library(data.table)
library(CMplot)

# 路径
wdir <- "/data/circ_ukb/AHua/project/GWAS/regenie/height/metaanalysis"
cad_file <- "/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas/categorical-20002-1075_hg38_withID.tsv"

# features 列表
features <- c("X45_LeastAxisLength", "X45_Sphericity", "X45_Elongation","X45_Flatness", "X45_Maximum2DDiameterColumn")

# === 1. 读入所有 height feature GWAS 文件并 merge ===
overall_df <- NULL
for (feat in features) {
  f <- paste0(wdir, "/height_", feat, "_all_metal1.txt")
  df <- fread(f, header=TRUE)
  
  # 假设文件里有列名: MarkerName, P-value
  df[, c("chr","pos","ref","alt") := tstrsplit(MarkerName, ":", fixed=TRUE)]
  df[, chr := as.integer(gsub("chr","", chr))]
  df[, pos := as.integer(pos)]
  
  df_sub <- df[, .(SNP = MarkerName, Chromosome = chr, Position = pos, P = `P-value`)]
  setnames(df_sub, "P", feat)
  
  if (is.null(overall_df)) {
    overall_df <- df_sub
  } else {
    overall_df <- merge(overall_df, df_sub, by=c("SNP","Chromosome","Position"))
  }
}

cat("Merged features:", ncol(overall_df)-3, "\n")



# CAD GWAS
cad <- fread(cad_file, header=TRUE)

# CAD construction
cad_sub <- cad[, .(
  SNP = ID,
  Chromosome = as.integer(chr),
  Position = pos,
  CAD = as.numeric(PVAL)
)]


cad_sub[, Chromosome := as.integer(Chromosome)]
overall_df[, Chromosome := as.integer(Chromosome)]


df_cad <- merge(overall_df, cad_sub, by=c("SNP","Chromosome","Position"))

cat("Merged columns:", colnames(df_cad), "\n")


setwd("/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas")
CMplot(df_cad, type="p", plot.type="c", col=c("#ADD8E6"), chr.labels=paste("Chr", 1:22, sep=""), r=0.4, outward=TRUE, cir.chr.h=1.3, chr.den.col="black", threshold=5e-7, ylab=c("X7_LeastAxisLength","X7_Maximum2DDiameterSlice","X7_SurfaceVolumeRatio","X7_Elongation","X7_MinorAxisLength","CAD"), file="pdf", dpi=300, file.output=TRUE, verbose=TRUE, width=10, height=10, amplify=TRUE, signal.col="red", signal.cex=1, file.name="height_features_cad_circle")
