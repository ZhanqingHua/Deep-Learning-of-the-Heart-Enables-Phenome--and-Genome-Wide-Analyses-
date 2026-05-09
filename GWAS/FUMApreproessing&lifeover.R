
library(data.table)

input_file <- "/data/circ_ukb/AHua/project/GWAS/regenie/height/metaanalysis/height_X7_Elongation_all_metal1.txt"
output_file <- "/data/circ_ukb/AHua/project/GWAS/regenie/height/FUMA/height_X7_Elongation_FUMA_ready_1.txt"

df <- fread(input_file)

df <- df[complete.cases(df), ]

valid_alleles <- c("A", "C", "G", "T", "a", "c", "g", "t")
df <- df[df$Allele1 %in% valid_alleles & df$Allele2 %in% valid_alleles, ]

df$A1 <- toupper(df$Allele1)
df$A2 <- toupper(df$Allele2)

df[, CHR := toupper(sub("^chr", "", tstrsplit(MarkerName, ":", fixed=TRUE)[[1]]))]
df[, BP := as.integer(tstrsplit(MarkerName, ":", fixed=TRUE)[[2]])]


df <- df[CHR %in% c(as.character(1:22), "X","Y")]

df_fuma <- data.frame(
  SNP = df$MarkerName,
  CHR = df$CHR,
  POS = df$BP,
  A1 = df$A1,
  A2 = df$A2,
  P = df$`P-value`,
  BETA = df$Effect,
  SE = df$StdErr,
  N = df$TotalSampleSize
)


df_fuma$P <- as.numeric(df_fuma$P)


library(data.table)
setDT(df_fuma)  

df_sig <- df_fuma[P < 0.05]
df_nonsig <- df_fuma[P >= 0.05]


set.seed(42)  
df_nonsig_sampled <- df_nonsig[sample(.N, size = floor(0.3 * .N))]


df_combined <- rbind(df_sig, df_nonsig_sampled)

）
df <- copy(df_combined)
df$CHR <- gsub("chr", "", df$CHR)                   ）
df$start <- as.integer(df$POS) - 1                  
df$end <- as.integer(df$POS)
df$CHR <- paste0("chr", df$CHR)                     


df_bed <- df[, .(CHR, start, end)]

df_bed[, start := format(start, scientific = FALSE)]
df_bed[, end := format(end, scientific = FALSE)]

）
fwrite(df_bed, "/data/circ_ukb/AHua/project/GWAS/regenie/height/FUMA/X7_SurfaceVolumeRatio.bed", sep = "\t", col.names = FALSE)




gwas <- df_combined

gwas[, SNP_key := paste0("chr", CHR, ":", POS)]

lifted <- fread("/data/circ_ukb/AHua/project/FUMA/X7_SurfaceVolumeRatio_changed.bed", 
                col.names = c("chr", "start", "end", "SNP_key", "status"))
lifted[, match_key := tstrsplit(SNP_key, "-", fixed = TRUE)[[1]]]
head(lifted$match_key)


merged <- merge(
  gwas,
  lifted[, .(match_key, new_chr = chr, new_pos = end)],
  by.x = "SNP_key", by.y = "match_key",
  all.x = FALSE
)

merged[, CHR := gsub("chr", "", new_chr)]
merged[, POS := new_pos]
merged[, c("new_chr", "new_pos") := NULL]


head(merged)
nrow(merged)

fwrite(merged, "/data/circ_ukb/AHua/project/FUMA/X7_SurfaceVolumeRatio_fuma_ready_hg19.txt", sep = "\t", quote = FALSE)
R.utils::gzip("/data/circ_ukb/AHua/project/FUMA/X7_SurfaceVolumeRatio_fuma_ready_hg19.txt", overwrite = TRUE)