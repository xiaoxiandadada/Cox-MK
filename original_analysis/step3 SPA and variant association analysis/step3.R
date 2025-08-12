
suppressPackageStartupMessages({
  library(gdsfmt)
  library(data.table)
  source("/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/KnockoffScreen.r")
  source("/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/SPACox.R")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript step3.R <GDS_FILE> <OUTDIR>")
}
gds_file <- args[1]
outdir <- args[2]

# 读取
obj.null <- readRDS("/lustre/home/acct-mashiyang1991/dhhhh020606/tzh2/null_hyper.rds")
gds <- openfn.gds(gds_file)

# 提取基础信息
samples <- read.gdsn(index.gdsn(gds, "sample.id"))
snps <- read.gdsn(index.gdsn(gds, "snp.pos"))
out_name <- sub("\\.gds$", "_p.csv", basename(gds$filename))

# 计算origin的MAF
geno_original <- read.gdsn(index.gdsn(gds, "original"))
rownames(geno_original) <- samples
colnames(geno_original) <- snps
maf_values <- SPACox(obj.null, geno_original)[, "MAF"]

# 处理所有数据集,保存p值
res_all <- matrix(nrow=length(snps), ncol=6)
colnames(res_all) <- c("pG0", paste0("pGk",1:5))
datasets <- c("original", paste0("knockoff",1:5))

for(i in seq_along(datasets)) {
  geno <- read.gdsn(index.gdsn(gds, datasets[i]))
  rownames(geno) <- samples
  colnames(geno) <- snps
  res_all[,i] <- SPACox(obj.null, geno)[, "p.value.spa"]
}

# 输出
fwrite(data.frame(SNP=snps, MAF=maf_values, res_all),
       file.path(outdir, out_name))

closefn.gds(gds)
