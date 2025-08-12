#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

## ================== 参数 ==================
TSS_FILE <- "/dssg/home/acct-mashiyang1991/dhhhh020606/.vep/homo_sapiens/105_GRCh38/transcript_gene_tss.txt"

IN_HYP_KNOCK <- "/dssg/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/hyper_counts/Hypertension_knockoff_sig.tsv"
IN_HYP_BONF  <- "/dssg/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/hyper_counts/Hypertension_bonferroni_sig.tsv"
IN_DM_KNOCK  <- "/dssg/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/dm_counts/Type_2_diabetes_knockoff_sig.tsv"
IN_DM_BONF   <- "/dssg/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/dm_counts/Type_2_diabetes_bonferroni_sig.tsv"

OUT_HYP_KNOCK <- sub("\\.tsv$", "_annot.tsv", IN_HYP_KNOCK)
OUT_HYP_BONF  <- sub("\\.tsv$", "_annot.tsv", IN_HYP_BONF)
OUT_DM_KNOCK  <- sub("\\.tsv$", "_annot.tsv", IN_DM_KNOCK)
OUT_DM_BONF   <- sub("\\.tsv$", "_annot.tsv", IN_DM_BONF)
## =========================================

## 1) 读 TSS 文件
read_tss <- function(tss_file) {
  tss <- fread(tss_file, header = FALSE,
               col.names = c("chr", "tss", "transcript_id", "gene_id", "gene_symbol"))
  # 转成正确类型，并去掉 NA
  tss[, chr := as.character(chr)]
  tss[, tss := as.numeric(tss)]
  tss <- tss[!is.na(chr) & !is.na(tss)]
  setkey(tss, chr, tss)
  tss[]
}

## 2) 对一个 SNP 表注释最近 TSS
annotate_nearest_tss <- function(snp_file, out_file, tss_dt) {
  if (!file.exists(snp_file)) {
    warning("File missing: ", snp_file)
    return(NULL)
  }
  snp <- fread(snp_file)
  # 类型转换
  snp[, chr := as.character(chr)]
  snp[, pos := as.numeric(pos)]
  
  # 按染色体分组，找到最近的 tss
  ann <- snp[, {
    this_tss <- tss_dt[chr == .BY$chr]
    if (nrow(this_tss)==0) {
      # 如果该 chr 没 TSS，就补 NA
      .(transcript_id=NA, gene_id=NA, gene_symbol=NA, tss=NA_real_, dist_bp=NA_integer_)
    } else {
      # 对每个 pos 找最近的 tss
      idx <- which.min(abs(this_tss$tss - pos))
      res <- this_tss[idx]
      .(transcript_id=res$transcript_id,
        gene_id=res$gene_id,
        gene_symbol=res$gene_symbol,
        tss=res$tss,
        dist_bp=abs(res$tss - pos))
    }
  }, by=.(chr, pos, pG0, mlp, W, kappa, tau, Qvalue)]
  
  # 输出：保留原来的 SNP 列 + 注释列
  out <- cbind(
    ann[, .(chr,pos,pG0,mlp,W,kappa,tau,Qvalue)],
    ann[, .(transcript_id,gene_id,gene_symbol,tss,dist_bp)]
  )
  fwrite(out, out_file, sep="\t")
  cat(sprintf("✓ %s -> %s  (n=%d)\n",
              basename(snp_file), basename(out_file), nrow(out)))
  invisible(out)
}

## ============ 主流程 ============
tss_dt <- read_tss(TSS_FILE)

annotate_nearest_tss(IN_HYP_KNOCK, OUT_HYP_KNOCK, tss_dt)
annotate_nearest_tss(IN_HYP_BONF,  OUT_HYP_BONF,  tss_dt)
annotate_nearest_tss(IN_DM_KNOCK,  OUT_DM_KNOCK,  tss_dt)
annotate_nearest_tss(IN_DM_BONF,   OUT_DM_BONF,   tss_dt)
