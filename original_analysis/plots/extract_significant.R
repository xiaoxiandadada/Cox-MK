#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

## --------- user params ----------
FDR_TARGET <- 0.10
BONF_P_CUT <- 5e-8

# 你的两个结果目录
IN_DIR_HYP <- "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/hyper_counts"
IN_DIR_DM  <- "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/dm_counts"
## --------------------------------

traits <- list(
  Hypertension      = file.path(IN_DIR_HYP, "Hypertension_global_snp_level.tsv"),
  "Type 2 diabetes" = file.path(IN_DIR_DM,  "Type 2 diabetes_global_snp_level.tsv")
)

export_sig <- function(trait, f){
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(invisible(NULL))
  }

  out_dir <- dirname(f)
  tag     <- gsub("[^A-Za-z0-9]+", "_", trait)  # 文件名安全

  dt <- fread(f)

  # 统一 chr 顺序并排序
  if (is.character(dt$chr)) {
    dt[, chr := as.integer(gsub("^chr", "", chr))]
  }
  dt <- dt[chr %in% 1:22]
  setorder(dt, chr, pos)

  # 选择要输出的列
  keep_cols <- c("chr","pos","pG0","mlp","W","kappa","tau","Qvalue")

  # Knockoff FDR
  knock <- dt[Qvalue <= FDR_TARGET, ..keep_cols]
  fwrite(knock,
         file.path(out_dir, sprintf("%s_knockoff_sig.tsv", tag)),
         sep = "\t")

  # Bonferroni
  bonf <- dt[pG0 <= BONF_P_CUT, ..keep_cols]
  fwrite(bonf,
         file.path(out_dir, sprintf("%s_bonferroni_sig.tsv", tag)),
         sep = "\t")

  cat(sprintf("[%s] knockoff=%d, bonferroni=%d  -> written to %s\n",
              trait, nrow(knock), nrow(bonf), out_dir))
}

invisible(mapply(export_sig, names(traits), traits))
