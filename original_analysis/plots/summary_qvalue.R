#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

## ========== User params ==========
FDR_TARGET   <- 0.10
BONF_P_CUT   <- 5e-8
P_FLOOR      <- 1e-300
M_KNOCK      <- 5
## =================================

## -------- q-value 计算函数 --------
MK.q.byStat <- function (kappa, tau, M, Rej.Bound = 10000){
  kappa <- as.numeric(kappa); tau <- as.numeric(tau)
  b  <- order(tau, decreasing = TRUE)
  c0 <- (kappa[b] == 0)

  ratio <- numeric(0); temp0 <- 0
  for(i in seq_along(b)){
    temp0 <- temp0 + c0[i]
    temp1 <- i - temp0
    temp_ratio <- (1/M + 1/M * temp1) / max(1, temp0)
    ratio <- c(ratio, temp_ratio)
    if(i > Rej.Bound) break
  }

  q <- rep(1, length(tau))
  index_bound <- suppressWarnings(max(which(tau[b] > 0)))
  if(!is.finite(index_bound)) index_bound <- 0

  for(i in seq_along(b)){
    temp.index <- i:min(length(b), Rej.Bound, index_bound)
    if(length(temp.index) == 0) next
    q[b[i]] <- min(ratio[temp.index]) * c0[i] + 1 - c0[i]
    if(i > Rej.Bound) break
  }
  return(q)
}
## ----------------------------------

read_all_q <- function(dir_q){
  fs <- list.files(dir_q, pattern = "_W_basic\\.csv$", full.names = TRUE)
  if(!length(fs)) stop("No *_W_basic.csv in ", dir_q)

  rbindlist(lapply(fs, function(f){
    dt <- fread(f)
    need <- c("chr","pos","pG0","W","kappa","tau")
    miss <- setdiff(need, names(dt))
    if(length(miss)){
      stop("Missing columns in ", basename(f), ": ", paste(miss, collapse = ", "))
    }
    dt <- dt[, ..need]
    dt[, `:=`(
      chr   = as.character(chr),
      pos   = as.numeric(pos),
      pG0   = as.numeric(pG0),
      W     = as.numeric(W),
      kappa = as.numeric(kappa),
      tau   = as.numeric(tau)
    )]
    dt[]
  }), fill = TRUE)
}

run_trait <- function(trait, dir_q, out_dir,
                      fdr = FDR_TARGET, p_bonf = BONF_P_CUT){

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ## 读取
  dt <- read_all_q(dir_q)

  ## 统一为 1..22 的整数，并排序
  dt[, chr := as.integer(gsub("^chr", "", chr))]
  dt <- dt[chr %in% 1:22]
  setorder(dt, chr, pos)

  ## 全基因组上重算 q
  dt[, Qvalue := MK.q.byStat(kappa, tau, M = M_KNOCK)]
  dt[, mlp := -log10(pmax(pG0, P_FLOOR))]

  ## Knockoff (全局 q) 显著
  sel_knock <- dt$Qvalue <= fdr
  W_thr     <- if(any(sel_knock)) min(dt$W[sel_knock], na.rm = TRUE) else NA_real_

  ## Bonferroni 显著（固定 5e-8）
  m_total      <- nrow(dt)
  p_cut        <- p_bonf
  p_cut_log10  <- -log10(p_cut)
  sel_bonf     <- dt$pG0 <= p_cut

  ## 标记并保存 SNP 级别表
  dt[, `:=`(
    is_knockoff = sel_knock,
    is_bonf     = sel_bonf
  )]

  fwrite(dt[, .(chr, pos, pG0, mlp, W, kappa, tau, Qvalue, is_knockoff, is_bonf)],
         file.path(out_dir, sprintf("%s_global_snp_level.tsv", trait)),
         sep = "\t")

  ## 计数与交集
  n_knock    <- sum(sel_knock)
  n_bonf     <- sum(sel_bonf)
  n_overlap  <- sum(sel_knock & sel_bonf)

  summary_dt <- data.table(
    trait        = trait,
    FDR_target   = fdr,
    W_thr        = W_thr,
    n_knockoff   = n_knock,
    n_bonf       = n_bonf,
    n_overlap    = n_overlap,
    m_total      = m_total,
    bonf_p_cut   = p_cut,
    bonf_log10   = p_cut_log10
  )

  fwrite(summary_dt,
         file.path(out_dir, sprintf("%s_counts_knockoff_vs_bonf.tsv", trait)),
         sep = "\t")

  cat(sprintf("\n[%s]\n  Knockoff (FDR=%.2f): %d SNPs (W_thr=%.3f)\n  Bonferroni (p<=5e-8): %d SNPs\n  Overlap: %d\n  m=%d, -log10(5e-8)=%.3f\n",
              trait, fdr, n_knock, W_thr, n_bonf, n_overlap, m_total, p_cut_log10))
}

## ================== RUN ==================
run_trait(
  trait   = "Hypertension",
  dir_q   = "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/score/hyper",
  out_dir = "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/hyper_counts"
)

run_trait(
  trait   = "Type 2 diabetes",
  dir_q   = "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/score/dm",
  out_dir = "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/dm_counts"
)
