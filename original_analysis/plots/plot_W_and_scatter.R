#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

## --------- params ----------
ALPHA_FDR    <- 0.10
ALPHA_BONF   <- 0.10
P_FLOOR      <- 1e-300
COL_RED      <- "#e41a1c"
COL_BLUE_D   <- "#1f78b4"   # 深蓝
COL_BLUE_L   <- "#a6cee3"   # 浅蓝
PT_SIZE_RED  <- 1.6
PT_SIZE_BLUE <- 1.0
M_KNOCK      <- 5
## ---------------------------

source("KnockoffScreen.r")   # MK.q.byStat()

read_all_q <- function(dir_q){
  fs <- list.files(dir_q, pattern = "_knockoff_q\\.csv$", full.names = TRUE)
  if(!length(fs)) stop("No *_knockoff_q.csv in ", dir_q)
  rbindlist(lapply(fs, function(f){
    m  <- str_match(basename(f), "chr(\\d+)_(\\d+)_(\\d+)_knockoff_q\\.csv")
    dt <- fread(f, select = c("SNP","pG0","W","Qvalue","kappa","tau"))
    setnames(dt, "SNP", "pos")
    dt[, `:=`(
      chr   = as.integer(m[2]),
      start = as.integer(m[3]),
      end   = as.integer(m[4]),
      pos    = as.numeric(pos),
      pG0    = as.numeric(pG0),
      W      = as.numeric(W),
      Qvalue = as.numeric(Qvalue),
      kappa  = as.numeric(kappa),
      tau    = as.numeric(tau)
    )]
    dt
  }), fill = TRUE)
}

plot_alpha010 <- function(trait, dir_q, out_dir,
                          alpha_fdr = ALPHA_FDR, alpha_bonf = ALPHA_BONF){
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  dt <- read_all_q(dir_q)
  dt[, Qvalue_global := MK.q.byStat(kappa, tau, M = M_KNOCK)]

  sel_q <- dt$Qvalue_global <= alpha_fdr
  W_thr <- if(any(sel_q)) min(dt$W[sel_q], na.rm = TRUE) else NA_real_
  dt[, sel_knock := if (is.finite(W_thr)) W >= W_thr else FALSE]

  m_total     <- nrow(dt)
  p_cut_log10 <- -log10(alpha_bonf / m_total)

  cat(sprintf("[%s] FDR=%.2f  W_thr=%.3f  #sel=%d/%d\n",
              trait, alpha_fdr, W_thr, sum(sel_q), m_total))

  dt[, mlp := -log10(pmax(pG0, P_FLOOR))]
  dt <- dt[is.finite(W) & is.finite(mlp) & is.finite(pos)]

  ## ---- build cumulative chromosome position
  chr_len <- dt[, .(chr_len = max(pos, na.rm = TRUE)), by = chr][order(chr)]
  chr_len[, offset := shift(cumsum(as.numeric(chr_len)), fill = 0)]
  dt <- merge(dt, chr_len[, .(chr, offset)], by = "chr")
  dt[, pos_cum := as.numeric(pos) + as.numeric(offset)]
  dt <- dt[is.finite(pos_cum)]
  chr_ticks <- chr_len[, .(chr, center = as.numeric(offset) + as.numeric(chr_len)/2)]

  ## 奇偶染色体交替蓝色（只用于非显著点）
  dt[, chr_parity := factor(chr %% 2, levels = c(0,1), labels = c("even","odd"))]

  ## 1) Manhattan: 底层用交替蓝色（不显著），阈值上方覆盖红色
  subtitle1 <- if (is.finite(W_thr)) sprintf("\u2015 red  W = %.3f", W_thr) else "no discoveries"

  p1 <- ggplot() +
    # non-significant points (parity colors)
    geom_point(data = dt[sel_knock == FALSE],
               aes(x = pos_cum, y = W, color = chr_parity),
               size = PT_SIZE_BLUE, alpha = 0.8) +
    scale_color_manual(values = c(even = COL_BLUE_L, odd = COL_BLUE_D), guide = "none") +
    # significant points (red)
    geom_point(data = dt[sel_knock == TRUE],
               aes(x = pos_cum, y = W),
               color = COL_RED, size = PT_SIZE_RED, alpha = 0.9) +
    scale_x_continuous(breaks = chr_ticks$center, labels = chr_ticks$chr) +
    labs(title    = sprintf("%s (FDR=%.2f, global q)", trait, alpha_fdr),
         subtitle = subtitle1,
         x = "Chromosome", y = "W statistic") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  if (is.finite(W_thr)) {
    p1 <- p1 + geom_hline(yintercept = W_thr, linetype = "dashed",
                          color = COL_RED, linewidth = 0.6)
  }

  ## 2) W vs -log10(p): 继续红/蓝（不显著统一蓝色）
  subtitle2 <- if (is.finite(W_thr)) {
    sprintf("\u2015 red  W = %.3f    \u2015 yellow  -log10(p) = %.3f (Bonferroni)",
            W_thr, p_cut_log10)
  } else {
    sprintf("\u2015 yellow  -log10(p) = %.3f (Bonferroni)", p_cut_log10)
  }

  p2 <- ggplot(dt, aes(x = mlp, y = W)) +
    geom_point(data = dt[sel_knock == FALSE],
               color = COL_BLUE_D, size = PT_SIZE_BLUE, alpha = 0.6) +
    geom_point(data = dt[sel_knock == TRUE],
               color = COL_RED, size = PT_SIZE_RED, alpha = 0.9) +
    labs(title    = "W vs. -log10(p)",
         subtitle = subtitle2,
         x = expression(-log[10](p)), y = "W statistic") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold")) +
    geom_vline(xintercept = p_cut_log10, linetype = "dashed",
               color = "goldenrod2", linewidth = 0.6)

  if (is.finite(W_thr)) {
    p2 <- p2 + geom_hline(yintercept = W_thr, linetype = "dashed",
                          color = COL_RED, linewidth = 0.6)
  }

  ## save
  ggsave(file.path(out_dir, sprintf("%s_FDR010_W_Manhattan_altBlue.png", trait)),
         plot = p1, width = 12, height = 4, dpi = 300)
  ggsave(file.path(out_dir, sprintf("%s_FDR010_W_vs_log10p.png", trait)),
         plot = p2, width = 6.5, height = 5, dpi = 300)
  ggsave(file.path(out_dir, sprintf("%s_FDR010_two_panels_altBlue.png", trait)),
         plot = p1 + p2 + plot_layout(widths = c(2,1)),
         width = 16, height = 5, dpi = 300)
}

## ---------------- run ----------------
plot_alpha010(
  trait   = "Hypertension",
  dir_q   = "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/score/hyper",
  out_dir = "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/hyper_global_alpha010_altBlue"
)

plot_alpha010(
  trait   = "Type 2 diabetes",
  dir_q   = "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/score/dm",
  out_dir = "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots/dm_global_alpha010_altBlue"
)
