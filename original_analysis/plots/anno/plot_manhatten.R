#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# 设置工作目录
base_dir <- "/dssg/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/plots"
plots <- list(
  list(
    trait     = "Hypertension",
    method    = "knockoff",
    snp_file  = file.path(base_dir, "hyper_counts/Hypertension_global_snp_level.tsv"),
    anno_file = file.path(base_dir, "anno/Hypertension_knockoff_sig_annot.tsv"),
    out_png   = "Hypertension_knockoff.png",
    title     = "Cox-Knockoff Hypertension",
    Pcol      = "W",
    threshold = 3.63073839346906
  ),
  list(
    trait     = "Hypertension",
    method    = "bonferroni",
    snp_file  = file.path(base_dir, "hyper_counts/Hypertension_global_snp_level.tsv"),
    anno_file = file.path(base_dir, "anno/Hypertension_bonferroni_sig_annot.tsv"),
    out_png   = "Hypertension_bonferroni.png",
    title     = "SPACox Hypertension",
    Pcol      = "mlp",
    threshold = -log10(5e-8)
  ),
  list(
    trait     = "Type 2 diabetes",
    method    = "knockoff",
    snp_file  = file.path(base_dir, "dm_counts/Type_2_diabetes_global_snp_level.tsv"),
    anno_file = file.path(base_dir, "anno/Type_2_diabetes_knockoff_sig_annot.tsv"),
    out_png   = "Type2Diabetes_knockoff.png",
    title     = "Cox-Knockoff Type 2 diabetes",
    Pcol      = "W",
    threshold = 4.07295494726266
  ),
  list(
    trait     = "Type 2 diabetes",
    method    = "bonferroni",
    snp_file  = file.path(base_dir, "dm_counts/Type_2_diabetes_global_snp_level.tsv"),
    anno_file = file.path(base_dir, "anno/Type_2_diabetes_bonferroni_sig_annot.tsv"),
    out_png   = "Type2Diabetes_bonferroni.png",
    title     = "SPACox Type 2 diabetes",
    Pcol      = "mlp",
    threshold = -log10(5e-8)
  )
)
## ———————————————————————————————— ##

for (cfg in plots) {
  # 1) 读入全 SNP
  dt <- fread(cfg$snp_file) %>%
    filter(!is.na(chr), chr %in% 1:22) %>%
    mutate(
      CHR = factor(chr, levels = 1:22, ordered = TRUE),
      BP  = pos
    ) %>%
    mutate(P = .data[[cfg$Pcol]], BPcum = NA_real_)

  # 2) 读注释并取前 10（保证 gene_symbol 唯一）
  ann <- fread(cfg$anno_file) %>%
    mutate(
      CHR = factor(chr, levels = 1:22, ordered = TRUE),
      BP  = pos,
      P   = .data[[cfg$Pcol]]
    ) %>%
    arrange(desc(P)) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    slice(1:10) %>%
    mutate(BPcum = NA_real_)

  # 3) 构建累积坐标
  chr_pos <- dt %>%
    group_by(CHR) %>%
    summarise(maxBP = max(BP, na.rm = TRUE)) %>%
    arrange(as.integer(as.character(CHR))) %>%
    mutate(offset = lag(cumsum(as.numeric(maxBP)), default = 0)) %>%
    select(CHR, offset)

  dt  <- left_join(dt,  chr_pos, by = "CHR") %>%
    mutate(BPcum = BP + offset)
  ann <- left_join(ann, chr_pos, by = "CHR") %>%
    mutate(BPcum = BP + offset)

  # 4) x 轴刻度
  x_ticks <- dt %>%
    group_by(CHR) %>%
    summarise(center = (min(BPcum) + max(BPcum)) / 2)

  # 5) 交替配色
  palette_chr <- rep(c("grey70","skyblue"), length.out = 22)
  names(palette_chr) <- levels(dt$CHR)

  # 6) 绘图
  p <- ggplot(dt, aes(x = BPcum, y = P)) +
    geom_point(aes(color = ifelse(P >= cfg$threshold, "Significant", as.factor(CHR))),
               alpha = 0.8, size = 1.3) +
    scale_color_manual(values = c("Significant" = "red", palette_chr),
                       guide  = "none") +
    geom_hline(yintercept = cfg$threshold,
               color = "red", linewidth = 1.2, linetype = "dashed") +
    geom_text_repel(
      data           = ann,
      aes(x = BPcum, y = P, label = gene_symbol),
      size         = 4,
      fontface     = "bold",
      color        = "red",
      box.padding  = 0.5,
      nudge_y      = max(dt$P, na.rm = TRUE) * 0.02,
      segment.color= "red",
      segment.size = 0.7,
      max.overlaps = 10,
      min.segment.length = 0.2
    ) +
    scale_x_continuous(breaks = x_ticks$center,
                       labels = x_ticks$CHR,
                       name   = "Chromosome") +
    scale_y_continuous(name   = ifelse(cfg$method=="bonferroni","-log10(p)","W statistic"),
                       expand = c(0,0),
                       limits = c(0, max(dt$P, na.rm=TRUE)*1.05)) +
    labs(title = cfg$title) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      axis.line.x = element_blank(),  
      axis.line.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title.x = element_text(face = "bold", size = 16),  
      axis.title.y = element_text(face = "bold", size = 16), 
      axis.text.x = element_text(size = 12),  
      axis.text.y = element_text(size = 12)   
    )

  # 7) 保存
  ggsave(cfg$out_png, p, width = 12, height = 5, dpi = 300)
  message("✔ Written: ", cfg$out_png)
}