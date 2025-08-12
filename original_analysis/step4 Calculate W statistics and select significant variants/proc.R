#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

## ================== 用户参数 ==================
ROOT_IN    <- "/lustre/home/acct-mashiyang1991/dhhhh020606/tzh2/step3_res"
COORD_FILE <- "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/coords.txt"
OUT_DIR    <- "/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/SPAcox/score/dm"
M_KNOCK    <- 5
## ============================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## -------- W 函数 --------
max.nth <- function(x, n){ sort(x, partial = length(x)-(n-1))[length(x)-(n-1)] }

MK.statistic <- function (T_0, T_k, method='median'){
  T_0 <- as.matrix(T_0); T_k <- as.matrix(T_k)
  T.temp <- cbind(T_0, T_k)
  T.temp[is.na(T.temp)] <- 0

  which.max.alt <- function(x){
    temp.index <- which(x == max(x))
    if(length(temp.index) != 1) temp.index[2] else temp.index[1]
  }
  kappa <- apply(T.temp, 1, which.max.alt) - 1

  if(method == 'max'){
    tau <- apply(T.temp, 1, max) - apply(T.temp, 1, max.nth, n = 2)
  } else {
    Get.OtherMedian <- function(x) median(x[-which.max(x)])
    tau <- apply(T.temp, 1, max) - apply(T.temp, 1, Get.OtherMedian)
  }
  cbind(kappa = kappa, tau = tau)
}
## --------------------------------

## 读 coords，拼输入/输出文件名
coords <- fread(COORD_FILE, col.names = c("chr","start","end"))
coords[, f_in  := file.path(ROOT_IN,
                            sprintf("chr%s_%s_%s_knockoff_p.csv", chr, start, end))]
coords[, f_out := file.path(OUT_DIR,
                            sprintf("chr%s_%s_%s_W_basic.csv", chr, start, end))]

cat(sprintf("[INFO] 共 %d 个块需要处理\n", nrow(coords)))

process_block <- function(row, M = M_KNOCK){
  f_in  <- row$f_in
  f_out <- row$f_out
  blk   <- sprintf("chr%s_%s_%s", row$chr, row$start, row$end)

  if(!file.exists(f_in)){
    message("[MISS] ", blk, " -> ", f_in)
    return(FALSE)
  }

  dt <- tryCatch(fread(f_in),
                 error = function(e){ message("[ERROR fread] ", blk, ": ", e$message); return(NULL) })
  if(is.null(dt)) return(FALSE)

  pgk_cols <- paste0("pGk", 1:M)
  need_cols <- c("SNP","pG0", pgk_cols)
  miss <- setdiff(need_cols, names(dt))
  if(length(miss)){
    message("[SKIP] ", blk, " 缺少列: ", paste(miss, collapse = ","))
    return(FALSE)
  }

  # 计算 W, kappa, tau
  Z.orig <- -log10(dt$pG0)
  Z.knk  <- as.matrix(-log10(dt[, ..pgk_cols]))
  MK     <- MK.statistic(Z.orig, Z.knk, method = "median")
  W      <- (Z.orig - apply(Z.knk, 1, median)) * (Z.orig >= apply(Z.knk, 1, max))

  # 构造 out
  out <- data.table(
    chr   = row$chr,
    start = row$start,
    end   = row$end,
    pos   = as.numeric(dt$SNP),
    pG0   = as.numeric(dt$pG0),
    dt[, ..pgk_cols],
    W     = as.numeric(W),
    kappa = as.numeric(MK[, "kappa"]),
    tau   = as.numeric(MK[, "tau"])
  )

  # 确保列顺序严格一致
  setcolorder(out, c("chr","start","end","pos","pG0", pgk_cols, "W","kappa","tau"))

  fwrite(out, f_out, sep = "\t")
  message("[OK]   ", blk, " -> ", f_out, "  (n=", nrow(out), ")")
  TRUE
}

res <- mclapply(split(coords, 1:nrow(coords)),
                function(r) process_block(r[1,]),
                mc.cores = max(1, detectCores() - 1))

cat(sprintf("[DONE] 成功写出 %d / %d 个块到 %s\n",
            sum(unlist(res)), length(res), OUT_DIR))
