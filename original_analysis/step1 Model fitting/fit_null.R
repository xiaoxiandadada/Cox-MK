args <- commandArgs(trailingOnly = TRUE)
events_csv <- args[1]
covars_csv <- args[2]
id_keep    <- args[3]
fam_file   <- args[4]
out_dir    <- args[5]

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(SPACox)
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- 1. 读数据 ----------
ev  <- fread(events_csv)
cv  <- fread(covars_csv)
ids <- fread(id_keep, header = FALSE)[[1]]
fam <- fread(fam_file,  header = FALSE)[, .(eid = V2)]

## ---------- 2. 合并 & 过滤 ----------
phe <- merge(ev, cv, by = "eid", all.x = TRUE)

# a) 只留指定白人无关个体
phe <- phe[eid %in% ids]

# b) 再按 .fam 顺序排好，确保 pIDs == gIDs
setkey(phe, eid)
phe <- phe[fam, nomatch = 0]

## ---------- 3. 表型信息 ----------
tr_info <- list(
  hyper = list(time = "survival_time_I10", event = "Phecode_I10"),
  dm    = list(time = "survival_time_T2D", event = "Phecode_T2D")
)

## ---------- 4. 循环拟合 NULL model ----------
for(tr in names(tr_info)){
  tm  <- tr_info[[tr]]$time
  evn <- tr_info[[tr]]$event

  obj <- SPACox_Null_Model(
          as.formula(paste0("Surv(", tm, ",", evn, ") ~ age + sex + ",
                            paste0("PC", 1:10, collapse = " + "))),
          data = phe,
          pIDs = phe$eid,
          gIDs = fam$eid
        )

  rds <- file.path(out_dir, paste0("null_", tr, ".rds"))
  saveRDS(obj, rds)
  cat("✅  保存: ", rds, "\n")
}
