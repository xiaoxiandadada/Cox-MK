suppressPackageStartupMessages({
  library(Matrix)
  library(BEDMatrix)
  library(gdsfmt)
})

## ---------- 参数 ----------
args <- commandArgs(TRUE)
if (length(args) < 3)
  stop("Usage: generate_knockoff_block_gds.R <chr> <start> <end>")
chr   <- as.integer(args[1]); start <- as.integer(args[2]); end <- as.integer(args[3])

in.dir  <- "split"
out.dir <- "generate"
dir.create(out.dir, FALSE, TRUE)
prefix <- sprintf("%s/ukb_split_%d_%d_%d", in.dir, chr, start, end)

## ---------- 读取 genotype ----------
G <- Matrix(as.matrix(BEDMatrix(paste0(prefix, ".bed"))), sparse = TRUE)

fam <- read.table(paste0(prefix, ".fam"), stringsAsFactors = FALSE)
sample.id <- fam$V2

bim <- read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE)
pos.all <- bim$V4
alleles.all <- paste(bim$V5, bim$V6, sep = "/")

# ## ---------- common 过滤 (MAF≥0.01) ----------
# maf  <- colMeans(G) / 2
# keep <- which(maf >= 0.01 & maf < 0.5)
# G <- G[, keep]; pos <- pos.all[keep]; alle <- alleles.all[keep]

# ---------- ① LD 过滤 (|r| ≥ 0.75, single-link) ----------
if (ncol(G) > 1) {
  sparse.cor <- function(x) {
    n <- nrow(x); cm <- colMeans(x)
    cov <- (as.matrix(crossprod(x)) - n * tcrossprod(cm)) / (n - 1)
    cov / tcrossprod(sqrt(diag(cov)))
  }
  cor.X <- sparse.cor(G)
  cl <- hclust(as.dist(1 - abs(cor.X)), "single")
  clusters <- cutree(cl, h = 1 - 0.75)
  rep.idx  <- match(unique(clusters), clusters)
  G <- G[, rep.idx]
  pos <- pos[rep.idx]
}

p <- ncol(G)                       # 最终 SNP 数

## ---------- 生成 5 套 knockoff (±100 kb 默认) ----------
source("KnockoffScreen.r")
G_k <- create.MK(
  G, pos, M = 5, bigmemory = FALSE,
  maxBP.neighbor = 1e5      # neighborhood ±100 kb
)

## ---------- 组织矩阵：样本 × 变异 (n × p) ----------
to_dense <- function(x) { m <- as.matrix(x); mode(m) <- "integer"; m }
geno.orig <- to_dense(G)              # n × p
geno.k    <- lapply(G_k, to_dense)    # 列表 5 个 n × p

## ---------- 写 GDS ----------
gds.fn <- file.path(out.dir,
  sprintf("chr%d_%d_%d_knockoff.gds", chr, start, end))
g <- createfn.gds(gds.fn)

add.gdsn(g, "sample.id",      sample.id, "string")
add.gdsn(g, "snp.pos",   pos, "int32")
add.gdsn(g, "original",  geno.orig,
         storage = "bit2", valdim = dim(geno.orig), compress = "LZMA_RA")
for (i in 1:5)
  add.gdsn(g, paste0("knockoff", i), geno.k[[i]],
           storage = "bit2", valdim = dim(geno.k[[i]]), compress = "LZMA_RA")

closefn.gds(g)
cat("✓ GDS written:", gds.fn, "\n")
