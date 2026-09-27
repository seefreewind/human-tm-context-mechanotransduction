#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("usage: 12_phase2R1_edger.R COUNTS DESIGN OUTPUT")
counts_path <- args[[1]]
design_path <- args[[2]]
out_path <- args[[3]]
suppressPackageStartupMessages(library(edgeR))
counts <- read.delim(counts_path, check.names=FALSE, stringsAsFactors=FALSE)
design_meta <- read.delim(design_path, check.names=FALSE, stringsAsFactors=FALSE)
sample_names <- design_meta$SRR
if (!all(sample_names %in% colnames(counts))) stop("design samples are absent from count matrix")
mat <- as.matrix(counts[, sample_names, drop=FALSE])
rownames(mat) <- counts$peak_id
storage.mode(mat) <- "integer"
design_meta$condition <- factor(design_meta$condition, levels=c("CONTROL", "DEX"))
design_meta$cell_line <- factor(design_meta$cell_line)
design <- model.matrix(~ cell_line + condition, data=design_meta)
y <- DGEList(counts=mat)
keep <- filterByExpr(y, design=design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
coef_name <- grep("conditionDEX", colnames(design), value=TRUE)
if (length(coef_name) != 1) stop("condition coefficient not uniquely identified")
test <- glmQLFTest(fit, coef=which(colnames(design) == coef_name))
tab <- topTags(test, n=Inf, sort.by="none")$table
contrast <- rep(0, ncol(design)); names(contrast) <- colnames(design); contrast[[coef_name]] <- 1
# Recent edgeR versions used in the validated rescue environment do not retain
# stdev.unscaled/var.post in glmQLFit objects. Preserve the tested QL statistics
# and leave SE explicitly unavailable instead of failing after the valid test.
se <- rep(NA_real_, nrow(mat))
all_ids <- counts$peak_id
out <- data.frame(peak_id=all_ids,
                  chr=counts$chr,
                  start=counts$start,
                  end=counts$end,
                  baseMean=rowMeans(mat),
                  log2FC=NA_real_,
                  SE=NA_real_,
                  statistic=NA_real_,
                  P=NA_real_,
                  FDR=NA_real_,
                  tested=FALSE,
                  stringsAsFactors=FALSE)
idx <- match(rownames(tab), out$peak_id)
out$log2FC[idx] <- tab$logFC
out$SE[idx] <- se[match(rownames(tab), rownames(out))]
out$statistic[idx] <- tab$F
out$P[idx] <- tab$PValue
out$FDR[idx] <- p.adjust(tab$PValue, method="BH")
out$tested[idx] <- TRUE
write.table(out, out_path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
