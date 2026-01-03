# Neuroblastoma-bone-marrow-niche-paper

setwd("E:/NB-Sal/S/Z")

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(grid)
library(pheatmap)

# Clean GEO cache
unlink(file.path(tempdir(), "GEO*"), recursive = TRUE, force = TRUE)
unlink(list.files(tempdir(), full.names = TRUE), recursive = TRUE, force = TRUE)


# LOAD GEO DATA
gset <- getGEO("GSE25624", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
#

#gset <- readRDS("E:/NB-Sal/S/GSE25624_gset.rds")
fvarLabels(gset) <- make.names(fvarLabels(gset))
# SAMPLE GROUPS

pheno <- pData(gset)
group <- factor(trimws(sub(".*-", "", pheno$title)))
stopifnot(all(levels(group) %in% c("H_BM", "NI_BM", "I_BM")))
gset$group <- group


# EXPRESSION NORMALIZATION

ex <- exprs(gset)

qx <- quantile(ex, c(0, .25, .5, .75, .99, 1), na.rm = TRUE)
if (qx[5] > 100 || (qx[6] - qx[1] > 50)) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

ex <- normalizeBetweenArrays(ex, method = "quantile")
exprs(gset) <- ex


# LIMMA DESIGN
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(gset, design)

cont.matrix <- makeContrasts(
  I_vs_H  = I_BM - H_BM,
  NI_vs_H = NI_BM - H_BM,
  I_vs_NI = I_BM - NI_BM,
  levels = design
)

fit2 <- eBayes(contrasts.fit(fit, cont.matrix))

#  EXTRACT ALL PROBES

tt_all <- lapply(colnames(cont.matrix), function(ct) {
  topTable(fit2, coef = ct, number = Inf, sort.by = "P")
})
names(tt_all) <- colnames(cont.matrix)


# PROBE-SAFE GENE COLLAPSING FUNCTION

get_gene_DEGs <- function(tt,
                          gene_col = "Gene.symbol",
                          fdr_cutoff = 0.05,
                          lfc_cutoff = NULL) {
  
  tt$ProbeID <- rownames(tt)   # for heatmap
  
  tt <- tt[tt$adj.P.Val < fdr_cutoff, ]
  tt <- tt[!is.na(tt[[gene_col]]) & tt[[gene_col]] != "", ]
  tt <- tt[order(tt$adj.P.Val), ]
  tt <- tt[!duplicated(tt[[gene_col]]), ]
  
  if (!is.null(lfc_cutoff)) {
    tt <- tt[abs(tt$logFC) >= lfc_cutoff, ]
  }
  
  tt
}

#  CHECK
gene_DEGs <- lapply(tt_all, get_gene_DEGs)

sapply(gene_DEGs, function(tt)
  length(unique(tt$Gene.symbol)) == nrow(tt)
)

# APPLY FILTERS
fdr_cutoff <- 0.05
lfc_cutoff <- 1

gene_DEGs_by_contrast <- lapply(
  tt_all,
  get_gene_DEGs,
  fdr_cutoff = fdr_cutoff,
  lfc_cutoff = lfc_cutoff
)


#  SAVE DEG TABLES (WITH ProbeID)
outdir <- "RESULTS/DEGs_BY_CONTRAST"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (ct in names(gene_DEGs_by_contrast)) {
  tt <- gene_DEGs_by_contrast[[ct]]
  
  write.table(
    data.frame(
      ProbeID     = tt$ProbeID,
      Gene.symbol = tt$Gene.symbol,
      log2FC      = tt$logFC,
      AveExpr     = tt$AveExpr,
      t           = tt$t,
      P.Value     = tt$P.Value,
      adj.P.Val   = tt$adj.P.Val
    ),
    file = file.path(outdir, paste0("GENE_DEGs_", ct, "_FDR05.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}


# SAVE PROCESSED GEO OBJECT

saveRDS(gset, "GSE25624_gset.rds")

#check
deg_I_vs_H <- read.delim("RESULTS/DEGs_BY_CONTRAST/GENE_DEGs_I_vs_H_FDR05.tsv")
"ProbeID" %in% colnames(deg_I_vs_H)


#read degs 
deg_dir <- "RESULTS/DEGs_BY_CONTRAST"

deg_I_vs_H <- read.delim(
  file.path(deg_dir, "GENE_DEGs_I_vs_H_FDR05.tsv"),
  stringsAsFactors = FALSE
)

deg_NI_vs_H <- read.delim(
  file.path(deg_dir, "GENE_DEGs_NI_vs_H_FDR05.tsv"),
  stringsAsFactors = FALSE
)

#subset I vs H from all samples
gset_IH <- gset[, gset$group %in% c("I_BM", "H_BM")]
gset_IH$group <- droplevels(gset_IH$group)

#NI vs H
gset_NIH <- gset[, gset$group %in% c("NI_BM", "H_BM")]
gset_NIH$group <- droplevels(gset_NIH$group)


make_heatmap <- function(deg_table,
                         title,
                         gset,
                         top_n = 50,
                         show_rownames = TRUE) {
  
  if (!is.null(top_n)) {
    deg_table <- deg_table[order(deg_table$adj.P.Val), ]
    deg_table <- head(deg_table, top_n)
  }
  
  expr <- exprs(gset)[deg_table$ProbeID, ]
  rownames(expr) <- deg_table$Gene.symbol
  
  expr <- t(scale(t(expr)))
  expr <- expr[apply(expr, 1, function(x) all(is.finite(x))), ]
  
  ann_col <- data.frame(Group = gset$group)
  rownames(ann_col) <- colnames(expr)
  
  ph <- pheatmap(
    expr,
    annotation_col = ann_col,
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    show_rownames  = show_rownames,
    show_colnames  = TRUE,
    fontsize_row   = 8,
    fontsize_col   = 10,
    border_color   = NA,
    main           = title
  )
  
  return(ph)
}
install.packages("pheatmap")
library(pheatmap)

hm_IH <- make_heatmap(
  deg_table = deg_I_vs_H,
  title     = "Metastatic vs Healthy (I vs H)",
  gset      = gset_IH,
  top_n     = 15
)

hm_NIH <- make_heatmap(
  deg_table     = deg_NI_vs_H,
  title         = "Localized vs Healthy (NI vs H)",
  gset          = gset_NIH,
  top_n         = 15,
  show_rownames = TRUE
)


#

# CORE GENE ANALYSIS

disease_core <- intersect(
  gene_DEGs_by_contrast[["I_vs_H"]]$Gene.symbol,
  gene_DEGs_by_contrast[["NI_vs_H"]]$Gene.symbol
)
metastasis_core <- intersect(
  gene_DEGs_by_contrast[["I_vs_H"]]$Gene.symbol,
  gene_DEGs_by_contrast[["I_vs_NI"]]$Gene.symbol
)

#volcano plot
make_volcano_df_gene <- function(gene_DEGs_by_contrast,
                                 fdr_cutoff = 0.05,
                                 lfc_cutoff = 1) {
  
  bind_rows(lapply(names(gene_DEGs_by_contrast), function(ct) {
    
    tt <- gene_DEGs_by_contrast[[ct]]
    
    tt %>%
      mutate(
        Gene = Gene.symbol,
        Contrast = ct,
        negLogFDR = -log10(adj.P.Val),
        Significant = adj.P.Val < fdr_cutoff & abs(logFC) > lfc_cutoff
      )
  }))
}
#volcano df
fdr_cutoff <- 0.05
lfc_cutoff <- 1

volcano_df <- make_volcano_df_gene(
  gene_DEGs_by_contrast,
  fdr_cutoff = fdr_cutoff,
  lfc_cutoff = lfc_cutoff
)

#construct labels
contrast_map <- c(
  "NI_vs_H" = "Localized vs Healthy",
  "I_vs_H"  = "Metastatic vs Healthy",
  "I_vs_NI" = "Metastatic vs Localized"
)

contrast_order <- c(
  "Localized vs Healthy",
  "Metastatic vs Healthy",
  "Metastatic vs Localized"
)

volcano_df <- volcano_df %>%
  mutate(
    Contrast_clean = factor(contrast_map[Contrast],
                            levels = contrast_order)
  )

#top 5 genes labeled
top5_labels <- volcano_df %>%
  filter(Significant) %>%
  group_by(Contrast_clean) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 5) %>%
  ungroup()
#plot 
p <- ggplot(volcano_df, aes(logFC, negLogFDR)) +
  
  geom_point(color = "grey75", alpha = 0.6, size = 1) +
  
  geom_point(
    data = subset(volcano_df, Significant),
    color = "red",
    size = 1
  ) +
  
  geom_text_repel(
    data = top5_labels,
    aes(label = Gene),
    size = 3,
    fontface = "bold",
    box.padding = 0.6,
    point.padding = 0.4,
    force = 5,
    max.overlaps = 10,
    segment.size = 0.4
  ) +
  
  facet_wrap(~ Contrast_clean, nrow = 1, drop = FALSE) +
  
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(fdr_cutoff),
    linetype = "dashed"
  ) +
  
  coord_cartesian(xlim = c(-6, 6), ylim = c(0, 14)) +
  
  labs(
    title = "Differential expression volcano plots (FDR < 0.05)",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )
ggsave("GSE25624_VOLCANO_GENELEVEL_FDR05.pdf",
       p, width = 14, height = 5)

ggsave("GSE25624_VOLCANO_GENELEVEL_FDR05.png",
       p, width = 14, height = 5, dpi = 300)

#Venn digram
library(VennDiagram)
library(grid)

gene_sets <- list(
  "Metastatic vs Healthy"   = gene_DEGs_by_contrast[["I_vs_H"]]$Gene.symbol,
  "Localized vs Healthy"    = gene_DEGs_by_contrast[["NI_vs_H"]]$Gene.symbol,
  "Metastatic vs Localized" = gene_DEGs_by_contrast[["I_vs_NI"]]$Gene.symbol
)
outdir <- "RESULTS/VENN_GENELEVEL_FDR05"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

venn.plot <- venn.diagram(
  x = gene_sets,
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.0,
  cat.dist = c(0.06, 0.06, 0.06),
  margin = 0.25,
  main = "GSE25624 – Gene-level DEGs (FDR < 0.05)",
  filename = NULL
)

pdf(file.path(outdir, "VENN_GENELEVEL_FDR05.pdf"),
    width = 7, height = 7)

png(file.path(outdir, "VENN_GENELEVEL_FDR05.png"),
    width = 7, height = 7, units = "in", res = 300)
grid.draw(venn.plot)
dev.off()
grid.newpage()
grid.draw(venn.plot)
dev.off()

sapply(gene_sets, length)

#Venn intersection
g <- lapply(gene_DEGs_by_contrast, `[[`, "Gene.symbol")

outdir <- "RESULTS/VENN_INTERSECTIONS"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write <- function(x, f)
  write.table(data.frame(Gene=x),
              file.path(outdir, f),
              sep="\t", row.names=FALSE, quote=FALSE)

write(Reduce(intersect, g), "CORE_ALL.tsv")
write(intersect(g$I_vs_H, g$NI_vs_H), "CORE_DISEASE.tsv")
write(intersect(g$I_vs_H, g$I_vs_NI), "CORE_METASTASIS.tsv")
write(setdiff(g$I_vs_H, union(g$NI_vs_H, g$I_vs_NI)), "UNIQUE_I_vs_H.tsv")
write(setdiff(g$NI_vs_H, union(g$I_vs_H, g$I_vs_NI)), "UNIQUE_NI_vs_H.tsv")
write(setdiff(g$I_vs_NI, union(g$I_vs_H, g$NI_vs_H)), "UNIQUE_I_vs_NI.tsv")

save.image("analysis_workspace.RData")
load("analysis_workspace.RData")

