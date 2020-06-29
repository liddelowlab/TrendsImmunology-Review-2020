#pipeline prepared and Mathys et al. (2019) data reanalyzed by Jessica S. Sadick
#majority of pipeline based on code originally deposited by https://github.com/HelenaLC and published in doi.org/10.1101/713412

#---------------------------------------------------------------------------------------------------
#Load libraries
library(cowplot)
library(ggplot2)
library(scater)
library(scds)
library(SingleCellExperiment)

#LOAD AND REFORMAT DATA
#Load raw counts
fastq_dirs <- list.dirs("file_path", recursive = FALSE, full.names = TRUE)
names(fastq_dirs) <- basename(fastq_dirs)
sce <- DropletUtils::read10xCounts(fastq_dirs)

#Rename row/colData colnames & SCE dimnames
names(rowData(sce)) <- c("ENSEMBL", "SYMBOL")
names(colData(sce)) <- c("sample_id", "barcode")
sce$sample_id <- factor(basename(sce$sample_id))
dimnames(sce) <- list(
  with(rowData(sce), paste(SYMBOL, sep = ".")),
  with(colData(sce), paste(barcode, sample_id, sep = ".")))

#Load metadata
md_dir <- file.path("file_path", "metadata_Mathys_2019.xlsx")
md <- readxl::read_excel(md_dir)
m <- match(sce$sample_id, md$`Sample ID`)

#Assign metadata variables
sce$group_id <- md$Characteristics[m]
sce$donor_id <- md$Donor[m]
sce$sex_id <- md$Sex[m]
sce$disease_id <- md$Disease[m]
sce$age_id <- md$Age[m]
sce$APOE_id <- md$APOE[m]
sce$braak_id <- md$braaksc[m]
sce$cerad_id <- md$ceradsc[m]
sce$MMSE_id <- md$MMSE[m]
sce$PMI_id <- md$PMI[m]

#Remove undetected genes
sce <- sce[Matrix::rowSums(counts(sce) > 0) > 0, ]
dim(sce)

#DOUBLET REMOVAL
#Split SCE by sample
cs_by_s <- split(colnames(sce), sce$sample_id)
sce_by_s <- lapply(cs_by_s, function(cs) sce[, cs])

#Run 'scds'
sce_by_s <- lapply(sce_by_s, function(u) 
  cxds_bcds_hybrid(bcds(cxds(u))))

#Remove doublets
sce_by_s <- lapply(sce_by_s, function(u) {
  #Compute expected number of doublets (10x)
  n_dbl <- ceiling(0.01 * ncol(u)^2 / 1e3)
  #Remove 'n_dbl' cells with highest doublet score
  o <- order(u$hybrid_score, decreasing = TRUE)
  u[, -o[seq_len(n_dbl)]]
})

#Merge back into single SCE
sce <- do.call(cbind, sce_by_s)

#CALCULATE QC METRICS
(mito <- grep("MT-", rownames(sce), value = TRUE))
sce <- calculateQCMetrics(sce, feature_controls = list(Mt = mito))

#FILTERING
#Get sample-specific outliers
cols <- c("total_counts", "total_features_by_counts", "pct_counts_Mt")
log <- c(TRUE, TRUE, FALSE)
type <- c("both", "both", "higher")

drop_cols <- paste0(cols, "_drop")
for (i in seq_along(cols))
  colData(sce)[[drop_cols[i]]] <- isOutlier(sce[[cols[i]]], 
                                            nmads = 2.5, type = type[i], log = log[i], batch = sce$sample_id)

sapply(drop_cols, function(i) 
  sapply(drop_cols, function(j)
    sum(sce[[i]] & sce[[j]])))

cd <- data.frame(colData(sce))
ps <- lapply(seq_along(cols), function (i) {
  p <- ggplot(cd, aes_string(x = cols[i], alpha = drop_cols[i])) +
    geom_histogram(bins = 100, show.legend = FALSE) +
    scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.4)) +
    facet_wrap(~sample_id, ncol = 1, scales = "free") + 
    theme_classic() + theme(strip.background = element_blank())
  if (log[i]) 
    p <- p + scale_x_log10()
  return(p)
})
plot_grid(plotlist = ps, ncol = 3)

layout(matrix(1:2, nrow = 1))
ol <- Matrix::rowSums(as.matrix(colData(sce)[drop_cols])) != 0
x <- sce$total_counts
y <- sce$total_features_by_counts
LSD::heatscatter(x, y, log="xy", main = "unfiltered", 
                 xlab = "Total counts", ylab = "Non-zero features")
LSD::heatscatter(x[!ol], y[!ol], log="xy", main = "filtered", 
                 xlab = "Total counts", ylab = "Non-zero features")

#Generate summary of cells kept
ns <- table(sce$sample_id)
ns_fil <- table(sce$sample_id[!ol])
print(rbind(
  unfiltered = ns, filtered = ns_fil, 
  "%" = ns_fil / ns * 100), digits = 0)

#Drop outlier cells
sce <- sce[, !ol]
dim(sce)

#Require count > 1 in at least 20 cells
sce <- sce[Matrix::rowSums(counts(sce) > 1) >= 20, ]
dim(sce)

#Save SCE
saveRDS(sce, file.path("file_path", "GRCh38_Mathys_SCE.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTERING
#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load data
sce <- readRDS(file.path("file_path", "GRCh38_Mathys_SCE.rds"))

#INTEGRATE
#Create SeuratObject
so <- CreateSeuratObject(
  counts = counts(sce),
  meta.data = data.frame(colData(sce)),
  project = "Mathys_10x_data")

#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so <- lapply(cells_by_sample, function(i)
  SubsetData(so, cells = i))

#Normalize, find variable genes, and scale
so <- lapply(so, NormalizeData, verbose = FALSE)
so <- lapply(so, FindVariableFeatures, nfeatures = 2e3,
             selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so <- lapply(so, ScaleData, verbose = FALSE)

#Save unintergrated object to feed into HPC for integration
saveRDS(so, file.path("file_path", "GRCh38_Mathys_SO_unint.rds"))

#---------------------------------------------------------------------------------------------------
#Integration completed on NYULMC HPC due to memory constraints
#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Set working directory
setwd("file_path")

#Load data
so <- readRDS(file.path("file_path", "GRCh38_Mathys_SO_unint.rds"))

#Change to reference-based integration because of memory constraints
#Choose 3 donors from each condition (i.e., non-symptomatic female, non-symptomatic male, AD female, and AD male) with the most nuclei captured
#Non-symptomatic female references
ROS1 <- which(names(so) == "ROS1")
ROS10 <- which(names(so) == "ROS10")
ROS12 <- which(names(so) == "ROS12")
#Non-symptomatic male  references
ROS15 <- which(names(so) == "ROS15")
ROS17 <- which(names(so) == "ROS17")
ROS24 <- which(names(so) == "ROS24")
#AD female references
ROS26 <- which(names(so) == "ROS26")
ROS32 <- which(names(so) == "ROS32")
ROS36 <- which(names(so) == "ROS36")
#AD male references
ROS37 <- which(names(so) == "ROS37")
ROS38 <- which(names(so) == "ROS38")
ROS42 <- which(names(so) == "ROS42")

as <- FindIntegrationAnchors(so, reference = c(ROS1,ROS10,ROS12,ROS15,ROS17,ROS24,ROS26,ROS32,ROS36,ROS37,ROS38,ROS42), verbose=TRUE)
so <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so) <- "integrated"
so <- ScaleData(so, display.progress = TRUE)

#DIMENSION REDUCTION
so <- RunPCA(so, npcs = 100, verbose = TRUE)
ElbowPlot(so, ndims = 100)
#Determine how many principal components to continue on with and update for tSNE and UMAP calculations: PC=30
so <- RunTSNE(so, reduction = "pca", dims = seq_len(30),
seed.use = 1, do.fast = TRUE, verbose = TRUE)
so <- RunUMAP(so, reduction = "pca", dims = seq_len(30),
seed.use = 1, verbose = TRUE)

#CLUSTERING
so <- FindNeighbors(so, reduction = "pca", dims = seq_len(30), verbose = TRUE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1))
so <- FindClusters(so, resolution = res, random.seed = 1, verbose = TRUE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
    p1 <- DimPlot(so, reduction = "tsne", group.by = u) + thm
    p2 <- DimPlot(so, reduction = "umap", group.by = u)
    lgd <- get_legend(p2)
    p2 <- p2 + thm
    list(p1, p2, lgd)
    plot_grid(p1, p2, lgd, nrow = 1,
    rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#Save seurat object
saveRDS(so, file.path("file_path", "GRCh38_Mathys_SO-ref_30PC.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION
#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(viridis)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so <- readRDS(file.path("file_path", "GRCh38_Mathys_SO-ref_30PC.rds"))
sce <- as.SingleCellExperiment(so, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Determine # clusters
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)
so <- SetIdent(so, value = "integrated_snn_res.0.1")

so@meta.data$cluster_id <- Idents(so)
sce$cluster_id <- Idents(so)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/Mathys_SO_cluster_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green"), gids)

#Generate relative cluster abundances
fqs <- prop.table(n_cells, margin = 2)
mat <- as.matrix(unclass(fqs))
Heatmap(mat,
        col = rev(brewer.pal(11, "RdGy")[-6]),
        name = "Frequency",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_title = "cluster_id",
        column_title = "sample_id",
        column_title_side = "bottom",
        rect_gp = gpar(col = "white"),
        cell_fun = function(i, j, x, y, width, height, fill)
          grid.text(round(mat[j, i] * 100, 2), x = x, y = y, 
                    gp = gpar(col = "white", fontsize = 8)))

#DR colored by cluster ID
cs <- sample(colnames(so), 5e3)
.plot_dr <- function(so, dr, id)
  DimPlot(so, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^MT-", x = rownames(so@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so@assays[["RNA"]])
so$percent.mito <- percent.mito

rb.genes <- grep(pattern = "^RP[SL]", x = rownames(so@assays[["RNA"]]), value = TRUE)
percent.rb <- Matrix::colSums(so@assays[["RNA"]][rb.genes, ])/Matrix::colSums(so@assays[["RNA"]])
so$percent.rb <- percent.rb

VlnPlot(object = so, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rb"), ncol = 4, pt.size = 0)

VlnPlot(object = so, features = c("nFeature_RNA"), pt.size = 0) + stat_summary(fun.y=median, geom="point", shape=23, size=2)

#Generate summary statistics per sample
library(data.table)
library(psych)

feature_by_sample <- as.data.frame(so$nFeature_RNA, row.names = so$sample_id)
feature_by_sample_table <- describeBy(feature_by_sample, group = so$sample_id, mat = TRUE)
write.csv(feature_by_sample_table, "file_path/Mathys_SO_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so$nCount_RNA, row.names = so$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "file_path/Mathys_SO_QC_count_by_sample.csv")

#Generate summary statistics per cluster
feature_by_cluster <- as.data.frame(so$nFeature_RNA, row.names = so$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "file_path/Mathys_SO_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so$nCount_RNA, row.names = so$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "file_path/Mathys_SO_QC_count_by_cluster.csv")

#DETERMINE WHAT CELL TYPES ARE PRESENT IN DATASET
#Find all markers
DefaultAssay(so) <- "integrated"
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so.markers, "file_path/Mathys_SO_genes.csv")

so.topmarkers <- so.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(so.topmarkers, "file_path/Mathys_SO_topgenes.csv")

#Visualize (examples)
DefaultAssay(so) <- "RNA"
DimPlot(so, reduction = "tsne") + theme(aspect.ratio = 1)

DotPlot(so, features = c("AQP4", "GFAP", "FGFR3", "CLDN10", "GJA1", "ALDH1L1", "SLC1A3", "SLC1A2", "HIST1H3E", "TPX2", "NUSAP1", "PPDPF"))
DotPlot(so, features = c("CLDN5", "NOSTRIN", "PDGFRB", "ANPEP", "FLT1", "ACTA2", "PECAM1", "VWF", "CSPG4"))
DotPlot(so, features = c("C1QB", "TYROBP"))
DotPlot(so, features = c("SNAP25", "STMN2", "SLC17A7", "GAD1"))
DotPlot(so, features = c("PDGFRA", "CSPG4", "SOX10"))
DotPlot(so, features = c("PLP1", "MBP", "MOG", "OLIG2"))

FeaturePlot(so, features = c("SOX9", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")

#---------------------------------------------------------------------------------------------------
#MAKE ASTROCYTE-SPECIFIC SEURAT OBJECT
#Assign cell type identity to clusters
so.renamed <- RenameIdents(so, `0` = "Oligo", `1` = "Neuro_A", `2` = "Neuro_B", `3` = "Neuro_C", `4` = "Neuro_D", `5`= "Astro", `6` = "Neuro_E", `7`= "Neuro_F", `8`= "OPC", `9`= "Neuro_G", `10`= "Micro", `11`= "Neuro_H", `12`= "Neuro_I", `13`= "Neuro_J", `14`= "Neuro_K", `15`= "Neuro_L", `16`= "Neuro_M")

#Remove donors with < 30 astrocytes captured (i.e., ROS13, ROS40, ROS43, ROS44, ROS45, ROS7, ROS8) to meet minimize threshold for future integration
so_1 <- subset(x = so.renamed, subset = sample_id == "ROS13", invert = TRUE)
so_2 <- subset(x = so_1, subset = sample_id == "ROS40", invert = TRUE)
so_3 <- subset(x = so_2, subset = sample_id == "ROS43", invert = TRUE)
so_4 <- subset(x = so_3, subset = sample_id == "ROS44", invert = TRUE)
so_5 <- subset(x = so_4, subset = sample_id == "ROS45", invert = TRUE)
so_6 <- subset(x = so_5, subset = sample_id == "ROS7", invert = TRUE)
so_7 <- subset(x = so_6, subset = sample_id == "ROS8", invert = TRUE)

#Subset object to only include astrocyte cluster
so_astro <- subset(x = so_7, idents = c("Astro"), invert = FALSE)

#---------------------------------------------------------------------------------------------------
#RECLUSTER NEW SEURAT OBJECT: ASTROCYTES
#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#INTEGRATE
#Make new SCE object because removed several donors from original SCE object due to low astrocyte captures
sce <- as.SingleCellExperiment(so_astro, assay = "RNA")

#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_astro <- lapply(cells_by_sample, function(i)
  SubsetData(so_astro, cells = i))

#Normalize, find variable genes, and scale
so_astro <- lapply(so_astro, NormalizeData, verbose = FALSE)
so_astro <- lapply(so_astro, FindVariableFeatures, nfeatures = 2e3,
                   selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so_astro <- lapply(so_astro, ScaleData, verbose = FALSE)

#Find anchors & integrate
#Decrease k.filter because of donors with low astrocyte captures
as <- FindIntegrationAnchors(so_astro, verbose = FALSE, k.filter = 30)
so_astro <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE)

#Scale integrated data
DefaultAssay(so_astro) <- "integrated"
so_astro <- ScaleData(so_astro, display.progress = FALSE)

#DIMENSION REDUCTION
so_astro <- RunPCA(so_astro, npcs = 50, verbose = FALSE)
ElbowPlot(so_astro, ndims = 50)

#Determine how many PCs to continue on with and update for tSNE and UMAP calculations: PC=15
so_astro <- RunTSNE(so_astro, reduction = "pca", dims = seq_len(15),
                    seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_astro <- RunUMAP(so_astro, reduction = "pca", dims = seq_len(15),
                    seed.use = 1, verbose = FALSE)

#CLUSTERING
so_astro <- FindNeighbors(so_astro, reduction = "pca", dims = seq_len(15), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_astro <- FindClusters(so_astro, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_astro, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_astro, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#Save seurat object
saveRDS(so_astro, file.path("file_path", "GRCh38_Mathys_SO_astro_15PC.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION FOR ASTROCYTE-SPECIFIC SEURAT OBJECT
#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(viridis)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so_astro <- readRDS(file.path("file_path", "GRCh38_Mathys_SO_astro_15PC.rds"))
sce <- as.SingleCellExperiment(so_astro, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Determine # clusters
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)
so_astro <- SetIdent(so_astro, value = "integrated_snn_res.0.3")

so_astro@meta.data$cluster_id <- Idents(so_astro)
sce$cluster_id <- Idents(so_astro)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/Mathys_SO_astro_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green"), gids)

#Generate relative cluster abundances
fqs <- prop.table(n_cells, margin = 2)
mat <- as.matrix(unclass(fqs))
Heatmap(mat,
        col = rev(brewer.pal(11, "RdGy")[-6]),
        name = "Frequency",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_title = "cluster_id",
        column_title = "sample_id",
        column_title_side = "bottom",
        rect_gp = gpar(col = "white"),
        cell_fun = function(i, j, x, y, width, height, fill)
          grid.text(round(mat[j, i] * 100, 2), x = x, y = y, 
                    gp = gpar(col = "white", fontsize = 8)))

#DR colored by cluster ID
cs <- sample(colnames(so_astro), 3e3)
.plot_dr <- function(so_astro, dr, id)
  DimPlot(so_astro, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so_astro, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so_astro, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#DETERMINE WHAT DEFINES EACH ASTROCYTE CLUSTER
#Find all markers
DefaultAssay(so_astro) <- "integrated"
so_astro.markers <- FindAllMarkers(so_astro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_astro.markers, "file_path/Mathys_SO_astro_genes.csv")

so_astro.topmarkers <- so_astro.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(so_astro.topmarkers, "file_path/Mathys_SO_astro_topgenes.csv")

#Visualize (examples)
DefaultAssay(so_astro) <- "RNA"
DimPlot(so_astro, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1)

DotPlot(so_astro, features = c("AQP4", "GFAP", "FGFR3", "CLDN10", "GJA1", "ALDH1L1", "SLC1A3", "SLC1A2", "HIST1H3E", "TPX2", "NUSAP1", "PPDPF"))
DotPlot(so_astro, features = c("CLDN5", "NOSTRIN", "PDGFRB", "ANPEP", "FLT1", "ACTA2", "PECAM1", "VWF", "CSPG4"))
DotPlot(so_astro, features = c("C1QB", "TYROBP"))
DotPlot(so_astro, features = c("SNAP25", "STMN2", "SLC17A7", "GAD1"))
DotPlot(so_astro, features = c("PDGFRA", "CSPG4", "SOX10"))
DotPlot(so_astro, features = c("PLP1", "MBP", "MOG", "OLIG2"))
DotPlot(so_astro, features = c("AQP4", "GFAP", "CLDN10", "GJA1", "ALDH1L1", "SLC1A3", "SLC1A2", "PDGFRA", "CSPG4", "SOX10", "PLP1", "MBP", "MOG", "C1QB", "TYROBP", "SNAP25", "STMN2", "SLC17A7", "GAD1", "CLDN5", "NOSTRIN", "PDGFRB")) + theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(so_astro, features = c("PDGFRB", "NOSTRIN", "CLDN5", "GAD1", "SLC17A7", "STMN2", "SNAP25", "TYROBP", "C1QB", "MOG", "PLP1", "SOX10", "CSPG4", "PDGFRA", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1", "CLDN10", "GFAP", "AQP4")) + theme(axis.text.x = element_text(angle = 45, hjust=1))

FeaturePlot(so_astro, features = c("SOX9", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")
