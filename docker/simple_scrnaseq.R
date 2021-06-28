#
#  Basic scRNAseq analysis pipeline based on the Seurat v4 tutorial.
#
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))

# Parse input data
option_list = list(
  make_option(c("-d", "--data_dir"), type="character", default=NULL, help="Path to a directory with a Cellranger counts matrix", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Path to a directory where output is written", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
data_dir = opt$data_dir
output_dir = opt$output_dir

############################################################################################################
# Read data and initialise Seurat object
############################################################################################################
cat(" >>> Running step 1 <<<\n")

# Load the given dataset
# pbmc.data <- Read10X(data.dir = "~/Documents/reference_data/human/10x_genomics/filtered_gene_bc_matrices/hg19//")
pbmc.data <- Read10X(data.dir = data_dir)

# Instantiate Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

############################################################################################################
# Basic QC
############################################################################################################
cat(" >>> Running step 2 <<<\n")

# Annotate MT expression
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
p = VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggplot2::ggsave(plot=p, filename=file.path(output_dir, "02_qc.png"), width=2496, height=1344, units="px")

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p = plot1 + plot2
ggplot2::ggsave(plot=p, filename=file.path(output_dir, "02_qc2.png"), width=2496, height=1344, units="px")

# remove outlier cells by QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Downsample - to keep memory down
pbmc <- subset(pbmc, cells = sample(Cells(pbmc), 1000))

############################################################################################################
# Normalisation and scaling
############################################################################################################
cat(" >>> Running step 3 <<<\n")

# Log normalise and find variable genes
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggplot2::ggsave(plot=plot1, filename=file.path(output_dir, "03_var_features.png"), width=1712, height=960, units="px")
ggplot2::ggsave(plot=plot2, filename=file.path(output_dir, "03_var_features_labelled.png"), width=1712, height=960, units="px")

# Perform uniform scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

############################################################################################################
# Clustering and marker finding
############################################################################################################
cat(" >>> Running step 4 <<<\n")

# Dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
p = DimPlot(pbmc, reduction = "pca")
ggplot2::ggsave(plot=p, filename=file.path(output_dir, "04_pca.png"), width=1400, height=865, units="px")

p = ElbowPlot(pbmc)
ggplot2::ggsave(plot=p, filename=file.path(output_dir, "04_elbow.png"), width=1920, height=1152, units="px")

# We're choosing 10 dimensions by default, find clusters and run UMAP dimension reduction
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
p = DimPlot(pbmc, reduction = "umap")
ggplot2::ggsave(plot=p, filename=file.path(output_dir, "04_umap.png"), width=1400, height=865, units="px")
p = LabelClusters(plot = p, id = 'ident', size=6)
ggplot2::ggsave(plot=p, filename=file.path(output_dir, "04_umap_labelled.png"), width=1400, height=865, units="px")

# Obtain marker genes
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers, file=file.path(output_dir, "04_markers.txt"), quote=F, sep="\t", row.names=F)

# Plot top markers
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
p = DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggplot2::ggsave(plot=p, filename=file.path(output_dir, "04_top_markers_heatmap.png"), width=2880, height=1936, units="px")

############################################################################################################
# Save results
############################################################################################################
cat(" >>> Running step 5 <<<\n")

saveRDS(pbmc, file = file.path(output_dir, "analysis_final.rds"))

cat(" >>> Done\n")





