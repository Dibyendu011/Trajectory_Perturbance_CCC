BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages('Seurat')
devtools::install_github("satijalab/seurat-wrappers")

files = list.files()
length(files)

# selecting files with "2i" or "Dox" in their name
files = grep("2i|Dox", files, value = TRUE)

# filtering selected files with .h5 extension
files = grep("h5", files, value = TRUE)
files
length(files)

# low RAM, keep only very significant files
files = files[c(1,7,13,21,33,35,45,57,65,79)]
files

# Load required libraries
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(hdf5r)
setwd("/home/ibab/Desktop/DIBYENDU/s3/RA/scrna/trajectory/")
getwd()

# Define function to process each input file
input_files = function(h5_path) {
  # Extract day information from filename using regex
  day = sub("^.*_D([^_]*)_Dox.*$", "\\1", h5_path)
  day = sub("^.*_D([^_]*)_2i.*$", "\\1", day)
  
  # Read 10X Genomics H5 file and create Seurat object
  data = Read10X_h5(paste0(h5_path))
  data = CreateSeuratObject(data, min.cells = 0, min.features = 200)
  
  # Calculate mitochondrial percentage
  data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "mt-")
  
  # Calculate quantiles for filtering
  lb = quantile(data[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02)
  ub = quantile(data[["nFeature_RNA"]]$nFeature_RNA, probs = 0.97)
  
  # Filter cells based on feature counts and mitochondrial percentage
  data = data[, data[["nFeature_RNA"]] > lb & data[["nFeature_RNA"]] < ub & data[["percent.mt"]] < 15]
  
  # Add day metadata
  data$day = day
  return(data)
}

# Apply input_files function to all files in the list
data_list = sapply(files, input_files)

# Merge all Seurat objects into one combined dataset
# First object is taken separately, then remaining objects are merged with it
data = merge(data_list[1]$GSM3195648_D0_Dox_C1_gene_bc_mat.h5, y = data_list[2:length(data_list)])

data

# Save entire Seurat object locally
saveRDS(data, file = "seurat_data.rds")

# analogous to adata.obs
data@meta.data

# analogous to adata.X(just transposed)
GetAssayData(data,layer="counts.1")[1:5, 1:5]

# make table if u want
df = as.data.frame(as.matrix(GetAssayData(data, layer = "counts.1")))

# get sum of all columns
colSums(GetAssayData(data,layer="counts.1"))

# check for NaN values
colSums(is.nan(GetAssayData(data,layer="counts.1")))

# to impute with mean
for (i in 1:ncol(GetAssayData(data,layer="counts.1"))) {
  col_mean <- mean(GetAssayData(data,layer="counts.1")[, i], na.rm = TRUE)  
  GetAssayData(data,layer="counts.1")[is.nan(GetAssayData(data,layer="counts.1")[, i]), i] <- col_mean
}

# Standard Seurat preprocessing pipeline
data = NormalizeData(object = data, verbose = FALSE)  # Normalize counts
data = FindVariableFeatures(object = data, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')  # Find variable genes
data = ScaleData(data, verbose = FALSE)  # Scale and center data
data = RunPCA(data, npcs = 30, verbose = FALSE)  # Run PCA dimensionality reduction
data = FindNeighbors(data, dims = 1:30)  # Find nearest neighbors for clustering

# Run UMAP for visualization
data = RunUMAP(data, reduction = "pca", dims = 1:30)

# Check available assays
names(data@assays)

# check reduction names
names(data@reductions)

# View scaled data matrix (cells x variable genes) -> ALWAYS use View
View(data[["RNA"]]$scale.data)

# View PCA cell embeddings (cells x PCs) -> Coordinates of cells in PCA space
View(Embeddings(data, reduction = "pca"))

# View PCA gene loadings (genes x PCs) -> Contribution of genes to each principal component
View(Loadings(data, reduction = "pca"))

# Visualize PCA
DimPlot(data, reduction = "pca")
ElbowPlot(data)  # Scree plot

# Set active assay to RNA
data@active.assay = 'RNA'

# check obs features
colnames(data@meta.data)

# Create UMAP plot colored by day information, can also plot for other features
DimPlot(data, group.by = c("day"))

# try View(data)

# Convert day to numeric for pseudotime analysis
data$dayint = data[[]]$day # adds dayint column to data@meta.data actually
# Convert "iPSC" to numeric value 20
data$dayint = ifelse(data$dayint == "iPSC", 20, data$dayint)
data$dayint = as.numeric(data$dayint)

# Visualize day information as a feature on UMAP. Almost the same as DimPlot(data, group.by = c("dayint")). Dimplot treats dayint as categorical values & Featureplot treats dayint as continuous values
FeaturePlot(data, "dayint")

saveRDS(data, file = "seurat_data.rds")

#data = readRDS("seurat_data.rds")

# converting seurat version 5 object to monocle object
# cds = SeuratWrappers::as.cell_data_set(data) -> doesnt work for seurat 5 yet, so we need to copy all contents of adata to a seurat version 4 object manually & then convert it to monocle object.
# Get all counts layers
counts_layers = Layers(data, search = "counts")
print(counts_layers)
# Combine all counts matrices (if they're from different samples)
combined_counts = do.call(cbind, lapply(counts_layers, function(layer) {
    LayerData(data, layer = layer)
}))
# Create v4 object
data_v4 = CreateSeuratObject(
    counts = combined_counts,
    meta.data = data[[]]
)
# Copy reductions
data_v4[["pca"]] = data[["pca"]]
data_v4[["umap"]] = data[["umap"]]
# Copy normalized data layer if needed
data_layers = Layers(data, search = "data")
if(length(data_layers) > 0) {
    combined_data = do.call(cbind, lapply(data_layers, function(layer) {
        LayerData(data, layer = layer)
    }))
    data_v4[["RNA"]]$data = combined_data
}
# Copy variable features
VariableFeatures(data_v4) = VariableFeatures(data)
# Convert to CellDataSet
cds = SeuratWrappers::as.cell_data_set(data_v4)

# Seurat v5 → Had separate layers: counts.1, counts.2, ..., counts.10. Monocle3 conversion → Combined all layers into single counts assay. We cannot see the individual count layers in the Monocle3 object.

# view the combined table
counts(cds)
dim(counts(cds))
colnames(counts(cds))
rownames(counts(cds))

# Perform clustering in Monocle3
cds = cluster_cells(cds)

# Plot cells colored by Monocle3 partitions. 
#Partition isn't same as clusters. Partitions are patches of data points that are significantly closer to each other on the dimention-reduced plot. It is fixed for a pair of PCs. Each patch can be a mosaic of more than one clusters with a cluster being shared between two patches. Clusters are groups of data points with similar expression. Clusters change with features considered for clustering. Partition is like a country & cluster ethnic/religious/political groups.
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")


# Steps of trajectory analysis : 
# 1. Euclidean distance is calculated between every data point & its k nearest neighbours. 
# 2. The MST (Minimum Spinning Tree) is then created to integrate the local distance informations & determine the most efficient path to connect all cells (in our case across all partitions). A trajectory passing through all partitions is formed. 
# 3. We determine a root data point. We can also let the algo choose a root based on which point around the centre has the max. no. of connections. 
# 4. Pseudotime (unitless quantity representing developmental path traveled by the cell so far) is calculated by (distance between cell A to cell B)+(distance between cell B to cell C) where cell A -> cell B -> cell C on the knn graph (biological relevance unknown). This is how we quantify the distances of all other cells from the root.

# Learn trajectory graph across all partitions. If use_partition = TRUE, learns trajectory graph within each partition.
cds = learn_graph(cds, use_partition = FALSE)

# Order cells along pseudotime trajectory
cds = order_cells(cds)

# Plot cells colored by pseudotime values
plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = FALSE)

# Add gene name information to Monocle3 object
rowData(cds)$gene_name = rownames(cds)
rowData(cds)$gene_short_name = rowData(cds)$gene_name

colData(cds)$cell_name = colnames(cds)
colData(cds)$gene_short_name = rowData(cds)$gene_name

# Plot expression of specific marker genes
plot_cells(cds,
           genes = c('Sox2', 'Nanog', 'Col6a2'),
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE, 
           min_expr = 3)
           
saveRDS(cds, file = "monocle_obj")
           
# Perform graph-based differential expression testing along pseudotime
cds_pt_res = graph_test(cds, neighbor_graph = "principal_graph", cores = 8)

saveRDS(cds_pt_res, file = "cds_pt_res")

# Load pre-computed pseudotime results (commented out - actual loading)
# cds_pt_res = readRDS("cds_pt_res.rds")

# Display pseudotime results
cds_pt_res

# Filter pseudotime results: remove NAs and keep significant results
cds_pt_res = na.omit(cds_pt_res)
cds_pt_res = cds_pt_res[cds_pt_res$p_value < 0.05 & cds_pt_res$status == "OK", ]
cds_pt_res

# Sort genes by Moran's I test statistic (spatial autocorrelation along trajectory)
cds_pt_res[order(-cds_pt_res$morans_test_statistic),] 

# Visualize top differentially expressed genes along trajectory
plot_cells(cds, genes = c("Col1a2", "Uba52", "Serpine1", "Dppa5a"),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)
           
# Genes that change early in trajectory (for detecting early markers of diseases)
early_genes = cds_pt_res[cds_pt_res$morans_I > 0.5 & 
                         cds_pt_res$q_value < 0.01, ]
early_genes = early_genes[order(early_genes$pseudotime_dependence), ]

# Genes expressed in advanced disease states
late_genes = cds_pt_res[cds_pt_res$morans_I > 0.5 & 
                        cds_pt_res$q_value < 0.01, ]
late_genes = late_genes[order(-late_genes$pseudotime_dependence), ]

# Find genes defining different disease outcomes
branch_genes = graph_test(cds, neighbor_graph = "principal_graph", cores = 8)
branch_specific = branch_genes[branch_genes$q_value < 0.01 & 
                              branch_genes$morans_I > 0.3, ]
           
# Interactively select subset of cells for focused analysis
cds_subset = choose_cells(cds)

# Display information about the selected subset
cds_subset

# Perform pseudotime analysis on the selected subset
cds_subset_pt_res = graph_test(cds_subset, neighbor_graph = "principal_graph", cores = 8)
cds_subset_pt_res = na.omit(cds_subset_pt_res)
cds_subset_pt_res = cds_subset_pt_res[cds_subset_pt_res$p_value < 0.05 & cds_subset_pt_res$status == "OK", ]
cds_subset_pt_res

# Sort subset results by Moran's I statistic
cds_subset_pt_res[order(-cds_subset_pt_res$morans_test_statistic),] 

# Visualize expression patterns in the subset
plot_cells(cds_subset, genes = c("Rpl7a", "Eef1a1", "Mgst1", "Lgals1"),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)
           
# Create a subset containing only specific genes of interest
cds_subset_subset = cds_subset[rowData(cds_subset)$gene_short_name %in% c("Rpl7a", "Eef1a1", "Mgst1", "Lgals1")]

# Plot gene expression dynamics along pseudotime
plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by = "dayint",
                         min_expr = 0.5)
                         
# Integration Example: Select features for integration across datasets
features = SelectIntegrationFeatures(object.list = data_list)

# Define function to scale data and run PCA for each dataset
scale_pca <- function(x) {
  x = ScaleData(x, features = features, verbose = FALSE)
  x = RunPCA(x, features = features, verbose = FALSE)
  return(x)
}

# Apply scaling and PCA to all datasets in the list
data_list = lapply(X = data_list, scale_pca)

# Find integration anchors using reciprocal PCA
anchors = FindIntegrationAnchors(object.list = data_list, anchor.features = features, reduction = "rpca")

# Save anchors for future use
saveRDS(anchors, file = "integration_anchors.rds")

# Integrate datasets using the found anchors
data = IntegrateData(anchorset = anchors)
