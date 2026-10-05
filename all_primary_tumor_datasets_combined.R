library(Seurat)
library(copykat)
library(dplyr)
library(ggplot2)
library(UCell)
library(patchwork)
library(stringr)
library(Matrix)
library(data.table)
library(parallel)
library(purrr)

tumor_df_list <- list()

################################################################################
# Breast cancer GSE176078 - 24489 tumor cells
################################################################################

print("#######################################################################")
print("Analysing GSE176078 - breast cancer")
print("#######################################################################")

#list_dirs <- list.dirs("/bigdata/projects_michal/Public_CTC_Tumor_PBMC_Data/breast_cancer_primary_tumor/GSE176078_RAW", full.names = TRUE)
list_dirs <- list.dirs("/home/michal_sieczczynski/GUMED/CTC_public_datasets/data/breast_cancer_primary_tumor/GSE176078_RAW", full.names = TRUE)
list_dirs <- list_dirs[grepl("/CID", list_dirs)]

process_sample <- function(dir) {
  metadata <- fread(file.path(dir, "metadata.csv"))
  print(View(metadata))
  tumor_barcodes <- metadata[celltype_major == "Cancer Epithelial", V1]

  if (length(tumor_barcodes) == 0) return(NULL)

  counts <- ReadMtx(
    cells = file.path(dir, "count_matrix_barcodes.tsv"),
    mtx = file.path(dir, "count_matrix_sparse.mtx"),
    features = file.path(dir, "count_matrix_genes.tsv"),
    feature.column = 1
  )
  counts <- counts[, colnames(counts) %in% tumor_barcodes]
  sample_name <- basename(dir)
  colnames(counts) <- paste0(sample_name, "_", colnames(counts))

  return(counts)
}

tumor_counts_list <- lapply(list_dirs, process_sample)
tumor_counts_list <- Filter(Negate(is.null), tumor_counts_list)
combined_counts <- do.call(cbind, tumor_counts_list)
rm(tumor_counts_list)

gse_176078_seurat <- CreateSeuratObject(counts = combined_counts)
rm(combined_counts)
gc()

gse_176078_seurat[["percent.mt"]] <- PercentageFeatureSet(gse_176078_seurat, pattern = "^MT-")
gse_176078_seurat <- subset(gse_176078_seurat,
                            subset = nFeature_RNA > 200 &
                              nCount_RNA > 1000 &
                              percent.mt < 20)

tumor_df_list[["gse_176078"]] <- t(gse_176078_seurat@assays$RNA$counts)
rm(gse_176078_seurat)
gc()

print("#######################################################################")
print("GSE176078 analyzed - breast cancer")
print("#######################################################################")

################################################################################
# Breast cancer OUR SPIKES - 1886 tumor cells
################################################################################

print("#######################################################################")
print("Analysing OUR SPIKES")
print("#######################################################################")

#list_dirs <- list.dirs("/bigdata/projects_michal/Public_CTC_Tumor_PBMC_Data/breast_cancer_primary_tumor/OUR_SPIKES", full.names = TRUE)
list_dirs <- list.dirs("/home/michal_sieczczynski/GUMED/CTC_public_datasets/data/breast_cancer_primary_tumor/OUR_SPIKES", full.names = TRUE)
list_dirs <- list_dirs[grepl("filtered_feature_bc_matrix", list_dirs)]

process_spike_sample <- function(dir) {

  data <- Read10X(dir)
  obj <- CreateSeuratObject(counts = data)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 20)

  ck <- copykat(
    rawmat = as.matrix(obj@assays$RNA$counts),
    id.type = "S",
    ngene.chr = 5,
    win.size = 25,
    n.cores = 10,
    output.seg = FALSE,
    plot.genes = FALSE
  )

  pred_df <- data.frame(
    cell.names = ck$prediction$cell.names,
    copykat_raw = ck$prediction$copykat.pred
  )

  status_map <- c("aneuploid" = "Tumor (Aneuploid)", "diploid" = "Normal (Diploid)")
  pred_df$final_status <- status_map[pred_df$copykat_raw]
  pred_df$final_status[is.na(pred_df$final_status)] <- "Unknown/Low Quality"

  rownames(pred_df) <- pred_df$cell.names
  obj <- AddMetaData(obj, metadata = pred_df)

  signatures <- list(Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19",
                                    "CDH1", "CLDN1", "DSP", "OCLN"))
  obj <- AddModuleScore_UCell(obj, features = signatures)

  meta <- obj@meta.data
  confirmed_cells <- rownames(meta[meta$Epithelial_UCell > 0.4 &
                                     meta$final_status == "Tumor (Aneuploid)", ])
 
  return(obj)# (obj@assays$RNA$counts[, confirmed_cells, drop = FALSE])
}

our_spikes_results <- lapply(list_dirs, process_spike_sample)

our_spikes_results <- Filter(Negate(is.null), our_spikes_results)

if(length(our_spikes_results) > 0) {
  our_spikes_data <- do.call(cbind, our_spikes_results)
  tumor_df_list[["our_spikes"]] <- t(our_spikes_data)
}

print("#######################################################################")
print("OUR SPIKES analyzed")
print("#######################################################################")

################################################################################
# Liver cancer GSE125449 - 1992 tumor cells
################################################################################

print("#######################################################################")
print("Analysing GSE125449 - liver cancer")
print("#######################################################################")

#liver_data_dir <- "/bigdata/projects_michal/Public_CTC_Tumor_PBMC_Data/liver_cancer_primary_tumor/GSE125449"
liver_data_dir <- "/home/michal_sieczczynski/GUMED/CTC_public_datasets/data/liver_cancer_primary_tumor/GSE125449"
list_dirs <- list.dirs(liver_data_dir, full.names = TRUE)
list_dirs <- list_dirs[grepl("Set", list_dirs)]

process_liver_set <- function(dir) {
  metadata <- fread(file.path(dir, "GSE125449_samples.txt"), sep = "\t")
  tumor_barcodes <- metadata[Type == "Malignant cell", `Cell Barcode`]

  if (length(tumor_barcodes) == 0) return(NULL)

  counts <- ReadMtx(
    cells = file.path(dir, "GSE125449_barcodes.tsv"),
    mtx = file.path(dir, "GSE125449_matrix.mtx"),
    features = file.path(dir, "GSE125449_genes.tsv"),
    feature.column = 2
  )
  counts <- counts[, colnames(counts) %in% tumor_barcodes, drop = FALSE]
  sample_name <- basename(dir)
  colnames(counts) <- paste0(sample_name, "_", colnames(counts))

  return(counts)
}

gse_125449_list <- lapply(list_dirs, process_liver_set)
gse_125449_list <- Filter(Negate(is.null), gse_125449_list)
all_gene_lists <- lapply(gse_125449_list, rownames)
common_genes <- Reduce(intersect, all_gene_lists)
gse_125449_list <- lapply(gse_125449_list, function(mat) {
  mat[common_genes, , drop = FALSE]
})

combined_counts <- do.call(cbind, gse_125449_list)
rm(gse_125449_list)
gc()

gse_125449_seurat <- CreateSeuratObject(counts = combined_counts)
rm(combined_counts)
gc()

gse_125449_seurat[["percent.mt"]] <- PercentageFeatureSet(gse_125449_seurat,
                                                          pattern = "^MT-")
gse_125449_seurat <- subset(gse_125449_seurat,
                            subset = nFeature_RNA > 200 &
                              nCount_RNA > 1000 &
                              percent.mt < 20)

tumor_df_list[["gse_125449"]] <- t(gse_125449_seurat@assays$RNA$counts)
rm(gse_125449_seurat)
gc()

print("#######################################################################")
print("GSE125449 analyzed - liver cancer")
print("#######################################################################")

################################################################################
# Prostate cancer GSE181294 - deleted normal cells and healthy donors and low stage patients
# 1,237 tumor cells 
################################################################################

print("#######################################################################")
print("Analysing GSE181294 - prostate cancer")
print("#######################################################################")

#prostate_tumor_data_path <- "/bigdata/projects_michal/Public_CTC_Tumor_PBMC_Data/prostate_cancer_primary_tumor/GSE181294_RAW"
prostate_tumor_data_path <- "/home/michal_sieczczynski/GUMED/CTC_public_datasets/data/prostate_cancer_primary_tumor/GSE181294_RAW"
prostate_tumor_dir_list <- list.dirs(prostate_tumor_data_path)
prostate_tumor_annotation_path = paste0(prostate_tumor_data_path, "/", "GSE181294_scRNAseq.ano.csv")
prostate_tumor_sample_info <- read.csv(prostate_tumor_annotation_path)
rownames(prostate_tumor_sample_info) <- prostate_tumor_sample_info$X
tumor_cells_indices <- rownames(prostate_tumor_sample_info[prostate_tumor_sample_info$cells == "Tumor",])

gse_181294_data <- list()

for(dir in prostate_tumor_dir_list)
{
  if(!is.na(str_extract_all(dir, "S\\d+")[[1]][1]))
  {
    file_list <- list.files(dir)
    barcode_file <- file_list[grep("barcode", file_list)]
    gene_file <- file_list[grep("genes", file_list)]
    matrix_file <- file_list[grep("count", file_list)]

    data <- ReadMtx(
      cells = paste0(dir, "/", barcode_file),
      mtx = paste0(dir, "/", matrix_file) ,
      features = paste0(dir, "/", gene_file),
      feature.column = 1,
    )
    data <- t(data)
  }
  else if(grepl("SCG-PCA", dir))
  {
    data_file_name <- list.files(dir)
    data_file_path <- paste0(dir, "/", data_file_name)
    data <- read.csv(data_file_path)
    rownames(data) <- data$X
    data$X <- NULL
    data <- t(data)
  }
  else
  {
    next
  }

  rownames(data) <- gsub("\\.", "-", rownames(data))
  common_cells <- intersect(rownames(data), tumor_cells_indices)

  if(length(common_cells) > 0)
  {
    sample_name <- tail(strsplit(dir, "/")[[1]], 1)
    gse_181294_data[[sample_name]] <- data[common_cells, , drop=FALSE]
  }
}

gse_181294_data <- do.call(rbind, gse_181294_data)

gse_181294_seurat <- CreateSeuratObject(counts=t(gse_181294_data))
gse_181294_seurat[["percent.mt"]] <- PercentageFeatureSet(gse_181294_seurat, pattern = "^MT-")
gse_181294_seurat <- subset(gse_181294_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 20)

tumor_df_list[["gse_181294"]] <- t(gse_181294_seurat@assays$RNA$counts)
rm(gse_181294_seurat)
gc()

print("#######################################################################")
print("GSE181294 analysed - prostate cancer")
print("#######################################################################")

################################################################################
# Lung cancer GSE131907 - 24784 tumor cells
################################################################################

print("#######################################################################")
print("Analysing GSE131907 - lung cancer")
print("#######################################################################")

#sample_info_path <- "/bigdata/projects_michal/Public_CTC_Tumor_PBMC_Data/lung_cancer_primary_tumor/GSE131907_Lung_Cancer_cell_annotation.txt"

sample_info_path <- "/home/michal_sieczczynski/GUMED/CTC_public_datasets/data/lung_cancer_primary_tumor/GSE131907_Lung_Cancer_cell_annotation.txt"
gse_131907_sample_info <- fread(sample_info_path)

tumor_cell_indices <- gse_131907_sample_info[Cell_subtype == "Malignant cells", Index]
rm(gse_131907_sample_info)

#data_path <- "/bigdata/projects_michal/Public_CTC_Tumor_PBMC_Data/lung_cancer_primary_tumor/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
data_path <- "/home/michal_sieczczynski/GUMED/CTC_public_datasets/data/lung_cancer_primary_tumor/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
full_data <- readRDS(data_path)

counts <- full_data[,colnames(full_data) %in% tumor_cell_indices]

rm(full_data)
gc()

nCounts <- colSums(counts)
nFeatures <- colSums(counts > 0)

mt_genes <- grep("^MT-", rownames(counts), value = TRUE, ignore.case = TRUE)
if(length(mt_genes) > 0) {
  pct_mt <- colSums(counts[mt_genes, , drop = FALSE]) / nCounts * 100
} else {
  pct_mt <- rep(0, ncol(counts))
}

keep_cells <- nFeatures > 200 & nCounts > 1000 & pct_mt < 20
filtered_counts <- counts[, keep_cells]
tumor_df_list[["gse_131907"]] <- t(filtered_counts)

rm(counts, filtered_counts)
gc()

print("#######################################################################")
print("GSE131907 analysed - lung cancer")
print("#######################################################################")

################################################################################
# Pancreatic cancer GSE194247 
################################################################################

print("#######################################################################")
print("Analysing GSE194247 - pancreatic cancer")
print("#######################################################################")

#pancreatic_cancer_data_path <- "/bigdata/projects_michal/Public_CTC_Tumor_PBMC_Data/pancreatic_cancer_primary_tumor/GSE194247_RAW"
pancreatic_cancer_data_path <- "/home/michal_sieczczynski/GUMED/CTC_public_datasets/data/pancreatic_cancer_primary_tumor/GSE194247_RAW"
pancreatic_cancer_dir_list <- list.dirs(pancreatic_cancer_data_path, full.names = TRUE)
pancreatic_cancer_dir_list <- pancreatic_cancer_dir_list[grepl("filtered_feature_bc_matrix", pancreatic_cancer_dir_list)]

process_pancreatic_sample <- function(dir) {

  data <- Read10X(dir)
  sample_id <- basename(dirname(dir))
  colnames(data) <- paste0(sample_id, "_", colnames(data))
  return(data)
}

gse_194247_list <- lapply(pancreatic_cancer_dir_list, process_pancreatic_sample)
combined_counts <- do.call(cbind, gse_194247_list)
rm(gse_194247_list)
gc()

gse_194247_seurat <- CreateSeuratObject(counts = combined_counts)
rm(combined_counts)
gc()

gse_194247_seurat[["percent.mt"]] <- PercentageFeatureSet(gse_194247_seurat, pattern = "^MT-")
gse_194247_seurat <- subset(gse_194247_seurat,
                            subset = nFeature_RNA > 200 &
                              nCount_RNA > 1000 &
                              percent.mt < 20)

tumor_df_list[["gse_194247"]] <- t(gse_194247_seurat@assays$RNA$counts)

rm(gse_194247_seurat)
gc()

print("#######################################################################")
print("GSE194247 analysed - pancreatic cancer")
print("#######################################################################")

################################################################################
# AGGREGATION
################################################################################

print("#######################################################################")
print("Start aggregation")
print("#######################################################################")

tumor_gene_list <- lapply(tumor_df_list, colnames)
tumor_common_all <- Reduce(intersect, tumor_gene_list)

tumor_df_list <- Map(function(mat, name) {

  mat <- mat[, tumor_common_all, drop = FALSE]
  rownames(mat) <- paste0(name, "_", rownames(mat))
  return(mat)
}, tumor_df_list, names(tumor_df_list))

combined_tumor_data <- do.call(rbind, tumor_df_list)
rm(tumor_df_list)
gc()

saveRDS(combined_tumor_data, file = "./data_to_combine/tumor_data.rds")

print("#######################################################################")
print("Aggregation finished - all data saved")
print("#######################################################################")
