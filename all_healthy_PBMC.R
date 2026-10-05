library(readxl)
library(dplyr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Seurat)
library(Matrix)

pbmc_df_list <- list()

################################################################################
# HD: All samples are healthy cells
################################################################################

pmbc_10k_counts <- Read10X("data/healthy_donors/SC3_v3_NextGem_SI_PBMC_10K_raw_feature_bc_matrix/raw_feature_bc_matrix/")
pbmc_10k_seurat <- CreateSeuratObject(counts=pmbc_10k_counts)
pbmc_10k_seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc_10k_seurat, pattern = "^MT-")
pbmc_10k_seurat <- subset(pbmc_10k_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

pmbc_10k_counts <- t(pbmc_10k_seurat@assays$RNA$counts)
rm(pbmc_10k_seurat)
pbmc_df_list[["pbmc_10k"]] <- pmbc_10k_counts

################################################################################
# HD: All samples are healthy cells
################################################################################

pbmc_8k_counts <- Read10X("data/healthy_donors/pbmc8k_raw_gene_bc_matrices/raw_gene_bc_matrices/GRCh38/")
pbmc_8k_seurat <- CreateSeuratObject(counts=pbmc_8k_counts)
pbmc_8k_seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc_8k_seurat, pattern = "^MT-")
pbmc_8k_seurat <- subset(pbmc_8k_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

pmbc_8k_counts <- t(pbmc_8k_seurat@assays$RNA$counts)
rm(pbmc_8k_seurat)
pbmc_df_list[["pbmc_8k"]] <- pmbc_8k_counts

################################################################################
# HD: All samples are healthy cells
################################################################################

pbmc_6k_counts <- Read10X("data/healthy_donors/pbmc6k_raw_gene_bc_matrices/matrices_mex/hg19/")
pbmc_6k_seurat <- CreateSeuratObject(counts=pbmc_6k_counts)
pbmc_6k_seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc_6k_seurat, pattern = "^MT-")
pbmc_6k_seurat <- subset(pbmc_6k_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

pbmc_6k_counts <- t(pbmc_6k_seurat@assays$RNA$counts)
rm(pbmc_6k_seurat)
pbmc_df_list[["pbmc_6k"]] <- pbmc_6k_counts

################################################################################
# HD: All samples are healthy cells
################################################################################

pbmc_4k_counts <- Read10X("data/healthy_donors/pbmc4k_raw_gene_bc_matrices/raw_gene_bc_matrices/GRCh38")
pbmc_4k_seurat <- CreateSeuratObject(counts=pbmc_4k_counts)
pbmc_4k_seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc_4k_seurat, pattern = "^MT-")
pbmc_4k_seurat <- subset(pbmc_4k_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

pbmc_4k_counts <- t(pbmc_4k_seurat@assays$RNA$counts)
rm(pbmc_4k_seurat)
pbmc_df_list[["pbmc_4k"]] <- pbmc_4k_counts

################################################################################
# HD: All samples are healthy cells
################################################################################

pbmc_3k_counts <- Read10X("data/healthy_donors/pbmc3k_raw_gene_bc_matrices/raw_gene_bc_matrices/hg19")
pbmc_3k_seurat <- CreateSeuratObject(counts=pbmc_3k_counts)
pbmc_3k_seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc_3k_seurat, pattern = "^MT-")
pbmc_3k_seurat <- subset(pbmc_3k_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

pbmc_3k_counts <- t(pbmc_3k_seurat@assays$RNA$counts)
rm(pbmc_3k_seurat)
pbmc_df_list[["pbmc_3k"]] <- pbmc_3k_counts

################################################################################
# HD: All samples are healthy cells
################################################################################

mat <- readMM("data/healthy_donors/E-HCAD-4-quantification-raw-files/E-HCAD-4.aggregated_filtered_counts.mtx")

genes <- read.table("data/healthy_donors/E-HCAD-4-quantification-raw-files/E-HCAD-4.aggregated_filtered_counts.mtx_rows",
                    stringsAsFactors = FALSE)

cells <- read.table("data/healthy_donors/E-HCAD-4-quantification-raw-files/E-HCAD-4.aggregated_filtered_counts.mtx_cols",
                    stringsAsFactors = FALSE)

gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = genes[,1], 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")
final_names <- ifelse(is.na(gene_symbols), genes, gene_symbols)

rownames(mat) <- final_names
colnames(mat) <- cells[,1]

mat <- mat[!duplicated(rownames(mat)), ]

E_HCAD_4_seurat <- CreateSeuratObject(counts=mat)
E_HCAD_4_seurat[["percent.mt"]] <- PercentageFeatureSet(E_HCAD_4_seurat, pattern = "^MT-")
E_HCAD_4_seurat <- subset(E_HCAD_4_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

mat <- t(E_HCAD_4_seurat@assays$RNA$counts)
rm(E_HCAD_4_seurat)

pbmc_df_list[["E_HCAD_4"]] <- mat

################################################################################
# HD: All samples are healthy cells
################################################################################

our_healthy_donors_path_list <- list.dirs("./data/healthy_donors/OUR_HEALTHY_DONORS/")
our_healthy_data <- list()

for(dir in our_healthy_donors_path_list)
{
  if(grepl("filtered_feature_bc_matrix", dir))
  {
    sample_name <- str_extract_all(dir, "No\\d+_CTRL")[[1]][1]
    data <- Read10X(data.dir = dir)
    data <- t(data)
    our_healthy_data[[sample_name]] <- data
  }
}

list_of_colnames <- lapply(our_healthy_data, colnames)
common_genes <- Reduce(intersect, list_of_colnames)

our_healthy_data <- lapply(our_healthy_data, function(mat) {
  mat[,common_genes, drop = FALSE]
})

our_healthy_data <- lapply(names(our_healthy_data), function(name) {
  df <- our_healthy_data[[name]]
  rownames(df) <- paste0(name, "_", rownames(df))
  return(df)
})
our_healthy_data <- do.call(rbind, our_healthy_data)

our_healthy_data_seurat <- CreateSeuratObject(counts=our_healthy_data)
our_healthy_data_seurat[["percent.mt"]] <- PercentageFeatureSet(our_healthy_data_seurat, pattern = "^MT-")
our_healthy_data_seurat <- subset(our_healthy_data_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

our_healthy_data <- our_healthy_data_seurat@assays$RNA$counts
rm(our_healthy_data_seurat)
pbmc_df_list[["our_healthy_donors"]] <- our_healthy_data

################################################################################
# HD: All samples are healthy cells
################################################################################

library(TENxPBMCData)

sce <- TENxPBMCData(dataset = "pbmc68k")
counts_mat <- counts(sce)

symbols <- rowData(sce)$Symbol
symbols[is.na(symbols)] <- rownames(sce)[is.na(symbols)]
rownames(counts_mat) <- symbols

colnames(counts_mat) <- paste0("Cell_", 1:ncol(counts_mat))
counts_mat <- counts_mat[!duplicated(rownames(counts_mat), fromLast = TRUE), ]

pbmc_68k_seurat <- CreateSeuratObject(counts=counts_mat)
pbmc_68k_seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc_68k_seurat, pattern = "^MT-")
pbmc_68k_seurat <- subset(pbmc_68k_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

pbmc_68k_counts <- t(pbmc_68k_seurat@assays$RNA$counts)
rm(pbmc_68k_seurat)
pbmc_df_list[["pbmc_68k"]] <- pbmc_68k_counts

################################################################################
# AGGREGATION
################################################################################

pbmc_gene_list <- lapply(pbmc_df_list, function(x) colnames(x))
pbmc_common_all <- Reduce(intersect, pbmc_gene_list)

output_file <- "./data_to_combine/pbmc_data.csv"
if (file.exists(output_file)) file.remove(output_file)

for (name in names(pbmc_df_list)) {
  
  temp_dt <- as.data.table(pbmc_df_list[[name]][, pbmc_common_all, drop = FALSE], 
                           keep.rownames = "row_names")
  
  temp_dt[, row_names := paste0(name, "_", row_names)]
  
  fwrite(temp_dt, 
         file = output_file, 
         append = TRUE, 
         col.names = (name == names(pbmc_df_list)[1]))
  
  rm(temp_dt)
}
print("Finished, data saved")
