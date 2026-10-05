library(Seurat)
library(data.table)
library(Matrix)
library(future)

plan(sequential)
options(future.globals.maxSize = Inf)

# Wczytanie danych
data <- fread("./data_to_combine/pbmc_data.csv")

# Nazwy genów
gene_names <- data$row_names

# Macierz counts
mat <- as.matrix(data[, !"row_names", with = FALSE])
rownames(mat) <- gene_names
mat <- as(mat, "dgCMatrix")

# Utworzenie obiektu Seurat
seurat_object <- CreateSeuratObject(counts = mat)

# Zapis raw counts
saveRDS(
  seurat_object[["RNA"]]$counts,
  "./ready_data/combined_healthy_data_raw_fixed.rds"
)

rm(data, mat)
gc()

# Normalizacja
seurat_object <- SCTransform(
  seurat_object,
  method = "glmGamPoi",
  vst.flavor = "v2",
  vars.to.regress = NULL,
  conserve.memory = TRUE,
  verbose = TRUE
)

# Zapis znormalizowanych danych
saveRDS(
  seurat_object[["SCT"]]$data,
  "./ready_data/combined_healthy_data_fixed.rds"
)

print("Data normalized and saved")

