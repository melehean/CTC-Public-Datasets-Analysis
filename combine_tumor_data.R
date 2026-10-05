library(Seurat)
library(data.table)
library(Matrix)
library(future)

plan(sequential)
options(future.globals.maxSize = Inf)

# Wczytanie danych
all_tumor_data <- as.data.table(
  readRDS("./data_to_combine/tumor_data.rds"),
  keep.rownames = "cell_id"
)

print("All tumor data read")

# Macierz counts: geny × komórki
mat <- t(as.matrix(all_tumor_data[, !"cell_id", with = FALSE]))
mat <- as(mat, "dgCMatrix")

# Nazwy komórek
colnames(mat) <- all_tumor_data$cell_id

rm(all_tumor_data)
gc()

print(dim(mat))

# Utworzenie obiektu Seurat
seurat_object <- CreateSeuratObject(counts = mat)

rm(mat)
gc()

# Zapis raw counts
saveRDS(
  seurat_object[["RNA"]]$counts,
  "./ready_data/combined_tumor_data_raw_fixed.rds"
)

print("Seurat object created")

# Normalizacja
seurat_object <- SCTransform(
  seurat_object,
  method = "glmGamPoi",
  vars.to.regress = NULL,
  verbose = TRUE,
  conserve.memory = TRUE,
  vst.flavor = "v2"
)

print("Data normalized")

# Zapis znormalizowanych danych
saveRDS(
  seurat_object[["SCT"]]$data,
  "./ready_data/combined_tumor_data_fixed.rds"
)

