library(Matrix)
library(SingleCellExperiment)
library(basilisk)

basilisk::setBasiliskShared(FALSE)

library(zellkonverter)

# ============================================================
# WCZYTAJ PONOWNIE ORYGINALNE RDS
# ============================================================

combined_healthy_data <- readRDS(
  "./ready_data/combined_healthy_data_fixed.rds"
)

combined_tumor_data <- readRDS(
  "./ready_data/combined_tumor_data_fixed.rds"
)
combined_tumor_data <- t(combined_tumor_data)

combined_ctc_data <- readRDS(
  "./ready_data/ctc_normalized_data_fixed.rds"
)
combined_ctc_data <- t(combined_ctc_data)


# ============================================================
# WSPÓLNE GENY
# Macierze są: cells x genes
# ============================================================

gn_h <- colnames(combined_healthy_data)
gn_t <- colnames(combined_tumor_data)
gn_c <- colnames(combined_ctc_data)

gene_lists <- list(gn_h, gn_t, gn_c)
gene_lists <- gene_lists[order(sapply(gene_lists, length))]

common_genes <- Reduce(intersect, gene_lists)

cat("Common genes:", length(common_genes), "\n")


# ============================================================
# USTAW IDENTYCZNĄ KOLEJNOŚĆ GENÓW
# ============================================================

combined_healthy_data <- combined_healthy_data[
  , common_genes, drop = FALSE
]

combined_tumor_data <- combined_tumor_data[
  , common_genes, drop = FALSE
]

combined_ctc_data <- combined_ctc_data[
  , common_genes, drop = FALSE
]


# kontrola
stopifnot(
  identical(
    colnames(combined_healthy_data),
    colnames(combined_tumor_data)
  )
)

stopifnot(
  identical(
    colnames(combined_healthy_data),
    colnames(combined_ctc_data)
  )
)

print(dim(combined_healthy_data))
print(dim(combined_tumor_data))
print(dim(combined_ctc_data))

head(rownames(combined_healthy_data))  # komórki
head(colnames(combined_healthy_data))  # geny


# ============================================================
# ZAPIS H5AD
# ============================================================

save_as_h5ad <- function(mat, output_file) {
  
  # mat jest cells x genes
  # SingleCellExperiment potrzebuje genes x cells
  sce <- SingleCellExperiment(
    assays = list(
      X = t(mat)
    )
  )
  
  zellkonverter::writeH5AD(
    sce,
    file = output_file,
    X_name = "X",
    version = "0.11.4"
  )
  
  rm(sce)
  gc()
}


save_as_h5ad(
  combined_healthy_data,
  "./ready_data/combined_healthy_data_fixed.h5ad"
)

save_as_h5ad(
  combined_tumor_data,
  "./ready_data/combined_tumor_data_fixed.h5ad"
)

save_as_h5ad(
  combined_ctc_data,
  "./ready_data/ctc_normalized_data_fixed.h5ad"
)