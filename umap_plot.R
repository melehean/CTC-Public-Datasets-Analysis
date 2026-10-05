library(Matrix)
library(Seurat)
library(ggplot2)

# 1. Pobierz tylko geny wspólne (oszczędzamy RAM nie trzymając wszystkiego)
# Wczytujemy, wyciągamy rownames i natychmiast usuwamy
m1 <- readRDS("./ready_data/combined_healthy_data.rds")
rn_h <- rownames(m1)
rm(m1); gc()

m2 <- readRDS("./ready_data/combined_tumor_data.rds")
rn_t <- rownames(m2)
rm(m2); gc()

m3 <- readRDS("./ready_data/ctc_normalized_data.rds")
rn_c <- rownames(m3)
rm(m3); gc()

common_genes <- Reduce(intersect, list(rn_h, rn_t, rn_c))

# 2. Funkcja do wczytywania i filtrowania macierzy "w locie"
load_and_subset <- function(path, genes, prefix) {
  mat <- readRDS(path)
  mat <- mat[genes, , drop = FALSE]
  # Dodajemy prefix do komórek od razu na macierzy
  colnames(mat) <- paste0(prefix, "_", colnames(mat))
  return(mat)
}

# 3. Łączymy macierze jedna po drugiej (najbardziej oszczędna metoda)
message("Combining matrices...")
full_mat <- load_and_subset("./ready_data/ctc_normalized_data.rds", common_genes, "CTC")

# Doczytujemy drugą i łączymy
tmp <- load_and_subset("./ready_data/combined_tumor_data.rds", common_genes, "Tumor")
full_mat <- cbind(full_mat, tmp)
rm(tmp); gc()

# Doczytujemy trzecią i łączymy
tmp <- load_and_subset("./ready_data/combined_healthy_data.rds", common_genes, "PBMC")
full_mat <- cbind(full_mat, tmp)
rm(tmp); gc()

# 4. Tworzymy JEDEN obiekt Seurat od razu na połączonej macierzy
# Dzięki temu nie ma "Layers" (counts.1, counts.2 itp.) i JoinLayers jest niepotrzebne
combined <- CreateSeuratObject(counts = full_mat, project = "Combined")
rm(full_mat); gc()

# 5. Odtwarzamy metadane na podstawie nazw komórek
# Skoro dodaliśmy prefiksy w kroku 2, teraz łatwo je wyciągnąć
cells <- Cells(combined)
combined$cell_class <- sub("^([^_]+)_.*", "\\1", cells)
# Batch: wyciągamy gse_123 z formatu PREFIX_gse_123_cellname
combined$batch <- sub("^[^_]+_([^_]+_[^_]+)_.*", "\\1", cells)

# 6. Procesowanie danych
# Jeśli dane były już znormalizowane (log-normalized), przepisujemy counts do data
# W Seurat v5 najlepiej zrobić to tak:
combined[["RNA"]]$data <- combined[["RNA"]]$counts

# FindVariableFeatures na 3000 genów
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)

# ScaleData - ograniczamy tylko do VariableFeatures, żeby zaoszczędzić RAM
combined <- ScaleData(combined, features = VariableFeatures(combined))

# Redukcje wymiarów
combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, dims = 1:30, verbose = FALSE)

# 7. Wykresy (użyj raster = TRUE dla dużej liczby komórek!)
p1 <- DimPlot(combined, group.by = "batch", label = FALSE, raster = TRUE) + 
  ggtitle("Batch Effect")

p2 <- DimPlot(combined, group.by = "cell_class", label = FALSE, raster = TRUE) + 
  ggtitle("Cell Identity")

ggsave("batch_effect.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("cell_identity.png", plot = p1, width = 8, height = 6, dpi = 300)