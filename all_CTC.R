library(readxl)
library(dplyr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Seurat)
library(Matrix)
library(ggplot2)

ctc_df_list <- list()
ctc_additional_control <- list()
mouse_samples <- list()

################################################################################
# BC: All samples are CTC
################################################################################
gse_51827_data <- read_xls("data/breast_cancer/GSE51827/GSE51827_readCounts.xls")
gse_51827_sample_info <- read_xls("data/breast_cancer/GSE51827/GSE51827_platform.xls")

gse_51827_data <- gse_51827_data %>%
  left_join(gse_51827_sample_info, by = c("...1" = "ID"))

gse_51827_data <- gse_51827_data %>% dplyr::select(-`...1`, -name, -uniGene, 
                                            -`Entrez GeneID`, -`hg19 knownGene ID`)
gse_51827_data <- gse_51827_data %>%
  group_by(symbol) %>%
  dplyr::slice(1) %>%
  ungroup()
gse_51827_data <- gse_51827_data %>%
  filter(!is.na(symbol))
gse_51827_data <- as.data.frame(gse_51827_data)
rownames(gse_51827_data) <- gse_51827_data$symbol
gse_51827_data <- gse_51827_data %>% dplyr::select(-symbol)
gse_51827_data <- t(gse_51827_data)
rm(gse_51827_sample_info)

ctc_df_list[["gse_51827"]] <- gse_51827_data

################################################################################
# BC: All samples are CTCs
################################################################################
gse_55807_data <- read.csv("data/breast_cancer/GSE55807_read_counts.txt", sep="\t")
gse_55807_data <- gse_55807_data %>% dplyr::select(-`X`, -name, -uniGene, 
                                            -`Entrez.GeneID`, -`hg19.knownGene.ID`)
gse_55807_data <- gse_55807_data %>%
  group_by(symbol) %>%
  dplyr::slice(1) %>%
  ungroup()
gse_55807_data <- gse_55807_data %>%
  filter(!is.na(symbol))
gse_55807_data <- as.data.frame(gse_55807_data)
rownames(gse_55807_data) <- gse_55807_data$symbol
gse_55807_data <- gse_55807_data %>% dplyr::select(-symbol)
gse_55807_data <- t(gse_55807_data)

ctc_df_list[["gse_55807"]] <- gse_55807_data

################################################################################
# BC: Class in row name
################################################################################
gse_67939_data <- read.csv("data/breast_cancer/GSE67939_readCounts.txt", sep="\t")
entrez_ids <- str_extract(gse_67939_data$X, "(?<=eg:)\\d+")
symbols <- mapIds(org.Hs.eg.db,
                  keys = entrez_ids,
                  column = "SYMBOL",
                  keytype = "ENTREZID",
                  multiVals = "first")
mapping_df <- data.frame(original_id = gse_67939_data$X,
                         entrez_id = entrez_ids,
                         symbol = symbols,
                         stringsAsFactors = FALSE)

gse_67939_data <- gse_67939_data %>%
  left_join(mapping_df, by = c("X" = "original_id" ))
gse_67939_data <- gse_67939_data %>% dplyr::select(-entrez_id, -X)
gse_67939_data$symbol <- as.character(gse_67939_data$symbol)
gse_67939_data <- gse_67939_data %>%
  filter(!is.na(symbol))

gse_67939_data <- gse_67939_data %>%
  group_by(symbol) %>%
  dplyr::slice(1) %>%
  ungroup()

gse_67939_data <- as.data.frame(gse_67939_data)
rownames(gse_67939_data) <- gse_67939_data$symbol
gse_67939_data <- gse_67939_data %>% dplyr::select(-symbol)
gse_67939_data <- t(gse_67939_data)

ctc_rownames <- rownames(gse_67939_data)[grepl("Brx", rownames(gse_67939_data))]
ctc_df_list[["gse_67939"]] <- gse_67939_data[ctc_rownames,]

################################################################################
# BC: All samples are CTCs
################################################################################

gse_75367_data <- read.csv("data/breast_cancer/GSE75367_readCounts.txt", sep="\t")
rownames(gse_75367_data) <- gse_75367_data$X
gse_75367_data <- t(gse_75367_data)

ctc_df_list[["gse_75367"]] <- gse_75367_data

################################################################################
# BC: All samples are CTCs
################################################################################

gse_86978_data <- read_xls("data/breast_cancer/GSE86978/GSE86978_readCounts.xls")
gse_86978_sample_info <- read_xls("data/breast_cancer/GSE86978/GSE86978_platform.xls")
gse_86978_data <- gse_86978_data %>%
  left_join(gse_86978_sample_info, by = c("...1" = "ID"))

gse_86978_data <- gse_86978_data %>% dplyr::select(-`...1`, -name, -uniGene, 
                                                   -`Entrez GeneID`, -`hg19 knownGene ID`)
gse_86978_data <- gse_86978_data %>%
  group_by(symbol) %>%
  dplyr::slice(1) %>%
  ungroup()
gse_86978_data <- gse_86978_data %>%
  filter(!is.na(symbol))
gse_86978_data <- as.data.frame(gse_86978_data)
rownames(gse_86978_data) <- gse_86978_data$symbol
gse_86978_data <- gse_86978_data %>% dplyr::select(-symbol)
gse_86978_data <- t(gse_86978_data)
rm(gse_86978_sample_info)

ctc_df_list[["gse_86978"]] <- gse_86978_data

################################################################################
# BC: Class in row name
################################################################################

gse_109761_data <- readRDS("data/breast_cancer/GSE109761_sce_hs.rds")
gse_109761_sample_info <- gse_109761_data@colData["cell_type"]
gse_109761_data <- gse_109761_data@assays@.xData$data$counts
gse_109761_data <- t(gse_109761_data)

classes <- gse_109761_sample_info[rownames(gse_109761_data), "cell_type"]
rownames(gse_109761_data) <- paste(rownames(gse_109761_data), classes, sep = "_")
rm(gse_109761_sample_info)

ctc_rownames <- rownames(gse_109761_data)[grepl("CTC", rownames(gse_109761_data))]
ctc_df_list[["gse_109761"]] <- gse_109761_data[ctc_rownames,]

################################################################################
# BC: Class in row name - MOUSE
################################################################################

gse_109761_data_mouse <- readRDS("data/breast_cancer/GSE109761_sce_mm.rds")
gse_109761_sample_info_mouse <- gse_109761_data_mouse@colData["cell_type"]
gse_109761_data_mouse <- gse_109761_data_mouse@assays@.xData$data$counts
gse_109761_data_mouse <- t(gse_109761_data_mouse)

classes <- gse_109761_sample_info_mouse[rownames(gse_109761_data_mouse), "cell_type"]
rownames(gse_109761_data_mouse) <- paste(rownames(gse_109761_data_mouse), classes, sep = "_")
rm(gse_109761_sample_info_mouse)

ctc_rownames <- rownames(gse_109761_data_mouse)[grepl("CTC", rownames(gse_109761_data_mouse))]
mouse_samples[["gse_109761"]] <- gse_109761_data_mouse[ctc_rownames,]


################################################################################
# BC: All samples are CTCs (deleted those with NA cell type)
################################################################################

folder_path <- "data/breast_cancer/GSE111065_RAW"
file_list <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)
data_list <- lapply(file_list, function(x) {
  read.csv(x, sep = "\t", row.names = 1) 
})
gse_111065_data <- do.call(cbind, data_list)
gse_111065_data <- gse_111065_data %>% 
  dplyr::select(matches("Br|LM"))
gse_111065_data <- t(gse_111065_data)

ctc_df_list[["gse_111065"]] <- gse_111065_data

################################################################################
# BC: Class in row name
################################################################################

folder_path <- "data/breast_cancer/GSE111842_RAW"
file_list <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)
data_list <- lapply(file_list, function(x) {
  df <- read.delim(x, sep = " ", row.names = 1)
  name <- gsub("\\.txt$", "", basename(x))
  colnames(df) <- name
  return(df)
})
gse_111842_data <- do.call(cbind, data_list)

ensembl_ids <- gsub("\\..*$", "", rownames(gse_111842_data))
ensembl_mapping <- select(org.Hs.eg.db, 
                    keys = ensembl_ids, 
                    column = "SYMBOL", 
                    keytype = "ENSEMBL")
ensembl_mapping <- ensembl_mapping[!duplicated(ensembl_mapping$ENSEMBL), ]

gse_111842_data$symbol <- ensembl_mapping$SYMBOL[match(ensembl_ids, ensembl_mapping$ENSEMBL)]
gse_111842_data <- na.omit(gse_111842_data)
gse_111842_data <- gse_111842_data[!duplicated(gse_111842_data$symbol), ]
rownames(gse_111842_data) <- gse_111842_data$symbol
gse_111842_data <- gse_111842_data %>% dplyr::select(-symbol)
gse_111842_data <- t(gse_111842_data)

ctc_rownames <- rownames(gse_111842_data)[grepl("CTC", rownames(gse_111842_data))]
ctc_df_list[["gse_111842"]] <- gse_111842_data[ctc_rownames,]

################################################################################
# BC: All samples are CTCs (deleted those treated with DTX)
################################################################################

gse_261194_data <- read.csv("data/breast_cancer/GSE261194_scRNAseq.counts.csv")
rownames(gse_261194_data) <- gse_261194_data$X
gse_261194_data <- gse_261194_data %>% dplyr::select(-X)

gse_261194_additional_control_data <- gse_261194_data %>% 
  dplyr::select(matches("DTX|Docetaxel"))
gse_261194_additional_control_data <- t(gse_261194_additional_control_data)
ctc_additional_control[["gse_261194"]] <- gse_261194_additional_control_data

gse_261194_data <- gse_261194_data %>% 
  dplyr::select(matches("Control|NT"))
gse_261194_data <- t(gse_261194_data)

ctc_df_list[["gse_261194"]] <- gse_261194_data

################################################################################
# BC: Class in row name
################################################################################

gse_268201_data <- read.csv("data/breast_cancer/GSE268201_Counts_genocode_v26_AMarkiewicz_240523.txt")
rownames(gse_268201_data) <- gse_268201_data$Geneid
gse_268201_data <- gse_268201_data %>% dplyr::select(-Geneid)
ensembl_ids <- gsub("\\..*$", "", rownames(gse_268201_data))
ensembl_mapping <- select(org.Hs.eg.db, 
                          keys = ensembl_ids, 
                          column = "SYMBOL", 
                          keytype = "ENSEMBL")
ensembl_mapping <- ensembl_mapping[!duplicated(ensembl_mapping$ENSEMBL), ]

gse_268201_data$symbol <- ensembl_mapping$SYMBOL[match(ensembl_ids, ensembl_mapping$ENSEMBL)]
gse_268201_data <- na.omit(gse_268201_data)
gse_268201_data <- gse_268201_data[!duplicated(gse_268201_data$symbol), ]
rownames(gse_268201_data) <- gse_268201_data$symbol
gse_268201_data <- gse_268201_data %>% dplyr::select(-symbol)
gse_268201_data <- t(gse_268201_data)

ctc_rownames <- rownames(gse_268201_data)[grepl("CTC|136|152|29|30|31|69|70", rownames(gse_268201_data))]
ctc_df_list[["gse_268201"]] <- gse_268201_data[ctc_rownames,]

################################################################################
# Lung Cancer: All cells are CTC/tumor; 041814 - primary tumor; rest is CTC
################################################################################

gse_74639_data <- read_xls("data/lung_cancer/GSE74639/GSE74639_readCounts.xls")
gse_74639_sample_info <- read_xls("data/lung_cancer/GSE74639/GSE74639_annotation_for_readCounts.xls")

gse_74639_data <- gse_74639_data %>%
  left_join(gse_74639_sample_info, by = c("...1" = "ID"))

gse_74639_data <- gse_74639_data %>% dplyr::select(-`...1`, -name, -uniGene, 
                                                   -`Entrez GeneID`)
gse_74639_data <- gse_74639_data %>%
  group_by(symbol) %>%
  dplyr::slice(1) %>%
  ungroup()
gse_74639_data <- gse_74639_data %>%
  filter(!is.na(symbol))
gse_74639_data <- as.data.frame(gse_74639_data)
rownames(gse_74639_data) <- gse_74639_data$symbol
gse_74639_data <- gse_74639_data %>% dplyr::select(-symbol)
gse_74639_data <- t(gse_74639_data)
rm(gse_74639_sample_info)

ctc_rownames <- rownames(gse_74639_data)[!grepl("041814", rownames(gse_74639_data))]
ctc_df_list[["gse_74639"]] <- gse_74639_data[ctc_rownames,]

################################################################################
# Liver Cancer: All cells are CTC
################################################################################

gse_117623_data <- read.csv("data/liver_cancer/GSE117623_RawCounts_Geo_CLDHCCWBC.csv")
rownames(gse_117623_data) <- gse_117623_data$X
gse_117623_data <- gse_117623_data %>% dplyr::select(-X)
gse_117623_data <- gse_117623_data %>% 
  dplyr::select(matches("HCC"))
gse_117623_data <- t(gse_117623_data)

ctc_df_list[["gse_117623"]] <- gse_117623_data

################################################################################
# Prostate Cancer: Cell type is in row names
################################################################################

gse_67980_data <- read.csv("data/prostate_cancer/GSE67980/GSE67980_readCounts.txt", sep="\t")
gse_67980_data <- na.omit(gse_67980_data)
gse_67980_data <- gse_67980_data[!duplicated(gse_67980_data$symbol), ]
rownames(gse_67980_data) <- gse_67980_data$symbol
gse_67980_data <- gse_67980_data %>% dplyr::select(-ID, -Entrez.GeneID, -uniGene,
                                                   -name, -symbol)
gse_67980_sample_info <- read.csv("data/prostate_cancer/GSE67980/GSE67980_sampleProperties.txt", sep="\t") 
rownames(gse_67980_sample_info) <- gse_67980_sample_info$title

gse_67980_sample_info_additional_data = gse_67980_sample_info[gse_67980_sample_info["source.name"] == "single cell from PCa cell line", ,drop=FALSE]
gse_67980_data_additional_data <- gse_67980_data[,rownames(gse_67980_sample_info_additional_data)]
gse_67980_data_additional_data <- t(gse_67980_data_additional_data)
rownames(gse_67980_data_additional_data) <- paste(rownames(gse_67980_data_additional_data), gse_67980_sample_info_additional_data$source.name, sep = "_")
ctc_additional_control[["gse_67980"]] <- gse_67980_data_additional_data

gse_67980_sample_info = gse_67980_sample_info[gse_67980_sample_info["source.name"] != "single cell from PCa cell line", ,drop=FALSE]
gse_67980_data <- gse_67980_data[,rownames(gse_67980_sample_info)]
gse_67980_data <- t(gse_67980_data)
rownames(gse_67980_data) <- paste(rownames(gse_67980_data), gse_67980_sample_info$source.name, sep = "_")

ctc_rownames <- rownames(gse_67980_data)[grepl("CTC", rownames(gse_67980_data))]
ctc_df_list[["gse_67980"]] <- gse_67980_data[ctc_rownames,]

################################################################################
# Pancreatic Cancer: All cells are CTC
################################################################################

gse_60407_data <- read_xls("data/pancreatic_cancer/GSE60407_readCounts.xls")
gse_60407_sample_info <- read_xls("data/pancreatic_cancer/GSE60407_platform.xls")
gse_60407_data <- gse_60407_data %>%
  left_join(gse_60407_sample_info, by = c("...1" = "ID"))

gse_60407_data <- gse_60407_data %>% dplyr::select(-`...1`, -name, -uniGene, 
                                                   -`Entrez GeneID`, -`hg19 knownGene ID`)
gse_60407_data <- gse_60407_data %>%
  group_by(symbol) %>%
  dplyr::slice(1) %>%
  ungroup()
gse_60407_data <- gse_60407_data %>%
  filter(!is.na(symbol))
gse_60407_data <- as.data.frame(gse_60407_data)
rownames(gse_60407_data) <- gse_60407_data$symbol
gse_60407_data <- gse_60407_data %>% dplyr::select(-symbol)
gse_60407_data <- t(gse_60407_data)
rm(gse_60407_sample_info)

ctc_df_list[["gse_60407"]] <- gse_60407_data

################################################################################
# AGGREGATION
################################################################################

ctc_gene_list <- lapply(ctc_df_list, function(x) colnames(x))
ctc_common_all <- Reduce(intersect, ctc_gene_list)

ctc_df_list <- lapply(ctc_df_list, function(x) {
  return(x[,ctc_common_all, drop = FALSE])
})
ctc_df_list <- lapply(names(ctc_df_list), function(name) {
  df <- ctc_df_list[[name]]
  rownames(df) <- paste0(name, "_", rownames(df))
  return(df)
})
ctc_combined_df <- do.call(rbind, ctc_df_list)

ctc_seurat_object <- CreateSeuratObject(counts=t(ctc_combined_df))
ctc_seurat_object <- subset(ctc_seurat_object, subset = !is.na(nFeature_RNA))
ctc_seurat_object[["percent.mt"]] <- PercentageFeatureSet(ctc_seurat_object, pattern = "^MT-")
ctc_seurat_object <- subset(ctc_seurat_object, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

cell_names <- colnames(ctc_seurat_object)
gse_ids <- sapply(strsplit(cell_names, "_"), function(x) x[2])
ctc_seurat_object$gse_id <- gse_ids

ctc_seurat_object <- SCTransform(ctc_seurat_object, 
                                 method = "glmGamPoi", 
                                 vars.to.regress = NULL, 
                                 verbose = TRUE, 
                                 conserve.memory = TRUE,
                                 vst.flavor = "v2")

saveRDS(ctc_seurat_object@assays$SCT$data, "./ready_data/ctc_normalized_data_fixed.rds")
saveRDS(ctc_seurat_object@assays$RNA$counts, "./ready_data/ctc_raw_data_fixed.rds")

################################################################################
# AGGREGATION - ADDITIONAL SAMPLES
################################################################################

ctc_gene_list <- lapply(ctc_additional_control, function(x) colnames(x))
ctc_common_all <- Reduce(intersect, ctc_gene_list)

ctc_additional_control <- lapply(ctc_additional_control, function(x) {
  return(x[,ctc_common_all, drop = FALSE])
})
ctc_additional_control <- lapply(names(ctc_additional_control), function(name) {
  df <- ctc_additional_control[[name]]
  rownames(df) <- paste0(name, "_", rownames(df))
  return(df)
})
ctc_combined_df <- do.call(rbind, ctc_additional_control)

ctc_seurat_object <- CreateSeuratObject(counts=t(ctc_combined_df))
ctc_seurat_object <- subset(ctc_seurat_object, subset = !is.na(nFeature_RNA))
ctc_seurat_object[["percent.mt"]] <- PercentageFeatureSet(ctc_seurat_object, pattern = "^MT-")
ctc_seurat_object <- subset(ctc_seurat_object, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)

cell_names <- colnames(ctc_seurat_object)
gse_ids <- sapply(strsplit(cell_names, "_"), function(x) x[2])
ctc_seurat_object$gse_id <- gse_ids

ctc_seurat_object <- SCTransform(ctc_seurat_object, 
                                 method = "glmGamPoi", 
                                 vars.to.regress = NULL, 
                                 verbose = TRUE, 
                                 conserve.memory = TRUE,
                                 vst.flavor = "v2")

saveRDS(ctc_seurat_object@assays$SCT$data, "./ready_data/ctc_additional_test_normalized_data.rds")
saveRDS(ctc_seurat_object@assays$RNA$counts, "./ready_data/ctc_additional_test_raw_data.rds")

################################################################################
# NORMALIZATION - MOUSE SAMPLE
################################################################################

mouse_gene_list <- lapply(mouse_samples, function(x) colnames(x))
mouse_common_all <- Reduce(intersect, mouse_gene_list)

mouse_samples <- lapply(mouse_samples, function(x) {
  return(x[, mouse_common_all, drop = FALSE])
})

mouse_samples <- lapply(names(mouse_samples), function(name) {
  df <- mouse_samples[[name]]
  rownames(df) <- paste0(name, "_", rownames(df))
  return(df)
})

mouse_combined_df <- do.call(rbind, mouse_samples)

mouse_seurat_object <- CreateSeuratObject(counts = t(mouse_combined_df))

mouse_seurat_object <- subset(
  mouse_seurat_object,
  subset = !is.na(nFeature_RNA)
)

# Mouse mitochondrial genes: mt-Nd1, mt-Co1, mt-Cytb, etc.
mouse_seurat_object[["percent.mt"]] <- PercentageFeatureSet(
  mouse_seurat_object,
  pattern = "^mt-"
)

mouse_seurat_object <- subset(
  mouse_seurat_object,
  subset = nFeature_RNA > 200 &
    nCount_RNA > 1000 &
    percent.mt < 10
)

# Jeden dataset w named list
mouse_seurat_object$gse_id <- names(mouse_gene_list)[1]

mouse_seurat_object <- SCTransform(
  mouse_seurat_object,
  method = "glmGamPoi",
  vars.to.regress = NULL,
  verbose = TRUE,
  conserve.memory = TRUE,
  vst.flavor = "v2"
)

saveRDS(
  mouse_seurat_object@assays$SCT$data,
  "./ready_data/mouse_additional_test_normalized_data.rds"
)

saveRDS(
  mouse_seurat_object@assays$RNA$counts,
  "./ready_data/mouse_additional_test_raw_data.rds"
)
