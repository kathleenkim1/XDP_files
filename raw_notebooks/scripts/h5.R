
library(Seurat)
library(qs)
library(tidyverse)
library(ggplot2)
library(zellkonverter)
library(SeuratDisk)


All_XDP_Cohorts_h5 = qread("final_sobj_wip2.qs")

Idents(All_XDP_Cohorts_h5) = "final_cell_class_RNA"

DefaultAssay(All_XDP_Cohorts_h5) <- "RNA"

# Keep only the RNA assay
# if seurat object version is 5 or later, counts is accessed with $, 4 or earlier with @

if (All_XDP_Cohorts_h5@version >= 5) {
  All_XDP_Cohorts_h5[["RNA"]] <- CreateAssayObject(counts = All_XDP_Cohorts_h5@assays$RNA$counts)
} else {
  All_XDP_Cohorts_h5[["RNA"]] <- CreateAssayObject(counts = All_XDP_Cohorts_h5@assays$RNA@counts)
}


# Remove all other assays
All_XDP_Cohorts_h5@assays <- list(RNA = All_XDP_Cohorts_h5[["RNA"]])

# Ensure metadata is properly aligned
if (!all(rownames(All_XDP_Cohorts_h5@meta.data) == colnames(All_XDP_Cohorts_h5))) {
  stop("Cell names in meta.data do not match counts matrix! Check alignment.")
}

# Save the cleaned Seurat object
SaveH5Seurat(All_XDP_Cohorts_h5, filename = "All_XDP_Cohorts_final.h5Seurat", overwrite = TRUE)

# Convert to H5AD
Convert("All_XDP_Cohorts_final.h5Seurat", dest = "h5ad", overwrite = TRUE)


