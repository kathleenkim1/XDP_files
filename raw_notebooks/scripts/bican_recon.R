
library(qs)
library(Seurat)
library(UCell)


bican_recon = qread("bican_recon_feb2025.qs")

options(future.globals.maxSize = 280 * 1024^3)

library(sctransform)
bican_recon = SCTransform(bican_recon, vars.to.regress = "pct_mito", verbose = FALSE)
DefaultAssay(bican_recon) = "SCT"
qsave(bican_recon, "bican_recon_feb2025_sct.qs")

# library(qs)
# library(Seurat)
# library(sctransform)
# library(future)
# 
# # Increase parallelization
# plan("multicore", workers = 8)  # Adjust based on available cores
# options(future.globals.maxSize = 300 * 1024^3)
# 
# # Load the large Seurat object
# bican_recon <- qread("bican_recon_feb2025.qs")
# 
# # Run SCTransform in chunks
# cell_subsets <- SplitObject(bican_recon, split.by = "orig.ident") 
# 
# subset_list <- lapply(cell_subsets, function(sub) {
#   SCTransform(sub, vars.to.regress = "pct_mito", verbose = FALSE, vst.flavor = "v2")
# })
# 
# # Merge subsets back
# bican_recon_sct <- merge(subset_list[[1]], y = subset_list[-1])
# DefaultAssay(bican_recon_sct) <- "SCT"
# 
# # Save output
# qsave(bican_recon_sct, "bican_recon_feb2025_sct.qs")