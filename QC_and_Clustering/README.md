Folder: kimkathl/XDP/QC_and_Clustering 

Last updated: 4/14/25

This folder should have the **map_my_cells_april2025.Rmd** notebook, which has the code for how the 3 current seurat objects (XDP_Cohort_1_2, BICAN_recon, XDP_recon) were converted into .h5ad files and ran with mapmycells script using hmba dataset. Script and dataset links provided by Bennett.
All 3 seurat objects were resaved with the new group, subclass, class identities. Also ran xdp cohort 1 and 2 neurons as a test. 

Something to keep in mind: I compared the labels between the mapmycells outputs for xdp cohort 1 and 2 neurons vs the full seurat object. There are a couple hundred mismatches between labels, but I'd say ~1% mismatch rate or less.

Within this folder there are 6 folders:
  -files_needed: contains all h5ad files used, precomputed stats and reference markers for HMBA_Human_BG reference
  -map_my_cells_xdp_cohort1_2_neurons: mapping output file + DimPlots
  -map_my_cells_xdp_cohort1_2: mapping output file + DimPlots
  -map_my_cells_xdp_recon: mapping output file + DimPlots
  -map_my_cells_xdp_bican_recon: mapping output file + DimPlots
  -cell_type_mapper: not important, just need to run mapmycells script and python