
library(harmony)
library(qs)
library(Seurat)
library(UCell)

All_XDP_Cohorts_CaHPut = qread("All_XDP_Cohorts_striatum_sct_final_051525.qs")


library(jsonlite)
json_data <- fromJSON("XDP/QC_and_Clustering/files_needed/query_markers.json")
head(json_data)
names(json_data)


D1_genes <- json_data[["Subclass_label/STR D1 MSN GABA"]]
length(D1_genes)
D2_genes <- json_data[["Subclass_label/STR D2 MSN GABA"]]
length(D2_genes)

combined_SPN = unique(c(D1_genes, D2_genes))
length(combined_SPN)

C_minus = qread("temp_disco/more_temp_1205/final_c_minus_genes.qs")
C_plus = qread("temp_disco/more_temp_1205/final_c_plus_genes.qs")

length(C_minus)
length(C_plus)

C_minus_intersect = intersect(C_minus, combined_SPN)
length(C_minus_intersect)
C_plus_intersect = intersect(C_plus, combined_SPN)
length(C_plus_intersect)

C_minus_intersect
C_plus_intersect
#39 shared genes, C minus
#13 shared genes, C plus

C_minus_new = setdiff(C_minus, combined_SPN)
C_plus_new = setdiff(C_plus, combined_SPN)

length(C_minus_new)
length(C_plus_new)

All_XDP_Cohorts_CaHPut@meta.data

library(UCell)
All_XDP_Cohorts_CaHPut <- AddModuleScore_UCell(All_XDP_Cohorts_CaHPut, features = list(Score = C_minus_new), name = 'C_minus_new')
All_XDP_Cohorts_CaHPut <- AddModuleScore_UCell(All_XDP_Cohorts_CaHPut, features = list(Score = C_plus_new), name = 'C_plus_new')

qsave(All_XDP_Cohorts_CaHPut@meta.data, "All_XDP_Cohorts_CaHPut_cminusplusnew_meta.qs")

rm(All_XDP_Cohorts_CaHPut)

xdp_recon = qread("xdp_recon_apr2025_sct_mask.qs")

xdp_recon <- AddModuleScore_UCell(xdp_recon, features = list(Score = C_minus_new), name = 'C_minus_new')
xdp_recon <- AddModuleScore_UCell(xdp_recon, features = list(Score = C_plus_new), name = 'C_plus_new')

qsave(xdp_recon@meta.data, "xdp_recon_cminusplusnew_meta.qs")
