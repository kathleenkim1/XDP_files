
library(harmony)
library(qs)
library(Seurat)
library(UCell)

All_XDP_Cohorts_CaHPut = qread("/broad/macosko/kimkathl/striatum_SCT.qs")


bican_recon_oligo_markers_top50_wm = qread("/broad/macosko/kimkathl/spatial_markers_wm_gm/bican_recon_oligo_markers_top50_wm.qs")
bican_recon_oligo_markers_top50_gm= qread("/broad/macosko/kimkathl/spatial_markers_wm_gm/bican_recon_oligo_markers_top50_gm.qs")

bican_recon_opc_markers_top50_wm= qread("/broad/macosko/kimkathl/spatial_markers_wm_gm/bican_recon_opc_markers_top50_wm.qs")
bican_recon_opc_markers_top50_gm= qread("/broad/macosko/kimkathl/spatial_markers_wm_gm/bican_recon_opc_markers_top50_gm.qs")

bican_recon_mg_markers_top50_wm= qread("/broad/macosko/kimkathl/spatial_markers_wm_gm/bican_recon_mg_markers_top50_wm.qs")
bican_recon_mg_markers_top50_gm= qread("/broad/macosko/kimkathl/spatial_markers_wm_gm/bican_recon_mg_markers_top50_gm.qs")

All_XDP_Cohorts_CaHPut <- AddModuleScore_UCell(All_XDP_Cohorts_CaHPut, features = list(Score = bican_recon_oligo_markers_top50_wm), name = 'bican_recon_oligo_markers_top50_wm')
All_XDP_Cohorts_CaHPut <- AddModuleScore_UCell(All_XDP_Cohorts_CaHPut, features = list(Score = bican_recon_oligo_markers_top50_gm), name = 'bican_recon_oligo_markers_top50_gm')

All_XDP_Cohorts_CaHPut <- AddModuleScore_UCell(All_XDP_Cohorts_CaHPut, features = list(Score = bican_recon_opc_markers_top50_wm), name = 'bican_recon_opc_markers_top50_wm')
All_XDP_Cohorts_CaHPut <- AddModuleScore_UCell(All_XDP_Cohorts_CaHPut, features = list(Score = bican_recon_opc_markers_top50_gm), name = 'bican_recon_opc_markers_top50_gm')

All_XDP_Cohorts_CaHPut <- AddModuleScore_UCell(All_XDP_Cohorts_CaHPut, features = list(Score = bican_recon_mg_markers_top50_wm), name = 'bican_recon_mg_markers_top50_wm')
All_XDP_Cohorts_CaHPut <- AddModuleScore_UCell(All_XDP_Cohorts_CaHPut, features = list(Score = bican_recon_mg_markers_top50_gm), name = 'bican_recon_mg_markers_top50_gm')

All_XDP_Cohort_meta = All_XDP_Cohorts_CaHPut@meta.data
All_XDP_Cohort_meta

softmax_cols = function(df, cols){
  z_cols = paste0(cols, "__z")
  softmax_cols = paste0(cols, "__softmax")
  
  # Z-score normalization of the selected columns
  df[z_cols] = scale(df[cols])
  
  # Apply softmax row-wise: exp of z-scores divided by row sum of exp(z-scores)
  df[softmax_cols] =  t(apply(df[z_cols], 1, function(x) exp(x) / sum(exp(x))))
  
  return(df)
}


columns = c("Scorebican_recon_oligo_markers_top50_wm", "Scorebican_recon_oligo_markers_top50_gm",
            "Scorebican_recon_opc_markers_top50_wm", "Scorebican_recon_opc_markers_top50_gm",
            "Scorebican_recon_mg_markers_top50_wm", "Scorebican_recon_mg_markers_top50_gm")


All_XDP_Cohort_meta = softmax_cols(All_XDP_Cohort_meta, columns)
All_XDP_Cohort_meta



All_XDP_Cohort_meta$WM_GM_oligo <- ifelse(
  All_XDP_Cohort_meta$Scorebican_recon_oligo_markers_top50_wm__softmax > All_XDP_Cohort_meta$Scorebican_recon_oligo_markers_top50_gm__softmax, "White Matter - Oligo", 
  ifelse(All_XDP_Cohort_meta$Scorebican_recon_oligo_markers_top50_gm__softmax > All_XDP_Cohort_meta$Scorebican_recon_oligo_markers_top50_wm__softmax, "Gray Matter - Oligo", "Uncertain")
)

All_XDP_Cohort_meta$WM_GM_opc <- ifelse(
  All_XDP_Cohort_meta$Scorebican_recon_opc_markers_top50_wm__softmax > All_XDP_Cohort_meta$Scorebican_recon_opc_markers_top50_gm__softmax, "White Matter - OPC", 
  ifelse(All_XDP_Cohort_meta$Scorebican_recon_opc_markers_top50_gm__softmax > All_XDP_Cohort_meta$Scorebican_recon_opc_markers_top50_wm__softmax, "Gray Matter - OPC", "Uncertain")
)

All_XDP_Cohort_meta$WM_GM_mg <- ifelse(
  All_XDP_Cohort_meta$Scorebican_recon_mg_markers_top50_wm__softmax > All_XDP_Cohort_meta$Scorebican_recon_mg_markers_top50_gm__softmax, "White Matter - OPC", 
  ifelse(All_XDP_Cohort_meta$Scorebican_recon_mg_markers_top50_gm__softmax > All_XDP_Cohort_meta$Scorebican_recon_mg_markers_top50_wm__softmax, "Gray Matter - OPC", "Uncertain")
)
All_XDP_Cohort_meta

qsave(All_XDP_Cohort_meta, "All_XDP_Cohort_meta_wm_gm_score_NEW.qs")
