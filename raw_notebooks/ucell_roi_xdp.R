
library(qs)
library(Seurat)
library(UCell)

xdp_cah_put = qread("xdp_cah_put_sct.qs")
astrocyte_roi_genes= qread("astrocyte_roi_genes.qs")
endo_roi_genes = qread("endo_roi_genes.qs")
oligo_roi_genes = qread("oligo_roi_genes.qs")
opc_roi_genes = qread("opc_roi_genes.qs")
neuron_roi_genes = qread("neuron_roi_genes.qs")
micro_roi_genes = qread("micro_roi_genes.qs")

xdp_cah_put <- AddModuleScore_UCell(xdp_cah_put,features = list(ROI = astrocyte_roi_genes),
                                                                 name = '_astrocyte'
)

xdp_cah_put <- AddModuleScore_UCell(xdp_cah_put,features = list(ROI = endo_roi_genes),
                                                                 name = '_endo'
)

xdp_cah_put <- AddModuleScore_UCell(xdp_cah_put,features = list(ROI = oligo_roi_genes),
                                                                 name = '_oligo'
)

xdp_cah_put <- AddModuleScore_UCell(xdp_cah_put,features = list(ROI = opc_roi_genes),
                                                                 name = '_opc'
)

xdp_cah_put <- AddModuleScore_UCell(xdp_cah_put,features = list(ROI = neuron_roi_genes),
                                                                 name = '_neuron'
)

xdp_cah_put <- AddModuleScore_UCell(xdp_cah_put, features = list(ROI = astrocyte_roi_genes),
                                                                 name = '_micro_roi_genes'
)

qsave(xdp_cah_put, "xdp_cah_put_sct_new.qs")