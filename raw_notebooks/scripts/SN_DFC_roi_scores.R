
library(qs)
library(Seurat)
library(UCell)

SN_sobj = qread("~/temphome/SN_DFC_sobjs_labeled/SN_sobj_all.qs")
DFC_sobj = qread("~/temphome/SN_DFC_sobjs_labeled/DFC_sobj_all.qs")


up_astrocyte_roi_genes= qread("up_astrocyte_roi_genes.qs")
up_endo_roi_genes = qread("up_endo_roi_genes.qs")
up_oligo_roi_genes = qread("up_oligo_roi_genes.qs")
up_opc_roi_genes = qread("up_opc_roi_genes.qs")
up_neuron_roi_genes = qread("up_neuron_roi_genes.qs")
up_micro_roi_genes = qread("up_micro_roi_genes.qs")

down_astrocyte_roi_genes= qread("down_astrocyte_roi_genes.qs")
down_oligo_roi_genes = qread("down_oligo_roi_genes.qs")
down_neuron_roi_genes = qread("down_neuron_roi_genes.qs")
down_micro_roi_genes = qread("down_micro_roi_genes.qs")

options(future.globals.maxSize = 300 * 1024^3)

library(sctransform)
SN_sobj = SCTransform(SN_sobj, vars.to.regress = "pct_mito", verbose = FALSE)
SN_sobj
DefaultAssay(SN_sobj) = "SCT"

SN_sobj <- AddModuleScore_UCell(SN_sobj,features = list(ROI_upreg = up_astrocyte_roi_genes),
                                                                 name = '_astrocyte'
)

SN_sobj <- AddModuleScore_UCell(SN_sobj, features = list(ROI_upreg = up_endo_roi_genes),
                                                                 name = '_endo'
)

SN_sobj <- AddModuleScore_UCell(SN_sobj, features = list(ROI_upreg = up_oligo_roi_genes),
                                                                 name = '_oligo'
)

SN_sobj <- AddModuleScore_UCell(SN_sobj, features = list(ROI_upreg = up_opc_roi_genes),
                                                                 name = '_opc'
)

SN_sobj <- AddModuleScore_UCell(SN_sobj,features = list(ROI_upreg = up_neuron_roi_genes),
                                                                 name = '_neuron'
)

SN_sobj <- AddModuleScore_UCell(SN_sobj,features = list(ROI_upreg = up_micro_roi_genes),
                                                                 name = '_mg'
)


#down
SN_sobj <- AddModuleScore_UCell(SN_sobj, features = list(ROI_downreg = down_astrocyte_roi_genes),
                                                                 name = '_astrocyte'
)

SN_sobj <- AddModuleScore_UCell(SN_sobj,features = list(ROI_downreg = down_oligo_roi_genes),
                                                                 name = '_oligo'
)


SN_sobj <- AddModuleScore_UCell(SN_sobj,features = list(ROI_downreg = down_neuron_roi_genes),
                                                                 name = '_neuron'
)

SN_sobj <- AddModuleScore_UCell(SN_sobj, features = list(ROI_downreg = down_micro_roi_genes),
                                                                 name = '_mg'
)

qsave(SN_sobj, "SN_sobj_sct_scores.qs")



DFC_sobj = SCTransform(DFC_sobj, vars.to.regress = "pct_mito", verbose = FALSE)
DFC_sobj
DefaultAssay(DFC_sobj) = "SCT"

DFC_sobj <- AddModuleScore_UCell(DFC_sobj,features = list(ROI_upreg = up_astrocyte_roi_genes),
                                name = '_astrocyte'
)

DFC_sobj <- AddModuleScore_UCell(DFC_sobj, features = list(ROI_upreg = up_endo_roi_genes),
                                name = '_endo'
)

DFC_sobj <- AddModuleScore_UCell(DFC_sobj, features = list(ROI_upreg = up_oligo_roi_genes),
                                name = '_oligo'
)

DFC_sobj <- AddModuleScore_UCell(DFC_sobj, features = list(ROI_upreg = up_opc_roi_genes),
                                name = '_opc'
)

DFC_sobj <- AddModuleScore_UCell(DFC_sobj,features = list(ROI_upreg = up_neuron_roi_genes),
                                name = '_neuron'
)

DFC_sobj <- AddModuleScore_UCell(DFC_sobj,features = list(ROI_upreg = up_micro_roi_genes),
                                name = '_mg'
)


#down
DFC_sobj <- AddModuleScore_UCell(DFC_sobj, features = list(ROI_downreg = down_astrocyte_roi_genes),
                                name = '_astrocyte'
)

DFC_sobj <- AddModuleScore_UCell(DFC_sobj,features = list(ROI_downreg = down_oligo_roi_genes),
                                name = '_oligo'
)


DFC_sobj <- AddModuleScore_UCell(DFC_sobj,features = list(ROI_downreg = down_neuron_roi_genes),
                                name = '_neuron'
)

DFC_sobj <- AddModuleScore_UCell(DFC_sobj, features = list(ROI_downreg = down_micro_roi_genes),
                                name = '_mg'
)

qsave(DFC_sobj, "DFC_sobj_sct_scores.qs")