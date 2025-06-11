
library(qs)
library(Seurat)
library(UCell)

clean_recon_sobj_with_neuron_subclusters = qread("clean_recon_sobj_with_neuron_subclusters_sct.qs")
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

clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_upreg = up_astrocyte_roi_genes),
                                                                 name = '_astrocyte'
)

clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_upreg = up_endo_roi_genes),
                                                                 name = '_endo'
)

clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_upreg = up_oligo_roi_genes),
                                                                 name = '_oligo'
)

clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_upreg = up_opc_roi_genes),
                                                                 name = '_opc'
)

clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_upreg = up_neuron_roi_genes),
                                                                 name = '_neuron'
)

clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_upreg = up_micro_roi_genes),
                                                                 name = '_mg'
)


#down
clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_downreg = down_astrocyte_roi_genes),
                                                                 name = '_astrocyte'
)

clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_downreg = down_oligo_roi_genes),
                                                                 name = '_oligo'
)


clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_downreg = down_neuron_roi_genes),
                                                                 name = '_neuron'
)

clean_recon_sobj_with_neuron_subclusters <- AddModuleScore_UCell(clean_recon_sobj_with_neuron_subclusters,
                                                                 features = list(ROI_downreg = down_micro_roi_genes),
                                                                 name = '_mg'
)


qsave(clean_recon_sobj_with_neuron_subclusters, "clean_recon_sobj_with_neuron_subclusters_sct_new.qs")