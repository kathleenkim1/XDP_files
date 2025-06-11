library(qs)
library(Seurat)
library(UCell)

df_subset = qread("temphome/RECON/Recon_microglia.qs")
Idents(df_subset) = "RNA_snn_res.0.2"
a = FindAllMarkers(df_subset, min.pct = 0.1, logfc.threshold = 0.5)
a_sig = subset(a, subset = p_val_adj < 0.05)
write.csv(a_sig, "microglia_interferon_genes.csv")


df_subset = qread("temphome/RECON/Recon_endo.qs")
Idents(df_subset) = "RNA_snn_res.0.7"
a = FindAllMarkers(df_subset, min.pct = 0.1, logfc.threshold = 0.5)
a_sig = subset(a, subset = p_val_adj < 0.05)
write.csv(a_sig, "endo_interferon_genes.csv")

df_subset = qread("temphome/RECON/Recon_oligo.qs")
Idents(df_subset) = "RNA_snn_res.0.4"
a = FindAllMarkers(df_subset, min.pct = 0.1, logfc.threshold = 0.5)
a_sig = subset(a, subset = p_val_adj < 0.05)
write.csv(a_sig, "oligo_interferon_genes.csv")

df_subset = qread("temphome/RECON/Recon_neuron.qs")
Idents(df_subset) = "RNA_snn_res.0.6"
a = FindAllMarkers(df_subset, min.pct = 0.1, logfc.threshold = 0.5)
a_sig = subset(a, subset = p_val_adj < 0.05)
write.csv(a_sig, "neuron_interferon_genes.csv")