

library(Seurat)
library(qs)
library(tidyverse)
library(sctransform)

final_sobj = qread("final_sobj_wip2.qs")

options(future.globals.maxSize = 400 * 1024^3)  
final_sobj = SCTransform(final_sobj, vars.to.regress = "pct_mito", verbose = FALSE)
DefaultAssay(final_sobj) = "SCT"
final_sobj

qsave(final_sobj, "final_sobj_wipsct.qs")
