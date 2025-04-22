Folder: XDP/Analysis/Gray_White_Matter_Scores

This folder should have 2 R notebooks:

**Recon_masks_April2025.Rmd** contains how the bican 3cm recon object and xdp 2cm recon object were drawn masks for caudate, putamen, internal capsule, nucelus accumbens, other (basically any other white matter or other non-CAP structures), and unplaced cells
Last updated: 4/10/2025

Original BICAN sobj --> Updated BICAN sobj with masks
"/broad/macosko/kimkathl/bican_recon_mar2025_sct.qs"" -- > "bican_recon_apr2025_sct_mask.qs" 

Original XDP sobj --> Updated XDP sobj with masks
"/broad/macosko/kimkathl/xdp_recon_full_mar2025.qs" --> "xdp_recon_apr2025_sct_mask.qs" 

The actual masks are saved in **mask_objects** folder for both recons

**Gray_White_Matter_Scores_April2025.Rmd** gets marker genes for WM vs GM for each cell class and the markers are all saved in **bican_marker_list.qs**
Last updated: 4/11/2025
