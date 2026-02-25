#### FIGURE CONFIGURATION 20-02-2026 MARNIX KOOPS


# Load packages and objects -----------------------------------------------


library(future)
plan("multisession", workers = 15) #not too much when working with very big objects
options(future.globals.maxSize = 100 * 1024 * 1024 * 1024)  # Set a different value (e.g., 1000 MiB)

setOldClass("ExpData")

options(bitmapType = "cairo")

# Vector of all packages you want to load
# Vector of all packages, duplicates removed
pkgs <- unique(c(
  "Seurat", "imcRtools", "SpatialExperiment", "ggplot2", "dplyr", "BiocStyle", "qs", 
  "BiocManager", "remotes", "cytomapper", "sf", "scuttle", "scran", "scater", "uwot",
  "rtracklayer", "patchwork", "pheatmap", "gridExtra", "Rphenograph", "igraph", "dittoSeq",
  "viridis", "bluster", "BiocParallel", "readxl", "readr", "hexbin", "openxlsx", "ggrastr",
  "tibble", "grid", "speckle", "limma", "edgeR", "tidyr", "cowplot", "scales", "RColorBrewer", "stringr", "ggpubr"
  ,"lmerTest"
))
# Load all packages
invisible(lapply(pkgs, library, character.only = TRUE))



source('/scripts/readWriteRDS_multiThreaded.R')

output_dir <- "~/Figures/"

Xenium_Object_Seurat <- mcreadRDS("~/Xenium_seurat_object.rds")
spe_Xenium <- mcreadRDS("~/Xenium_as_spatial_experiment.rds" )
spe_IMC <- mcreadRDS("~/IMC_as_spatial_object.rds")



# Colors and orders ------------------------------------------------------

cluster_colors_umap <- c(  #fir annot even less broad
  "Other" = "#999999", 
  "Erythroid" = "#E78AC3",  
  "Endothelium" = "#81C784", 
  "Stroma" = "#3CB44B",                   # Stromal / Myofibroblast
  "Osteolineage" = "#808000",   
  "Adipocytes" = "#A65628", 
  "MKC lineage" = "#D4A5A5", 
  "VSMC" = "#FFD92F" ,# New: brown
  "Interferon stimulated cells" = "#984EA3",  
  "Neutrophil Progenitors" = "#D1E5F0",
  "Neutrophils" = "#46F0F0",
  "Mature neutrophils" = "#91C3E6", # Myelocytes
  "Basophils Eosinophils Mast cells" = "#FB8072",  # Ba Eo Ma
  "Macrophages" = "#4575B4",              # Macrophages
  "Monocytes" = "#6699CC",     
  "Dendritic cells" = "#F781BF",          # New: bright pink
  "T NK cells" = "#FF7F00", 
  "B/Plasma cells" = "#D73027" ,
  "Adipocyte"                      = "#A65628",  # brown
  "MKC Lineage"                    = "#D4A5A5",  # pale red
  "Neutrophil Progenitor"          = "#D1E5F0",  # light blue
  "Neutrophil"                     = "#46F0F0",
  "Mature Neutrophil"              = "#46F0F0",  #same colour for UMAP
  #"Mature Neutrophil"              = "#91C3E6",  # steel blue
  "Basophil/Eosinophil/Mast Cell"  = "#FB8072",  # salmon
  "Macrophage"                     = "#4575B4",  # dark blue
  "Monocyte"                       = "#6699CC",  # medium blue
  "Dendritic Cell"                 = "#F781BF",  # magenta
  "T/NK Cell"                      = "#FF7F00",  # orange
  "B/Plasma Cell"                  = "#D73027"   # red# New: yellow
)

cluster_colors <- c(  #fir annot even less broad
  "Other" = "#999999", 
  "Erythroid" = "#E78AC3",  
  "Endothelium" = "#81C784", 
  "Stroma" = "#3CB44B",                   # Stromal / Myofibroblast
  "Osteolineage" = "#808000",   
  "Adipocytes" = "#A65628", 
  "MKC lineage" = "#D4A5A5", 
  "VSMC" = "#FFD92F" ,# New: brown
  "Interferon stimulated cells" = "#984EA3",  
  "Neutrophil Progenitors" = "#D1E5F0",
  "Neutrophils" = "#46F0F0",
  "Mature neutrophils" = "#91C3E6", # Myelocytes
  "Basophils Eosinophils Mast cells" = "#FB8072",  # Ba Eo Ma
  "Macrophages" = "#4575B4",              # Macrophages
  "Monocytes" = "#6699CC",     
  "Dendritic cells" = "#F781BF",          # New: bright pink
  "T NK cells" = "#FF7F00", 
  "B/Plasma cells" = "#D73027" ,
  "Adipocyte"                      = "#A65628",  # brown
  "MKC Lineage"                    = "#D4A5A5",  # pale red
  "Neutrophil Progenitor"          = "#D1E5F0",  # light blue
  "Neutrophil"                     = "#46F0F0",
  #"Mature Neutrophil"              = "#46F0F0",  #same colour for UMAP
  "Mature Neutrophil"              = "#91C3E6",  # steel blue
  "Basophil/Eosinophil/Mast Cell"  = "#FB8072",  # salmon
  "Macrophage"                     = "#4575B4",  # dark blue
  "Monocyte"                       = "#6699CC",  # medium blue
  "Dendritic Cell"                 = "#F781BF",  # magenta
  "T/NK Cell"                      = "#FF7F00",  # orange
  "B/Plasma Cell"                  = "#D73027"   # red# New: yellow
)


cluster_colors_with_IMC <- c(  #fir annot even less broad
  "Other" = "#999999", 
  "Erythroid" = "#E78AC3",  
  "Endothelium" = "#81C784", 
  "Stroma" = "#3CB44B",                   # Stromal / Myofibroblast
  "Osteolineage" = "#808000",   
  "Adipocytes" = "#A65628", 
  "MKC lineage" = "#D4A5A5", 
  "VSMC" = "#FFD92F" ,# New: brown
  #"Interferon stimulated cells" = "#984EA3",  
  #"Neutrophils" = "#D1E5F0",
  "Neutrophils" = "#46F0F0",
  #"Mature neutrophils" = "#91C3E6", # Myelocytes
  #"Basophils Eosinophils Mast cells" = "#FB8072",  # Ba Eo Ma
  "Myeloid" = "#4575B4",              # Macrophages
  #"Monocytes" = "#6699CC",     
  "Dendritic cells" = "#F781BF",          # New: bright pink
  "T NK cells" = "#FF7F00", 
  "B/Plasma cells" = "#D73027"    # New: yellow
)

celltype_colors_IMC <- c(
  "Undefined_Cell"                       = "lightgray",# Undefined Cell
  "Fibrotic_Tissue"              = "darkgray",  # Fibrogen Undefined Cell
  "Ki67_likely_erythroid"          = "pink",
  "MKC_Endothelial"     = "#90ee90", 
  "CD105_Structural_Cell"       =  "#1C750C",#"#710193",
  "Neutrophil_Progenitor"           =  "#00FFFF",     # Neutrophil Progenitor# Neutrophil
  "Neutrophil"                  = "#00BFFF", 
  "Macrophage"            = "#005BFF",   # Macrophage
  "CD14_Monocyte"               = "#c51b8a",   # Antigen Presenting Cell
  "HLA-DR_Cell"               = "purple",   # Antigen Presenting Cell
  "CD8_T_Cell"                  = "#FA7921",  
  "CD4_T_Cell"                  = "#FFE66D",   # CD4 T cell
  "Plasma_Cell"                 = "#FF6666"      # Plasma Cell
)



status_colors <- c(  #fir annot even less broad
  
  "CBM" = "#46F0F0",
  "IMC_CBM" = "#46F0F0",
  "SMM" = "#4575B4",              # Macrophages
  "NDMM" = "#FF7F00", 
  "IMC_NDMM" = "#46F0F0",
  "PCL" = "#D73027",
  "Relapse" = "#984EA3" #New: yellow
)

patient_order <-c(
  "CBM_1", "CBM_2", "CBM_3", "CBM_4", "CBM_5",
  "SMM_1", "SMM_2", "SMM_3", "SMM_4",
  "NDMM_1", "NDMM_2", "NDMM_3", "NDMM_4", "NDMM_5",
  "PCL_1", "PCL_2", "PCL_3",
  "Relapse_1", "Relapse_2", "Relapse_3"
)



stromal_cluster_colors <- c(
  # General Stroma and MSCs (greens/teals, like before)
  "Mix"                      = "#999999",
  "CXCL14 MSC"                   = "#3CB44B",
  "Inflamed_structural_cell" = "#984EA3",# Bright green
  # Sea green
  "THY1 MSC"                 = "#66C2A5",  # Soft teal
  "APOD MSC"                 = "#8DD3C7",  # Pale turquoise
  "CXCL14 MSC"                 = "#A6D854",  # Light green
  "Fibro-like MSC"                = "#B3DE69",  # Lime green
  "LEPR MSC"                = "#2E8B57",  # Yellow
  "Adipocytes"             = "#FFD92F",
  "VSMC"            = "#D4A5A5",# Pale yellow
  # Endothelial (shades of red)
  "EC"                       = "#D73027",  # Vivid red
  "SEC"                      = "#F46D43",  # Orange-red
  "AEC"                      = "#FDAE61",  # Peach/light coral
  # Osteolineage (whitish to light brown)
  "Osteoblasts"        = "#D9D9D9",  # Light gray (bone tone)
  "Osteo-fibroblastic MSC"      = "#A65628"  # Slightly darker gray
  # Inflammatory/activated
  # Strong red
  # Purple
  # Misc
  # Mid gray
)

broad_cluster_colors <- c(
  # General Stroma and MSCs (greens/teals, like before)
  
  "Structural"                   = "#3CB44B",
  "Myeloid"                 = "#4575B4",  # Pale turquoise
  "Other"             = "black",
  "Lymphoid"                      = "#D73027")  # Orange-red


broad_cluster_colors_2 <- c(
  "MKC Lineage" = "black", "Mature Neutrophil" = "#4575B4", 
  "LTF Immature Neutrophil 3" = "#4575B4", "MMP9 Immature Neutrophil 4" = "#4575B4",
  "CEACAM6 Immature Neutrophil 1" = "#4575B4", "OLFM4 Immature Neutrophil 2" = "#4575B4",
  "Neutrophil Progenitor" = "#4575B4", "EC" = "#3CB44B",
  "AEC" = "#3CB44B", "SEC" = "#3CB44B",
  "Erythroid" = "black", "LEPR MSC" = "#3CB44B",
  "THY1 MSC" = "#3CB44B", "CXCL14 MSC" = "#3CB44B",
  "Fibro MSC" = "#3CB44B", "APOD MSC" = "#3CB44B",
  "Osteo-Fibroblastic MSC" = "#3CB44B", "VSMC" = "#3CB44B",
  "Conventional DC" = "#4575B4", "Plasmacytoid DC" = "#4575B4",
  "Macrophage" = "#4575B4", "Classical Monocyte" = "#4575B4",
  "Non Classical Monocyte" = "#4575B4", "CD8 T Cell" = "#D73027",
  "Naive or CM CD4 T Cell" = "#D73027", "NK or Cytotoxic T Cell" = "#D73027",
  "Activated or Exhausted T Cell" = "#D73027", "Regulatory T Cell" = "#D73027",
  "B or Plasma Cell" = "#D73027", "Plasma Cell" = "#D73027",
  "Cycling B or Plasma Cell" = "#D73027", "Interferon Stimulated Cell" = "black",
  "Interleukin Producing Cell" = "black", "Basophil Eosinophil" = "#4575B4",
  "Mast Cell" = "#4575B4", "Adipocyte" = "#3CB44B",
  "Osteoblast" = "#3CB44B"
)



niche_colors <- c(
  "Erythroid niche"              = "#9837db",  # pinkish
  "Erythroid Niche"              = "#9837db",  # pinkish
  "Megakaryocyte T cell niche"   = "#FFD92F",  # yellow
  "Megakaryocyte Niche"   = "#FFD92F",  # yellow
  "Neutrophil Progenitor niche"  = "#91C3E6",  # light blue
  "Neutrophil Progenitor Niche"  = "#91C3E6",  # light blue
  "Myeloid niche"                = "#4575B4",  # dark blue
  "Myeloid Niche"                = "#4575B4",  # dark blue
  "Mature Myeloid niche"         = "#46F0F0",  #
  "Immune-Endosteal Niche"         = "#46F0F0",  # cyan
  "Structural niche"             = "#3CB44B",  # green
  "Stromal Niche"             = "#3CB44B", 
  "Immune-Stromal Niche" = "#808000"  ,         # Darker brown  
  "Phagocytic niche"             = "#FF7F00",
  "T cell–Myeloid Niche"             = "#FF7F00",  # orange-red (plasma/tumor family)
  "Tumor-associated niche"       = "#D73027", 
  "Tumor-Associated Niche"       = "#D73027",  # medium red
  "Proliferative Tumor niche"    = "#8B0000",  # very dark red
  "Tumor niche"                  = "#FC4E2A",
  "Tumor Niche"                  = "#FC4E2A",
  "Erythroid"              = "#9837db",  # pinkish
  "Megakaryocyte T cell"   = "#FFD92F",  # yellow
  "Megakaryocyte"   = "#FFD92F",  # yellow
  "Neutrophil Progenitor"  = "#91C3E6",  # light blue
  "Myeloid"                = "#4575B4",  # dark blue
  "Mature Myeloid"         = "#46F0F0",  #
  "Immune-Endosteal"         = "#46F0F0",  # cyan
  "Structural"             = "#3CB44B",  # green
  "Stromal"             = "#3CB44B", 
  "Immune-Stromal" = "#808000"  ,         # Darker brown  
  "Phagocytic"             = "#FF7F00",
  "T cell–Myeloid"             = "#FF7F00",  # orange-red (plasma/tumor family)
  "Tumor-associated"       = "#D73027", 
  "Tumor-Associated"       = "#D73027",  # medium red
  "Proliferative Tumor"    = "#8B0000",  # very dark red
  "Tumor"                  = "#FC4E2A"
)

niche_colors_2 <- c(
  "Erythroid niche"              = "#dd1c77",  # pinkish
  "Erythroid Niche"              = "#dd1c77",  # pinkish
  "Megakaryocyte T cell niche"   = "#FFD92F",  # yellow
  "Megakaryocyte Niche"   = "#FFD92F",  # yellow
  "Neutrophil Progenitor niche"  = "#91C3E6",  # light blue
  "Neutrophil Progenitor Niche"  = "#91C3E6",  # light blue
  "Myeloid niche"                = "#4575B4",  # dark blue
  "Myeloid Niche"                = "#4575B4",  # dark blue
  "Mature Myeloid niche"         = "#46F0F0",  #
  "Immune-Endosteal Niche"         = "#46F0F0",  # cyan
  "Structural niche"             = "#3CB44B",  # green
  "Stromal Niche"             = "#3CB44B", 
  "Immune-Stromal Niche" = "#808000"  ,         # Darker brown  
  "Phagocytic niche"             = "#FF7F00",
  "T cell–Myeloid Niche"             = "#FF7F00",  # orange-red (plasma/tumor family)
  "Tumor-associated niche"       = "#9837db", 
  "Tumor-Associated Niche"       = "#9837db",  # medium red
  "Proliferative Tumor niche"    = "#8B0000",  # very dark red
  "Tumor niche"                  = "#FC4E2A",
  "Tumor Niche"                  = "#FC4E2A",
  "Erythroid"              = "#dd1c77",  # pinkish
  "Megakaryocyte T cell"   = "#FFD92F",  # yellow
  "Megakaryocyte"   = "#FFD92F",  # yellow
  "Neutrophil Progenitor"  = "#91C3E6",  # light blue
  "Myeloid"                = "#4575B4",  # dark blue
  "Mature Myeloid"         = "#46F0F0",  #
  "Immune-Endosteal"         = "#46F0F0",  # cyan
  "Structural"             = "#3CB44B",  # green
  "Stromal"             = "#3CB44B", 
  "Immune-Stromal" = "#808000"  ,         # Darker brown  
  "Phagocytic"             = "#FF7F00",
  "T cell–Myeloid"             = "#FF7F00",  # orange-red (plasma/tumor family)
  "Tumor-associated"       = "#9837db", 
  "Tumor-Associated"       = "#9837db",  # medium red
  "Proliferative Tumor"    = "#8B0000",  # very dark red
  "Tumor"                  = "#FC4E2A"
)

cluster_order_dense_vs_sparse_hetero <- c(
  "Osteo-Fibroblastic MSC", "MMP9 Immature Neutrophil 4", "Neutrophil Progenitor",
  "LTF Immature Neutrophil 3", "Basophil Eosinophil", "CEACAM6 Immature Neutrophil 1",
  "OLFM4 Immature Neutrophil 2", "Erythroid", "Classical Monocyte",
  "Non Classical Monocyte", "Interleukin Producing Cell", "Regulatory T Cell",
  "Osteoblast", "Mast Cell", "Mature Neutrophil",
  "MKC Lineage", "APOD MSC", "CXCL14 MSC",
  "LEPR MSC", "EC", "Interferon Stimulated Cell",
  "NK or Cytotoxic T Cell", "VSMC", "Naive or CM CD4 T Cell",
  "Adipocyte", "Plasmacytoid DC", "CD8 T Cell",
  "Fibro MSC", "AEC", "SEC", "Macrophage",
  "Conventional DC", "THY1 MSC", "Activated or Exhausted T Cell"
)

cluster_order_NDMM_dense_vs_CBM <- c(
  "MMP9 Immature Neutrophil 4", "Neutrophil Progenitor",
  "OLFM4 Immature Neutrophil 2", "CEACAM6 Immature Neutrophil 1", "LTF Immature Neutrophil 3", 
  "Non Classical Monocyte","Osteo-Fibroblastic MSC","Basophil Eosinophil",  "Interleukin Producing Cell", 
  "Erythroid", "Mature Neutrophil", "Naive or CM CD4 T Cell","Classical Monocyte",
  "Adipocyte","Osteoblast","Regulatory T Cell","Mast Cell", "NK or Cytotoxic T Cell","MKC Lineage" ,"APOD MSC",  "EC", "LEPR MSC","CXCL14 MSC","VSMC", "Interferon Stimulated Cell","Plasmacytoid DC", "CD8 T Cell",
  "Fibro MSC", "AEC",  
  "Activated or Exhausted T Cell","THY1 MSC","Macrophage","SEC","Conventional DC"
)

cluster_order_NDMM_vs_CBM <- c(
  "OLFM4 Immature Neutrophil 2","Non Classical Monocyte","MMP9 Immature Neutrophil 4", "Neutrophil Progenitor","Regulatory T Cell", 
  "CEACAM6 Immature Neutrophil 1","Basophil Eosinophil",  "Interleukin Producing Cell", 
  "MKC Lineage" , "LTF Immature Neutrophil 3", 
  "Adipocyte", "NK or Cytotoxic T Cell","Naive or CM CD4 T Cell",
  "EC","Osteoblast","Erythroid", "CXCL14 MSC",
  "Osteo-Fibroblastic MSC","Fibro MSC", "Mature Neutrophil",
  "Activated or Exhausted T Cell","Classical Monocyte","LEPR MSC", "Plasmacytoid DC",
  "CD8 T Cell","Interferon Stimulated Cell","VSMC","Mast Cell","APOD MSC", 
  "AEC",  "THY1 MSC","SEC","Macrophage","Conventional DC", "Cycling B or Plasma Cell", "B or Plasma Cell", "Plasma Cell"
)

cluster_order_NDMM_vs_CBM_ordered_per_broad_celltype <- c(
  "Interleukin Producing Cell", "MKC Lineage", "Erythroid", "Interferon Stimulated Cell",
  "Regulatory T Cell", "NK or Cytotoxic T Cell", "Naive or CM CD4 T Cell", "Activated or Exhausted T Cell",
  "CD8 T Cell", "Cycling B or Plasma Cell", "B or Plasma Cell", "Plasma Cell",
  "Adipocyte", "EC", "Osteoblast", "CXCL14 MSC",
  "Osteo-Fibroblastic MSC", "Fibro MSC", "LEPR MSC", "VSMC",
  "APOD MSC", "AEC", "THY1 MSC", "SEC",
  "OLFM4 Immature Neutrophil 2", "Non Classical Monocyte", "MMP9 Immature Neutrophil 4", "Neutrophil Progenitor",
  "CEACAM6 Immature Neutrophil 1", "Basophil Eosinophil", "LTF Immature Neutrophil 3", "Mature Neutrophil",
  "Classical Monocyte", "Plasmacytoid DC", "Mast Cell", "Macrophage",
  "Conventional DC"
)


cluster_ordered_for_s_curves <- c(
  
  "EC","AEC","SEC","Osteoblast","Osteo-Fibroblastic MSC",
  "APOD MSC","Fibro MSC","LEPR MSC",  
  "CXCL14 MSC","THY1 MSC","Adipocyte","VSMC", 
  "Regulatory T Cell","Naive or CM CD4 T Cell","NK or Cytotoxic T Cell", 
  "Activated or Exhausted T Cell","CD8 T Cell","Interleukin Producing Cell", "Undefined Cell",
  "Interferon Stimulated Cell","Erythroid","Erythroid Lineage","MKC Lineage","Basophil Eosinophil","Mast Cell","Neutrophil Progenitor","Early Immature Neutrophil","Late Immature Neutrophil","CEACAM6 Immature Neutrophil 1",
  "OLFM4 Immature Neutrophil 2", "LTF Immature Neutrophil 3", "MMP9 Immature Neutrophil 4","Immature Neutrophil",
  "Mature Neutrophil","Classical Monocyte", "Non Classical Monocyte",
  "Macrophage", "Conventional DC", "Plasmacytoid DC",
  "B or Plasma Cell","Cycling B or Plasma Cell","Plasma Cell"
)

cluster_ordered_for_s_curves_broad <- c(
  
  "Endothelium","Osteolineage","Stroma",
  "Adipocyte","VSMC", 
  "T/NK Cell","Other", "Erythroid","MKC Lineage","Basophil/Eosinophil/Mast Cell","Neutrophil Progenitor","Neutrophil",
  "Mature Neutrophil","Monocyte", "Macrophage", "Dendritic Cell", "B/Plasma Cell"
)

common_group_colors <- c(
  "Other"                 = "darkgray",       # from Fibrotic_Tissue / Undefined_Cell
  "Erythroid"             = "pink",           # from Ki67_likely_erythroid
  "MKC or Endothelial"    = "#90ee90",        # from MKC_Endothelial
  "Stroma"                = "#1C750C",         # from CD105_Structural_Cell
  "Neutrophil Progenitor" = "#00FFFF",        # from Neutrophil_Progenitor
  "Neutrophil"            = "#00BFFF",        # from Neutrophil
  "Macrophage"            = "#005BFF",        # from Macrophage
  "Monocyte"              = "#c51b8a",        # from CD14_Monocyte
  "APC"                  = "purple",    
  "CD8 T Cell"            = "#FA7921",        # from CD8_T_Cell
  "CD4 T Cell"            = "#FFE66D",        # from CD4_T_Cell
  "Plasma Cell"           = "#FF6666"       # from Plasma_Cell
)



region_colors <- c("CBM" = "cyan", "unqualified" = "#4CAF50", "normal PC percentage" =  "#FFE066", #, "normal PC percentage" =  "#FFFFB3"
                   "sparse" ="#377EB8","dispersed" ="#B080FF", "dense"= "#D72638") #, "dense_otheroption_toopink"= "#FB8072")  # Replace with your preferred order


status_order <- c("CBM", "SMM", "NDMM", "PCL","Relapse")

# Define your custom cluster order

ordered_sample_groups <- c(
  # CBM
  "CBM_1", "CBM_1B", "CBM_2", "CBM_2B", "CBM_3", "CBM_3B", "CBM_4", "CBM_4B", "CBM_5",
  # SMM
  "SMM_1", "SMM_2", "SMM_3", "SMM_4",
  # NDMM
  "NDMM_1", "NDMM_2", "NDMM_3", "NDMM_3B", "NDMM_4", "NDMM_4B", "NDMM_5", "NDMM_5B",
  # PCL
  "PCL_1", "PCL_2", "PCL_2B", "PCL_3",
  # Relapse
  "Relapse_1", "Relapse_2", "Relapse_3")



# Functions ---------------------------------------------------------------



prepare_metadata_frame <- function(meta, exclude_tumor = FALSE, cluster_col, region_col, tumor_col, levels_tumor, celltype_levels = NULL) {
  meta <- meta
  if (exclude_tumor) {
    #tumor_cells <- c("Plasma cells")
    tumor_cells <-c("Plasma cells", "B cells / Plasma cells", "Cycling B/Plasma cells","Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell", "B/Plasma Cell")
    meta <- meta[!meta[[cluster_col]] %in% tumor_cells, ]
  }
  
  #meta <- seurat_obj@meta.data
  df <- as.data.frame(table(meta[[cluster_col]], meta[[region_col]]))
  names(df) <- c("cluster_celltype", "small_region", "Freq")
  
  df$total_cells <- ave(df$Freq, df$small_region, FUN = sum)
  df$percentages <- df$Freq / df$total_cells * 100
  
  # Ensure character for safe matching
  region_match <- match(as.character(df$small_region), as.character(meta[[region_col]]))
  
  # Add annotations from metadata
  df$patient_id <- meta$Tissue_ext[region_match]
  df$ObjectName_anonymous <- factor(meta$ObjectName_anonymous[region_match], levels = unique(meta$ObjectName_anonymous))
  df$tumor_percentage_grouped <- factor(meta[[tumor_col]][region_match], levels = levels_tumor)
  df$Status <- factor(meta$Status[region_match], levels = unique(meta$Status))
  df$Status_simp <- factor(meta$Status_simp[region_match], levels = unique(meta$Status_simp))
  df$Status_simp_2 <- factor(meta$Status_simp_2[region_match], levels = unique(meta$Status_simp_2))
  df$ObjectName_anonymous_grouped <- factor(meta$ObjectName_anonymous_grouped[region_match], levels = unique(meta$ObjectName_anonymous_grouped))
  df$quality_group_nFeature <- factor(meta$quality_group_nFeature[region_match], levels = unique(meta$quality_group_nFeature))
  #df$bone_zone <- factor(meta$bone_zone[region_match], levels = c("Endosteal", "Peri-endosteal", "Intermediate Marrow", "Central Marrow"))
  
  if (!is.null(celltype_levels)) {
    df$cluster_celltype <- factor(df$cluster_celltype, levels = celltype_levels)
  }
  
  return(df)
}



aggregate_percentages_dynamic <- function(df, group_col) {
  group_col_sym <- rlang::sym(group_col)
  
  df <- df %>%
    group_by(ObjectName_anonymous_grouped, !!group_col_sym) %>% # change this if you want to group tissues
    mutate(Freq_per_patient_per_group = sum(Freq)) %>%
    group_by(ObjectName_anonymous_grouped, cluster_celltype, !!group_col_sym) %>%
    mutate(Freq_per_celltype_per_patient_per_group = sum(Freq)) %>%
    ungroup()
  
  # Create a new percentage column name dynamically
  percent_col_name <- paste0("percentages_per_patient_per_", group_col)
  
  df[[percent_col_name]] <- df$Freq_per_celltype_per_patient_per_group / df$Freq_per_patient_per_group * 100
  
  return(df)
}


plot_patient_deltas_boxplot_with_stat_same_order_capped_test_archi_lmer <- function(
    df, group1, group2, method = c("delta", "ratio"),
    cluster_order_defined,
    cluster_colors = NULL,
    tumor = FALSE,
    adjust_p = TRUE,
    save_stats_path = "~/R_analyses/scRNAseq/Xenium_run3/p_values_LMER_test2.xlsx",
    comparison_name = NULL
) {
  
  method <- match.arg(method)
  eps <- 1e-5
  
  group2_means <- df %>%
    filter(tumor_percentage_grouped035ext_hex_average_method2 == group2) %>%
    group_by(clusters) %>%
    summarise(mean_group2 = mean(Freq), .groups = "drop")
  
  result_df <- df %>%
    filter(tumor_percentage_grouped035ext_hex_average_method2 == group1) %>%
    inner_join(group2_means, by = "clusters") %>%
    mutate(
      delta = Freq - mean_group2,
      ratio = log2((Freq + eps) / (mean_group2 + eps))
    ) %>%
    mutate(value = if (method == "delta") delta else ratio)
  
  result_df$clusters <- factor(result_df$clusters, levels = cluster_order_defined)
  
  df <- df %>% mutate(patient_id = sample)
  
  df_agg <- df %>%
    group_by(patient_id, clusters,ObjectName_anonymous_grouped, tumor_percentage_grouped035ext_hex_average_method2) %>%
    summarise(Freq = mean(Freq), .groups = "drop")
  
  #stat_df <- df_agg %>%
  # filter(tumor_percentage_grouped035ext_hex_average_method2 %in% c(group1, group2)) %>%
  #group_by(clusters) %>%
  #summarise(
  # p_val = tryCatch({
  #  t.test(
  #   Freq[tumor_percentage_grouped035ext_hex_average_method2 == group1],
  #  Freq[tumor_percentage_grouped035ext_hex_average_method2 == group2],
  # var.equal = TRUE,     # <-- Student's t-test (assumes equal variance)
  #paired = FALSE
  #)$p.value
  #}, error = function(e) NA),
  #.groups = "drop"
  # ) %>%
  # mutate(
  #  p_used = if (adjust_p) p.adjust(p_val, method = "BH") else p_val,
  #  stars = case_when(
  #   is.na(p_used) ~ "",
  #   p_used < 0.001 ~ "***",
  #  p_used < 0.01  ~ "**",
  ##   p_used < 0.05  ~ "*",
  #   TRUE ~ ""
  #  )
  #  )
  
  # Get all clusters
  clusters <- unique(df_agg$clusters)
  
  stat_df <- purrr::map_dfr(clusters, function(cluster_of_interest) {
    
    tmp <- df_agg %>%
      filter(
        clusters == cluster_of_interest,
        tumor_percentage_grouped035ext_hex_average_method2 %in% c(group1, group2)
      )%>%
      mutate(
        tumor_percentage_grouped035ext_hex_average_method2 =
          factor(tumor_percentage_grouped035ext_hex_average_method2,
                 levels = c(group2, group1))   # <-- IMPORTANT
      )
    
    
    p_val <- tryCatch({
      model <- lmerTest::lmer(
        Freq ~ tumor_percentage_grouped035ext_hex_average_method2 +
          (1 | ObjectName_anonymous_grouped),
        data = tmp,
        REML = FALSE
      )
      coef_name <- paste0(
        "tumor_percentage_grouped035ext_hex_average_method2",
        group1
      )
      
      summary(model)$coefficients[coef_name, "Pr(>|t|)"]
      #p_val <- summary(model)$coefficients[2, "Pr(>|t|)"]
    }, error = function(e) NA)
    
    tibble(
      clusters = cluster_of_interest,
      p_val = p_val
    )
  })
  
  # Add multiple testing correction and dynamic naming
  stat_df <- stat_df %>%
    mutate(
      p_used = p.adjust(p_val, method = "BH"),
      stars = case_when(
        is.na(p_used) ~ "",
        p_used < 0.001 ~ "***",
        p_used < 0.01  ~ "**",
        p_used < 0.05  ~ "*",
        TRUE ~ ""
      ),
      comparison = paste(group1, "vs", group2)
    )
  
  result_df <- result_df %>%
    left_join(stat_df, by = "clusters") %>%
    mutate(value_capped = pmin(pmax(value, -4.9), 4.9))
  
  comparison_name <- paste(group1, "vs", group2)
  
  safe_sheet <- gsub("[^A-Za-z0-9_]", "_", comparison_name)
  safe_sheet <- substr(safe_sheet, 1, 31)
  
  # Save to Excel
  if (!is.null(save_stats_path)) {
    if (!file.exists(save_stats_path)) {
      wb <- createWorkbook()
    } else {
      wb <- loadWorkbook(save_stats_path)
    }
    
    addWorksheet(wb, safe_sheet)
    writeData(wb, safe_sheet, result_df)
    saveWorkbook(wb, save_stats_path, overwrite = TRUE)
  }
  
  
  star_x_pos <- result_df %>%
    group_by(clusters) %>%
    summarise(
      median_val = median(value, na.rm = TRUE),
      xpos = ifelse(median_val >= 0, max(value_capped, na.rm = TRUE), min(value_capped, na.rm = TRUE)),
      label = first(stars),
      .groups = "drop"
    )
  
  if (!is.null(cluster_colors)) {
    axis_text_color <- function(x) cluster_colors[as.character(x)]
  } else {
    axis_text_color <- function(x) rep("black", length(x))
  }
  
  result_df$clusters <- droplevels(result_df$clusters)
  cluster_levels <- levels(result_df$clusters)
  
  cluster_groups1 <- tibble::tribble(
    ~group,        ~start,        ~end,
    "Myeloid",     "Plasmacytoid DC",      "Basophil Eosinophil",
    "Other",     "MKC Lineage",        "Undefined Cell",
    "T/NK Cell",       "CD8 T Cell",     "Regulatory T Cell")
  
  cluster_groups2 <- tibble::tribble(
    ~group,        ~start,        ~end,
    "B Lin",     "Plasma Cell",      "B or Plasma Cell",
    "Myeloid",     "Plasmacytoid DC",      "Basophil Eosinophil",
    "Other",     "MKC Lineage",        "Undefined Cell",
    "T/NK Cell",       "CD8 T Cell",     "Regulatory T Cell")
  
  # Select based on tumor flag 
  cluster_groups <- if (isTRUE(tumor)) cluster_groups2 else cluster_groups1
  
  
  bg_rects <- cluster_groups %>%
    mutate(
      ymin = pmin(match(start, cluster_levels),
                  match(end,   cluster_levels)) - 0.5,
      ymax = pmax(match(start, cluster_levels),
                  match(end,   cluster_levels)) + 0.5
    )
  
  hline_df <- bg_rects %>%
    transmute(
      yintercept = ymin - 0
    )
  
  
  ggplot(result_df, aes(x = value_capped, y = clusters)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(
      data = hline_df,
      aes(yintercept = yintercept),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.6
    ) +
    geom_boxplot(aes(group = clusters), fill = "gray90", color = "black", 
                 width = 0.6, outlier.shape = NA, size = 0.9) +
    geom_jitter(aes(color = sample), height = 0.2, size = 2, alpha = 0.8) +
    geom_text(data = star_x_pos,
              aes(x = xpos * 1.05, y = clusters, label = label),
              inherit.aes = FALSE,
              hjust = ifelse(star_x_pos$xpos > 0, 0, 1),
              size = 5) +
    coord_cartesian(xlim = c(-5, 5)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(
        size = 12,
        color = axis_text_color(levels(result_df$clusters))  # <-- NEW
      ),
      axis.text.x = element_text(size = 12)
    ) +
    labs(
      x = ifelse(method == "delta",
                 paste("Delta from mean of", group2),
                 paste("log2 ratio vs mean of", group2)),
      y = "Cluster",
      title = paste("Per-patient cluster", method, ":", group1, "vs", group2)
    )
}




plot_patient_deltas_boxplot_with_stat_same_order_capped_normal_T_test <- function(
    df, group1, group2, method = c("delta", "ratio"),
    cluster_order_defined,
    cluster_colors = NULL,
    tumor = FALSE,
    adjust_p = TRUE,
    save_stats_path = NULL,
    comparison_name = NULL
) {
  
  method <- match.arg(method)
  eps <- 1e-5
  
  group2_means <- df %>%
    filter(Status_simp_2 == group2) %>%
    group_by(clusters) %>%
    summarise(mean_group2 = mean(Freq), .groups = "drop")
  
  result_df <- df %>%
    filter(Status_simp_2 == group1) %>%
    inner_join(group2_means, by = "clusters") %>%
    mutate(
      delta = Freq - mean_group2,
      ratio = log2((Freq + eps) / (mean_group2 + eps))
    ) %>%
    mutate(value = if (method == "delta") delta else ratio)
  
  result_df$clusters <- factor(result_df$clusters, levels = cluster_order_defined)
  
  df <- df %>% mutate(patient_id = sample)
  
  df_agg <- df %>%
    group_by(patient_id, clusters, Status_simp_2) %>%
    summarise(Freq = mean(Freq), .groups = "drop")
  
  stat_df <- df_agg %>%
    filter(Status_simp_2 %in% c(group1, group2)) %>%
    group_by(clusters) %>%
    summarise(
      p_val = tryCatch({
        t.test(
          Freq[Status_simp_2 == group1],
          Freq[Status_simp_2 == group2],
          var.equal = TRUE
        )$p.value
      }, error = function(e) NA),
      mean_diff = mean(Freq[Status_simp_2 == group1], na.rm = TRUE) -
        mean(Freq[Status_simp_2 == group2], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      p_used = if (adjust_p) p.adjust(p_val, method = "BH") else p_val,
      stars = case_when(
        is.na(p_used) ~ "",
        p_used < 0.001 ~ "***",
        p_used < 0.01  ~ "**",
        p_used < 0.05  ~ "*",
        TRUE ~ ""
      ),
      comparison = ifelse(is.null(comparison_name), paste(group1, "vs", group2), comparison_name)
    )
  
  
  comparison_name <- paste(group1, "vs", group2)
  # Make it safe for Excel sheet name (max 31 chars, no illegal characters)
  safe_sheet <- gsub("[^A-Za-z0-9_]", "_", comparison_name)
  safe_sheet <- substr(safe_sheet, 1, 31)
  
  # Save to Excel
  if (!is.null(save_stats_path)) {
    if (!file.exists(save_stats_path)) {
      wb <- createWorkbook()
    } else {
      wb <- loadWorkbook(save_stats_path)
    }
    
    addWorksheet(wb, safe_sheet)
    writeData(wb, safe_sheet, stat_df)
    saveWorkbook(wb, save_stats_path, overwrite = TRUE)
  }
  
  
  result_df <- result_df %>%
    left_join(stat_df, by = "clusters") %>%
    mutate(value_capped = pmin(pmax(value, -4.9), 4.9))
  
  star_x_pos <- result_df %>%
    group_by(clusters) %>%
    summarise(
      median_val = median(value, na.rm = TRUE),
      xpos = ifelse(median_val >= 0, max(value_capped, na.rm = TRUE), min(value_capped, na.rm = TRUE)),
      label = unique(stars),
      .groups = "drop"
    )
  
  if (!is.null(cluster_colors)) {
    axis_text_color <- function(x) cluster_colors[as.character(x)]
  } else {
    axis_text_color <- function(x) rep("black", length(x))
  }
  
  result_df$clusters <- droplevels(result_df$clusters)
  cluster_levels <- levels(result_df$clusters)
  
  cluster_groups1 <- tibble::tribble(
    ~group,        ~start,        ~end,
    "Myeloid",     "Plasmacytoid DC",      "Basophil Eosinophil",
    "Other",     "MKC Lineage",        "Undefined Cell",
    "T/NK Cell",       "CD8 T Cell",     "Regulatory T Cell")
  
  cluster_groups2 <- tibble::tribble(
    ~group,        ~start,        ~end,
    "B Lin",     "Plasma Cell",      "B or Plasma Cell",
    "Myeloid",     "Plasmacytoid DC",      "Basophil Eosinophil",
    "Other",     "MKC Lineage",        "Undefined Cell",
    "T/NK Cell",       "CD8 T Cell",     "Regulatory T Cell")
  
  # Select based on tumor flag 
  cluster_groups <- if (isTRUE(tumor)) cluster_groups2 else cluster_groups1
  
  
  bg_rects <- cluster_groups %>%
    mutate(
      ymin = pmin(match(start, cluster_levels),
                  match(end,   cluster_levels)) - 0.5,
      ymax = pmax(match(start, cluster_levels),
                  match(end,   cluster_levels)) + 0.5
    )
  
  hline_df <- bg_rects %>%
    transmute(
      yintercept = ymin - 0
    )
  
  
  ggplot(result_df, aes(x = value_capped, y = clusters)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(
      data = hline_df,
      aes(yintercept = yintercept),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.6
    ) +
    geom_boxplot(aes(group = clusters), fill = "gray90", color = "black", 
                 width = 0.6, outlier.shape = NA, size = 0.9) +
    geom_jitter(
      aes(x = value_capped, color = sample),
      height = 0.2,
      size = 2,
      alpha = 0.8
    )+
    geom_text(data = star_x_pos,
              aes(x = xpos * 1.05, y = clusters, label = label),
              inherit.aes = FALSE,
              hjust = ifelse(star_x_pos$xpos > 0, 0, 1),
              size = 5) +
    coord_cartesian(xlim = c(-5, 5)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(
        size = 12,
        color = axis_text_color(levels(result_df$clusters))  # <-- NEW
      ),
      axis.text.x = element_text(size = 12)
    ) +
    labs(
      x = ifelse(method == "delta",
                 paste("Delta from mean of", group2),
                 paste("log2 ratio vs mean of", group2)),
      y = "Cluster",
      title = paste("Per-patient cluster", method, ":", group1, "vs", group2)
    )
}





plot_patient_deltas_boxplot_with_stat_same_order_capped_paired <- function(
    df,  cluster_order_defined,
    cluster_colors = NULL,method = c("delta", "ratio"),
    tumor = FALSE,
    adjust_p = TRUE 
) {
  
  eps <- 1e-5
  
  df <- df %>%
    mutate(
      patient_id = sub("_(dense|sparse)$", "", sample)
    )
  
  df_agg <- df %>%
    group_by(patient_id, clusters, arch_grouped) %>%
    summarise(Freq = mean(Freq), .groups = "drop")
  
  # Pivot to wide format so dense and sparse are paired per patient
  paired_df <- df_agg %>%
    pivot_wider(
      id_cols = c(patient_id , clusters),
      names_from = arch_grouped,
      values_from = Freq
    )
  
  result_df <- paired_df %>%
    mutate(
      delta = dense - sparse   ,
      ratio = log2((dense + eps) / (sparse + eps))
    ) %>%
    mutate(value = if (method == "delta") delta else ratio)
  
  # Then compute paired t-test per cluster
  stat_df <- paired_df %>%
    group_by(clusters) %>%
    summarise(
      p_val = tryCatch({
        t.test(dense, sparse, paired = TRUE)$p.value
      }, error = function(e) NA),
      .groups = "drop"
    ) %>%
    mutate(
      p_used = if (adjust_p) p.adjust(p_val, method = "BH") else p_val,
      stars = case_when(
        is.na(p_used) ~ "",
        p_used < 0.001 ~ "***",
        p_used < 0.01  ~ "**",
        p_used < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  # Prepare plotting df
  
  
  result_df <- result_df %>%
    left_join(stat_df, by = "clusters") %>%
    mutate(value_capped = pmin(pmax(value, -4.9), 4.9))
  
  star_x_pos <- result_df %>%
    group_by(clusters) %>%
    summarise(
      median_val = median(value, na.rm = TRUE),
      xpos = ifelse(median_val >= 0, max(value_capped, na.rm = TRUE), min(value_capped, na.rm = TRUE)),
      label = unique(stars),
      .groups = "drop"
    )
  
  if (!is.null(cluster_colors)) {
    axis_text_color <- function(x) cluster_colors[as.character(x)]
  } else {
    axis_text_color <- function(x) rep("black", length(x))
  }
  
  result_df$clusters <- droplevels(result_df$clusters)
  cluster_levels <- levels(result_df$clusters)
  
  cluster_groups1 <- tibble::tribble(
    ~group,        ~start,        ~end,
    "Myeloid",     "Plasmacytoid DC",      "Basophil Eosinophil",
    "Other",     "MKC Lineage",        "Undefined Cell",
    "T/NK Cell",       "CD8 T Cell",     "Regulatory T Cell")
  
  cluster_groups2 <- tibble::tribble(
    ~group,        ~start,        ~end,
    "B Lin",     "Plasma Cell",      "B or Plasma Cell",
    "Myeloid",     "Plasmacytoid DC",      "Basophil Eosinophil",
    "Other",     "MKC Lineage",        "Undefined Cell",
    "T/NK Cell",       "CD8 T Cell",     "Regulatory T Cell")
  
  # Select based on tumor flag 
  cluster_groups <- if (isTRUE(tumor)) cluster_groups2 else cluster_groups1
  
  
  bg_rects <- cluster_groups %>%
    mutate(
      ymin = pmin(match(start, cluster_levels),
                  match(end,   cluster_levels)) - 0.5,
      ymax = pmax(match(start, cluster_levels),
                  match(end,   cluster_levels)) + 0.5
    )
  
  hline_df <- bg_rects %>%
    transmute(
      yintercept = ymin - 0
    )
  
  
  df_plot <- result_df
  df_plot$clusters <- factor(df_plot$clusters, levels = cluster_order_defined)
  
  ggplot(df_plot, aes(x = value_capped, y = clusters)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(
      data = hline_df,
      aes(yintercept = yintercept),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.6
    ) +
    geom_boxplot(aes(group = clusters), fill = "gray90", color = "black", 
                 width = 0.6, outlier.shape = NA, size = 0.9) +
    geom_jitter(aes(color = patient_id), height = 0.2, size = 2, alpha = 0.8) +
    geom_text(data = star_x_pos,
              aes(x = xpos * 1.05, y = clusters, label = label),
              inherit.aes = FALSE,
              hjust = ifelse(star_x_pos$xpos > 0, 0, 1),
              size = 5) +
    coord_cartesian(xlim = c(-5, 5)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(
        size = 12
      ),
      axis.text.x = element_text(size = 12)
    ) +
    labs(
      x = ifelse(method == "delta",
                 paste("Delta"),
                 paste("log2 ratio vs mean of")),
      y = "Cluster",
      title = paste("Per-patient cluster", method)
    )
}



plot_patient_deltas_boxplot_with_stat_same_order_capped_paired_broad <- function(
    df,  cluster_order_defined,
    cluster_colors = NULL,method = c("delta", "ratio"),
    tumor = FALSE,
    adjust_p = TRUE 
) {
  
  df <- df %>%
    mutate(
      patient_id = sub("_(dense|sparse)$", "", sample)
    )
  
  df_agg <- df %>%
    group_by(patient_id, clusters, arch_grouped) %>%
    summarise(Freq = mean(Freq), .groups = "drop")
  
  # Pivot to wide format so dense and sparse are paired per patient
  paired_df <- df_agg %>%
    pivot_wider(
      id_cols = c(patient_id , clusters),
      names_from = arch_grouped,
      values_from = Freq
    )
  
  result_df <- paired_df %>%
    mutate(
      delta = dense - sparse   ,
      ratio = log2((dense + eps) / (sparse + eps))
    ) %>%
    mutate(value = if (method == "delta") delta else ratio)
  result_df <- paired_df %>%
    mutate(
      delta = dense - sparse   ,
      ratio = log2((dense + eps) / (sparse + eps))
    ) %>%
    mutate(value = if (method == "delta") delta else ratio)
  
  # Then compute paired t-test per cluster
  stat_df <- paired_df %>%
    group_by(clusters) %>%
    summarise(
      p_val = tryCatch({
        t.test(dense, sparse, paired = TRUE)$p.value
      }, error = function(e) NA),
      .groups = "drop"
    ) %>%
    mutate(
      p_used = if (adjust_p) p.adjust(p_val, method = "BH") else p_val,
      stars = case_when(
        is.na(p_used) ~ "",
        p_used < 0.001 ~ "***",
        p_used < 0.01  ~ "**",
        p_used < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  # Prepare plotting df
  
  
  result_df <- result_df %>%
    left_join(stat_df, by = "clusters") %>%
    mutate(value_capped = pmin(pmax(value, -2.4), 2.4))
  
  star_x_pos <- result_df %>%
    group_by(clusters) %>%
    summarise(
      median_val = median(value, na.rm = TRUE),
      xpos = ifelse(median_val >= 0, max(value_capped, na.rm = TRUE), min(value_capped, na.rm = TRUE)),
      label = unique(stars),
      .groups = "drop"
    )
  
  if (!is.null(cluster_colors)) {
    axis_text_color <- function(x) cluster_colors[as.character(x)]
  } else {
    axis_text_color <- function(x) rep("black", length(x))
  }
  
  result_df$clusters <- droplevels(result_df$clusters)
  cluster_levels <- levels(result_df$clusters)
  
  cluster_groups1 <- tibble::tribble(
    ~group,        ~start,        ~end,
    "Myeloid",     "Plasmacytoid DC",      "Basophil Eosinophil",
    "Other",     "MKC Lineage",        "Undefined Cell",
    "T/NK Cell",       "CD8 T Cell",     "Regulatory T Cell")
  
  cluster_groups2 <- tibble::tribble(
    ~group,        ~start,        ~end,
    "B Lin",     "Plasma Cell",      "B or Plasma Cell",
    "Myeloid",     "Plasmacytoid DC",      "Basophil Eosinophil",
    "Other",     "MKC Lineage",        "Undefined Cell",
    "T/NK Cell",       "CD8 T Cell",     "Regulatory T Cell")
  
  # Select based on tumor flag 
  cluster_groups <- if (isTRUE(tumor)) cluster_groups2 else cluster_groups1
  
  
  bg_rects <- cluster_groups %>%
    mutate(
      ymin = pmin(match(start, cluster_levels),
                  match(end,   cluster_levels)) - 0.5,
      ymax = pmax(match(start, cluster_levels),
                  match(end,   cluster_levels)) + 0.5
    )
  
  hline_df <- bg_rects %>%
    transmute(
      yintercept = ymin - 0
    )
  
  
  df_plot <- result_df
  df_plot$clusters <- factor(df_plot$clusters, levels = cluster_order_defined)
  
  df_plot <- df_plot %>% mutate(patient_id = patient_id )
  
  
  ggplot(df_plot, aes(x = value_capped, y = clusters)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(
      data = hline_df,
      aes(yintercept = yintercept),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.6
    ) +
    geom_boxplot(aes(group = clusters), fill = "gray90", color = "black", 
                 width = 0.6, outlier.shape = NA, size = 0.9) +
    geom_jitter(aes(color = patient_id), height = 0.2, size = 2, alpha = 0.8) +
    geom_text(data = star_x_pos,
              aes(x = xpos * 1.05, y = clusters, label = label),
              inherit.aes = FALSE,
              hjust = ifelse(star_x_pos$xpos > 0, 0, 1),
              size = 5) +
    coord_cartesian(xlim = c(-2.5, 2.5)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(
        size = 12,
        color = axis_text_color(levels(df_plot$clusters))  # <-- NEW
      ),
      axis.text.x = element_text(size = 12)
    ) +
    labs(
      x = ifelse(method == "delta",
                 paste("Delta"),
                 paste("log2 ratio vs mean of")),
      y = "Cluster",
      title = paste("Per-patient cluster", method)
    )
}




transfer_metadata_columns_flex <- function(from_obj, to_obj, columns_to_transfer, link_from = NULL, link_to = NULL) {
  # Extract metadata from Seurat, SpatialExperiment, or assume data.frame
  if (inherits(from_obj, "Seurat")) {
    from_meta <- from_obj@meta.data
    from_keys <- if (!is.null(link_from)) from_meta[[link_from]] else rownames(from_meta)
  } else if (inherits(from_obj, "SpatialExperiment")) {
    from_meta <- as.data.frame(colData(from_obj))
    from_keys <- if (!is.null(link_from)) from_meta[[link_from]] else rownames(from_meta)
  } else {
    from_meta <- from_obj
    from_keys <- if (!is.null(link_from)) from_meta[[link_from]] else rownames(from_meta)
  }
  
  # Handle to_obj class
  if (inherits(to_obj, "Seurat")) {
    to_meta <- to_obj@meta.data
    target_class <- "Seurat"
  } else if (inherits(to_obj, "SpatialExperiment")) {
    to_meta <- as.data.frame(colData(to_obj))
    target_class <- "SPE"
  } else {
    to_meta <- to_obj
    target_class <- "data.frame"
  }
  
  to_keys <- if (!is.null(link_to)) to_meta[[link_to]] else rownames(to_meta)
  original_rownames <- rownames(to_meta)
  
  # Prepare metadata to join
  metadata_df <- from_meta[, columns_to_transfer, drop = FALSE]
  metadata_df$.__link__. <- from_keys
  
  # Remove overlapping columns from target
  cols_to_remove <- intersect(columns_to_transfer, colnames(to_meta))
  if (length(cols_to_remove)) {
    to_meta <- to_meta[, !(colnames(to_meta) %in% cols_to_remove), drop = FALSE]
  }
  to_meta$.__link__. <- to_keys
  
  # Join
  joined_meta <- dplyr::left_join(to_meta, metadata_df, by = ".__link__.")
  joined_meta$.__link__. <- NULL
  rownames(joined_meta) <- original_rownames
  
  # Return in original class
  if (target_class == "Seurat") {
    to_obj@meta.data <- joined_meta
    return(to_obj)
  } else if (target_class == "SPE") {
    colData(to_obj) <- S4Vectors::DataFrame(joined_meta)
    return(to_obj)
  } else {
    return(joined_meta)
  }
}


# Function
make_volcano_plot <- function(
    data,
    log2fc_col   = "log2FoldChange",
    padj_col     = "padj",
    gene_col     = NULL,                 # if NULL, uses rownames
    up_label     = "CTC-high",
    down_label   = "CTC-low",
    other_label  = "Not significant",
    padj_cut     = 0.05,
    fc_cut       = 1.5,
    color_up     = "#6A0DAD",
    color_down   = "#996699",
    color_other  = "#D9D9D9",
    title        = NULL,
    subtitle     = NULL,
    label_all_sig = TRUE,                # label all significant DEGs
    label_top_n   = NULL,                # or label top N by padj
    label_genes   = NULL,                # or label only this vector
    label_size    = 4,
    max_overlaps  = 100,
    seed          = 123,
    legend_position="right",
    # Y-axis handling
    y_upper_mode  = c("max", "quantile"),# choose "max" (default) or "quantile"
    y_quantile    = 0.999,               # used only if y_upper_mode = "quantile"
    y_extra       = 0.5,                 # extra headroom added to upper limit
    expand_y      = c(0, 0.02),          # extra expansion for y (bottom, top)
    # Padj clamping
    clamp_floor    = NULL                # NULL: half of min positive padj; else fixed (e.g. 1e-6)
) {
  # --- Dependencies ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required.")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Package 'rlang' is required.")
  
  y_upper_mode <- match.arg(y_upper_mode)
  
  # --- Extract columns safely ---
  if (!all(c(log2fc_col, padj_col) %in% colnames(data))) {
    stop("Input data must contain columns: ", paste(c(log2fc_col, padj_col), collapse = ", "))
  }
  df <- data[!is.na(data[[log2fc_col]]) & !is.na(data[[padj_col]]), , drop = FALSE]
  
  # Gene names
  if (is.null(gene_col)) {
    if (!is.null(rownames(df))) {
      df$gene <- rownames(df)
    } else {
      stop("No gene_col provided and rownames are NULL. Please specify 'gene_col'.")
    }
  } else {
    if (!gene_col %in% colnames(df)) stop("gene_col '", gene_col, "' not found in data.")
    df$gene <- as.character(df[[gene_col]])
  }
  
  # --- Padj clamping: avoid -log10(0) inflation ---
  padj_vec <- df[[padj_col]]
  if (is.null(clamp_floor)) {
    min_pos <- suppressWarnings(min(padj_vec[padj_vec > 0], na.rm = TRUE))
    floor_padj <- if (is.finite(min_pos)) min_pos / 2 else 1e-6
  } else {
    floor_padj <- clamp_floor
  }
  padj_plot <- padj_vec
  padj_plot[is.na(padj_plot)] <- 1   # NA treated as non-significant
  padj_plot[padj_plot <= 0]   <- floor_padj
  
  df$padj_clamped <- padj_plot
  df$neglog10padj <- -log10(df$padj_clamped)
  lfc <- df[[log2fc_col]]
  
  # --- Classification (Up/Down/Other) ---
  is_sig  <- (df$padj_clamped <= padj_cut) & (abs(lfc) >= fc_cut)
  is_up   <- is_sig & (lfc >= fc_cut)
  is_down <- is_sig & (lfc <= -fc_cut)
  
  df$category <- ifelse(is_up, up_label,
                        ifelse(is_down, down_label, other_label))
  df$category <- factor(df$category, levels = c(up_label, down_label, other_label))
  
  # --- Compute generous Y upper bound ---
  if (y_upper_mode == "max") {
    y_upper <- max(df$neglog10padj, na.rm = TRUE) + y_extra
  } else {
    y_upper <- max(stats::quantile(df$neglog10padj, y_quantile, na.rm = TRUE),
                   -log10(padj_cut)) + y_extra
  }
  y_lower <- 0
  
  # --- Base plot (tidy-eval, no aes_string) ---
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = !!rlang::sym(log2fc_col), y = neglog10padj, color = category)
  ) +
    ggplot2::geom_point(size = 2, alpha = 1) +
    ggplot2::geom_vline(xintercept = c(-fc_cut, fc_cut),
                        linetype = "dashed", color = "grey50", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut),
                        linetype = "dashed", color = "grey50", linewidth = 0.4) +
    ggplot2::scale_color_manual(
      values = setNames(c(color_up, color_down, color_other),
                        c(up_label, down_label, other_label))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = expression(Log[2]~"fold change"),
      y = expression(-Log[10]~"(padj)"),
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.position = legend_position
    ) +
    ggplot2::coord_cartesian(ylim = c(y_lower, y_upper)) +
    ggplot2::scale_y_continuous(expand = expand_y)
  
  # --- Labels selection ---
  df_sig <- df[is_sig, , drop = FALSE]
  
  if (!is.null(label_genes)) {
    df_lab <- df_sig[df_sig$gene %in% label_genes, , drop = FALSE]
  } else if (isTRUE(label_all_sig)) {
    df_lab <- df_sig
  } else {
    # Top 20 per direction (up/down)
    df_lab_up <- df_sig[df_sig$category == up_label, , drop = FALSE]
    df_lab_down <- df_sig[df_sig$category == down_label, , drop = FALSE]
    
    df_lab_up <- df_lab_up[order(df_lab_up$padj_clamped, decreasing = FALSE), , drop = FALSE]
    df_lab_down <- df_lab_down[order(df_lab_down$padj_clamped, decreasing = FALSE), , drop = FALSE]
    
    df_lab <- rbind(
      head(df_lab_up, label_top_n),
      head(df_lab_down, label_top_n)
    )
  }
  
  #--- Exclude genes whose names start with "ENSG0" from labeling 
  # Robust to leading/trailing spaces and case (e.g., "ensg0000...")
  df_lab <- df_lab[!grepl("^\\s*ENSG0", df_lab$gene, ignore.case = TRUE), , drop = FALSE]
  
  # Split labels by up/down to fix colors and keep legend clean
  df_lab_up   <- df_lab[df_lab$category == up_label, , drop = FALSE]
  df_lab_down <- df_lab[df_lab$category == down_label, , drop = FALSE]
  
  if (nrow(df_lab_up) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = df_lab_up,
      ggplot2::aes(x = !!rlang::sym(log2fc_col), y = neglog10padj, label = gene),
      fontface = "italic",
      size = label_size,
      max.overlaps = max_overlaps,
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.color = "grey40",
      segment.size = 0.3,
      color = color_up,          # fixed; no legend entry
      show.legend = FALSE,
      seed = seed
    )
  }
  if (nrow(df_lab_down) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = df_lab_down,
      ggplot2::aes(x = !!rlang::sym(log2fc_col), y = neglog10padj, label = gene),
      fontface = "italic",
      size = label_size,
      max.overlaps = max_overlaps,
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.color = "grey40",
      segment.size = 0.3,
      color = color_down,        # fixed; no legend entry
      show.legend = FALSE,
      seed = seed
    )
  }
  
  return(p)
}






# Figure 1A umap ---------------------------------------------------------------

Idents(Xenium_Object_Seurat  ) <- Xenium_Object_Seurat@meta.data $annot_broad_run3_correct_ifndiv

n_cells <- ncol(Xenium_Object_Seurat)
plot <- DimPlot(Xenium_Object_Seurat, label = T, repel = T, raster = T, cols = cluster_colors)+
  ggtitle(paste0("ALl (n = ", format(n_cells, big.mark = ","), ")"))

plot

pdf(file= file.path(output_dir,"UMAP_broad.pdf"), 
    width = 15, height = 15, useDingbats = FALSE, onefile = FALSE)

print(plot)

dev.off()

# Figure 1B/D UMAP dot plot structural ---------------------------------------------------------------

structural_cells <- qread("~/structural_cells.qs" )

# Assign Idents first
Idents(structural_cells) <- structural_cells$annot_structure_second
# Rename Idents with singular names and consistent formatting
structural_cells <- RenameIdents(object = structural_cells,
                                                   "Osteoblasts" = "Osteoblast",
                                                   "Adipocytes" = "Adipocyte",
                                                   "Fibro-like MSC" = "Fibro MSC",
                                                   "Osteo-fibroblastic MSC" = "Osteo-Fibroblastic MSC"
)
structural_cells$annot_structure_second_correct  <- Idents(structural_cells)


stromal_cluster_colors <- c(
  "Mix"                      = "#999999",
  "CXCL14 MSC"                   = "#3CB44B",
  "Inflamed_structural_cell" = "#984EA3",# Bright green
  "THY1 MSC"                 = "#66C2A5",  # Soft teal
  "APOD MSC"                 = "#8DD3C7",  # Pale turquoise
  "CXCL14 MSC"                 = "#A6D854",  # Light green
  "Fibro MSC"                = "#B3DE69",  # Lime green
  "LEPR MSC"                = "#2E8B57",  # Yellow
  "Adipocyte"             = "#FFD92F",
  "VSMC"            = "#D4A5A5",# Pale yellow
  "EC"                       = "#D73027",  # Vivid red
  "SEC"                      = "#F46D43",  # Orange-red
  "AEC"                      = "#FDAE61",  # Peach/light coral
  "Osteoblast"        = "#D9D9D9",  # Light gray (bone tone)
  "Osteo-Fibroblastic MSC"      = "#A65628"  # Slightly darker gray
)

Idents(structural_cells  ) <- structural_cells@meta.data $annot_structure_second_correct
n_cells <- ncol(structural_cells)
plot <- DimPlot(structural_cells , label = T, repel = T, raster = T, reduction = "harmony_umap_structural2", cols =  stromal_cluster_colors ) #, cols = cluster_colors) #+  

plot 
pdf(file= file.path(output_dir,"Figure_1","Umap_stroma.pdf"), 
    width = 13, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


ordered_precise_levels_adjusted_corrected <- c(
  "Interleukin Producing Cell",
  "Interferon Stimulated Cell",
  "Erythroid",
  "EC",
  "SEC",
  "AEC",
  "Osteoblast",
  "Osteo-Fibroblastic MSC",
  "Fibro MSC",
  "LEPR MSC",
  "THY1 MSC",
  "APOD MSC",
  "CXCL14 MSC",
  "Adipocyte",
  "VSMC",
  "MKC Lineage",
  "Basophil Eosinophil",
  "Mast Cell",
  "Neutrophil Progenitor",
  "CEACAM6 Immature Neutrophil 1",
  "OLFM4 Immature Neutrophil 2",
  "LTF Immature Neutrophil 3",
  "MMP9 Immature Neutrophil 4",
  "Mature Neutrophil",
  "Classical Monocyte",
  "Non Classical Monocyte",
  "Macrophage",
  "Conventional DC",
  "Plasmacytoid DC",
  "CD8 T Cell",
  "Activated or Exhausted T Cell",
  "NK or Cytotoxic T Cell",
  "Naive or CM CD4 T Cell",
  "Regulatory T Cell",
  "B or Plasma Cell",
  "Cycling B or Plasma Cell",
  "Plasma Cell"
)


ordered_genes_structure <- c(
  
  "BTNL9", "DNASE1L3", "ENG",  "VWF", "PLVAP", "SEMA3G", "SPP1", "BGLAP", "COL1A1",
  "TNC","COMP", "LEPR","NGFR", "THY1", "NOTCH3","STEAP4",  "APOD", "MGP", "CXCL14",
  "PLIN4", "ADIPOQ", "LPL", "TIMP4", "MYH11", "ACTA2")

Idents(structural_cells) <- structural_cells@meta.data$annot_structure_second_correct
dotplot_df <- DotPlot(structural_cells  ,features = ordered_genes_structure, assay = "SCT")
dot_data <- dotplot_df$data
#dot_data$id <- factor(dot_data$id, levels = cluster_colors) #choose!
dot_data$id_factor <- factor(dot_data$id, levels = ordered_precise_levels_adjusted_corrected) #choose!

library(ggplot2)



plot <- ggplot(dot_data, aes(x = id_factor, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "lightblue", high = "firebrick") +
  theme_minimal() +
  coord_flip()+
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Cluster names
    axis.text.y = element_text(size = 14, face = "italic")         # Gene names
  )+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = "Gene", color = "Avg Expr (scaled)", size = "% expressed")


plot 
pdf(file= file.path(output_dir,"Figure_1","Dotplot_stroma.pdf"), 
    width = 13, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Figure 1D Barplot ---------------------------------------------------------------

meta_df <- Xenium_Object_Seurat@meta.data
meta_df_NDMM_SMM_CBM_PCL <- meta_df[meta_df$Status_simp_2 %in% c("NDMM", "SMM", "CBM", "PCL"), ]

df_plot <- meta_df %>%
  dplyr::select(
    ObjectName = ObjectName_anonymous_grouped,
    annot_evenless_broad_run3 = annot_broad_run3_correct_ifndiv 
  )

# Step 3: Calculate proportions
df_plot <- df_plot %>%
  group_by(ObjectName, annot_evenless_broad_run3) %>%
  tally() %>%  # Count occurrences of each combination
  ungroup() %>%
  group_by(ObjectName) %>%
  mutate(proportion = n / sum(n))  # Calculate proportion within each ObjectName

# Ensure ObjectName is a factor with custom order
df_plot$ObjectName_fac <- factor(df_plot$ObjectName, levels = patient_order)
#df_plot$ObjectName <- factor(df_plot$ObjectName, levels = architecture_order)

df_plot <- df_plot %>% filter(!(ObjectName %in% c(NA, "unqualified")))

df_plot$cluster_factor <- factor(
  df_plot$annot_evenless_broad_run3,        # replace with your actual cluster column
  levels = names(cluster_colors)  # order according to the color vector
)

# join IMC data
library(dplyr)

plot <- ggplot(df_plot, aes(x = ObjectName_fac, y = proportion, fill = cluster_factor       )) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot
  scale_fill_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x-axis labels and increase size
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove background
        axis.title.x = element_blank(),  # Optional: Remove x-axis title
        axis.title.y = element_blank(),  # Optional: Remove y-axis title
        text = element_text(size = 15))#+
#coord_flip()  # Flip x and y axes  # Customize text size if needed
plot


pdf(file= file.path(output_dir,"Figure_1/barplot_patients_annot_broad_with_relapse3.pdf"), 
    width = 8.5, height = 9, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

# Figure 1E-H images only tumor cells ---------------------------------------------------------------

meta_df <- Xenium_Object_Seurat@meta.data
meta_df$y_flipped <- -meta_df$y_corrected_together
bounds_df <-  read.xlsx("patient_coordinates_x_and_y.xlsx")


pdf(file= file.path(output_dir,"Figure_1","Patient_only_tumor_cellen.pdf"))

patients <- unique(meta_df$Tissue_ext)

for (pat in patients) {
  df_sub <- meta_df %>% filter(Tissue_ext == pat)
  bounds <- bounds_df %>% filter(Patient == pat)

  df_sub <- df_sub %>%
    mutate(color_group = ifelse(annot_evenless_broad_run3 == "B/Plasma cells",
                                "Plasma", "Other"))
  
  p <- ggplot(df_sub, aes(x = x_corrected_together, 
                          y = y_flipped)) +
    geom_point(aes(color = color_group), size = 0.5) +
    labs(title = pat) +
    scale_color_manual(values = c("Other" = "gray60", "Plasma" = "red")) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    annotate("segment", 
             x = bounds$xlim_min + 200, 
             xend = bounds$xlim_min + 1200, 
             y = -bounds$ylim_min - 200, 
             yend = -bounds$ylim_min - 200, 
             colour = "black", size = 0.8) +
    annotate("text", 
             x = bounds$xlim_min + 700, 
             y = -bounds$ylim_min - 350, 
             label = "1 mm", vjust = 0, hjust = 0.5) +
    coord_fixed(xlim = c(bounds$xlim_min, bounds$xlim_max),
                ylim = c(-bounds$ylim_max, -bounds$ylim_min))
  
  print(p)
}

dev.off()


library(ggplot2)
library(dplyr)
library(hexbin)  # needed for geom_hex
library(openxlsx)


meta_df <- Xenium_Object_Seurat@meta.data
meta_df$y_flipped <- -meta_df$y_corrected_together
bounds_df <-  read.xlsx("patient_coordinates_x_and_y.xlsx")

pdf(file= file.path(output_dir,"Figure_1","Patient_only_tumor_cellen_hexbin.pdf"), width = 12, height = 9)

library(ggplot2)
library(dplyr)
library(ggrastr)  # rasterized geom_point

# Open PDF with Cairo for better raster support
cairo_pdf("Figures/Figure_1/Patient_only_tumor_cellen_rast.pdf", width = 12, height = 9)

patients <- unique(meta_df$Tissue_ext)

for (pat in patients) {
  df_sub <- meta_df %>% filter(Tissue_ext == pat)
  bounds <- bounds_df %>% filter(Patient == pat)
  
  # Define color group: Plasma vs Other
  df_sub <- df_sub %>%
    mutate(color_group = ifelse(annot_evenless_broad_run3 == "B/Plasma cells",
                                "Plasma", "Other"))
  
  p <- ggplot(df_sub, aes(x = x_corrected_together, y = y_flipped)) +
    geom_point_rast(aes(color = color_group), size = 0.5) +  # rasterized points
    scale_color_manual(values = c("Other" = "gray60", "Plasma" = "red")) +
    labs(title = pat) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    annotate("segment", 
             x = bounds$xlim_min + 200, 
             xend = bounds$xlim_min + 1200, 
             y = -bounds$ylim_min - 200, 
             yend = -bounds$ylim_min - 200, 
             colour = "black", size = 0.8) +
    annotate("text", 
             x = bounds$xlim_min + 700, 
             y = -bounds$ylim_min - 350, 
             label = "1 mm", vjust = 0, hjust = 0.5) +
    coord_fixed(xlim = c(bounds$xlim_min, bounds$xlim_max),
                ylim = c(-bounds$ylim_max, -bounds$ylim_min))
  
  print(p)
}

dev.off()


# Figure 1I Barplot interacting PC - PC  ---------------------------------------------------------------


folder <- "/interactions_PC_PC/Per_status_per_archi/"
files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)

library(ggplot2)
library(data.table)



all_summary_patient <- rbindlist(lapply(files, function(f) {
  dt <- readRDS(f)
  dt[, .(n_neighbors, n_cells, patient)]
}))

# Aggregate by patient
agg_patient <- all_summary_patient[, .(n_cells = sum(n_cells)), by = .(patient, n_neighbors)]

# Total PCs per patient
agg_patient[, total_pcs := sum(n_cells), by = patient]

# Percentages
agg_patient[, percent := n_cells / total_pcs * 100]

# Cumulative percentages
agg_patient[, pct_at_least1 := sum(n_cells[n_neighbors >= 1]) / sum(n_cells) * 100, by = patient]
agg_patient[, pct_at_least2 := sum(n_cells[n_neighbors >= 2]) / sum(n_cells) * 100, by = patient]


# Optional: reshape wide
pct_patient_wide <- dcast(
  agg_patient, 
  patient + pct_at_least1 + pct_at_least2  ~ n_neighbors, 
  value.var = "percent", fill = 0
)
setnames(pct_patient_wide, old = as.character(0:max(all_summary_patient$n_neighbors)), 
         new = paste0("pct_", 0:max(all_summary_patient$n_neighbors)))



# Copy table
dt <- copy(pct_patient_wide)

# Extract disease stage from patient name
dt[, disease_stage := sub("_.*", "", patient)]

# Compute mean and SD per disease stage
# Compute mean and SD per disease stage
summary_dt <- dt[, .(
  mean_at_least1 = mean(pct_at_least1),
  sd_at_least1   = sd(pct_at_least1),
  mean_at_least2 = mean(pct_at_least2),
  sd_at_least2   = sd(pct_at_least2),
), by = disease_stage]

summary_dt$disease_stage_fac <- factor(summary_dt$disease_stage, 
                                       levels = c("CBM", "SMM", "NDMM", "PCL", "Relapse"))

status_colors <- c(  #fir annot even less broad
  
  "CBM" = "#46F0F0",
  "IMC_CBM" = "#46F0F0",
  "SMM" = "#4575B4",              # Macrophages
  "NDMM" = "#FF7F00", 
  "IMC_NDMM" = "#46F0F0",
  "PCL" = "#D73027",
  "Relapse" = "#984EA3" #New: yellow
)

dt[, disease_stage_fac := factor(
  disease_stage,
  levels = c("CBM", "SMM", "NDMM", "PCL", "Relapse")
)]

dt2 <- dt[!disease_stage %in% "Relapse"]



welch_1 <- compare_means(
  pct_at_least1 ~ disease_stage_fac,
  data = dt2,
  method = "t.test",
  var.equal = FALSE,
  p.adjust.method = "BH"
)

welch_2 <- compare_means(
  pct_at_least2 ~ disease_stage_fac,
  data = dt2,
  method = "t.test",
  var.equal = FALSE,
  p.adjust.method = "bonferroni"
)



welch_2 <-  welch_2 %>%
  mutate(
    p.adj.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

welch_1 <-  welch_1 %>%
  mutate(
    p.adj.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )


p1 <- ggplot(summary_dt, aes(x = disease_stage_fac, y = mean_at_least1, fill = disease_stage_fac)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_at_least1 - sd_at_least1, ymax = mean_at_least1 + sd_at_least1), width = 0.2) +
  labs(title = "Percentage of PCs interacting with at least 1 PC", x = "Disease stage", y = "Mean Percentage ± SD") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = status_colors) +
  theme(legend.position = "none")

p2 <- ggplot(summary_dt, aes(x = disease_stage_fac, y = mean_at_least2, fill = disease_stage_fac)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_at_least2 - sd_at_least2, ymax = mean_at_least2 + sd_at_least2), width = 0.2) +
  labs(title = "Percentage of PCs interacting with at least 2 PCs", x = "Disease stage", y = "Mean Percentage ± SD") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = status_colors) +
  theme(legend.position = "none")

p1
p2

pdf(file= file.path(output_dir,"Figure_1/barplot_pc_pc_atleast1.pdf"), 
    width = 6, height = 9, useDingbats = FALSE, onefile = FALSE)
print(p1)
dev.off()

pdf(file= file.path(output_dir,"Figure_1/barplot_pc_pc_atleast2.pdf"), 
    width = 6, height = 9, useDingbats = FALSE, onefile = FALSE)
print(p2)
dev.off()


summary_dt_no_relapse <- summary_dt[summary_dt$disease_stage != "Relapse",]

dt_no_relapse <- dt[dt$disease_stage != "Relapse",]

ymax <- max(
  summary_dt_no_relapse$mean_at_least1 + summary_dt_no_relapse$sd_at_least1,
  summary_dt_no_relapse$mean_at_least2 + summary_dt_no_relapse$sd_at_least2
)



# Helper to compute a safe y max per plot (based on the raw jitter data)
ymax_atleast1 <- max(dt_no_relapse$pct_at_least1, na.rm = TRUE)
ymax_atleast2 <- max(dt_no_relapse$pct_at_least2, na.rm = TRUE)
ymax_atleast3 <- max(dt_no_relapse$pct_at_least3, na.rm = TRUE)

# A function to add a y.position per comparison row
add_y_positions <- function(df, base_y, step = 4) {
  # Order pairs in a stable way and assign increasing y levels
  df <- df[order(df$group1, df$group2), ]
  df$y.position <- base_y + seq_len(nrow(df)) * step
  df
}

# Build annotation tables with y.position
anno_T1 <- add_y_positions(welch_1, base_y = ymax_atleast1 - 20)  # for p1
anno_T2 <- add_y_positions(welch_2, base_y = ymax_atleast2 - 20)  # for p2

anno_T1_sig <- subset(anno_T1, p.adj.signif != "ns")
anno_T2_sig <- subset(anno_T2, p.adj.signif != "ns")



p1 <- ggplot(summary_dt_no_relapse, aes(x = disease_stage_fac, y = mean_at_least1, fill = disease_stage_fac)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(
    aes(ymin = mean_at_least1 - sd_at_least1,
        ymax = mean_at_least1 + sd_at_least1),
    width = 0.2
  ) +
  geom_jitter(
    data = dt_no_relapse,
    aes(x = disease_stage_fac, y = pct_at_least1),
    inherit.aes = FALSE,
    width = 0.2,
    size = 2.5,
    alpha = 0.8,
    color = "black"
  ) +
  labs(
    title = "Percentage of PCs interacting with at least 1 PC",
    x = "Disease stage",
    y = "Percentage per patient"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = status_colors) +
  theme(legend.position = "none")


p1_star <- p1 +
  stat_pvalue_manual(
    anno_T1_sig,
    label = "p.adj.signif",
    y.position = "y.position",
    xmin = "group1", xmax = "group2",
    tip.length = 0.01, size = 5, bracket.size = 0.6
  ) +
  coord_cartesian(ylim = c(0, max(anno_T1_sig$y.position) + 2))



p2 <- ggplot(summary_dt_no_relapse, aes(x = disease_stage_fac, y = mean_at_least2, fill = disease_stage_fac)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(
    aes(ymin = mean_at_least2 - sd_at_least2,
        ymax = mean_at_least2 + sd_at_least2),
    width = 0.2
  ) +
  geom_jitter(
    data = dt_no_relapse,
    aes(x = disease_stage_fac, y = pct_at_least2),
    inherit.aes = FALSE,
    width = 0.2,
    size = 2.5,
    alpha = 0.8,
    color = "black"
  ) +
  labs(
    title = "Percentage of PCs interacting with at least 2 PCs",
    x = "Disease stage",
    y = "Percentage per patient"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = status_colors) +
  theme(legend.position = "none")

p2_star <- p2 +
  stat_pvalue_manual(
    anno_T2_sig,
    label = "p.adj.signif",
    y.position = "y.position",
    xmin = "group1", xmax = "group2",
    tip.length = 0.01, size = 5, bracket.size = 0.6
  ) +
  coord_cartesian(ylim = c(0, max(anno_T1_sig$y.position) + 2))



p1 <- p1 + ylim(0, 100)
p2 <- p2 + ylim(0, 100)


library(patchwork)

p2_mod <- p2_star +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

combined <- p1_star + p2_mod + 
  plot_layout(ncol = 3)

combined


pdf(file= file.path(output_dir,"Figure_1/barplot_pc_pc_atleast1and2_no_relapsewithdots_WELCH_ttest.pdf"), 
    width = 8, height = 9, useDingbats = FALSE, onefile = FALSE)
print(combined)
dev.off()



# Figure 2B barplot tumor grouped ---------------------------------------------------------------

meta_df <- Xenium_Object_Seurat@meta.data
meta_df_NDMM_SMM_CBM_PCL <- meta_df[meta_df$Status_simp_2 %in% c("NDMM", "SMM", "CBM", "PCL", "Relapse"), ]


df_plot <- meta_df_NDMM_SMM_CBM_PCL %>%
  dplyr::select(
    ObjectName = ObjectName_anonymous_grouped,
    annot_evenless_broad_run3 = tumor_percentage_grouped035ext_hex_average_method2
  )



# Step 3: Calculate proportions
df_plot <- df_plot %>%
  group_by(ObjectName, annot_evenless_broad_run3) %>%
  tally() %>%  # Count occurrences of each combination
  ungroup() %>%
  group_by(ObjectName) %>%
  mutate(proportion = n / sum(n))  # Calculate proportion within each ObjectName

# Ensure ObjectName is a factor with custom order

df_plot$ObjectName_fac <- factor(df_plot$ObjectName, levels = patient_order)


#df_plot$ObjectName <- factor(df_plot$ObjectName, levels = architecture_order)

df_plot <- df_plot %>% filter(!(ObjectName %in% c(NA, "unqualified")))

df_plot$cluster_factor <- factor(
  df_plot$annot_evenless_broad_run3,        # replace with your actual cluster column
  levels = names(region_colors)  # order according to the color vector
)



plot <- ggplot(df_plot, aes(x = ObjectName_fac, y = proportion, fill = cluster_factor       )) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot
  scale_fill_manual(values = region_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x-axis labels and increase size
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove background
        axis.title.x = element_blank(),  # Optional: Remove x-axis title
        axis.title.y = element_blank(),  # Optional: Remove y-axis title
        text = element_text(size = 15))#+
#coord_flip()  # Flip x and y axes  # Customize text size if needed
plot

#plot <- dittoBarPlot(seurat_HBM_MM_ownSCT , var = "annot_MM_HBM_more_simple", group.by = "ObjectName", color.panel = cluster_colors) #,x.reorder = x_reorder_indices) #, colors = cluster_colors) #, scale = "count")  

pdf(file= file.path(output_dir,"Figure_2/barplot_patients_tumorgrouped_horizontal.pdf"), 
    width = 11, height = 4, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

# Figure 2D  ---------------------------------------------------------------

df_all_combined <-  readRDS("S_curves_combined_new_order.rds")

df_summary <- df_all_combined %>%
  group_by(clusters, group_comparison) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    mean_value   = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

library(ggplot2)

df_all_filtered <- df_all_combined[df_all_combined$group_comparison %in% c("all_SMM vs CBM", "all_PCL vs NDMM",  "all_PCL vs CBM", 
                                                                           "all_NDMM vs CBM","normal_NDMM vs CBM", "dense_vs_all_NDMM vs CBM",
                                                                           "all_NDMM vs SMM", "sparse_NDMM vs SMM","dense_vs_all_NDMM vs SMM"  ),]
df_all_filtered <- df_all_combined[df_all_combined$group_comparison %in% c(
  "normal_NDMM vs CBM", "dense_vs_all_NDMM vs CBM" ),]

df_all_filtered <- df_all_combined[df_all_combined$group_comparison %in% c(
  "all_NDMM vs CBM","normal_NDMM vs CBM", "dense_vs_all_NDMM vs CBM" ),]



comparison_order <- c("all_SMM vs CBM", "all_PCL vs NDMM",  "all_PCL vs CBM", 
                      "all_NDMM vs CBM","normal_NDMM vs CBM", "dense_vs_all_NDMM vs CBM",
                      "all_NDMM vs SMM", "sparse_NDMM vs SMM","dense_vs_all_NDMM vs SMM"  )
df_all_filtered$group_comparison <- factor(df_all_filtered$group_comparison, levels = comparison_order)



clip_value <- 5  # maximum |value| to display

df_all_filtered <- df_all_filtered %>%
  mutate(value_clipped = pmin(pmax(value, -clip_value), clip_value))



# Filter only the two comparisons used for ordering
df_order_calc <- df_all_filtered %>%
  filter(group_comparison %in% c("normal_NDMM vs CBM", "dense_vs_all_NDMM vs CBM")) %>%
  group_by(clusters, group_comparison) %>%
  summarise(mean_value = mean(value_clipped, na.rm = TRUE), .groups = "drop")

# Reshape wider so each cluster has two columns
df_order_wide <- df_order_calc %>%
  tidyr::pivot_wider(
    names_from = group_comparison,
    values_from = mean_value
  )

# Compute the difference
df_order_wide <- df_order_wide %>%
  mutate(diff_order = `normal_NDMM vs CBM` - `dense_vs_all_NDMM vs CBM`)

# Create the ordered factor for clusters
cluster_order <- df_order_wide %>%
  arrange(diff_order) %>%
  pull(clusters)

# Apply to your dataframe
#df_all_filtered$clusters <- factor(df_all_filtered$clusters, levels = cluster_order)

# or the same order
#df_all_filtered$clusters <- factor(df_all_filtered$clusters, levels = rev(cluster_order_NDMM_vs_CBM))

df_all_filtered$clusters <- factor(df_all_filtered$clusters, levels = rev(cluster_ordered_for_s_curves))


# Optional: compute positions for vertical lines between cell types
celltype_positions <- 1:length(unique(df_all_filtered$clusters))
vline_positions <- celltype_positions[-length(celltype_positions)] + 0.5

colors <- c( "#FFE066","#D72638") #, "dense_otheroption_toopink"= "#FB8072")  # Replace with your preferred order
colors <- c( "#A7C7E7", "#FFE066","#D72638") #, "dense_otheroption_toopink"= "#FB8072")  # Replace with your preferred order



plot <- ggplot(df_all_filtered, aes(x = clusters, y = value_clipped, fill = group_comparison)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
               outlier.shape = NA ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = vline_positions, linetype = "dotted", color = "gray70") +
  # Optional: add jitter
  # geom_jitter(aes(color = group_comparison),
  #             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
  #             size = 2, alpha = 0.6) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(-5, 4)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = "Log2 ratio (clipped at ±5)",
    x = "Cell type",
    fill = "Comparison"
  )



plot

pdf(file= file.path(output_dir,"Sup_Figure_3","S_curves_linked_CBM_normalanddensesame_order3plots.pdf"), width = 16, height = 6)
print(plot)
dev.off()




# Figure 3A dense vs sparse linked ---------------------------------------------------------------

meta_df <- Xenium_Object_Seurat@meta.data

metadata_onlynormalsparseanddense <- meta_df[!meta_df$tumor_percentage_grouped035ext_hex_average_method2 %in% 
                                               c("unqualified", "dispersed"), ]
metadata_onlynormalsparseanddense_heterogenous <- metadata_onlynormalsparseanddense %>%
  filter(ObjectName_anonymous_grouped %in% c("PCL_2", "NDMM_3", "NDMM_4", "NDMM_5","PCL_1"))

#metadata_onlynormalsparseanddense_heterogenous <- metadata_onlynormalsparseanddense %>%
# filter(ObjectName_anonymous_grouped %in% c("PCL_2", "NDMM_3", "NDMM_4", "NDMM_5"))


metadata_onlynormalsparseanddense_heterogenous$arch_grouped <- ifelse(
  metadata_onlynormalsparseanddense_heterogenous$tumor_percentage_grouped035ext_hex_average_method2_CBM %in% c("dense", "dispersed"),
  "dense",
  "sparse"
)

metadata_onlynormalsparseanddense_heterogenous <- metadata_onlynormalsparseanddense_heterogenous[!metadata_onlynormalsparseanddense_heterogenous$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves%in% 
                                                                                                   c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell"), ]

metadata_onlynormalsparseanddense_heterogenous$sample_cond <- paste(metadata_onlynormalsparseanddense_heterogenous$ObjectName_anonymous_grouped, metadata_onlynormalsparseanddense_heterogenous$arch_grouped, sep = "_")

metadata_onlynormalsparseanddense_heterogenous$sample_cond <- as.character(metadata_onlynormalsparseanddense_heterogenous$sample_cond)


# 3. Compute proportions per patient (normal+sparse combined)
prop_df <- metadata_onlynormalsparseanddense_heterogenous %>%
  group_by(sample_cond,
           annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(sample_cond) %>%
  mutate(
    Freq = cell_count / sum(cell_count) * 100
  ) %>%
  ungroup()



# 4. create wide table like propeller output

prop_trans_res_sparse <- prop_df %>%
  rename(
    sample = sample_cond,
    clusters = annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves
  ) %>%
  select(sample, clusters, Freq)


# 5. join status

metadata_unique <- metadata_onlynormalsparseanddense_heterogenous %>%
  distinct(sample_cond, Status_simp_2,arch_grouped)

prop_trans_res_sparse <- prop_trans_res_sparse %>%
  left_join(metadata_unique,
            by = c("sample" = "sample_cond"))


# 6. clean nonsense values

prop_trans_res_sparse <- prop_trans_res_sparse %>%
  filter(!is.na(Freq),
         !is.nan(Freq))

prop_trans_res_sparse



pdf(file= file.path(output_dir,"Figure_2","dense_vs_sparse_not_normalized.pdf"), width = 8, height = 6)


p <- plot_patient_deltas_boxplot_with_stat_same_order_capped_paired(prop_trans_res_sparse, cluster_order_defined = cluster_ordered_for_s_curves, tumor = FALSE, method = "ratio", adjust_p = TRUE) +
  ggtitle(paste0("Ratio Proportions Dense vs Sparse Paired"))

print(p)

dev.off()




meta_df <- Xenium_Object_Seurat@meta.data
metadata_onlynormalsparseanddense <- meta_df[!meta_df$tumor_percentage_grouped035ext_hex_average_method2 %in% 
                                               c("dispersed", "unqualified"), ]
metadata_onlynormalsparseanddense_heterogenous <- metadata_onlynormalsparseanddense %>%
  filter(ObjectName_anonymous_grouped %in% c("PCL_2", "NDMM_3", "NDMM_4", "NDMM_5"))


dataframe_no_tumor <- prepare_metadata_frame(
  meta = metadata_onlynormalsparseanddense_heterogenous,
  exclude_tumor = TRUE,
  cluster_col = "annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv",
  region_col = "small_region_hex1",
  tumor_col = "tumor_percentage_grouped035ext_hex_average_method2",
  levels_tumor = unique(Xenium_Object_Seurat$tumor_percentage_grouped035ext_hex_average_method2)
)

dataframe_no_tumor <- aggregate_percentages_dynamic(dataframe_no_tumor, group_col = "tumor_percentage_grouped")


### check if dataframe_no_tumor is still full!!
mean_freqs <- dataframe_no_tumor %>%
  group_by(ObjectName_anonymous_grouped , tumor_percentage_grouped,cluster_celltype  ) %>%
  summarise(mean_freq = mean(Freq), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tumor_percentage_grouped, values_from = mean_freq)

library(dplyr)
library(tidyr)
library(ggplot2)

# Convert wide to long for plotting
mean_freqs_long <- mean_freqs %>%
  pivot_longer(
    cols = c(sparse, dense),           # only the tumor density columns you care about
    names_to = "tumor_group",
    values_to = "mean_freq"
  )

# Optional: order the tumor group factor
mean_freqs_long$tumor_group <- factor(mean_freqs_long$tumor_group, levels = c("sparse", "dense"))

# Step 1: Compute delta or ratio for each patient and cluster
delta_df <- mean_freqs %>%
  rowwise() %>%
  mutate(
    sparse_avg = mean(c(sparse, `normal PC percentage`), na.rm = TRUE),  # average of sparse + normal
    log2_ratio_dense_sparse = log2(dense / sparse_avg)                   # log2 ratio
  ) %>%
  ungroup()

#delta_df <- delta_df %>%
#mutate(log2_ratio_dense_sparse = log2((dense + 0.001) / (sparse + 0.001)))
delta_df <- delta_df %>%
  filter(!is.na(log2_ratio_dense_sparse), is.finite(log2_ratio_dense_sparse))
# Step 2: Order clusters by mean log2 ratio across patients (largest at top)
cluster_order <- delta_df %>%
  group_by(cluster_celltype) %>%
  summarise(median_log2_ratio = median(log2_ratio_dense_sparse, na.rm = TRUE)) %>%
  arrange(median_log2_ratio) %>%
  pull(cluster_celltype)

delta_df <- delta_df %>%
  mutate(cluster_celltype_fac = factor(cluster_celltype, levels = cluster_order))



delta_df$cluster_celltype_fac <- factor(delta_df$cluster_celltype_fac, levels = cluster_order)
delta_df$patient_id <- as.factor(delta_df$ObjectName_anonymous_grouped )

# Step 1: calculate significance per cluster
# Step 1: keep only relevant columns and rows
paired_df <- delta_df %>%
  dplyr::select(patient_id, cluster_celltype_fac, sparse_avg, dense) %>%
  filter(!is.na(sparse_avg) & !is.na(dense))  # remove rows missing one group

# Step 2: compute paired t-test per cluster
sig_df <- paired_df %>%
  group_by(cluster_celltype_fac) %>%
  summarise(
    p_value = tryCatch(
      t.test(dense, sparse_avg, paired = TRUE)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(star = ifelse(!is.na(p_value) & p_value < 0.05, "*", ""))

# Compute star positions dynamically per cluster
star_positions <- delta_df %>%
  group_by(cluster_celltype_fac) %>%
  summarise(
    median_val = median(log2_ratio_dense_sparse, na.rm = TRUE),
    max_val = max(log2_ratio_dense_sparse, na.rm = TRUE),
    min_val = min(log2_ratio_dense_sparse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(sig_df, by = "cluster_celltype_fac") %>%
  mutate(
    x_pos = ifelse(median_val >= 0, max_val * 1.2, min_val * 1.2)
  )

# Plot with dynamic star positions
plot <- ggplot(delta_df, aes(x = log2_ratio_dense_sparse, y = cluster_celltype_fac)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(width = 0.5, size = 0.9 ,alpha = 0.3, outlier.shape = NA, fill = "gray90", color = "black", width = 0.6)+
  geom_jitter(aes(color = ObjectName_anonymous_grouped), height = 0.15, size = 2) +
  geom_text(
    data = star_positions,
    aes(x = x_pos, y = cluster_celltype_fac, label = star),
    inherit.aes = FALSE,
    size = 5
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )+
  #theme_minimal(base_size = 12) +
  labs(
    x = "log2(Dense / Sparse)",
    y = "Cell type (ordered by dense/sparseratio)",
    title = "Change in cell frequencies: Dense vs Sparse"
  )
plot

pdf(file= file.path(output_dir,"Figure_2","densevssparse_hetero.pdf"), 
    width = 9, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Figure 3B Horizontal all combined per archi  ---------------------------------------------------------------


meta_df <- Xenium_Object_Seurat@meta.data
metadata_onlynormalsparseanddense <- meta_df[!meta_df$tumor_percentage_grouped035ext_hex_average_method2 %in% 
                                               c("unqualified"), ]


low_groups <- metadata_onlynormalsparseanddense %>%
  group_by(ObjectName_anonymous_grouped, tumor_percentage_grouped035ext_hex_average_method2) %>%
  summarise(
    n_regions = n(),
    .groups = "drop"
  ) %>%
  group_by(ObjectName_anonymous_grouped) %>%
  mutate(
    total_regions = sum(n_regions),
    prop = n_regions / total_regions
  ) %>%
  ungroup() %>%
  filter(prop < 0.03) %>% ## if a region is less than 3%, it is removed
  dplyr::select(ObjectName_anonymous_grouped, tumor_percentage_grouped035ext_hex_average_method2)

low_groups_df <- low_groups %>%
  rename(ObjectName_anonymous_grouped = ObjectName_anonymous_grouped,
         tumor_percentage_grouped035ext_hex_average_method2 = tumor_percentage_grouped035ext_hex_average_method2)


# Remove rows in dataframe_no_tumor that match these patient-tumor combinations

metadata_onlynormalsparseanddense <- metadata_onlynormalsparseanddense %>%
  anti_join(low_groups_df, by = c("ObjectName_anonymous_grouped", "tumor_percentage_grouped035ext_hex_average_method2"))

metadata_onlynormalsparseanddense_heterogenous <- metadata_onlynormalsparseanddense

metadata_onlynormalsparseanddense_heterogenous <- metadata_onlynormalsparseanddense_heterogenous[!metadata_onlynormalsparseanddense_heterogenous$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% 
                                                                                                   c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell"), ]

metadata_onlynormalsparseanddense_heterogenous$sample_cond <- paste(metadata_onlynormalsparseanddense_heterogenous$ObjectName_anonymous_grouped, metadata_onlynormalsparseanddense_heterogenous$tumor_percentage_grouped035ext_hex_average_method2, sep = "_")

metadata_onlynormalsparseanddense_heterogenous$sample_cond <- as.character(metadata_onlynormalsparseanddense_heterogenous$sample_cond)




# 3. Compute proportions per patient (normal+sparse combined)
prop_df <- metadata_onlynormalsparseanddense_heterogenous %>%
  group_by(sample_cond,
           annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(sample_cond) %>%
  mutate(
    Freq = cell_count / sum(cell_count) * 100
  ) %>%
  ungroup()



# 4. create wide table like propeller output

prop_trans_res_horizontal <- prop_df %>%
  rename(
    sample = sample_cond,
    clusters = annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv
  ) %>%
  select(sample, clusters, Freq)


# 5. join status

metadata_unique <- metadata_onlynormalsparseanddense_heterogenous %>%
  distinct(sample_cond, Status_simp_2, tumor_percentage_grouped035ext_hex_average_method2, ObjectName_anonymous_grouped )

prop_trans_res_horizontal <- prop_trans_res_horizontal %>%
  left_join(metadata_unique,
            by = c("sample" = "sample_cond"))


# 6. clean nonsense values

prop_trans_res_horizontal <- prop_trans_res_horizontal %>%
  filter(!is.na(Freq),
         !is.nan(Freq))

prop_trans_res_horizontal$clusters <- as.factor(prop_trans_res_horizontal$clusters)


statuses <- unique(prop_trans_res_horizontal$tumor_percentage_grouped035ext_hex_average_method2)
combs_ordered <- expand.grid(group1 = statuses, group2 = statuses, stringsAsFactors = FALSE)
combs_ordered <- combs_ordered[combs_ordered$group1 != combs_ordered$group2, ]

# Order by group1, then group2 alphabetically
combs_ordered <- combs_ordered[order(combs_ordered$group1, combs_ordered$group2), ]



pdf(file= file.path(output_dir,"Figure_2","ratio_boxplots_allArchiesCombined_adjusted_normalTtest_multisample_adjusted9.pdf"), width = 8, height = 6)

for (i in seq_len(nrow(combs_ordered))) {
  group1 <- combs_ordered$group1[i]
  group2 <- combs_ordered$group2[i]
  
  p <- plot_patient_deltas_boxplot_with_stat_same_order_capped_test_archi_lmer(prop_trans_res_horizontal , group1 = group1, group2 = group2, method = "ratio",cluster_order_defined = cluster_ordered_for_s_curves, tumor = FALSE) +
    ggtitle(paste0("Ratio Proportions: ", group1, " vs ", group2))
  
  print(p)
}

dev.off()


mcsaveRDS(stat_df, "~/stat_df_normalpcvsdensehorizontalarchiplot.rds")


# Figure 3C  ----------------------------------------------------------------

# see Figure 4I

# Figure 3D proportion per broad annotation per architecture ---------------------------

meta_df <- Xenium_Object_Seurat@meta.data


## choose 1 or 2
##1
# Aggregate proportions per broad cell type, subtype, and status
composition_df <- meta_df %>%
  group_by(tumor_percentage_grouped035ext_hex_average_method2_CBM, annot_even_less_broad_correct_ifndiv, annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(tumor_percentage_grouped035ext_hex_average_method2_CBM, annot_even_less_broad_correct_ifndiv) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%  # fraction within broad type
  ungroup()


##2
composition_df <- meta_df %>%
  group_by(tumor_percentage_grouped035ext_hex_average_method2_CBM, annot_most_broad_correct_ifndiv, annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(tumor_percentage_grouped035ext_hex_average_method2_CBM, annot_most_broad_correct_ifndiv) %>%
  mutate(prop = n_cells / sum(n_cells),
         annot_even_less_broad_correct_ifndiv =   annot_most_broad_correct_ifndiv    ) %>%  # fraction within broad type
  ungroup()


broad_types <- unique(composition_df$annot_even_less_broad_correct_ifndiv)


composition_df <- composition_df %>%
  filter(tumor_percentage_grouped035ext_hex_average_method2_CBM != "unqualified")

composition_df$tumor_percentage_grouped035ext_hex_average_method2_CBM <- 
  factor(composition_df$tumor_percentage_grouped035ext_hex_average_method2_CBM,
         levels = c("CBM", "normal PC percentage","sparse", "dispersed", "dense"))

# Plot: one facet per broad cell type, bars show subtypes per status
ggplot(composition_df, aes(x = tumor_percentage_grouped035ext_hex_average_method2_CBM, y = prop, fill = annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves)) +
  geom_bar(stat = "identity", position = "stack") +  # stacked bar
  facet_wrap(~ annot_even_less_broad_correct_ifndiv, scales = "free_y") +  # one facet per broad type
  scale_y_continuous(labels = scales::percent_format()) +
  #scale_fill_manual(values = cluster_colors) +  # your predefined colors per subtype
  theme_minimal() +
  labs(
    x = "Status",
    y = "Fraction within broad cell type",
    fill = "Subtype"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Output PDF
pdf(file= file.path(output_dir,"Sup_Figure_5","composition_per_broadcelltype_per_region_test_withlegend.pdf"), width = 4, height = 6)

for (bt in broad_types) {
  
  # Filter data for this broad type
  df_sub <- composition_df %>% filter(annot_even_less_broad_correct_ifndiv == bt)
  
  # Keep only subtypes present in this subset for legend
  present_subtypes <- unique(df_sub$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves)
  
  # Optional: subset cluster_colors to only the present subtypes
  cluster_colors_sub <- cluster_colors[present_subtypes]
  
  n_subtypes <- length(unique(df_sub$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves))
  palette_colors <- brewer.pal(min(12, n_subtypes), "Set3")  # nice pastel palette
  names(palette_colors) <- unique(df_sub$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves)
  
  p <- ggplot(df_sub, aes(
    x = tumor_percentage_grouped035ext_hex_average_method2_CBM,
    y = prop,
    fill = annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves
  )) +
    geom_bar(stat = "identity", position = "stack", width = 0.6, color = "black", size = 0.2) +  # thinner bars, black outlines
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = palette_colors) +  # use nicer palette
    theme_minimal(base_size = 14) +
    labs(
      x = "Tumor Architecture",
      y = "Fraction within broad cell type",
      fill = "Subtype",
      title = paste("Composition of", bt)
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),  # remove vertical grid
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.3), # light horizontal lines
      #legend.position = "none",
      legend.key.size = unit(0.5, "cm"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Print plot to PDF (each print is a new page)
  print(p)
}
dev.off()


#  Figure 3F ratio dense vs sparse IMC -------------------------------------------------------------------

meta_df <- colData(spe_IMC)
#meta_df_NDMM_SMM_CBM_PCL <- meta_df[meta_df$Status_simp_2 %in% c("NDMM", "SMM", "CBM", "PCL"), ]

metadata_onlynormalsparseanddense <- meta_df[!meta_df$tumor_percentage_grouped0.35 %in% 
                                               c("unqualified", "dispersed"), ]


metadata_onlynormalsparseanddense <- metadata_onlynormalsparseanddense[metadata_onlynormalsparseanddense$status %in% 
                                                                         c("NDMM", "PCL"), ]


#at least 2000 dense cells
metadata_onlynormalsparseanddense <- metadata_onlynormalsparseanddense[!metadata_onlynormalsparseanddense$ObjectName_anonymous_grouped %in% 
                                                                         c("PCL_3", "PCL_1", "NDMM_1"), ]
#metadata_onlynormalsparseanddense_heterogenous <- metadata_onlynormalsparseanddense %>%
# filter(ObjectName_anonymous_grouped %in% c("PCL_2", "NDMM_3", "NDMM_4", "NDMM_5"))


# For IMC metadata (meta_df2)
metadata_onlynormalsparseanddense$common_group <- dplyr::case_when(
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Ki67_likely_erythroid") ~ "Erythroid",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Plasma_Cell") ~ "Plasma Cell",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD105_Structural_Cell", "Stroma") ~ "Stroma",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("MKC_Endothelial") ~ "MKC or Endothelial",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Macrophage") ~ "Macrophage",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Neutrophil") ~ "Neutrophil",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Neutrophil_Progenitor") ~ "Neutrophil Progenitor",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD14_Monocyte") ~ "Monocyte",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD8_T_Cell") ~ "CD8 T Cell",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD4_T_Cell") ~ "CD4 T Cell",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("HLA-DR_Cell") ~ "APC",
  metadata_onlynormalsparseanddense$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Fibrotic_Tissue", "Undefined_Cell") ~ "Other",
  TRUE ~ "Other"
)



metadata_onlynormalsparseanddense$arch_grouped <- ifelse(
  metadata_onlynormalsparseanddense$tumor_percentage_grouped0.35 %in% c("dense", "dispersed"),
  "dense",
  "sparse"
)

metadata_onlynormalsparseanddense_heterogenous <- metadata_onlynormalsparseanddense

metadata_onlynormalsparseanddense_heterogenous <- metadata_onlynormalsparseanddense_heterogenous[!metadata_onlynormalsparseanddense_heterogenous$common_group %in% 
                                                                                                   c("Plasma Cell"), ]

metadata_onlynormalsparseanddense_heterogenous$sample_cond <- paste(metadata_onlynormalsparseanddense_heterogenous$ObjectName_anonymous_grouped, metadata_onlynormalsparseanddense_heterogenous$arch_grouped, sep = "_")

metadata_onlynormalsparseanddense_heterogenous$sample_cond <- as.character(metadata_onlynormalsparseanddense_heterogenous$sample_cond)


props_densevssparse_linked     <- getTransformedProps(clusters = metadata_onlynormalsparseanddense_heterogenous $common_group, 
                                                      sample = metadata_onlynormalsparseanddense_heterogenous$sample_cond,
                                                      transform = "logit")


metadata_onlynormalsparseanddense_heterogenous
metadata_unique <- metadata_onlynormalsparseanddense_heterogenous[!duplicated(metadata_onlynormalsparseanddense_heterogenous$sample_cond), ]
rownames(metadata_unique) <- metadata_unique$sample_cond

#metadata_unique$source <- factor(c(rep("control", 5), rep("myeloma", 19)), levels=c("control","myeloma"))
model  <- model.matrix(~0 + arch_grouped , data = metadata_unique)



props_densevssparse_linked$cluster <- rownames(props_densevssparse_linked)


prop_trans_res_densevssparse_linked <- as.data.frame(props_densevssparse_linked$Proportions)



metadata_join <- as_tibble(metadata_unique) %>%
  rename(sample = sample_cond) %>%   # now dplyr::rename will work
  select(sample, status, arch_grouped, ObjectName_anonymous_grouped)
# Or Source, if it's uppercase
prop_trans_res_densevssparse_linked <- prop_trans_res_densevssparse_linked %>%
  left_join(metadata_join, by = "sample")

### very important, remove non sense samples!

prop_trans_res_densevssparse_linked <- prop_trans_res_densevssparse_linked[!is.nan(prop_trans_res_densevssparse_linked$Freq) & !is.na(prop_trans_res_densevssparse_linked$Freq), ]


prop_trans_res_densevssparse_linked <- prop_trans_res_densevssparse_linked[
  !prop_trans_res_densevssparse_linked$clusters %in% c("Plasma_Cell"),]


prop_trans_res_densevssparse_linked <- prop_trans_res_densevssparse_linked %>%
  mutate(Status_simp_2 = status)



pdf(file= file.path(output_dir,"Sup_Figure_5_IMC","ratio_boxplots_ordered_combinations_annot_precise_rem_low_quali_IFNdiv_stats__densevssparse_linked_same_order_capped_s_curve_order_withlines_adjusted2.pdf"), width = 6, height = 3)


p <- plot_patient_deltas_boxplot_with_stat_same_order_capped_paired_broad(prop_trans_res_densevssparse_linked, cluster_order_defined = names(common_group_colors), tumor = FALSE, method = "ratio") +
  ggtitle(paste0("Ratio Proportions Dense vs Sparse Paired"))

print(p)


dev.off()


# Figure 3G neutrophil Madelon BAFF ---------------------------------------------------------------
Neutrophils <- mcreadRDS("~/neutrophils_SCTtransformed.rds")

Idents(Neutrophils) <- Neutrophils$annot
# Rename Idents with singular names and consistent formatting
Neutrophils <- RenameIdents(object = Neutrophils,
                            "MatNeu1"               = "Late Immature Neutrophil",
                            "MatNeu2"              = "Mature Neutrophil",
                            "ImmNeu"            = "Early Immature Neutrophil",
                            "PreNeu2"             = "Early Immature Neutrophil",
                            "PreNeu1"             = "Early Immature Neutrophil",
                            "Myelocytes"             = "Neutrophil Progenitor",
                            "MKI67+ myelocytes"             = "Neutrophil Progenitor",
                            "Erythrocytes"             = "Erythroid")

Neutrophils$annot_Marnix  <- Idents(Neutrophils)

Neutrophils <- subset(Neutrophils, subset =annot_Marnix != "Erythroid" )

# Specify the order you want
Neutrophils$annot_Marnix <- factor(Neutrophils$annot_Marnix,
                                   levels = c("Neutrophil Progenitor",
                                              "Early Immature Neutrophil",
                                              "Late Immature Neutrophil",
                                              "Mature Neutrophil"))

Idents(Neutrophils) <- Neutrophils$annot_Marnix




VlnPlot(Neutrophils, features = c("TNFSF13B"))
VlnPlot(Neutrophils,
        features = c("TNFSF13B"),
        pt.size = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # optional for better x-axis label


dotplot_df <- DotPlot(Neutrophils  ,features =  c("TNFSF13B"), assay = "SCT")
dot_data <- dotplot_df$data
#dot_data$id <- factor(dot_data$id, levels = cluster_colors) #choose!
#dot_data$id_factor <- factor(dot_data$id, levels = names(cluster_colors_umap )) #choose!

#library(ggplot2)

#broad_annot_map <- meta_df %>%
# select(id = annot_merged_final_with_macro_noslash, broader = annot_broad_run3) %>%
#distinct()
#dot_data <- dot_data %>%
# left_join(broad_annot_map, by = "id")

#dot_data_sub <- dot_data %>%
# filter(broader %in% c("Endothelium", "Structural cells"))

plot <-  ggplot(dot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "lightblue", high = "firebrick") +
  #scale_color_viridis_c(option = "C")  # or "C", "B", "A"
  #facet_grid(. ~ broader, scales = "free_x", space = "free_x") +
  theme_minimal() +
  coord_flip()+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Cluster names
    axis.text.y = element_text(size = 14, face = "italic")         # Gene names
  )+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = "Gene", color = "Avg Expr (scaled)", size = "% expressed")
plot

pdf(file= file.path(output_dir,"Sup_Figure_3","dotplotBAFF_Madelon.pdf"), 
    width = 7, height = 3, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

dotplot_df <- DotPlot(
  Neutrophils,
  features = "TNFSF13B",
  split.by = "source"
)

dot_data <- dotplot_df$data

dot_data$features.plot <- factor(dot_data$features.plot,
                                 levels = "TNFSF13B")

dot_data <- dot_data %>%
  tidyr::separate(
    col = id,
    into = c("celltype", "source"),
    sep = "_",
    extra = "merge"
  )


dot_data$celltype <- factor(dot_data$celltype, levels = levels(Neutrophils$annot_Marnix))


library(ggplot2)

plot <- ggplot(dot_data, aes(x = source, y = celltype)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  #coord_flip() +
  scale_color_gradient(low = "lightblue", high = "firebrick") +
  theme_minimal(base_size = 14) +
  labs(
    x = "",
    y = "Cell type",
    color = "Avg expr (scaled)",
    size = "% expressed"
  ) +
  scale_x_discrete(limits = c("NDMM", "control"))
plot

pdf(file= file.path(output_dir,"Sup_Figure_3","dotplotBAFF_Madelon_splitbystatus.pdf"), 
    width = 7, height = 3, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()




# Figure 3H TACIpos in PCs --------------------------------------------

# Make sure the gene name matches exactly (case sensitive)
gene <- "TNFRSF13B"

Xenium_Object_Seurat$TNFRSF13B_status <- ifelse(
  Seurat::GetAssayData(Xenium_Object_Seurat, assay = "RNA", layer = "counts")[gene, ] > 0,
  "pos",
  "neg"
)

library(dplyr)

# 1) Subset to Plasma Cells
pc_cells <- subset(Xenium_Object_Seurat, subset = annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves == "Plasma Cell")

pc_cells$TNFRSF13B_status_2 <- ifelse(
  Seurat::GetAssayData(pc_cells, assay = "RNA", layer = "counts")[gene, ] > 1,
  "pos",
  "neg"
)

# 2) Calculate per patient
pc_pct <- pc_cells@meta.data %>%
  group_by(Status_simp_2) %>%
  summarise(
    n_pc = n(),
    n_pc_TNFRSF13B_pos = sum(TNFRSF13B_status == "pos"),
    pct_TNFRSF13B_pos = n_pc_TNFRSF13B_pos / n_pc * 100
  ) %>%
  arrange(Status_simp_2)

pc_pct

pc_pct_2 <- pc_cells@meta.data %>%
  group_by(ObjectName_anonymous_grouped) %>%
  summarise(
    n_pc = n(),
    n_pc_TNFRSF13B_pos = sum(TNFRSF13B_status == "pos"),
    pct_TNFRSF13B_pos = n_pc_TNFRSF13B_pos / n_pc * 100
  ) %>%
  arrange(ObjectName_anonymous_grouped)

pc_pct_2

pc_pct_by_patient_group <- pc_cells@meta.data %>%
  group_by(ObjectName_anonymous_grouped, tumor_percentage_grouped035ext_hex_average_method2) %>%
  summarise(
    n_pc = n(),
    n_pc_TNFRSF13B_pos = sum(TNFRSF13B_status == "pos"),
    pct_TNFRSF13B_pos = n_pc_TNFRSF13B_pos / n_pc * 100
  ) %>%
  arrange(ObjectName_anonymous_grouped, tumor_percentage_grouped035ext_hex_average_method2)

pc_pct_by_patient_group



pc_pct_by_patient_group_filtered <- pc_pct_by_patient_group %>%
  filter(
    tumor_percentage_grouped035ext_hex_average_method2 != "unqualified",
    ObjectName_anonymous_grouped %in% c( "PCL_1","PCL_2","NDMM_3", "NDMM_4", "NDMM_5")
  )

library(ggplot2)

ggplot(pc_pct_by_patient_group, 
       aes(x = tumor_percentage_grouped035ext_hex_average_method2, y = pct_TNFRSF13B_pos)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.8) +
  facet_wrap(~ ObjectName_anonymous_grouped, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Tumor percentage group",
    y = "% TNFRSF13B+ Plasma Cells",
    title = "TNFRSF13B positivity in plasma cells per patient and tumor group"
  )


pc_pct_by_patient_group_filtered <- pc_pct_by_patient_group %>%
  filter(
    tumor_percentage_grouped035ext_hex_average_method2 != "unqualified",
    ObjectName_anonymous_grouped %in% c( "PCL_2","NDMM_3", "NDMM_4", "NDMM_5")
  )

ggplot(pc_pct_by_patient_group_filtered,
       aes(
         x = tumor_percentage_grouped035ext_hex_average_method2,
         y = pct_TNFRSF13B_pos,
         group = ObjectName_anonymous_grouped,      # connect by patient
         color = ObjectName_anonymous_grouped  # NDMM vs PCL color
       )) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    x = "Tumor percentage group",
    y = "% TNFRSF13B+ Plasma Cells",
    color = "Disease stage",
    title = "TNFRSF13B positivity in plasma cells"
  ) 

library(dplyr)
library(ggplot2)

# Compute normalized TNFRSF13B per patient
pc_pct_by_patient_group_norm <- pc_pct_by_patient_group_filtered %>%
  # Filter unqualified regions first if needed
  filter(tumor_percentage_grouped035ext_hex_average_method2 != "unqualified") %>%
  group_by(ObjectName_anonymous_grouped) %>%
  filter(any(tumor_percentage_grouped035ext_hex_average_method2 == "normal PC percentage")) %>%
  mutate(
    # Define baseline = lowest tumor percentage (assume "scarce" / "normal PC")
    baseline_pct = pct_TNFRSF13B_pos[tumor_percentage_grouped035ext_hex_average_method2 == "normal PC percentage"],
    baseline_pct = ifelse(length(baseline_pct) == 0, 1, baseline_pct),  # avoid empty baseline
    pct_TNFRSF13B_norm = pct_TNFRSF13B_pos / baseline_pct
  ) %>%
  ungroup()

plot <- ggplot(pc_pct_by_patient_group_norm,
               aes(
                 x = tumor_percentage_grouped035ext_hex_average_method2,
                 y = pct_TNFRSF13B_norm,
                 group = ObjectName_anonymous_grouped,
                 color = ObjectName_anonymous_grouped
               )) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    x = "Tumor percentage group",
    y = "Normalized % TNFRSF13B+ PCs (relative to normal)",
    color = "Patient",
    title = "Normalized TNFRSF13B positivity in plasma cells"
  ) +
  scale_y_continuous(limits = c(0, NA))  # start at 0


plot

pdf(file= file.path(output_dir,"Sup_Figure_6","TACIpos_PCS_perarchi.pdf"), 
    width = 8, height = 6, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()



library(dplyr)

paired_df <- pc_pct_by_patient_group_filtered %>%
  filter(
    tumor_percentage_grouped035ext_hex_average_method2 %in% c("normal PC percentage", "dense")
  ) %>%
  select(
    ObjectName_anonymous_grouped,
    tumor_percentage_grouped035ext_hex_average_method2,
    pct_TNFRSF13B_pos
  ) %>%
  tidyr::pivot_wider(
    names_from = tumor_percentage_grouped035ext_hex_average_method2,
    values_from = pct_TNFRSF13B_pos
  )

t.test(
  paired_df$`normal PC percentage`,
  paired_df$dense,
  paired = TRUE
)


plot <- ggplot(
  pc_pct_by_patient_group_filtered,
  aes(
    x = tumor_percentage_grouped035ext_hex_average_method2,
    y = pct_TNFRSF13B_pos
  )
) +
  ## Boxplot (main signal)
  geom_boxplot(
    outlier.shape = NA,
    fill = "#F2F2F2",
    color = "black",
    width = 0.55,
    linewidth = 0.4
  ) +
  
  ## Paired patient trajectories (context)
  geom_line(
    aes(group = ObjectName_anonymous_grouped),
    color = "grey65",
    alpha = 0.6,
    linewidth = 0.4
  ) +
  
  ## Patient points
  geom_point(
    aes(color = ObjectName_anonymous_grouped),
    size = 2.4
  ) +
  
  ## Optional: nicer discrete colors
  scale_color_brewer(palette = "Dark2") +
  
  ## Clean theme
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(linewidth = 0.5)
  ) +
  
  labs(
    x = "Tumor density group",
    y = "% TNFRSF13B⁺ plasma cells",
    color = "Patient",
    title = "TNFRSF13B expression in plasma cells across tumor density"
  )

pdf(file= file.path(output_dir,"Sup_Figure_6","TACIpos_PCS_perarchi.pdf"), 
    width = 6, height = 6, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Figure 3 I-L Volcanoplots -----------------------------------------------




combined_DEG <- readRDS("~/DEG.rds")

combined_DEG_sub <- combined_DEG[combined_DEG$celltype == "Stroma",]

# Code to use to plot
volcano_ctcvsbm <- make_volcano_plot(
  data        = combined_DEG_sub,
  log2fc_col  = "avg_log2FC",
  padj_col    = "p_val_adj",
  gene_col    = "gene",  # use rownames as gene names, or set to your gene column name
  up_label    = "Dense",
  down_label  = "Sparse & Normal",
  other_label = "Not significant",
  padj_cut    = 0.005,
  fc_cut      = 0.025,
  color_up    = "#D82638",
  color_down  = "#397FB9", 
  color_other = "#D9D9D9",
  title       = "analysis",
  subtitle    = "specific_analysis",
  label_all_sig = FALSE,
  # label_genes   = labels_filtered,   # only label these if significant
  label_top_n   = 40,
  max_overlaps  = 20,
  # or FALSE + label_top_n = 50
  seed        = 123
)
volcano_ctcvsbm

pdf(file= file.path(dir, file), 
    width = 8, height = 6, onefile = FALSE)
print(volcano_ctcvsbm)
dev.off()






# Figure 4A Cellular Neighborhood ---------------------------------------------------------------

library(pheatmap)

for_plot <- prop.table(table(spe_Xenium$cn_celltypes_12_named, as.character(spe_Xenium$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv)), margin = 1)

niche_counts <- table(spe_Xenium$cn_celltypes_12_named)
neighborhood_percentage <- niche_counts / sum(niche_counts) * 100
# Step 2: Create a data frame for row annotation (with the percentage of the neighborhood size)
row_annotation <- data.frame(Neighborhood_Percentage = neighborhood_percentage)
# Convert row_annotation to a named vector
neighborhood_percentage_vector <- setNames(row_annotation$Neighborhood_Percentage.Freq, row_annotation$Neighborhood_Percentage.Var1)
row_annotation_df <- data.frame(Neighborhood_Percentage = neighborhood_percentage_vector)

plot <- pheatmap( for_plot,
                  color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
                  scale = "column",
                  annotation_row = row_annotation_df,  # tilt column labels
                  border_color = NA,
                  angle_col = 45  # tilt column labels
)
plot

pdf(file= file.path(output_dir,"Figure_3","CN_alltogether.pdf"), width = 12, height = 6)
print(plot)
dev.off()



# Figure 4B CN per tumor region  ---------------------------------------------------------------

meta_df <- Xenium_Object_Seurat@meta.data
#meta_df_NDMM_SMM_CBM_PCL <- meta_df[meta_df$Status_simp_2 %in% c("NDMM", "SMM", "CBM", "PCL"), ]


df_plot <- meta_df %>%
  dplyr::select(
    ObjectName = tumor_percentage_grouped035ext_hex_average_method2_CBM,
    annot_evenless_broad_run3 = cn_celltypes_12_named
  )



# Step 2: Reorder 'annot_MM_HBM_more_simple' based on your custom order
# Ensure meta_df$cluster_name is a factor ordered by your color vector


# Step 3: Calculate proportions
df_plot <- df_plot %>%
  group_by(ObjectName, annot_evenless_broad_run3) %>%
  tally() %>%  # Count occurrences of each combination
  ungroup() %>%
  group_by(ObjectName) %>%
  mutate(proportion = n / sum(n))  # Calculate proportion within each ObjectName

# Ensure ObjectName is a factor with custom order
df_plot$ObjectName_fac <- factor(df_plot$ObjectName, levels = names(region_colors))



#df_plot$ObjectName <- factor(df_plot$ObjectName, levels = architecture_order)

df_plot <- df_plot %>% filter(!(ObjectName %in% c(NA, "unqualified")))

df_plot$cluster_factor <- factor(
  df_plot$annot_evenless_broad_run3,        # replace with your actual cluster column
  levels = names(niche_colors)  # order according to the color vector
)

# join IMC data
library(dplyr)

#df_plot_combined <- bind_rows(df_plot, df_cbm_combined_renamed)
#df_plot_combined$ObjectName <- factor(df_plot_combined$ObjectName, levels = names(status_colors))

#df_plot_combined$cluster_factor2 <- factor(df_plot_combined$cluster_factor,        # replace with your actual cluster column
# levels = names(cluster_colors_with_IMC)  # order according to the color vector
#)
# Prepare data frame
#df_plot <- df_plot %>%
#  mutate(ObjectName_fac = factor(ObjectName, levels = patient_order))

plot <- ggplot(df_plot, aes(x = ObjectName_fac, y = proportion, fill = cluster_factor       )) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot
  scale_fill_manual(values = niche_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x-axis labels and increase size
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove background
        axis.title.x = element_blank(),  # Optional: Remove x-axis title
        axis.title.y = element_blank(),  # Optional: Remove y-axis title
        text = element_text(size = 15))
#coord_flip()  # Flip x and y axes  # Customize text size if needed
plot

#plot <- dittoBarPlot(seurat_HBM_MM_ownSCT , var = "annot_MM_HBM_more_simple", group.by = "ObjectName", color.panel = cluster_colors) #,x.reorder = x_reorder_indices) #, colors = cluster_colors) #, scale = "count")  

pdf(file= file.path(output_dir,"Figure_3","barplot_tumor_grouped_CNs.pdf"), 
    width = 8, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

# Figure 4C CN per status  ---------------------------------------------------------------

meta_df <- Xenium_Object_Seurat@meta.data
meta_df_NDMM_SMM_CBM_PCL <- meta_df[meta_df$Status_simp_2 %in% c("NDMM", "SMM", "CBM", "PCL", "Relapse"), ]


df_plot <- meta_df_NDMM_SMM_CBM_PCL %>%
  dplyr::select(
    ObjectName = Status_simp_2,
    annot_evenless_broad_run3 = cn_celltypes_12_named
  )



# Step 2: Reorder 'annot_MM_HBM_more_simple' based on your custom order
# Ensure meta_df$cluster_name is a factor ordered by your color vector


# Step 3: Calculate proportions
df_plot <- df_plot %>%
  group_by(ObjectName, annot_evenless_broad_run3) %>%
  tally() %>%  # Count occurrences of each combination
  ungroup() %>%
  group_by(ObjectName) %>%
  mutate(proportion = n / sum(n))  # Calculate proportion within each ObjectName

# Ensure ObjectName is a factor with custom order
df_plot$ObjectName_fac <- factor(df_plot$ObjectName, levels = status_order)



#df_plot$ObjectName <- factor(df_plot$ObjectName, levels = architecture_order)

df_plot <- df_plot %>% filter(!(ObjectName %in% c(NA, "unqualified")))

df_plot$cluster_factor <- factor(
  df_plot$annot_evenless_broad_run3,        # replace with your actual cluster column
  levels = names(niche_colors)  # order according to the color vector
)

# join IMC data
library(dplyr)

#df_plot_combined <- bind_rows(df_plot, df_cbm_combined_renamed)
#df_plot_combined$ObjectName <- factor(df_plot_combined$ObjectName, levels = names(status_colors))

#df_plot_combined$cluster_factor2 <- factor(df_plot_combined$cluster_factor,        # replace with your actual cluster column
# levels = names(cluster_colors_with_IMC)  # order according to the color vector
#)
# Prepare data frame
#df_plot <- df_plot %>%
#  mutate(ObjectName_fac = factor(ObjectName, levels = patient_order))

plot <- ggplot(df_plot, aes(x = ObjectName_fac, y = proportion, fill = cluster_factor       )) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot
  scale_fill_manual(values = niche_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x-axis labels and increase size
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove background
        axis.title.x = element_blank(),  # Optional: Remove x-axis title
        axis.title.y = element_blank(),  # Optional: Remove y-axis title
        text = element_text(size = 15))
#coord_flip()  # Flip x and y axes  # Customize text size if needed
plot

#plot <- dittoBarPlot(seurat_HBM_MM_ownSCT , var = "annot_MM_HBM_more_simple", group.by = "ObjectName", color.panel = cluster_colors) #,x.reorder = x_reorder_indices) #, colors = cluster_colors) #, scale = "count")  

pdf(file= file.path(output_dir,"Figure_3","barplot_statuses_CNs.pdf"), 
    width = 8, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Figure 4J absolute interactions ---------------------------------------------------------------

plot_data <- readRDS("~/combined.rds")


status_colors <- c(  #fir annot even less broad
  
  "CBM" = "#46F0F0",
  "IMC_CBM" = "#46F0F0",
  "SMM" = "#4575B4",              # Macrophages
  "NDMM" = "#FF7F00", 
  "IMC_NDMM" = "#46F0F0",
  "PCL" = "#D73027",
  "Relapse" = "#984EA3" #New: yellow
)


# 3️⃣ Compute average % per celltype per architecture across all disease stages
avg_dt <- plot_data %>%
  group_by(architecture, celltype) %>%
  summarise(mean_percent_all_stages = mean(percent, na.rm = TRUE), .groups = "drop")

# 4️⃣ Determine ordering per architecture
arch_ordering <- avg_dt %>%
  group_by(architecture) %>%
  arrange(desc(mean_percent_all_stages)) %>%
  summarise(celltype_order = list(celltype)) %>%
  tibble::deframe()  # named list: architecture -> ordered celltypes

# 5️⃣ Plot PDF
pdf(file= file.path(output_dir,"Figure_4/PlasmaCell_Interactions_Ordered_by_Avg_CBMNDMM_Stages_broad.pdf"), width = 9, height = 7)

architectures <- unique(plot_data$architecture)
plot_data <- plot_data[plot_data$DiseaseStage %in% c("CBM", "NDMM"),]

for (arch in architectures) {
  data_arch <- plot_data %>% filter(architecture == arch)
  
  # Order celltypes by average percent across all stages
  celltype_levels <- arch_ordering[[arch]]
  data_arch <- data_arch %>%
    mutate(celltype = factor(celltype, levels = celltype_levels))
  
  p <- ggplot(data_arch, aes(x = celltype, y = percent)) +
    geom_boxplot(aes(fill = DiseaseStage), outlier.shape = NA) +
    #geom_jitter(aes(color = DiseaseStage), width = 0.2, size = 2) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Plasma Cell Interactions -", arch),
         x = "Cell Type",
         y = "Percentage of Plasma Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  library(ggplot2)
  
  p <- ggplot(data_arch, aes(x = celltype, y = percent, fill = DiseaseStage)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    # Optional: add jittered points
    # geom_jitter(aes(color = DiseaseStage), width = 0.2, size = 2, alpha = 0.7) +
    scale_fill_manual(values = status_colors) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Plasma Cell Interactions -", arch),
      x = "Cell Type",
      y = "Percentage of Plasma Cells",
      fill = "Disease Stage"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12)
    )
  
  print(p)
  
}

dev.off()


#Figure 4J Absolute interactions architectures in 1 plot per diseasestage ----------



pdf("Figures/Sup_Figure_6/PlasmaCell_Interactions_Celltypes_x_Architectures_grouped_perDiseaseStagex.pdf",
    width = 10, height = 7)

arch_order <- c("normal PC percentage", "sparse", "dispersed", "dense")
region_colors <- c("CBM" = "cyan", "unqualified" = "#4CAF50", "normal PC percentage" =  "#FFE066", #, "normal PC percentage" =  "#FFFFB3"
                   "sparse" ="#377EB8","dispersed" ="#B080FF", "dense"= "#D72638") #, "dense_otheroption_toopink"= "#FB8072")  # Replace with your preferred order

# Compute global ordering of celltypes (most → least abundant) from plot_data
global_celltype_order <- plot_data %>%
  group_by(celltype) %>%
  summarise(global_mean = mean(percent, na.rm = TRUE)) %>%
  arrange(desc(global_mean)) %>%
  pull(celltype)

# Ensure consistent ordering in the plot
plot_data2 <- plot_data %>%
  mutate(
    celltype = factor(celltype, levels = global_celltype_order),
    architecture = factor(architecture, levels = c("normal PC percentage", "sparse", "dispersed", "dense"))
  )
# Loop over disease stages → one page per stage
# Loop over disease stages → one page per stage
for (stage in unique(plot_data2$DiseaseStage)) {
  
  df_stage <- plot_data2 %>% filter(DiseaseStage == stage)
  
  p <- ggplot(df_stage, aes(
    x = celltype,
    y = percent,
    fill = architecture
  )) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      outlier.shape = NA,
      linewidth = 0.3
    ) +
    scale_fill_manual(values = region_colors) +
    labs(
      title = paste("Plasma Cell Interactions –", stage),
      x = "Cell Type",
      y = "Percentage of Plasma Cells",
      fill = "Architecture"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 13),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 12, face = "bold")
    )
  
  print(p)
}

dev.off()

## statistics

plot_data2_NDMM <- plot_data2[plot_data2$DiseaseStage == "NDMM",]

library(dplyr)
library(broom)

normal_level <- "normal PC percentage"
dense_level  <- "dense"

plot_data2_NDMM$architecture <- factor(
  plot_data2_NDMM$architecture,
  levels = c(normal_level, dense_level)
)

tt_short <- plot_data2_NDMM %>%
  group_by(celltype) %>%
  do({
    tt <- t.test(percent ~ architecture, data = ., var.equal = TRUE)
    tibble(
      mean_normal = mean(.$percent[.$architecture == normal_level]),
      mean_dense  = mean(.$percent[.$architecture == dense_level]),
      p_value     = tt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),   # <-- Benjamini–Hochberg
    p_signif = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001  ~ "***",
      p_adj < 0.01   ~ "**",
      p_adj < 0.05   ~ "*",
      TRUE           ~ "ns"
    ))


tt_short


# Figure 4G 3 spatial combined ---------------------------------------------------------------

## Interaction Analysis
plot_data <- readRDS("~/combined_per_patient.rds")

out_all_no_duplicates <- plot_data[!(plot_data$group_by %in% c( "02390_2_NDMM", "12928_2_PCL", "52998_2_NDMM","17034_2_NDMM","14719_3_CBM", "8667_3_CBM",  
                                                                "9414_3_CBM", "17336_3_CBM")), ]
out_all_noduplicates_noCBM <- out_all_no_duplicates[!(out_all_no_duplicates$group_by %in% c( "14719_2_CBM",  "14719_3_CBM", "8667_2_CBM","8667_3_CBM","9141_2_CBM","9414_2_CBM", "9414_3_CBM" ,  "17336_3_CBM", "17336_2_CBM", "14719_2_CBM"  ,   "14719_3_CBM")), ]

out_all_tib <- as_tibble(out_all_no_duplicates)

# Step 2 — summarize interactions
df <- out_all_tib %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE), .groups = "drop")


# Filter for plasma cell interactions
df_pc <- df %>%
  filter(from_label == "Plasma Cell")

# Arrange by sum_sigval (descending)
df_pc <- df_pc %>%
  arrange(sum_sigval)

# Reorder interaction labels so bars are sorted
df_pc <- df_pc %>%
  mutate(interaction = paste(from_label, "←→", to_label)) %>%
  mutate(interaction = factor(interaction, levels = interaction))

df_pc <- df_pc %>%
  mutate(other_cell = str_replace(interaction, "Plasma Cell ←→ ", ""))


df_pc_filt <- df_pc[
  !df_pc$other_cell %in% c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell"),
]

df_pc_filt <- df_pc_filt %>%
  mutate(
    norm_sigval = case_when(
      sum_sigval > 0 ~ sum_sigval / max(sum_sigval),
      sum_sigval < 0 ~ sum_sigval / abs(min(sum_sigval)),
      TRUE ~ 0
    )
  )

library(forcats)

df_pc_filt <- df_pc_filt %>%
  # rank by sum_sigval descending, tie-break alphabetically
  mutate(
    tie_rank = dense_rank(desc(sum_sigval)),         # rank 1 = largest sum_sigval
    alpha_order = row_number()                        # will re-order ties alphabetically
  ) %>%
  arrange(tie_rank, other_cell)   # first by rank, then by alphabetical order

# Plot
plot <- ggplot(df_pc_filt, aes(x = reorder(other_cell, norm_sigval), y = norm_sigval, fill = norm_sigval)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#3B3E98", mid = "gray95", high = "#842525", midpoint = 0,
    name = "Sum sigval"
  ) +
  labs(
    x = "Interacting cell type",
    y = "Sum of significance value",
    title = "Plasma Cell Interactions"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
plot
pdf(file= file.path(output_dir,"Figure_4","Interactions_all_withPCranked_filt.pdf"), 
    width = 8, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

interaction_df <- df_pc[!(df_pc$to_label %in% c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell")), ]
interaction_df <- interaction_df %>%
  arrange(desc(sum_sigval)) %>%
  mutate(order_rank = row_number())  # new column with ranking


### neighborhood

for_plot <- prop.table(table(spe_Xenium$cn_celltypes_12_named, as.character(spe_Xenium$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv)), margin = 1)

niche_counts <- table(spe_Xenium$cn_celltypes_12_named)
neighborhood_percentage <- niche_counts / sum(niche_counts) * 100
# Step 2: Create a data frame for row annotation (with the percentage of the neighborhood size)
row_annotation <- data.frame(Neighborhood_Percentage = neighborhood_percentage)
# Convert row_annotation to a named vector
neighborhood_percentage_vector <- setNames(row_annotation$Neighborhood_Percentage.Freq, row_annotation$Neighborhood_Percentage.Var1)
row_annotation_df <- data.frame(Neighborhood_Percentage = neighborhood_percentage_vector)


# Convert for_plot to numeric matrix if not already
for_plot_matrix <- as.matrix(for_plot)
storage.mode(for_plot_matrix) <- "numeric"

# Scale columns (z-score per cell type)
for_plot_scaled <- scale(for_plot_matrix)  # scale by column (default)
for_plot_df_neighborhood <- as.data.frame(for_plot_scaled)
# Rename columns for clarity
colnames(for_plot_df_neighborhood) <- c("Neighborhood", "Celltype", "Enrichment_Score")

# Optionally filter to a single neighborhood
selected_neighborhood1 <- "Tumor-Associated"  # change this as needed, or NULL for all
selected_neighborhood2 <- "Tumor"  # change this as needed, or NULL for all

for_plot_df_neighborhood_filt <- for_plot_df_neighborhood
if (!is.null(selected_neighborhood1)) {
  for_plot_df_neighborhood_filt <- for_plot_df_neighborhood %>% filter(Neighborhood == selected_neighborhood1)
}

# Order cell types by decreasing proportion (within the filtered neighborhood)
for_plot_df_neighborhood_filt <- for_plot_df_neighborhood_filt %>%
  arrange(Enrichment_Score) %>%
  mutate(Celltype = factor(Celltype, levels = unique(Celltype)))

# Plot
library(ggplot2)
library(dplyr)

# Define neighborhoods to combine
selected_neighborhoods <- c("Tumor-Associated", "Tumor")  # change as needed

# Filter and sum enrichment scores per cell type
for_plot_df_neighborhood_combined <- for_plot_df_neighborhood %>%
  filter(Neighborhood %in% selected_neighborhoods) %>%
  group_by(Celltype) %>%
  summarise(Enrichment_Score = sum(Enrichment_Score, na.rm = TRUE), .groups = "drop") %>%
  arrange(Enrichment_Score) %>%
  mutate(Celltype = factor(Celltype, levels = unique(Celltype)))



# Reorder Celltype by Enrichment_Score for better visualization
for_plot_df_neighborhood_combined <- for_plot_df_neighborhood_combined %>%
  mutate(Celltype = reorder(Celltype, Enrichment_Score))


for_plot_df_neighborhood_PC2_filt <- for_plot_df_neighborhood_filt[
  !for_plot_df_neighborhood_filt$Celltype %in% c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell"),
]

plot <- ggplot(
  for_plot_df_neighborhood_PC2_filt,
  aes(x = Celltype, y = Enrichment_Score, fill = Enrichment_Score)
) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradientn(
    colors = c(
      "#2A2D7C",   # deep red  (−2)
      "#F4F4F4",   # off-white (0)
      "#8A181A"    # deep blue (+2)
    ),
    values = scales::rescale(c(-2, 0, 2)),
    limits = c(-2, 2),     # 🔴 FIXED scale
    oob = scales::squish,  # clip values outside ±2
    name = "Enrichment"
  ) +
  labs(
    x = "Cell type",
    y = "Enrichment score",
    title = ifelse(
      is.null(selected_neighborhoods),
      "All neighborhoods",
      selected_neighborhoods
    )
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 17),
    axis.text.x = element_text(size = 17),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 17, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

plot


pdf(file= file.path(output_dir,"Figure_4","Interactions_tumor_associated_niche_ranked_filt3.pdf"), 
    width = 9, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# dense vs sparse all patients without CBM, run 2F
# Perform one-sided Wilcoxon tests (dense > normal)

library(openxlsx)
library(dplyr)
library(purrr)

# 1. Get all sheet names
sheets <- getSheetNames("p_values_LMER_test2.xlsx")

# 2. Read all sheets + add sheet name column
dense_vs_normal_archiLMER <- map_dfr(sheets, function(sh) {
  read.xlsx("p_values_LMER_test2.xlsx", sheet = sh) %>%
    mutate(comparison = sh)
})

dense_vs_normal_archiLMER <- dense_vs_normal_archiLMER[dense_vs_normal_archiLMER$comparison =="dense_vs_normal_PC_percentage", ]
dense_vs_normal_archiLMER <- dense_vs_normal_archiLMER[dense_vs_normal_archiLMER$tumor_percentage_grouped035ext_hex_average_method2  =="dense", ]

significance_df_plot <- dense_vs_normal_archiLMER
significance_df_plot$logp_value <- -log10(significance_df_plot$p_val)

significance_df_plot2 <- significance_df_plot  %>% 
  distinct(clusters, comparison, .keep_all = TRUE)       # keep first row per  

# Arrange by sum_sigval (descending)
# Reorder factor levels by logp_value
significance_df_plot2 <- significance_df_plot2 %>%
  arrange(logp_value) %>%
  mutate(cluster_celltype = factor(clusters, levels = unique(clusters)))


library(dplyr)
library(ggplot2)

# Prepare the data
significance_df_plot2 <- significance_df_plot2 %>%
  mutate(logp_value = -log10(p_val)) %>%
  arrange(desc(ratio)) %>%  # order by log fold change ascending
  mutate(cluster_celltype = factor(clusters, levels = rev(clusters)))  # largest on top

# Plot logFC bars, colored by -log10(p-value)
plot <- ggplot(significance_df_plot2, aes(x = cluster_celltype, y = ratio, fill = logp_value)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#2166ac",    # dark blue
    mid = "gray95",     # soft gray
    high = "#b2182b",   # dark red
    midpoint = median(significance_df_plot2$ratio),
    name = "-log10(p-value)"
  ) +
  labs(
    x = "Cell type",
    y = "Log2 fold change (Dense / Normal PC)",
    title = "Log fold change of cell types in Dense vs Normal PC regions"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

plot

pdf(file= file.path(output_dir,"Figure_4","Dense_vs_normal_ranked_LMER_filt.pdf"), 
    width = 8, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


significance_df_plot2 <- significance_df_plot2 %>%
  mutate(
    logp_value = -log10(p_val),
    score = ratio * logp_value  # combined effect size and significance
  )

# Prepare the data
significance_df_plot2 <- significance_df_plot2 %>%
  arrange(desc(score)) %>%  # order by log fold change ascending
  mutate(cluster_celltype = factor(clusters, levels = rev(clusters)))  # largest on top


significance_df_plot2 <- significance_df_plot2 %>%
  mutate(
    sig_label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

plot <- ggplot(significance_df_plot2, aes(x = cluster_celltype, y = score, fill = score)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sig_label), hjust = ifelse(significance_df_plot2$score > 0, -0.1, 1.1), size = 5) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "gray95",
    high = "#b2182b",
    midpoint = 0,
    name = "Combined score"
  ) +
  labs(
    x = "Cell type",
    y = "Combined score (logFC × -log10(p-value))",
    title = "Cell type enrichment in Dense vs Normal PC regions"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

plot

pdf(file= file.path(output_dir,"Figure_4","Dense_vs_normal_ranked_LMER_filt_combined_score.pdf"), 
    width = 8, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


signed_competition_rank_with_zero <- function(x) {
  out <- numeric(length(x))
  
  # Zeros stay zero
  out[x == 0] <- 0
  
  # Positive values
  pos <- x > 0
  if (any(pos)) {
    # rank positives ascending so largest positive = largest rank
    out[pos] <- rank(x[pos], ties.method = "min")
  }
  
  # Negative values
  neg <- x < 0
  if (any(neg)) {
    # rank negatives descending so most negative = most negative rank
    out[neg] <- -rank(-x[neg], ties.method = "min")
  }
  
  return(out)
}



x <- c(5, 3, 4,4,4, 0, -3, -2)
signed_dense_rank(x)


for_plot_df_neighborhood_ranked <- for_plot_df_neighborhood_filt[!(for_plot_df_neighborhood_filt$Celltype %in% c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell")), ]
for_plot_df_neighborhood_ranked1 <- for_plot_df_neighborhood_ranked %>%
  arrange(desc(Enrichment_Score)) %>%
  mutate(order_rank = row_number())  # new column with ranking
for_plot_df_neighborhood_ranked2 <- for_plot_df_neighborhood_ranked %>%
  mutate(rank1 = signed_competition_rank_with_zero(Enrichment_Score))

interaction_df <- df_pc[!(df_pc$to_label %in% c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell")), ]
interaction_df1 <- interaction_df %>%
  arrange(desc(sum_sigval)) %>%
  mutate(order_rank = row_number())  # new column with ranking

interaction_df2 <- interaction_df %>%
  mutate(rank1 = signed_competition_rank_with_zero(sum_sigval))


significance_df_plot_densevsnormal <- significance_df_plot2
significance_df_plot_densevsnormal1 <- significance_df_plot_densevsnormal %>%
  arrange(desc(ratio)) %>%
  mutate(order_rank = row_number())  # new column with ranking
significance_df_plot_densevsnormal2 <- significance_df_plot_densevsnormal %>%
  mutate(rank1 = signed_competition_rank_with_zero(score))


# Example: df1, df2, df3 with same rownames
# Each has a column "rank"

combined1 <- interaction_df %>%
  select(celltype = to_label, rank1 = order_rank) %>%
  inner_join(for_plot_df_neighborhood_ranked %>% 
               select(celltype = Celltype, rank2 = order_rank),
             by = "celltype") %>%
  inner_join(significance_df_plot_densevsnormal %>%
               select(celltype = cluster_celltype, rank3 = order_rank),
             by = "celltype") %>%
  mutate(rank_sum = rank1 + rank2 + rank3) %>%
  arrange(rank_sum)


combined2 <- interaction_df2 %>%
  select(celltype = to_label, rank1 = rank1) %>%
  inner_join(for_plot_df_neighborhood_ranked2 %>% 
               select(celltype = Celltype, rank2 = rank1),
             by = "celltype") %>%
  inner_join(significance_df_plot_densevsnormal2 %>%
               select(celltype = cluster_celltype, rank3 = rank1),
             by = "celltype") %>%
  mutate(rank_sum = rank1 + rank2 + rank3) %>%
  arrange(rank_sum)


# Plot stacked barplot
combined2 <- combined2 %>%
  arrange(desc(rank_sum)) %>%
  mutate(celltype = factor(celltype, levels = celltype))

ggplot(combined2, aes(x = celltype, y = rank_sum, fill = rank_sum)) +
  geom_col() +
  coord_flip() +  
  scale_fill_gradient2(low ="darkred", mid = "gray", high =  "darkblue", midpoint = 40) +
  labs(x = "Interaction", y = "Sum sigval", 
       title = "Plasma Cell Interactions (ranked by 3 spatial methods)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Invert rank_sum
combined2 <- combined2 %>%
  mutate(rank_sum_inverted = max(rank_sum) - rank_sum + min(rank_sum))  # ensures positive values

library(ggplot2)
library(dplyr)

# Prepare data
combined2 <- combined2 %>%
  arrange(rank_sum_inverted) %>%  # largest on top
  mutate(celltype = factor(celltype, levels = rev(unique(celltype))))

# Plot
plot <- ggplot(combined2, aes(x = celltype, y = rank_sum, fill = rank_sum)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#2166ac",    # dark blue
    mid = "gray95",     # soft gray
    high = "#b2182b",   # dark red
    midpoint = 0,
    name = "Inverted Sum sigval"
  ) +
  labs(
    x = "Cell type",
    y = "Sum ranks",
    title = "Plasma Cell Interactions (ranked by 3 spatial methods, inverted)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

plot

pdf(file= file.path(output_dir,"Figure_4","3_spatial_methods_ranked_withLMER_posneg2.pdf"), 
    width = 8, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

# Plot
library(viridis)
plot <- ggplot(combined2, aes(x = celltype, y = rank_sum, fill = rank_sum)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradient2(
    low = "purple",    # teal
    mid = "#f7f7f7",    # gray/white
    high = "#1b9e77",   # coral
    midpoint = 0,
    name = "Inverted Sum sigval"
  )+
  labs(
    x = "Cell type",
    y = "Sum ranks",
    title = "Plasma Cell Interactions (ranked by 3 spatial methods, inverted)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

plot
pdf(file= file.path(output_dir,"Figure_4","3_spatial_methods_ranked_withLMER_posneg2_newColour.pdf"), 
    width = 8, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


#### ranked dotplot

library(dplyr)

prep_rank_df <- function(
    df,
    cell_col,
    score_col,
    pval_col = NULL,
    method_name,
    rank_by = c("pval", "abs_score")
) {
  
  rank_by <- match.arg(rank_by)
  
  out <- df %>%
    dplyr::rename(
      Celltype = {{ cell_col }},
      score    = {{ score_col }}
    ) %>%
    mutate(
      score = as.numeric(score),
      Method = method_name
    )
  
  if (!is.null(pval_col) && rank_by == "pval") {
    out <- out %>%
      dplyr::rename(pval = {{ pval_col }}) %>%
      mutate(
        pval = as.numeric(pval),
        rank = rank(pval, ties.method = "average"),
        rank_pct = 1 - (rank - 1) / (max(rank) - 1)
      )
  } else {
    out <- out %>%
      mutate(
        pval = NA_real_,
        rank = rank(-(score), ties.method = "average"),
        rank_pct = 1 - (rank - 1) / (max(rank) - 1)
      )
  }
  
  out %>%
    select(Celltype, Method, score, pval, rank, rank_pct)
}

df_dense_ranked <- prep_rank_df(
  significance_df_plot2,
  cell_col  = clusters,
  score_col = score ,
  method_name = "Dense vs Sparse"
)


df_neigh_ranked <- prep_rank_df(
  for_plot_df_neighborhood_ranked,
  cell_col  = Celltype,
  score_col = Enrichment_Score,
  method_name = "PC Neighborhood"
)

df_interact_ranked <- prep_rank_df(
  interaction_df ,
  cell_col  = to_label,
  score_col = sum_sigval,
  method_name = "PC Interaction"
)

plot_df <- bind_rows(
  df_dense_ranked,
  df_neigh_ranked#,
  #df_interact_ranked
)

cell_order_tom <- plot_df %>%
  filter(Method == "Dense vs Sparse")%>%
  group_by(Celltype) %>%
  summarise(score) %>%
  arrange(desc(score)) %>%
  pull(Celltype)

cell_order <- plot_df %>%
  group_by(Celltype) %>%
  summarise(mean_rank = mean(rank_pct, na.rm = TRUE)) %>%
  arrange(desc(mean_rank)) %>%
  pull(Celltype)


plot_df$Celltype <- factor(plot_df$Celltype, levels = rev(cell_order_tom))

plot_df <- plot_df %>%
  group_by(Method) %>%
  mutate(score_scaled = score / max(abs(score), na.rm = TRUE)) %>%
  ungroup()

plot_df$Method  <- factor(plot_df$Method, levels = c("Dense vs Sparse", "PC Neighborhood", "PC Interaction"))


plot <- ggplot(
  plot_df,
  aes(
    x = Method,
    y = Celltype,
    size = rank_pct,
    fill = score_scaled
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  scale_size_continuous(
    range = c(1.5, 6),
    name = "Rank percentile"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Relative effect"
  ) +
  theme_minimal(base_size = 17) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


plot
pdf(file= file.path(output_dir,"Figure_4","2_spatial_methods_ranked_dotplot_otherorder.pdf"), 
    width = 8, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

# Get list of cell types positive in all methods
positive_all_methods <- plot_df %>%
  group_by(Celltype) %>%
  summarize(all_positive = all(score > 0), .groups = "drop") %>%
  filter(all_positive) %>%
  pull(Celltype)

positive_all_methods


# Figure 4I interactions combined architecture and disease stage normalized-----------------------------------------


out_all <- readRDS("~/interactions.rds")
out_all_no_duplicates <- out_all[!(out_all$group_by %in% c( "02390_2_NDMM", "12928_2_PCL", "52998_2_NDMM","17034_2_NDMM","14719_3_CBM", "8667_3_CBM",  
                                                            "9414_3_CBM", "17336_3_CBM")), ]
rename_vector <- c(
  "7459_3_PCL_MM" = "7459_3_PCL",
  "04000_3_PCL_MM" = "04000_3_PCL"
)

# Apply renaming to the group_by column
out_all_no_duplicates <- as_tibble(out_all_no_duplicates) %>%
  mutate(group_by = ifelse(group_by %in% names(rename_vector),
                           rename_vector[group_by],
                           group_by))

out_all_no_duplicates <- as_tibble(out_all_no_duplicates) %>%
  mutate(status = str_extract(group_by, "(?<=_)[^_]+$"))



out_all_CBM_no_duplicates  <- out_all_no_duplicates[out_all_no_duplicates$group_by %in% c( "14719_2_CBM",   "8667_2_CBM","9141_2_CBM","9414_2_CBM",  "17336_3_CBM", "14719_2_CBM"), ]

out_all_per_archi_stage <- readRDS("~/interactions_combined_wo_duplicated.rds")
out_all_per_archi_stage_no_duplicates <- out_all_per_archi_stage[!(out_all_per_archi_stage$group_by %in% c( "02390_2_NDMM", "12928_2_PCL", "52998_2_NDMM","17034_2_NDMM","14719_3_CBM", "8667_3_CBM",  
                                                                                                            "9414_3_CBM", "17336_3_CBM")), ]
out_all_per_archi_stage_no_duplicates_no_CBM <- out_all_per_archi_stage_no_duplicates[!(out_all_per_archi_stage_no_duplicates$group_by %in%  c("14719_2_CBM", "14719_3_CBM", "8667_2_CBM", "8667_3_CBM", 
                                                                                                                                               "9141_2_CBM", "9414_2_CBM", "9414_3_CBM", 
                                                                                                                                               "9414_3_CBM", "17336_3_CBM")), ]



# Split by architecture
split_architecture_wo_CBM <- split(out_all_per_archi_stage_no_duplicates_no_CBM, out_all_per_archi_stage_no_duplicates_no_CBM$architecture)
split_architecture_with_CBM <- split(out_all_per_archi_stage_no_duplicates, out_all_per_archi_stage_no_duplicates$architecture)
split_disease_stage <- split(out_all_no_duplicates, out_all_no_duplicates$status)
combined <- c(split_disease_stage, split_architecture_wo_CBM)
combined$all <- out_all_no_duplicates
#combined$all_without_CBM <- out_all_no_duplicates_no_CBM

combined$unqualified <- NULL
# Extract the plasma cell interactions for each architecture
# Step 2 — summarize interactions


# Step 1: Summarize interactions per dataset
plasma_interactions <- lapply(combined, function(df) {
  df %>%
    group_by(group_by, from_label, to_label) %>%   # group_by = patient/sample
    summarize(sum_sigval = sum(sigval, na.rm = TRUE), .groups = "drop") %>%
    filter(from_label == "Plasma Cell")
})

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Assume you have a list of dataframes, one per interaction analysis
# each with columns: from_label, to_label, sigval
# Example: pc_list <- list(analysis1 = df1, analysis2 = df2, analysis3 = df3)


# or
# Convert to tibble
out_all_tib <- as_tibble(combined$all)
df_pc <- out_all_tib %>% filter(from_label == "Plasma Cell")

# Summarize interactions per 'to_label'
df_pc <- df_pc %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE), .groups = "drop")

# Keep only Plasma Cell interactions

# Compute average interaction per 'to_label'
avg_sig <- df_pc %>%
  group_by(to_label) %>%
  summarize(avg_sigval = mean(sum_sigval, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(avg_sigval))   # highest → lowest
# Extract ordering
col_order <- avg_sig$to_label
col_order <- rev(avg_sig$to_label)


# 2. Add a 'dataset' column to each dataframe and combine
combined_df <- bind_rows(
  lapply(names(plasma_interactions), function(name) {
    df <- plasma_interactions[[name]]
    df$dataset <- name
    df
  })
)

# Step 3: Compute number of patients per dataset for normalization
samples_per_dataset <- combined_df %>%
  group_by(dataset) %>%
  summarise(n_samples = n_distinct(group_by), .groups = "drop")

# Step 4: Average interaction per patient
combined_df <- combined_df %>%
  group_by(dataset, from_label, to_label) %>%
  summarise(sum_sigval = sum(sum_sigval), .groups = "drop") %>%
  left_join(samples_per_dataset, by = "dataset") %>%
  mutate(avg_sigval = sum_sigval / n_samples)

# 4. Reorder factors
combined_df <- combined_df %>%
  mutate(
    to_label = factor(to_label, levels = col_order),          # x-axis: plasma cell
    from_label = factor(from_label, levels = unique(from_label)),  # y-axis stacked
    dataset = factor(dataset, levels = names(plasma_interactions))      # for grouping in y-axis
  )

# Compute average sum_sigval per dataset (y-axis)
dataset_order <- combined_df %>%
  group_by(dataset) %>%
  summarise(avg_sig_all = mean(avg_sigval, na.rm = TRUE)) %>%
  arrange(desc(avg_sig_all)) %>%   # highest avg first
  pull(dataset)

dataset_order <- c( "CBM","normal PC percentage", "sparse",  "dispersed"  , "SMM",  "Relapse","NDMM", "dense","PCL", "all_without_CBM", "all" )

# Reorder dataset factor for y-axis stacking
combined_df <- combined_df %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

# Remove Plasma Cell from y-axis labels (since all are from Plasma Cell)
combined_df <- combined_df %>%
  mutate(from_label = "")  # blank out y-axis labels


dataset_order <- c( "all_without_CBM", "all", "CBM" , "SMM",  "Relapse","NDMM","PCL","normal PC percentage", "sparse",  "dispersed",  "dense" )

# Reorder dataset factor for y-axis stacking
combined_df <- combined_df %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

combined_df <- combined_df %>%
  filter(dataset != "Relapse")

# Plot
plot <- ggplot(combined_df, aes(x = from_label, y = to_label, fill = avg_sigval)) +
  geom_tile(color = "white") +  # optional: white borders between tiles
  facet_grid(. ~ dataset, scales = "free_x", space = "free_x", switch = "x") +  # datasets on top
  scale_fill_gradient2(low = muted("blue"), mid = "#f7f7f7", high = muted("red"), midpoint = 0) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_blank(),  # hide x-axis text (all from Plasma Cell)
    axis.text.y = element_text(size = 14, face = "plain"),  # show cell types
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "plain"),
    strip.text.x.top = element_text(size = 14, face = "plain", angle = 45),
    strip.text.x = element_text(size = 14, face = "plain", angle = 45),  # dataset labels on top
    panel.grid = element_blank(),
    legend.title = element_text(size = 14, face = "plain"),
    legend.text = element_text(size = 12)
  ) +
  labs(
    y = "Other Cell Types",  # y-axis shows cell types
    fill = "Interaction Score",
    title = "Interactions of Plasma Cells with Other Cell Types"
  )

plot

pdf(file= file.path(output_dir,"Figure_4","interaction_combined_regionsZeroiswhite_combinedvertical.pdf"), 
    width = 13, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()




combined_df_archi <- combined_df %>%
  filter(dataset %in% c("normal PC percentage", "sparse", "dispersed", "dense", "all"))

combined_df_archi_no_CBM_with_CBM <- combined_df %>%
  filter(dataset %in% c("CBM", "normal PC percentage", "sparse", "dispersed", "dense", "all"))


combined_df_state <- combined_df %>%
  filter(dataset %in% c("CBM", "SMM", "Relapse", "NDMM", "PCL", "all"))


# Plot
plot <- ggplot(combined_df_archi_no_CBM_with_CBM, aes(x = to_label, y = from_label, fill = avg_sigval)) +
  geom_tile(color = "white") +  # optional: white borders between tiles
  #coord_flip() +# datasets on left
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y", switch = "y") +  
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_blank(),  # hide y-axis labels
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    strip.text.y.left = element_text(size = 14, face = "bold", angle = 0),  # dataset labels on left
    panel.grid = element_blank(),  # remove internal grid lines
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  
  labs(
    x = "Other Cell Types",  # y-axis removed; x-axis shows the other cells
    fill = "Interaction Score",
    title = "Interactions of Plasma Cells with Other Cell Types"
  )

plot

pdf(file= file.path(output_dir,"Figure_3","interaction_combined_regionsZeroiswhite.pdf"), 
    width = 13, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

combined_df_state <- combined_df_state[combined_df_state$dataset != "Relapse",]

# Plot
plot <- ggplot(combined_df_state, aes(x = to_label, y = from_label, fill = avg_sigval)) +
  geom_tile(color = "white") +  # optional: white borders between tiles
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y", switch = "y") +  # datasets on left
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_blank(),  # hide y-axis labels
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    strip.text.y.left = element_text(size = 14, face = "bold", angle = 0),  # dataset labels on left
    panel.grid = element_blank(),  # remove internal grid lines
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  labs(
    x = "Other Cell Types",  # y-axis removed; x-axis shows the other cells
    fill = "Interaction Score",
    title = "Interactions of Plasma Cells with Other Cell Types"
  )

plot
pdf(file= file.path(output_dir,"Figure_3","interaction_combined_disease_stage_zeroiswhite.pdf"), 
    width = 13, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()



# Supplementary -----------------------------------------------------------


# Sup Figure 1A dotplot broad ---------------------------------------------------------------------

cluster_colors_umap <- c(  #fir annot even less broad
  "Other" = "#999999", 
  "Erythroid" = "#E78AC3",  
  "Endothelium" = "#81C784", 
  "Stroma" = "#3CB44B",                   # Stromal / Myofibroblast
  "Osteolineage" = "#808000",   
  "Adipocytes" = "#A65628", 
  "MKC lineage" = "#D4A5A5", 
  "VSMC" = "#FFD92F" ,# New: brown
  "Interferon stimulated cells" = "#984EA3",  
  "Neutrophil Progenitors" = "#D1E5F0",
  "Neutrophils" = "#46F0F0",
  "Mature neutrophils" = "#91C3E6", # Myelocytes
  "Basophils Eosinophils Mast cells" = "#FB8072",  # Ba Eo Ma
  "Macrophages" = "#4575B4",              # Macrophages
  "Monocytes" = "#6699CC",     
  "Dendritic cells" = "#F781BF",          # New: bright pink
  "T NK cells" = "#FF7F00", 
  "B/Plasma cells" = "#D73027" ,
  "Adipocyte"                      = "#A65628",  # brown
  "MKC Lineage"                    = "#D4A5A5",  # pale red
  "Neutrophil Progenitor"          = "#D1E5F0",  # light blue
  "Neutrophil"                     = "#46F0F0",
  "Mature Neutrophil"              = "#46F0F0",  #same colour for UMAP
  #"Mature Neutrophil"              = "#91C3E6",  # steel blue
  "Basophil/Eosinophil/Mast Cell"  = "#FB8072",  # salmon
  "Macrophage"                     = "#4575B4",  # dark blue
  "Monocyte"                       = "#6699CC",  # medium blue
  "Dendritic Cell"                 = "#F781BF",  # magenta
  "T/NK Cell"                      = "#FF7F00",  # orange
  "B/Plasma Cell"                  = "#D73027"   # red# New: yellow
)


ordered_genes_broad <- c(
  
  "ALAS2", "GYPA", "BTNL9","ENG",  "VWF", "CD34",
  "SPARC", "PCOLCE", "LEPR", "COL3A1","THY1",  "SPP1", "BGLAP", "COL1A1", "MYH11", "ACTA2",
  "PLIN4", "ADIPOQ", "THBS1", "MMRN1", "SELP","CTSG", "ELANE", "MPO", "LTF",
  "ITGAM", "MMP9", "LYZ", "OLFM4",
  "FCGR3B", "MME", "SELL", "FCGR2A","CLC", "GATA2", "MS4A2", "CPA3", "KIT", "IL1RL1",
  "CD163", "C1QC", "CD5L", "FCGR3A",  "CD14","VCAN", "VSIG4", "CLEC10A", "HLA-DRA", "CD74", "IRF8", "CD8A", "TRAC", "CCL5",  "SSR4", "SDC1", "PRDM1","SLAMF7")

Idents(Xenium_Object_Seurat) <- Xenium_Object_Seurat@meta.data$annot_broad_run3_correct_ifndiv
dotplot_df <- DotPlot(Xenium_Object_Seurat  ,features = ordered_genes_broad, assay = "SCT")
dot_data <- dotplot_df$data
#dot_data$id <- factor(dot_data$id, levels = cluster_colors) #choose!
dot_data$id_factor <- factor(dot_data$id, levels = names(cluster_colors_umap )) #choose!

library(ggplot2)

#broad_annot_map <- meta_df %>%
# select(id = annot_merged_final_with_macro_noslash, broader = annot_broad_run3) %>%
#distinct()
#dot_data <- dot_data %>%
# left_join(broad_annot_map, by = "id")

#dot_data_sub <- dot_data %>%
# filter(broader %in% c("Endothelium", "Structural cells"))

plot <-  ggplot(dot_data, aes(x = id_factor, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "lightblue", high = "firebrick") +
  #scale_color_viridis_c(option = "C")  # or "C", "B", "A"
  #facet_grid(. ~ broader, scales = "free_x", space = "free_x") +
  theme_minimal() +
  coord_flip()+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Cluster names
    axis.text.y = element_text(size = 14, face = "italic")         # Gene names
  )+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = "Gene", color = "Avg Expr (scaled)", size = "% expressed")


pdf(file= file.path(output_dir,"Sup_Figure_1","dotplotBroad_more_room.pdf"), 
    width = 16, height = 5, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

# Sup Figure 1B QC ---------------------------------------------------------------------

ordered_sample_groups
status_colors <- c(  #fir annot even less broad
  
  "CBM" = "#46F0F0",
  "SMM" = "#4575B4",              # Macrophages
  "NDMM" = "#FF7F00", 
  "PCL" = "#D73027",
  "Relapse" = "#984EA3" #New: yellow
)

library(Seurat)
library(ggplot2)
library(dplyr)

# Ensure the ordering
Xenium_Object_Seurat@meta.data$ObjectName_anonymous <- factor(
  Xenium_Object_Seurat@meta.data$ObjectName_anonymous,
  levels = ordered_sample_groups
)

# Build a data frame for plotting
df <- Xenium_Object_Seurat@meta.data %>%
  group_by(ObjectName_anonymous, Status_simp_2) %>%
  summarise(count = n(), .groups = "drop")

# Plot
ggplot(df, aes(x = ObjectName_anonymous, y = count, fill = Status_simp_2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = status_colors) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Sample group", y = "Cell count", fill = "Status")

# Example: violin plot for a gene
plot <- VlnPlot(
  Xenium_Object_Seurat,
  features = "nFeature_RNA",  # change this to your feature of interest
  group.by = "ObjectName_anonymous",
  split.by = "Status_simp_2",
  cols = status_colors,
  pt.size = 0
) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "right"
  )

plot
pdf(file= file.path(output_dir,"Sup_Figure_1","QC_nFeature.pdf"), 
    width = 7, height = 4, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Sup Figure 1C patient overview---------------------------------------------------------------

library(ggplot2)
library(dplyr)

# Example: assuming your Seurat object is 'Xenium_Object_Seurat'
meta_df <- Xenium_Object_Seurat@meta.data %>%
  distinct(ObjectName, Status_simp_2, Tissue_ext)
# Get patient-level disease stage (assuming each patient has one stage)
# If multiple cells per patient, we collapse to unique per patient
patient_stage <- meta_df %>%
  group_by(Tissue_ext) %>%    # replace with your patient ID column name
  summarise(Disease_stage = unique(Status_simp_2)) %>%
  ungroup()
# Count patients per stage
stage_counts <- patient_stage %>%
  count(Disease_stage) %>%
  mutate(Percentage = n / sum(n) * 100)
meta_df <- Xenium_Object_Seurat@meta.data %>%
  distinct(ObjectName, Status_simp_2, Tissue_ext)
# Mark patients with duplicates
patient_counts <- meta_df %>%
  group_by(ObjectName, Status_simp_2) %>%
  summarise(n_samples = n(), .groups = "drop")
counts <- meta_df %>%
  count(Status_simp_2)
# Plot
# Define the order
patient_counts$Status_simp_2 <- factor(patient_counts$Status_simp_2, levels = status_order)

plot <- ggplot(patient_counts, aes(x = Status_simp_2, fill = n_samples > 1)) +
  geom_bar(color = "gray") +
  scale_fill_manual(values = c("FALSE" = "#D1E5F0", "TRUE" ="#35978F"),
                    labels = c("Single Run", "Double Run")) +
  labs(x = "Disease stage", y = "Number of patients", fill = "Duplicate runs",
       title = "Cohort composition by patient") +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot

pdf(file= file.path(output_dir,"Patient_Overview.pdf"), 
    width = 7, height = 10, useDingbats = FALSE, onefile = FALSE)

print(plot)

dev.off()


# Sup Figure 1D umap split per patient ---------------------------------------------------------------------


Idents(Xenium_Object_Seurat) <- Xenium_Object_Seurat@meta.data$annot_even_less_broad_correct_ifndiv
DimPlot(Xenium_Object_Seurat)


df_umap <- as.data.frame(Embeddings(Xenium_Object_Seurat, "harmony_umap")) %>%
  mutate(
    patient = Xenium_Object_Seurat$ObjectName_anonymous_grouped,   # adjust column name
    celltype = Xenium_Object_Seurat$annot_broad_run3_correct_ifndiv     # optional coloring
  )

df_umap$patient <- factor(df_umap$patient, levels = ordered_sample_groups)

# Install ggrastr if needed
# install.packages("ggrastr")
library(ggplot2)
library(ggrastr)

p <- ggplot(df_umap, aes(x = harmonyumap_1, y = harmonyumap_2)) +
  geom_point_rast(aes(color = celltype), size = 0.15, alpha = 0.4, shape = 16) +
  scale_color_manual(values = cluster_colors)+   # define consistent colors
  facet_wrap(~ patient, ncol = 5) +   # adjust ncol to fit journal layout
  coord_equal() +
  theme_void() +
  theme(
    strip.text = element_text(size = 6),      # small patient labels
    strip.background = element_blank()
  )

p
pdf(file= file.path(output_dir,"Sup_Figure_1","UMAP_broad_splitted_rasterized2.pdf"), 
    width = 15, height = 15, useDingbats = FALSE, onefile = FALSE)

print(p)

dev.off()

# Sup Figure 1E dotplot all myeloid -----------------------------------------------------

myeloid_all <- mcreadRDS("~/myeloid_all.rds" )


ordered_precise_levels_adjusted_corrected <- c(
  "Interleukin Producing Cell",
  "Interferon Stimulated Cell",
  "Erythroid",
  "EC",
  "SEC",
  "AEC",
  "Osteoblast",
  "Osteo-Fibroblastic MSC",
  "Fibro MSC",
  "LEPR MSC",
  "THY1 MSC",
  "APOD MSC",
  "CXCL14 MSC",
  "Adipocyte",
  "VSMC",
  "MKC Lineage",
  "Basophil Eosinophil",
  "Mast Cell",
  "Neutrophil Progenitor",
  "CEACAM6 Immature Neutrophil 1",
  "OLFM4 Immature Neutrophil 2",
  "LTF Immature Neutrophil 3",
  "MMP9 Immature Neutrophil 4",
  "Mature Neutrophil",
  "Classical Monocyte",
  "Non Classical Monocyte",
  "Macrophage",
  "Conventional DC",
  "Plasmacytoid DC",
  "CD8 T Cell",
  "Activated or Exhausted T Cell",
  "NK or Cytotoxic T Cell",
  "Naive or CM CD4 T Cell",
  "Regulatory T Cell",
  "B or Plasma Cell",
  "Cycling B or Plasma Cell",
  "Plasma Cell"
)


ordered_genes_myeloid_only <- c(
  "IL27", "IL12B", "CLC","PGLYRP1", "CPA3", "KIT", "IL1RL1", "CTSG", "ELANE", "MPO","LYZ","CEACAM6",
  "OLFM4","LTF", "CHIT1", "PADI4", "MMP9",  
  "FCGR3B","MME", "FCGR2A", "CD14","VCAN", "FCGR3A", "CXCR1", "CD163", "C1QC", "CD5L", "HLA-DRA", "CD74", "CD44","IRF8","GZMB",
  "CLEC4C")

Idents(myeloid_all) <- myeloid_all@meta.data$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv
dotplot_df <- DotPlot(myeloid_all  ,features = ordered_genes_myeloid_only, assay = "SCT")
dot_data <- dotplot_df$data
#dot_data$id <- factor(dot_data$id, levels = cluster_colors) #choose!
dot_data$id_factor <- factor(dot_data$id, levels = ordered_precise_levels_adjusted_corrected) #choose!

library(ggplot2)


ggplot(dot_data, aes(x = id_factor, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "lightblue", high = "firebrick") +
  #scale_color_viridis_c(option = "C")  # or "C", "B", "A"
  #facet_grid(. ~ broader, scales = "free_x", space = "free_x") +
  theme_minimal() +
  coord_flip()+
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Cluster names
    axis.text.y = element_text(size = 14, face = "italic")         # Gene names
  )+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = "Gene", color = "Avg Expr (scaled)", size = "% expressed")

plot <- ggplot(dot_data, aes(x = id_factor, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "lightblue", high = "firebrick") +
  theme_minimal() +
  coord_flip()+
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Cluster names
    axis.text.y = element_text(size = 14, face = "italic")         # Gene names
  )+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = "Gene", color = "Avg Expr (scaled)", size = "% expressed")



pdf(file= file.path(output_dir,"Sup_Figure_1","myeloid_all.pdf"), 
    width = 15, height = 5, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Sup Figure 1F myeloid cells  ---------------------------------------------------------------------

myeloid <- qread("~/myeloid.qs" )


# Assign Idents first
Idents(myeloid) <- myeloid$annot_myeloid_only_just_for_paper
# Rename Idents with singular names and consistent formatting
myeloid <- RenameIdents(object = myeloid,
                                  "MKC_lineage" = "MKC Lineage",
                                  "FCGR3B mature neutrophils" = "FCGR3B Mature Neutrophil",
                                  "SEC" = "SEC",
                                  "Erythroid" = "Erythroid",
                                  "Neutrophil Progenitors" = "Neutrophil Progenitor",
                                  "LEPR MSC" = "LEPR MSC",
                                  "LTF_neutrophils" = "LTF Neutrophil",
                                  "MMP9_neutrophils" = "MMP9 Neutrophil",
                                  "VSMC" = "VSMC",
                                  "Conventional dendritic cells" = "Conventional DC",
                                  "EC" = "EC",
                                  "CD8+ T cells" = "CD8 T Cell",
                                  "B cells or Plasma cells" = "B or Plasma Cell",
                                  "Interferon stimulated cells" = "Interferon Stimulated Cell",
                                  "Interleukin producing cells" = "Interleukin Producing Cell",
                                  "Basophil Eosinophil" = "Basophil Eosinophil",
                                  "Adipocytes" = "Adipocyte",
                                  "Osteoclast-like" = "Osteoblast",
                                  "Plasma cells" = "Plasma Cell",
                                  "Cycling B or Plasma cells" = "Cycling B or Plasma Cell",
                                  "THY1 MSC" = "THY1 MSC",
                                  "Macrophages" = "Macrophage",
                                  "Naive or Central Memory T cells" = "Naive or CM CD4 T Cell",
                                  "CEACAM_neutrophils" = "CEACAM Neutrophil",
                                  "NK cells or Cytotoxic T cells" = "NK or Cytotoxic T Cell",
                                  "OLFM4_neutrophils" = "OLFM4 Neutrophil",
                                  "Non classical monocytes" = "Non Classical Monocyte",
                                  "CXCL14 MSC" = "CXCL14 MSC",
                                  "Stroma" = "Fibro MSC",
                                  "Classical Monocytes" = "Classical Monocyte",
                                  "Plasmacytoid DCs" = "Plasmacytoid DC",
                                  "AEC" = "AEC",
                                  "APOD MSC" = "APOD MSC",
                                  "COL1A1 Osteolineage " = "Osteo-Fibroblastic MSC",
                                  "OLFM4_neutrophil" = "OLFM4 Neutrophil",
                                  "Activated T cells or Exhausted T cells" = "Activated or Exhausted T Cell",
                                  "Regulatory T cells" = "Regulatory T Cell",
                                  "Mast cells" = "Mast Cell"
)
myeloid$annot_myeloid_only_just_for_paper_correct  <- Idents(myeloid)

Idents(myeloid) <- myeloid$annot_myeloid_only_just_for_paper_correct
myeloid <- RenameIdents(object = myeloid,
                                  "MKC lineage" = "MKC Lineage",
                                  "FCGR3B mature neutrophil" = "Mature Neutrophil",
                                  "Neutrophil progenitor" = "Neutrophil Progenitor",
                                  "LTF neutrophil" = "LTF Immature Neutrophil 3",
                                  "MMP9 neutrophil" = "MMP9 Immature Neutrophil 4",
                                  "Interleukin producing cell" = "Interleukin Producing Cell",
                                  "Basophil Eosinophil" = "Basophil Eosinophil",
                                  "Macrophage" = "Macrophage",
                                  "CEACAM neutrophil" = "CEACAM6 Immature Neutrophil 1",
                                  "OLFM4 neutrophil" = "OLFM4 Immature Neutrophil 2",
                                  "Non classical monocyte" = "Non Classical Monocyte",
                                  "Classical monocyte" = "Classical Monocyte",
                                  "Plasmacytoid DC" = "Plasmacytoid DC",
                                  "Mast cell" = "Mast Cell"
)
myeloid$annot_myeloid_only_just_for_paper_correct_correct  <- Idents(myeloid)

Idents(myeloid) <- myeloid$annot_myeloid_only_just_for_paper_correct
myeloid <- RenameIdents(object = myeloid,
                                  "MKC lineage" = "MKC Lineage",
                                  "FCGR3B mature neutrophil" = "Mature Neutrophil",
                                  "Neutrophil progenitor" = "Neutrophil Progenitor",
                                  "LTF neutrophil" = "LTF Immature Neutrophil 3",
                                  "MMP9 neutrophil" = "MMP9 Immature Neutrophil 4",
                                  "Interleukin producing cell" = "Interleukin Producing Cell",
                                  "Basophil Eosinophil" = "Basophil Eosinophil",
                                  "Macrophage" = "Mononuclear Phagocyte",
                                  "CEACAM neutrophil" = "CEACAM6 Immature Neutrophil 1",
                                  "OLFM4 neutrophil" = "OLFM4 Immature Neutrophil 2",
                                  "Non classical monocyte" = "Non Classical Monocyte",
                                  "Classical monocyte" = "Mononuclear Phagocyte",
                                  "Plasmacytoid DC" = "Mononuclear Phagocyte",
                                  "Mast cell" = "Mast Cell"
)
myeloid$annot_myeloid_only_just_for_paper_correct_correct_mononuclear_together  <- Idents(myeloid)


myeloid_colors_paper_1 <- c(
  "Neutrophil Progenitor" = "#E6F2FF",             # very pale blue
  "CEACAM6 Immature Neutrophil 1" = "#B3D9FF",     # light pastel blue
  "Classical Monocyte" = "#99CCFF",                # medium pastel blue
  "LTF Immature Neutrophil 3" = "#80BFFF",         # soft blue
  "Non Classical Monocyte" = "#66FFCC",           # tealish pastel, more unique
  "Basophil Eosinophil" = "#FFCCE5",              # soft pink
  "MMP9 Immature Neutrophil 4" = "#3399FF",        # medium blue
  "Mature Neutrophil" = "#0073CC",                # darker pastel blue
  "Macrophage" = "#336699",                        # medium-dark pastel blue, softer than before
  "Plasmacytoid DC" = "#CC99FF", 
  "Conventional DC" = "#FF99CC",                         # soft pink
  # pastel purple
  "Interleukin Producing Cell" = "#FFB266" ,      # lighter purple
  "Mast Cell" = "#FF99CC",                         # soft pink
  "OLFM4 Immature Neutrophil 2" = "#66CCFF" ,
  "Mononuclear Phagocyte" = "#336699"# soft cyan-blue
)



myeloid   @reductions$harmony_umap_myeloid_only
Idents(myeloid       ) <- myeloid@meta.data$annot_myeloid_only_just_for_paper_correct_correct_mononuclear_together
n_cells <- ncol(myeloid)
plot <- DimPlot(myeloid, label = T, reduction = "harmony_umap_myeloid_only",  cols = myeloid_colors_paper_1) + NoLegend()+
  ggtitle(paste0("Mononuclear Phagocytes (n = ", format(n_cells, big.mark = ","), ")"))

plot 
pdf(file= file.path(output_dir,"Sup_Figure_1","Umap_myeloid.pdf"), 
    width = 13, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()



# Sup Figure 1G mononuclear_phagocytes ------------------------------------------------

mononuclear_phagocytes <- qread("~//mononuclear_phagocytes.qs" )

# Assign Idents first
Idents(mononuclear_phagocytes) <- mononuclear_phagocytes$annot_mononuclear_phagocytes_simple
# Rename Idents with singular names and consistent formatting
mononuclear_phagocytes <- RenameIdents(object = mononuclear_phagocytes,
                                                "Macrophages" = "Macrophage",
                                                "Classical Monocytes" = "Classical Monocyte",
                                                "Conventional dendritic cells" = "Conventional DC",
                                                "Plasmacytoid DCs" = "Plasmacytoid DC"
)
mononuclear_phagocytes$annot_mononuclear_phagocytes_simple_correct  <- Idents(mononuclear_phagocytes)


myeloid_colors_paper <- c(
  # light pastel blue
  "Classical Monocyte" = "#FFB266" ,                 # medium pastel blue
  "LTF Immature Neutrophil 3" = "#80BFFF",         # soft blue
  "Non Classical Monocyte" = "#66FFCC",           # tealish pastel, more unique
  "Basophil Eosinophil" = "#FFCCE5",              # soft pink
  "MMP9 Immature Neutrophil 4" = "#3399FF",        # medium blue
  "Mature Neutrophil" = "#0073CC",                # darker pastel blue
  "Macrophage" = "#99CCFF",                 # medium-dark pastel blue, softer than before
  "Plasmacytoid DC" = "#CC99FF",                  # pastel purple
  "Interleukin Producing Cell" = "#FFB266" ,      # lighter purple
  "Conventional DC" = "#FF99CC",                         # soft pink
  "OLFM4 Immature Neutrophil 2" = "#66CCFF" ,
  "Mononuclear Phagocyte" = "#336699"# soft cyan-blue
)

mononuclear_phagocytes   @reductions$harmony_umap_mononuclear_phagocytes
Idents(mononuclear_phagocytes       ) <- mononuclear_phagocytes@meta.data$annot_mononuclear_phagocytes_simple_correct
n_cells <- ncol(mononuclear_phagocytes)
plot <- DimPlot(mononuclear_phagocytes, label = T, reduction = "harmony_umap_mononuclear_phagocytes",  cols = myeloid_colors_paper) + NoLegend()+
  ggtitle(paste0("Mononuclear Phagocytes (n = ", format(n_cells, big.mark = ","), ")"))

plot 
pdf(file= file.path(output_dir,"Sup_Figure_1","Umap_mononuclear_phago.pdf"), 
    width = 13, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Sup Figure 1H percentage Macrophage scRNA vs Xenium -----------------------------------------------------------


#run this and skip to end
combined_percent_true_mac <- readRDS("combined_percent_true_mac.rds")

# or run entire pipeline
total_myeloid_madelon_2023 <- mcreadRDS("~/myeloid_scRNAseq.rds" )


# Copy metadata
meta_df <- total_myeloid_madelon_2023@meta.data

# Remove neutrophils
meta_no_neut <- meta_df[!meta_df$annotation_combined_simple %in% c("Mature neutrophils", "ImmNeu", "PreNeu", "Proliferating myeloblasts", " Myeloblasts"), ]

# Calculate percentage of True_macrophages per patient
percent_true_mac_scRNAseq <- meta_no_neut %>%
  group_by(source) %>%
  summarise(
    total_myeloid = n(),
    true_mac = sum(annotation_combined_simple == "True_macrophages"),
    percent_true_mac = true_mac / total_myeloid * 100
  ) %>%
  ungroup()

# Calculate percentage of True_macrophages per patient
percent_true_mac_scRNAseq <- meta_no_neut %>%
  group_by(patient, source) %>%
  summarise(
    total_myeloid = n(),
    true_mac = sum(annotation_combined_simple == "True_macrophages"),
    percent_true_mac = true_mac / total_myeloid * 100
  ) %>%
  ungroup()


# Plot
ggplot(percent_true_mac_scRNAseq, aes(x = source, y = percent_true_mac)) +
  geom_col(fill = "#FB8072") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    x = "Patient",
    y = "Percentage of True Macrophages (%)",
    title = "Fraction of True Macrophages in Myeloid Cells in aspirates after Ficol"
  )


# Copy metadata
meta_df <- Xenium_Object_Seurat@meta.data

# Remove neutrophils
meta_no_neut <- meta_df[meta_df$annot_broad_run3 %in% c("Myeloid"), ]


# Calculate percentage of True_macrophages per patient
percent_true_mac_xenium <- meta_no_neut %>%
  group_by(Status_simp_2, ObjectName_anonymous_grouped) %>%
  summarise(
    total_myeloid = n(),
    true_mac = sum(annot_merged_final_with_macro_noslash == "Macrophages"),
    percent_true_mac = true_mac / total_myeloid * 100
  ) %>%
  ungroup()

# Plot
ggplot(percent_true_mac_xenium, aes(x = Status_simp_2, y = percent_true_mac)) +
  geom_col(fill = "#FB8072") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    x = "Patient",
    y = "Percentage of True Macrophages (%)",
    title = "Fraction of True Macrophages in Myeloid Cells in aspirates after Ficol"
  )

# Add method column to each dataset
percent_true_mac_xenium2 <- percent_true_mac_xenium %>%
  mutate(method = "Xenium",
         source = as.character(Status_simp_2))  # optional, unify column name

percent_true_mac_scRNAseq2 <- percent_true_mac_scRNAseq %>%
  mutate(
    method = "scRNA-seq",
    Status_simp_2 = ifelse(source == "control", "CBM", as.character(source)),
    ObjectName_anonymous_grouped = patient)
# Combine datasets
combined_percent_true_mac <- bind_rows(percent_true_mac_xenium2 %>% select(Status_simp_2, total_myeloid, true_mac, percent_true_mac, method, ObjectName_anonymous_grouped),
                                       percent_true_mac_scRNAseq2 %>% select(Status_simp_2, total_myeloid, true_mac, percent_true_mac, method,ObjectName_anonymous_grouped))


saveRDS(combined_percent_true_mac, "combined_percent_true_mac.rds")

df_filt <- combined_percent_true_mac %>% 
  filter(!(Status_simp_2 %in% c("PCL", "Relapse", "SMM")))

# Publication-ready plot
plot <- ggplot(df_filt, aes(x = Status_simp_2, y = percent_true_mac, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +  # black border for clarity
  scale_fill_manual(values = c(
    "scRNA-seq" = "#A6CEE3",  # very soft blue
    "Xenium"   = "#fb8072"   # soft orange
  ))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "True Macrophages (% of Myeloid, no Neutrophils)",
    title = "Comparison of True Macrophages between scRNA-seq and Xenium",
    fill = "Method"
  )

# Compute summary stats per condition
summary_df <- df_filt %>%
  group_by(Status_simp_2, method) %>%
  summarise(
    mean_percent = mean(percent_true_mac),
    sd_percent = sd(percent_true_mac),
    .groups = "drop"
  )

# Compute p-values for CBM and NDMM
pvals <- df_filt %>%
  group_by(Status_simp_2) %>%
  summarise(
    p_value = wilcox.test(
      percent_true_mac[method == "scRNA-seq"],
      percent_true_mac[method == "Xenium"]
    )$p.value,
    .groups = "drop"
  )

# For labeling significance
pvals <- pvals %>%
  mutate(stars = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

# Plot
plot <- ggplot(df_filt, aes(x = Status_simp_2, y = percent_true_mac, fill = method)) +
  
  # Mean bars
  stat_summary(fun = mean, geom = "col",
               position = position_dodge(width = 0.8),
               width = 0.7, color = "black") +
  
  # Error bars
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar",
               position = position_dodge(width = 0.8),
               width = 0.2, color = "black") +
  
  # Individual patient dots
  geom_point(aes(color = method),
             position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.8),
             size = 2.5, alpha = 0.9) +
  
  # p-value stars
  geom_text(data = pvals,
            aes(x = Status_simp_2,
                y = max(df_filt$percent_true_mac) + 5,
                label = stars),
            inherit.aes = FALSE,
            size = 5) +
  
  scale_fill_manual(values = c("scRNA-seq" = "#A6CEE3", "Xenium" = "#FB8072")) +
  scale_color_manual(values = c("scRNA-seq" = "#1F78B4", "Xenium" = "#E31A1C")) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 12),
    panel.grid = element_blank()
  ) +
  
  labs(
    x = "Condition",
    y = "True Macrophages (% of Myeloid, no Neutrophils)",
    title = "Comparison of True Macrophages between scRNA-seq and Xenium",
    fill = "Method",
    color = "Method"
  )

plot

pdf(file= file.path(output_dir,"Sup_Figure_1/Macrophage_perc_compare_with_stat_with_dots.pdf"), 
    width = 7, height = 6, useDingbats = FALSE, onefile = FALSE)

print(plot)

dev.off()


# Sup Figure 1I UMAP T cells  ---------------------------------------------------------------------

Tcells <- readRDS("~/Tcells.rds")



tcell_colors_warm <- c(
  "CD8 T Cell" = "#FF9999",                     # soft orange
  "NK or Cytotoxic T Cell" = "#FF704D",         # coral-orange
  "Activated or Exhausted T Cell" = "#FFA500", # soft red
  "Naive or CM CD4 T Cell" = "#FF6666",         # medium red
  "Regulatory T Cell" =  "#FFB266"              # light pink-red
)



Tcells   @reductions$harmony_umap_Tcells
Idents(Tcells       ) <- Tcells@meta.data$annot_tcell_limited_correct
n_cells <- ncol(Tcells)
plot <- DimPlot(Tcells, label = T, reduction = "harmony_umap_Tcells",  cols = tcell_colors_warm) + NoLegend()+
  ggtitle(paste0("T cells (n = ", format(n_cells, big.mark = ","), ")"))

pdf(file= file.path(output_dir,"Sup_Figure_1","UMAP_Tcells.pdf"), 
    width = 6, height = 5, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

# Sup Figure 1J dotplot T cells  ---------------------------------------------------------------------


ordered_precise_levels_adjusted_corrected <- c(
  "Interleukin Producing Cell",
  "Interferon Stimulated Cell",
  "Erythroid",
  "EC",
  "SEC",
  "AEC",
  "Osteoblast",
  "Osteo-Fibroblastic MSC",
  "Fibro MSC",
  "LEPR MSC",
  "THY1 MSC",
  "APOD MSC",
  "CXCL14 MSC",
  "Adipocyte",
  "VSMC",
  "MKC Lineage",
  "Basophil Eosinophil",
  "Mast Cell",
  "Neutrophil Progenitor",
  "CEACAM6 Immature Neutrophil 1",
  "OLFM4 Immature Neutrophil 2",
  "LTF Immature Neutrophil 3",
  "MMP9 Immature Neutrophil 4",
  "Mature Neutrophil",
  "Classical Monocyte",
  "Non Classical Monocyte",
  "Macrophage",
  "Conventional DC",
  "Plasmacytoid DC",
  "CD8 T Cell",
  "Activated or Exhausted T Cell",
  "NK or Cytotoxic T Cell",
  "Naive or CM CD4 T Cell",
  "Regulatory T Cell",
  "B or Plasma Cell",
  "Cycling B or Plasma Cell",
  "Plasma Cell"
)


ordered_genes_T_cells <- c(
  "CCL4", "CCL5", "CD8A", "CCL3", "CCL3L1", "RGS1","TOX","CD69", "PRF1", "KLRD1", "IL2RB",
  "KLRB1", "CD247",  "IL7R", "NELL2", "TCF7", "CD4", "CCR4","CTLA4",
  "FOXP3")

Idents(Tcells) <- Tcells@meta.data$annot_tcell_limited_correct
dotplot_df <- DotPlot(Tcells  ,features = ordered_genes_T_cells, assay = "SCT")
dot_data <- dotplot_df$data
#dot_data$id <- factor(dot_data$id, levels = cluster_colors) #choose!
dot_data$id_factor <- factor(dot_data$id, levels = ordered_precise_levels_adjusted_corrected) #choose!

library(ggplot2)

plot <- ggplot(dot_data, aes(x = id_factor, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "lightblue", high = "firebrick") +
  theme_minimal() +
  coord_flip()+
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Cluster names
    axis.text.y = element_text(size = 14, face = "italic")         # Gene names
  )+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = "Gene", color = "Avg Expr (scaled)", size = "% expressed")

plot

pdf(file= file.path(output_dir,"Sup_Figure_1","dotplot_Tcells.pdf"), 
    width = 10, height = 3.5, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Sup Figure 1K dotplot all ------------------------------------------------





ordered_precise_levels_adjusted_corrected <- c(
  "Interleukin Producing Cell",
  "Interferon Stimulated Cell",
  "Erythroid",
  "EC",
  "SEC",
  "AEC",
  "Osteoblast",
  "Osteo-Fibroblastic MSC",
  "Fibro MSC",
  "LEPR MSC",
  "THY1 MSC",
  "APOD MSC",
  "CXCL14 MSC",
  "Adipocyte",
  "VSMC",
  "MKC Lineage",
  "Basophil Eosinophil",
  "Mast Cell",
  "Neutrophil Progenitor",
  "CEACAM6 Immature Neutrophil 1",
  "OLFM4 Immature Neutrophil 2",
  "LTF Immature Neutrophil 3",
  "MMP9 Immature Neutrophil 4",
  "Mature Neutrophil",
  "Classical Monocyte",
  "Non Classical Monocyte",
  "Macrophage",
  "Conventional DC",
  "Plasmacytoid DC",
  "CD8 T Cell",
  "Activated or Exhausted T Cell",
  "NK or Cytotoxic T Cell",
  "Naive or CM CD4 T Cell",
  "Regulatory T Cell",
  "B or Plasma Cell",
  "Cycling B or Plasma Cell",
  "Plasma Cell"
)


ordered_genes_adjusted <- c(
  "IL27", "IL12B", "IFIT3", "MX1", "ISG15", "ALAS2", "GYPA",
  "BTNL9", "DNASE1L3", "ENG",  "VWF", "PLVAP", "SEMA3G", "SPP1", "BGLAP", "COL1A1",
  "TNC","COMP", "LEPR","NGFR", "THY1", "NOTCH3","STEAP4",  "APOD", "MGP", "CXCL14",
  "PLIN4", "ADIPOQ", "LPL", "TIMP4", "MYH11", "ACTA2",
  "MMRN1", "SELP",
  "CLC","PGLYRP1", "CPA3", "KIT", "IL1RL1", "CTSG", "ELANE", "MPO","LYZ","CEACAM6",
  "OLFM4","LTF", "CHIT1", "PADI4", "MMP9",  
  "FCGR3B","MME", "FCGR2A","CXCR2", "CD14","VCAN", "TNFSF13", "FCGR3A", "CXCR1", "CD163", "C1QC", "CD5L", "HLA-DRA", "CD74", "CD44","IRF8","GZMB",
  "CLEC4C",
  "CCL4", "CCL5", "CD8A", "CCL3", "CCL3L1", "RGS1","TOX","CD69", "PRF1", "KLRD1", "IL2RB",
  "KLRB1", "CD247",  "IL7R", "NELL2", "TCF7", "CD4", "CCR4","CTLA4",
  "FOXP3",
  "MS4A1", "CD19", "STMN1", "MKI67", "SSR4", "SDC1", "PRDM1","SLAMF7", "IRF4"
)

Idents(Xenium_Object_Seurat) <- Xenium_Object_Seurat@meta.data$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv
dotplot_df <- DotPlot(Xenium_Object_Seurat  ,features = ordered_genes_adjusted, assay = "SCT")
dot_data <- dotplot_df$data
#dot_data$id <- factor(dot_data$id, levels = cluster_colors) #choose!
dot_data$id_factor <- factor(dot_data$id, levels = ordered_precise_levels_adjusted_corrected) #choose!

library(ggplot2)



plot <- ggplot(dot_data, aes(x = id_factor, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "lightblue", high = "firebrick") +
  theme_minimal() +
  coord_flip()+
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  # Cluster names
    axis.text.y = element_text(size = 16, face = "italic")         # Gene names
  )+
  labs(x = "Cluster", y = "Gene", color = "Avg Expr (scaled)", size = "% expressed")

plot

pdf(file= file.path(output_dir,"Sup_Figure_1","dotplot_all_withAPRILandCXCR2.pdf"), 
    width = 30, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Sup Figure 2E Hexagonal grid on image ---------------------------------------------------------------


meta_df <- Xenium_Object_Seurat@meta.data
meta_df_02390 <- meta_df[meta_df$Tissue_ext == "NDMM_02390_3",]

# --- Step 1: Create the fixed bounding box for the study area ---
bbox <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 200000, ymax = 200000), crs = NA))

# --- Step 2: Define offsets for 6 different grid starting positions ---
offsets <- list(
  c(0, 0),    # no shift
  c(20, 0),   # shift 20 units right
  c(0, 20),   # shift 20 units up
  c(20, 20),  # shift 20 units right and 20 units up
  c(40, 0),   # shift 40 units right
  c(0, 40)    # shift 40 units up
)



library(ggplot2)
library(sf)
library(dplyr)

# --- Step 0: select grid offset ---
i <- 1
dx <- offsets[[i]][1]
dy <- offsets[[i]][2]

# --- Step 1: convert metadata coordinates to sf ---
centroids_df <- meta_df_02390[, c("x_corrected_together", "y_corrected_together")]
cell_points <- st_as_sf(centroids_df, coords = c("x_corrected_together", "y_corrected_together"), crs = NA)

# --- Step 2: get tissue bounding box ---
tissue_bbox <- st_as_sfc(st_bbox(cell_points))

# --- Step 3: create the hex grid and shift by offset ---
hex_grid <- st_make_grid(bbox, cellsize = 125, square = FALSE, what = "polygons")
hex_grid_shifted <- st_geometry(hex_grid) + c(dx, dy)
hex_grid_sf <- st_sf(geometry = hex_grid_shifted)
hex_grid_sf$region_name <- paste0("region_", seq_len(nrow(hex_grid_sf)))

# --- Step 4: crop hex grid to tissue bounding box ---
hex_grid_cropped <- st_intersection(hex_grid_sf, tissue_bbox)

# --- Base layer: cell density or tissue outline ---
p <- ggplot() +
  geom_sf(data = cell_points, color = "#1f78b4", size = 0.2, alpha = 0.5) +  # cells
  geom_sf(data = hex_grid_cropped, fill = NA, color = "black", size = 0.4, alpha = 1) +  # grid
  coord_sf() +  # use coord_sf() for sf objects
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  ggtitle(paste("Hexagonal grid overlay (offset", i, ")"))

print(p)

pdf(file= file.path(output_dir,"Sup_Figure_3","Grid_overlay.pdf"), width = 15, height = 6)
print(p)
dev.off()

library(ggplot2)
library(ggrastr)
library(sf)


# Extract coordinates and assign names
cell_coords_df <- st_coordinates(cell_points) %>% 
  as.data.frame()
colnames(cell_coords_df) <- c("x", "y")  # assign names manually

# Rasterized plot
p_raster <- ggplot() +
  geom_point_rast(
    data = cell_coords_df,
    aes(x = x, y = y),
    color = "#1f78b4",
    size = 0.2,
    alpha = 0.5
  ) +
  geom_sf(data = hex_grid_cropped, fill = NA, color = "black", size = 0.4, alpha = 1) +
  coord_sf() +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  ggtitle(paste("Hexagonal grid overlay (offset", i, ")"))

# Save to PDF
pdf(file = file.path(output_dir, "Sup_Figure_3", "Grid_overlay_rasterized.pdf"),
    width = 15, height = 6)
print(p_raster)
dev.off()



## legend


library(ggplot2)
library(tibble)
library(grid)
celltypes <- c(
  "Other",
  "Erythroid",
  "Endothelium",
  "Stroma",
  "Osteolineage",
  "Adipocytes",
  "MKC lineage",
  "VSMC",
  "Neutrophil",
  "Basophil/Eosinophil/Mast Cell",
  "Macrophage/Monocyte",
  "Dendritic Cell",
  "T/NK Cell",
  "B/Plasma Cell"
)

colors <- c(
  "#999999",   # Other
  "#E78AC3",   # Erythroid
  "#81C784",   # Endothelium
  "#3CB44B",   # Stroma
  "#808000",   # Osteolineage
  "#A65628",   # Adipocytes
  "#D4A5A5",   # MKC lineage
  "#FFD92F",   # VSMC
  "#46F0F0",   # Neutrophil
  "#FB8072",   # Basophil/Eosinophil/Mast Cell
  "#4575B4",   # Macrophage/Monocyte
  "#F781BF",   # Dendritic Cell
  "#FF7F00",   # T/NK Cell
  "#D73027"    # B/Plasma Cell
)
legend_df <- tibble(celltype = celltypes, color = colors)
legend_df$celltype_fac <- factor(legend_df$celltype, levels = celltypes)

# --- Create plot to extract legend ---
p <- ggplot(legend_df, aes(x = 1, y = celltype_fac, color = celltype_fac)) +
  geom_point(size = 5) +
  scale_color_manual(values = setNames(colors, celltypes)) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "lines"),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# --- Extract legend only ---
legend_only <- cowplot::get_legend(p)

# --- Save as PDF ---
pdf(file= file.path(output_dir,"Sup_Figure_3/legend_02390_broad.pdf"), width = 3, height = 4, useDingbats = FALSE, onefile = FALSE)
grid.newpage()
grid.draw(legend_only)
dev.off()

# --- 2️⃣ Horizontal legend (side by side) ---
p <- p +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 5), direction = "horizontal"))

legend_horizontal <- cowplot::get_legend(p)

pdf(file = file.path(output_dir, "Figure_2/legend_HandE_tumorgrouped2A_horizontal.pdf"),
    width = 6, height = 2, useDingbats = FALSE, onefile = FALSE)
grid.newpage()
grid.draw(legend_horizontal)
dev.off()



# Sup Figure 3A  only normal Vertical box plot comparing stages  ----------------------------------------------------------------

# 1. subset meta
meta_df <- Xenium_Object_Seurat@meta.data

meta_df_no_tumor <- meta_df %>%
  filter(!annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves %in%
           c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell"))

metadata_normal_sparse <- meta_df_no_tumor %>%
  filter(tumor_percentage_grouped035ext_hex_average_method2 %in% c("normal PC percentage"))

# 2. (Optional) remove low-sample patients
low_groups <- metadata_normal_sparse %>%
  group_by(ObjectName_anonymous_grouped, tumor_percentage_grouped035ext_hex_average_method2) %>%
  summarise(n_regions = n(), .groups = "drop") %>%
  group_by(ObjectName_anonymous_grouped) %>%
  mutate(total_regions = sum(n_regions),
         prop = n_regions / total_regions) %>%
  ungroup() %>%
  filter(tumor_percentage_grouped035ext_hex_average_method2 %in% c("normal PC percentage")) %>%
  group_by(ObjectName_anonymous_grouped) %>%
  summarise(prop_sum = sum(prop)) %>%
  filter(prop_sum < 0.03) %>%
  pull(ObjectName_anonymous_grouped)

metadata_normal_sparse <- metadata_normal_sparse %>%
  filter(!ObjectName_anonymous_grouped %in% low_groups)

# 3. Compute proportions per patient (normal+sparse combined)
prop_df <- metadata_normal_sparse %>%
  group_by(ObjectName_anonymous_grouped,
           annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(ObjectName_anonymous_grouped) %>%
  mutate(
    Freq = cell_count / sum(cell_count) * 100
  ) %>%
  ungroup()

# 4. create wide table like propeller output

prop_trans_res_normal <- prop_df %>%
  rename(
    sample = ObjectName_anonymous_grouped,
    clusters = annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves
  ) %>%
  select(sample, clusters, Freq)


# 5. join status

metadata_unique <- metadata_normal_sparse %>%
  distinct(ObjectName_anonymous_grouped, Status_simp_2)

prop_trans_res_normal <- prop_trans_res_normal %>%
  left_join(metadata_unique,
            by = c("sample" = "ObjectName_anonymous_grouped"))


# 6. clean nonsense values

prop_trans_res_normal <- prop_trans_res_normal %>%
  filter(!is.na(Freq),
         !is.nan(Freq))

prop_trans_res_normal


statuses <- unique(prop_trans_res_normal$Status_simp_2)
combs_ordered <- expand.grid(group1 = statuses, group2 = statuses, stringsAsFactors = FALSE)
combs_ordered <- combs_ordered[combs_ordered$group1 != combs_ordered$group2, ]

# Order by group1, then group2 alphabetically
combs_ordered <- combs_ordered[order(combs_ordered$group1, combs_ordered$group2), ]




pdf(file= file.path(output_dir,"Figure_2","only_normal_normalTtest.pdf"), width = 8, height = 6)

for (i in seq_len(nrow(combs_ordered))) {
  group1 <- combs_ordered$group1[i]
  group2 <- combs_ordered$group2[i]
  
  p <- plot_patient_deltas_boxplot_with_stat_same_order_capped_normal_T_test(metadata_normal_sparse , group1 = group1, group2 = group2, method = "ratio", cluster_order_defined = cluster_ordered_for_s_curves) +
    ggtitle(paste0("Ratio Proportions: ", group1, " vs ", group2))
  
  print(p)
}

dev.off()


# Sup Figure 3B only normal and sparse Vertical box plot comparing stages ----------------------------------------------------------------

# 1. subset meta
meta_df <- Xenium_Object_Seurat@meta.data

meta_df_no_tumor <- meta_df %>%
  filter(!annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves %in%
           c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell"))

metadata_normal_sparse <- meta_df_no_tumor %>%
  filter(tumor_percentage_grouped035ext_hex_average_method2 %in% c("normal PC percentage", "sparse"))

# 2. (Optional) remove low-sample patients
low_groups <- metadata_normal_sparse %>%
  group_by(ObjectName_anonymous_grouped, tumor_percentage_grouped035ext_hex_average_method2) %>%
  summarise(n_regions = n(), .groups = "drop") %>%
  group_by(ObjectName_anonymous_grouped) %>%
  mutate(total_regions = sum(n_regions),
         prop = n_regions / total_regions) %>%
  ungroup() %>%
  filter(tumor_percentage_grouped035ext_hex_average_method2 %in% c("normal PC percentage", "sparse")) %>%
  group_by(ObjectName_anonymous_grouped) %>%
  summarise(prop_sum = sum(prop)) %>%
  filter(prop_sum < 0.03) %>%
  pull(ObjectName_anonymous_grouped)

metadata_normal_sparse <- metadata_normal_sparse %>%
  filter(!ObjectName_anonymous_grouped %in% low_groups)

# 3. Compute proportions per patient (normal+sparse combined)
prop_df <- metadata_normal_sparse %>%
  group_by(ObjectName_anonymous_grouped,
           annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(ObjectName_anonymous_grouped) %>%
  mutate(
    Freq = cell_count / sum(cell_count) * 100
  ) %>%
  ungroup()



# 4. create wide table like propeller output

prop_trans_res_sparse <- prop_df %>%
  rename(
    sample = ObjectName_anonymous_grouped,
    clusters = annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves
  ) %>%
  select(sample, clusters, Freq)


# 5. join status

metadata_unique <- metadata_normal_sparse %>%
  distinct(ObjectName_anonymous_grouped, Status_simp_2)

prop_trans_res_sparse <- prop_trans_res_sparse %>%
  left_join(metadata_unique,
            by = c("sample" = "ObjectName_anonymous_grouped"))


# 6. clean nonsense values

prop_trans_res_sparse <- prop_trans_res_sparse %>%
  filter(!is.na(Freq),
         !is.nan(Freq))

prop_trans_res_sparse
prop_trans_res_sparse_notnormalized <- prop_trans_res_sparse


statuses <- unique(prop_trans_res_sparse$Status_simp_2)
combs_ordered <- expand.grid(group1 = statuses, group2 = statuses, stringsAsFactors = FALSE)
combs_ordered <- combs_ordered[combs_ordered$group1 != combs_ordered$group2, ]

# Order by group1, then group2 alphabetically
combs_ordered <- combs_ordered[order(combs_ordered$group1, combs_ordered$group2), ]




pdf(file= file.path(output_dir,"Figure_2","only_normal_and_sparse_normalTtest.pdf"), width = 8, height = 6)

for (i in seq_len(nrow(combs_ordered))) {
  group1 <- combs_ordered$group1[i]
  group2 <- combs_ordered$group2[i]
  
  p <- plot_patient_deltas_boxplot_with_stat_same_order_capped_normal_T_test(prop_trans_res_sparse , group1 = group1, group2 = group2, method = "ratio", cluster_order_defined = cluster_ordered_for_s_curves) +
    ggtitle(paste0("Ratio Proportions: ", group1, " vs ", group2))
  
  print(p)
}
dev.off()


# Sup Figure 3C ratio all vs all Vertical box plot comparing stages ------------------------------------------

# 1. subset meta
meta_df <- Xenium_Object_Seurat@meta.data

meta_df_no_tumor <- meta_df %>%
  filter(!annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves %in%
           c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell"))



# 3. Compute proportions per patient (normal+sparse combined)
prop_df <- meta_df_no_tumor %>%
  group_by(ObjectName_anonymous_grouped,
           annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(ObjectName_anonymous_grouped) %>%
  mutate(
    Freq = cell_count / sum(cell_count) * 100
  ) %>%
  ungroup()



# 4. create wide table like propeller output

prop_trans_res_all <- prop_df %>%
  rename(
    sample = ObjectName_anonymous_grouped,
    clusters = annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves
  ) %>%
  select(sample, clusters, Freq)


# 5. join status

metadata_unique <- meta_df_no_tumor %>%
  distinct(ObjectName_anonymous_grouped, Status_simp_2)

prop_trans_res_all <- prop_trans_res_all %>%
  left_join(metadata_unique,
            by = c("sample" = "ObjectName_anonymous_grouped"))


# 6. clean nonsense values

prop_trans_res_all <- prop_trans_res_all %>%
  filter(!is.na(Freq),
         !is.nan(Freq))

prop_trans_res_all


statuses <- unique(prop_trans_res_all$Status_simp_2)
combs_ordered <- expand.grid(group1 = statuses, group2 = statuses, stringsAsFactors = FALSE)
combs_ordered <- combs_ordered[combs_ordered$group1 != combs_ordered$group2, ]

# Order by group1, then group2 alphabetically
combs_ordered <- combs_ordered[order(combs_ordered$group1, combs_ordered$group2), ]




pdf(file= file.path(output_dir,"Figure_2","all_vs_all_normal_Ttest.pdf"), width = 8, height = 6)

for (i in seq_len(nrow(combs_ordered))) {
  group1 <- combs_ordered$group1[i]
  group2 <- combs_ordered$group2[i]
  
  p <- plot_patient_deltas_boxplot_with_stat_same_order_capped_normal_T_test(prop_trans_res_all , group1 = group1, group2 = group2, method = "ratio", cluster_order_defined = cluster_ordered_for_s_curves) +
    ggtitle(paste0("Ratio Proportions: ", group1, " vs ", group2))
  
  print(p)
}

dev.off()

# Sup Figure 3D-E dense vs all Vertical box plot comparing stages ----------------------------------------------------------------

# 1. subset meta
meta_df <- Xenium_Object_Seurat@meta.data

meta_df_no_tumor <- meta_df %>%
  filter(!annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves %in%
           c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell"))

metadata_dense_NDMM_against_whole_rest <- meta_df_no_tumor %>%
  filter(!(Status_simp_2 == "NDMM" & 
             tumor_percentage_grouped035ext_hex_average_method2 %in% 
             c("normal PC percentage", "sparse", "unqualified", "dispersed")))

#metadata_dense_NDMM_against_whole_rest <- metadata_dense_NDMM_against_whole_rest %>%
# filter(!ObjectName_anonymous_grouped %in%
#         c("NDMM_1"))



# 3. Compute proportions per patient (normal+sparse combined)
prop_df <- metadata_dense_NDMM_against_whole_rest %>%
  group_by(ObjectName_anonymous_grouped,
           annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(ObjectName_anonymous_grouped) %>%
  mutate(
    Freq = cell_count / sum(cell_count) * 100
  ) %>%
  ungroup()



# 4. create wide table like propeller output

prop_trans_res_densevsall <- prop_df %>%
  rename(
    sample = ObjectName_anonymous_grouped,
    clusters = annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves
  ) %>%
  select(sample, clusters, Freq)


# 5. join status

metadata_unique <- metadata_dense_NDMM_against_whole_rest %>%
  distinct(ObjectName_anonymous_grouped, Status_simp_2)

prop_trans_res_densevsall <- prop_trans_res_densevsall %>%
  left_join(metadata_unique,
            by = c("sample" = "ObjectName_anonymous_grouped"))


# 6. clean nonsense values

prop_trans_res_densevsall <- prop_trans_res_densevsall %>%
  filter(!is.na(Freq),
         !is.nan(Freq))

prop_trans_res_densevsall


statuses <- unique(prop_trans_res_densevsall$Status_simp_2)
combs_ordered <- expand.grid(group1 = statuses, group2 = statuses, stringsAsFactors = FALSE)
combs_ordered <- combs_ordered[combs_ordered$group1 != combs_ordered$group2, ]

# Order by group1, then group2 alphabetically
combs_ordered <- combs_ordered[order(combs_ordered$group1, combs_ordered$group2), ]




pdf(file= file.path(output_dir,"Figure_2","dense_vs_all_normalTtest_with_NDMM_1.pdf"), width = 8, height = 6)

for (i in seq_len(nrow(combs_ordered))) {
  group1 <- combs_ordered$group1[i]
  group2 <- combs_ordered$group2[i]
  
  p <- plot_patient_deltas_boxplot_with_stat_same_order_capped_normal_T_test(prop_trans_res_densevsall , group1 = group1, group2 = group2, method = "ratio", cluster_order_defined = cluster_ordered_for_s_curves) +
    ggtitle(paste0("Ratio Proportions: ", group1, " vs ", group2))
  
  print(p)
}

dev.off()








# Sup Figure 3G compare S curves linked SMM ---------------------------------------------------------------

df_all_combined <-  readRDS("S_curves_combined_new_order.rds")

df_summary <- df_all_combined %>%
  group_by(clusters, group_comparison) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    mean_value   = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

library(ggplot2)

df_all_filtered <- df_all_combined[df_all_combined$group_comparison %in% c("all_SMM vs CBM", "all_PCL vs NDMM",  "all_PCL vs CBM", 
                                                                           "all_NDMM vs CBM","normal_NDMM vs CBM", "dense_vs_all_NDMM vs CBM",
                                                                           "all_NDMM vs SMM", "sparse_NDMM vs SMM","dense_vs_all_NDMM vs SMM"  ),]



df_all_filtered <- df_all_combined[df_all_combined$group_comparison %in% c("all_NDMM vs SMM", "sparse_NDMM vs SMM","dense_vs_all_NDMM vs SMM"  ),]
#df_all_filtered <- df_all_combined[df_all_combined$group_comparison %in% c( "sparse_NDMM vs SMM","dense_vs_all_NDMM vs SMM"  ),]


comparison_order <- c("all_SMM vs CBM", "all_PCL vs NDMM",  "all_PCL vs CBM", 
                      "all_NDMM vs CBM","normal_NDMM vs CBM", "dense_vs_all_NDMM vs CBM",
                      "all_NDMM vs SMM", "sparse_NDMM vs SMM","dense_vs_all_NDMM vs SMM"  )
df_all_filtered$group_comparison <- factor(df_all_filtered$group_comparison, levels = comparison_order)



clip_value <- 5  # maximum |value| to display

df_all_filtered <- df_all_filtered %>%
  mutate(value_clipped = pmin(pmax(value, -clip_value), clip_value))



# Filter only the two comparisons used for ordering
df_order_calc <- df_all_filtered %>%
  filter(group_comparison %in% c("sparse_NDMM vs SMM", "dense_vs_all_NDMM vs SMM")) %>%
  group_by(clusters, group_comparison) %>%
  summarise(mean_value = mean(value_clipped, na.rm = TRUE), .groups = "drop")

# Reshape wider so each cluster has two columns
df_order_wide <- df_order_calc %>%
  tidyr::pivot_wider(
    names_from = group_comparison,
    values_from = mean_value
  )

# Compute the difference
df_order_wide <- df_order_wide %>%
  mutate(diff_order = `sparse_NDMM vs SMM` - `dense_vs_all_NDMM vs SMM`)

# Create the ordered factor for clusters
cluster_order <- df_order_wide %>%
  arrange(diff_order) %>%
  pull(clusters)

# Apply to your dataframe
df_all_filtered$clusters <- factor(df_all_filtered$clusters, levels = cluster_order)

#or same order
df_all_filtered$clusters <- factor(df_all_filtered$clusters, levels = rev(cluster_ordered_for_s_curves))


# Optional: compute positions for vertical lines between cell types
celltype_positions <- 1:length(unique(df_all_filtered$clusters))
vline_positions <- celltype_positions[-length(celltype_positions)] + 0.5

colors <- c("cyan","#377EB8","#D72638") #, "dense_otheroption_toopink"= "#FB8072")  # Replace with your preferred order
colors <- c( "#A7C7E7", "#FFE066","#D72638") #, "dense_otheroption_toopink"= "#FB8072")  # Replace with your preferred order



plot <- ggplot(df_all_filtered, aes(x = clusters, y = value_clipped, fill = group_comparison)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = vline_positions, linetype = "dotted", color = "gray70") +
  # Optional: add jitter
  # geom_jitter(aes(color = group_comparison),
  #             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
  #             size = 2, alpha = 0.6) +
  coord_cartesian(ylim = c(-5, 4)) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = "Log2 ratio (clipped at ±5)",
    x = "Cell type",
    fill = "Comparison"
  )

plot

pdf(file= file.path(output_dir,"Sup_Figure_3","S_curves_linked_SMM_sparseanddensesame_order3plots.pdf"), width = 16, height = 6)
print(plot)
dev.off()


# 1️⃣ Aggregate data: compute mean and SD per cluster and comparison
df_summary <- df_all_filtered %>%
  group_by(clusters, group_comparison) %>%
  summarise(
    mean_value = mean(value_clipped, na.rm = TRUE),
    sd_value   = sd(value_clipped, na.rm = TRUE),
    .groups = "drop"
  )

# Optional: order clusters if needed
df_summary$clusters <- factor(df_summary$clusters, levels = cluster_order)  

# 2️⃣ Plot bars with SD
plot <- ggplot(df_summary, aes(x = clusters, y = mean_value, fill = group_comparison)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = "Log2 ratio (clipped at ±5)",
    x = "Cell type",
    fill = "Comparison"
  )

plot


# Sup Figure 3I Horizontal all combined per status wo relapse ---------------------------------------------------------------

meta_df <- Xenium_Object_Seurat@meta.data


dataframe_no_tumor <- prepare_metadata_frame(
  meta = meta_df,
  exclude_tumor = TRUE,
  cluster_col = "annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv",
  region_col = "small_region_hex1",
  tumor_col = "Status_simp_2",
  levels_tumor = unique(Xenium_Object_Seurat$tumor_percentage_grouped035ext_hex_average_method2)
)

dataframe_no_tumor <- aggregate_percentages_dynamic(dataframe_no_tumor, group_col = "Status_simp_2")


dataframe_no_tumor_no_relapse <- dataframe_no_tumor[dataframe_no_tumor$Status_simp_2 != "Relapse",] 

# Step 1: Define the biological order
status_order <- c("CBM", "SMM", "NDMM", "PCL")


ordered_celltypes <- dataframe_no_tumor_no_relapse %>%
  group_by(Status_simp_2, ObjectName_anonymous_grouped, cluster_celltype) %>%
  summarise(mean_freq_patient = median(percentages_per_patient_per_Status_simp_2, na.rm = TRUE), .groups = "drop") %>%
  group_by(Status_simp_2, cluster_celltype) %>%
  summarise(mean_freq = median(mean_freq_patient, na.rm = TRUE), .groups = "drop") %>%
  mutate(Status_simp_2 = factor(Status_simp_2, levels = status_order)) %>%
  arrange(Status_simp_2, desc(mean_freq)) %>%
  group_by(cluster_celltype) %>%
  summarise(
    peak_status = Status_simp_2[which.max(mean_freq)],
    peak_freq = max(mean_freq)
  ) %>%
  mutate(peak_status = factor(peak_status, levels = status_order)) %>%
  arrange(peak_status, desc(peak_freq)) %>%
  pull(cluster_celltype)


ordered_celltypes_written <- c(
  "Neutrophil Progenitor", "MMP9 Immature Neutrophil 4", "OLFM4 Immature Neutrophil 2", 
  "LTF Immature Neutrophil 3", "Non Classical Monocyte", "CEACAM6 Immature Neutrophil 1", 
  "Osteoblast", "Interleukin Producing Cell", "Mast Cell", "Adipocyte", 
  "VSMC", "CXCL14 MSC", "AEC", "APOD MSC", "Activated or Exhausted T Cell", 
  "Interferon Stimulated Cell", "Regulatory T Cell", "LEPR MSC", "Erythroid", 
  "Macrophage", "MKC Lineage", "Fibro MSC", "CD8 T Cell", "EC", 
  "SEC", "Conventional DC", "Plasmacytoid DC", "Osteo-Fibroblastic MSC", 
  "THY1 MSC", "Mature Neutrophil", "Basophil Eosinophil", "Classical Monocyte", 
  "NK or Cytotoxic T Cell", "Naive or CM CD4 T Cell"
)



dataframe_plot <- dataframe_no_tumor_no_relapse

# Reorder factor levels for plotting
dataframe_plot$cluster_celltype <- factor(
  dataframe_plot$cluster_celltype,
  levels = rev(cluster_order_NDMM_vs_CBM)
)

# Reorder factor levels for plotting
dataframe_plot$cluster_celltype <- factor(
  dataframe_plot$cluster_celltype,
  levels = rev(cluster_ordered_for_s_curves)
)


# Step 1: Get the max observed value per celltype
max_values <- dataframe_plot %>%
  group_by(cluster_celltype) %>%
  summarize(max_pct = max(percentages_per_patient_per_Status_simp_2 , na.rm = TRUE)) %>%
  arrange(desc(max_pct))

# Step 2: Define high vs low abundance clusters based on threshold
abundance_threshold <- 10 # e.g., 10%
high_abundance_clusters <- max_values %>%
  filter(max_pct > abundance_threshold) %>%
  pull(cluster_celltype)

# Step 3: Add abundance class to your dataframe
dataframe_plot <- dataframe_plot %>%
  mutate(abundance_class = ifelse(cluster_celltype %in% high_abundance_clusters, "High Abundant Cell Types", "Low Abundant Cell Types"))


# Step 4: Plot with facets
library(dplyr)

dataframe_plot <- dataframe_plot %>%
  distinct(
    cluster_celltype,
    ObjectName_anonymous_grouped,
    Status_simp_2,
    .keep_all = TRUE
  )


# Step 1: make sure tumor_percentage_grouped is ordered
dataframe_plot$Status_simp_2_fac <- factor(
  dataframe_plot$Status_simp_2,
  levels = c("CBM", "SMM", "NDMM", "PCL") # <-- adjust if you prefer
)

dataframe_tumor_unique <- dataframe_plot %>%
  distinct(cluster_celltype, ObjectName_anonymous_grouped, Status_simp_2, .keep_all = TRUE)

clusters_to_remove <- c("B cells or Plasma cells", "Cycling B or Plasma cells", "Plasma Cell","B or Plasma Cell", "Cycling B or Plasma Cell", "Plasma cells") # change to your unwanted clusters

dataframe_tumor_unique <- dataframe_tumor_unique %>%
  filter(!cluster_celltype %in% clusters_to_remove & !is.na(cluster_celltype) &!is.na(Status_simp_2_fac))

#dataframe_tumor_unique <- dataframe_tumor %>%
#filter(!is.na(tumor_percentage_grouped_fac) & !is.na(cluster_celltype))
niche_colors_adjusted <- c(
  "CBM" = "#0072B2",               # dark blue instead of cyan
  "unqualified" = "#4CAF50",       # keep green
  "normal PC percentage" = "#E69F00", # darker orange instead of pale yellow
  "sparse" = "#377EB8",            # keep blue
  "dispersed" = "#6A0DAD",         # keep purple
  "dense" = "darkred",              # keep red
  "dense_otheroption_toopink" = "#FB8072" # keep salmon
)

## add statistics

library(dplyr)




## place stars yourself in Illustrator

plot <- ggplot(dataframe_tumor_unique, aes(x = cluster_celltype, y = percentages_per_patient_per_Status_simp_2)) +
  geom_boxplot(
    aes(fill = Status_simp_2_fac),
    position = position_dodge(width = 0.8),
    alpha = 0.3,
    outlier.shape = NA,
    width = 0.6,
  ) +
  geom_jitter(
    aes(color = Status_simp_2_fac, group = Status_simp_2_fac),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, jitter.height = 0),
    size = 1,
    alpha = 0.6   # transparency to avoid overplotting
  ) +
  
  facet_wrap(~ abundance_class, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 13, face = "bold"
    ),
    axis.text.y = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = status_colors) +
  scale_color_manual(values = status_colors) +
  labs(
    title = "Cell proportions per cluster across tumor architecture",
    subtitle = "Faceted by abundance grouping",
    x = "Cluster",
    y = "Proportion (Freq)"
  )

plot

pdf(file= file.path(output_dir,"Sup_Figure_2","Status_simp_barplotallhorizontal_celltypebasedorder.pdf"), 
    width = 20, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Sup Figure 4A IMC data UMAP ---------------------------

# Load necessary libraries
library(ggplot2)
library(dplyr)
## activating packages
library(BiocStyle)
library(imcRtools)
library(BiocManager)
library(remotes)
library(imcRtools)
library(cytomapper) 
library(dplyr)
# Load necessary libraries
library(sf)

library(dplyr)
library(dittoSeq)
library(scater)
library(patchwork)
library(cowplot)
library(viridis)



plot <- dittoDimPlot(spe_IMC, 
                     var = "cluster_celltype_integrated_202511_good_MKC_CD14_seperate", 
                     reduction.use = "UMAP_harmony", 
                     size = 0.2,
                     do.label = TRUE) +
  scale_color_manual(values = celltype_colors_IMC) 
#theme(legend.title = element_blank()) +
plot

pdf(file= file.path(output_dir,"Sup_Figure_5","IMC_UMAP_2.pdf"), 
    width = 15, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

plot2 <- dittoDimPlot(spe_IMC, 
                      var = "cluster_celltype_integrated_202511_good_MKC_CD14_seperate", 
                      reduction.use = "UMAP_harmony", 
                      size = 0.2,
                      do.label = TRUE,
                      do.raster = TRUE) +
  scale_color_manual(values = celltype_colors_IMC) 
#theme(legend.title = element_blank()) +
plot2


#Create the plot and rasterize points
# Rasterize all point layers in the plot
plot2 <- ggrastr::rasterise(plot, layers = "geom_point", dpi = 300)
plot2

# Save to PDF
pdf(file = file.path(output_dir, "Sup_Figure_5", "IMC_UMAP_lowerres.pdf"), 
    width = 15, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot2)
dev.off()


# Sup Figure 4B IMC data Heatmap ---------------------------


# 1️⃣ Channels to use (exclude certain markers)
rowData(spe_IMC)$use_channel2 <- !grepl(
  "Background|DNA|Radius|Eccentricity|Area|CD63|FOXP3|CLEC9A|GZMB", 
  rownames(spe_IMC)
)


## by cell type mean
celltype_mean <- aggregateAcrossCells(as(spe_IMC, "SingleCellExperiment"),  
                                      ids = spe_IMC$cluster_celltype_integrated_202511_good_MKC_CD14_seperate, 
                                      statistics = "mean",
                                      use.assay.type = "exprs_asinh", 
                                      subset.row = rownames(spe_IMC)[rowData(spe_IMC)$use_channel2])

# No scaling
plot <- dittoHeatmap(celltype_mean,
                     genes = rownames(spe_IMC)[rowData(spe_IMC)$use_channel2],
                     assay = "exprs_asinh", cluster_cols = FALSE, 
                     scale = "none",
                     heatmap.colors = viridis(100),
                     scaled.to.max = FALSE,
                     annot.by = c("cluster_celltype_integrated_202511_good_MKC_CD14_seperate"))

plot

# 1️⃣ Extract clusters present in your data
present_clusters <- unique(celltype_mean$cluster_celltype_integrated_202511_good_MKC_CD14_seperate)

# 2️⃣ Define desired order for the heatmap (use your original vector order)
desired_order <- names(celltype_colors_IMC[celltype_colors_IMC %in% celltype_colors_IMC[present_clusters]])

# 3️⃣ Make the annotation column a factor with levels in desired order
celltype_mean$cluster_celltype_integrated_202511_good_MKC_CD14_seperate <- 
  factor(celltype_mean$cluster_celltype_integrated_202511_good_MKC_CD14_seperate,
         levels = desired_order)

# 4️⃣ Subset and order colors to match factor levels
ordered_colors <- celltype_colors_IMC[desired_order]

# 5️⃣ Generate heatmap
plot2 <- dittoHeatmap(
  celltype_mean,
  genes = rownames(spe_IMC)[rowData(spe_IMC)$use_channel2],
  assay = "exprs_asinh",
  cluster_cols = FALSE,
  scale = "none",
  heatmap.colors = viridis(100),
  scaled.to.max = FALSE,
  annot.by = "cluster_celltype_integrated_202511_good_MKC_CD14_seperate",
  annot.colors = ordered_colors  # vector directly
)

plot2


pdf(file= file.path(output_dir,"Sup_Figure_5_IMC","IMC_heatmap_ordered_noGZMB.pdf"), 
    width = 15, height = 8, useDingbats = FALSE, onefile = FALSE)
print(plot2)
dev.off()


# Sup Figure 4C patient overview -------------------------------------------


library(dplyr)
library(ggplot2)

## run first the code above

#meta_df2 <- as.data.frame(meta_df2)

plot_df <- meta_df2 %>%
  distinct(patient_id, Status_simp_2, shared_patient) %>%   # count patients uniquely
  count(Status_simp_2, shared_patient, name = "n")

plot <- ggplot(plot_df, aes(x = Status_simp_2, y = n, fill = shared_patient)) +
  geom_col(position = "stack") +
  labs(
    x = "Status",
    y = "Number of patients",
    fill = "Run on ST Xenium platform"
  ) +
  theme_classic()

plot

pdf(file= file.path(output_dir,"Sup_Figure_5","IMC_patient_overview.pdf"), 
    width = 6, height = 7, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()



# Sup Figure 4D IMC data barplot per tumor region ---------------------------


meta_df <- colData(spe_IMC)
#meta_df_NDMM_SMM_CBM_PCL <- meta_df[meta_df$Status_simp_2 %in% c("NDMM", "SMM", "CBM", "PCL"), ]

df_plot <- as.data.frame(meta_df) %>%
  dplyr::select(
    ObjectName = tumor_percentage_grouped0.35_CBM,
    annot_evenless_broad_run3 = cluster_celltype_integrated_202511_good_MKC_CD14_seperate 
  )

# Step 2: Reorder 'annot_MM_HBM_more_simple' based on your custom order
# Ensure meta_df$cluster_name is a factor ordered by your color vector

# Step 3: Calculate proportions
df_plot <- df_plot %>%
  group_by(ObjectName, annot_evenless_broad_run3) %>%
  tally() %>%  # Count occurrences of each combination
  ungroup() %>%
  group_by(ObjectName) %>%
  mutate(proportion = n / sum(n))  # Calculate proportion within each ObjectName

# Ensure ObjectName is a factor with custom order
df_plot$ObjectName_fac <- factor(df_plot$ObjectName, levels = names(region_colors))
#df_plot$ObjectName <- factor(df_plot$ObjectName, levels = architecture_order)

df_plot <- df_plot %>% filter(!(ObjectName %in% c(NA, "unqualified")))

df_plot$cluster_factor <- factor(
  df_plot$annot_evenless_broad_run3,        # replace with your actual cluster column
  levels = names(celltype_colors_IMC)  # order according to the color vector
)

# join IMC data
library(dplyr)

plot <- ggplot(df_plot, aes(x = ObjectName_fac, y = proportion, fill = cluster_factor       )) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot
  scale_fill_manual(values = celltype_colors_IMC) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x-axis labels and increase size
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove background
        axis.title.x = element_blank(),  # Optional: Remove x-axis title
        axis.title.y = element_blank(),  # Optional: Remove y-axis title
        text = element_text(size = 15))
plot



pdf(file= file.path(output_dir,"Sup_Figure_4","IMC_barplot_tumor_regions.pdf"), 
    width = 8, height = 12, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Sup Figure 4E IMC data barplot combined ---------------------------


meta_df <- Xenium_Object_Seurat @meta.data

meta_df2 <- colData(spe_IMC)




# For IMC metadata (meta_df2)
meta_df2$common_group <- dplyr::case_when(
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Ki67_likely_erythroid") ~ "Erythroid",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Plasma_Cell") ~ "Plasma Cell",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD105_Structural_Cell", "Stroma") ~ "Stroma",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("MKC_Endothelial") ~ "MKC or Endothelial",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Macrophage") ~ "Macrophage",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Neutrophil") ~ "Neutrophil",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Neutrophil_Progenitor") ~ "Neutrophil Progenitor",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD14_Monocyte") ~ "Monocyte",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD8_T_Cell") ~ "CD8 T Cell",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD4_T_Cell") ~ "CD4 T Cell",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("HLA-DR_Cell") ~ "APC",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Fibrotic_Tissue", "Undefined_Cell") ~ "Other",
  TRUE ~ "Other"
)

# For the other metadata (meta_df)
meta_df$common_group <- dplyr::case_when(
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Erythroid") ~ "Erythroid",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell") ~ "Plasma Cell",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Osteoblast" ,"LEPR MSC","THY1 MSC","APOD MSC","Osteo-Fibroblastic MSC","Fibro MSC","CXCL14 MSC","VSMC", "Adipocyte") ~ "Stroma",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("MKC Lineage","SEC","EC","AEC") ~ "MKC or Endothelial",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Macrophage") ~ "Macrophage",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Mature Neutrophil","LTF Immature Neutrophil 3","MMP9 Immature Neutrophil 4","OLFM4 Immature Neutrophil 2","CEACAM6 Immature Neutrophil 1") ~ "Neutrophil",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Neutrophil Progenitor") ~ "Neutrophil Progenitor",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Classical Monocyte","Non Classical Monocyte") ~ "Monocyte",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("CD8 T Cell","NK or Cytotoxic T Cell","Activated or Exhausted T Cell") ~ "CD8 T Cell",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Regulatory T Cell" ,"Naive or CM CD4 T Cell") ~ "CD4 T Cell",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Plasmacytoid DC","Conventional DC","Plasmacytoid DC") ~ "APC",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Interleukin Producing Cell","Basophil Eosinophil","Mast Cell") ~ "Other",
  TRUE ~ "Other"
)


common_group_colors <- c(
  "Other"                 = "darkgray",       # from Fibrotic_Tissue / Undefined_Cell
  "Erythroid"             = "pink",           # from Ki67_likely_erythroid
  "MKC or Endothelial"    = "#90ee90",        # from MKC_Endothelial
  "Stroma"                = "#1C750C",         # from CD105_Structural_Cell
  "Neutrophil Progenitor" = "#00FFFF",        # from Neutrophil_Progenitor
  "Neutrophil"            = "#00BFFF",        # from Neutrophil
  "Macrophage"            = "#005BFF",        # from Macrophage
  "Monocyte"              = "#c51b8a",        # from CD14_Monocyte
  "APC"                  = "purple",    
  "CD8 T Cell"            = "#FA7921",        # from CD8_T_Cell
  "CD4 T Cell"            = "#FFE66D",        # from CD4_T_Cell
  "Plasma Cell"           = "#FF6666"       # from Plasma_Cell
)



meta_df$Technique <- "Xenium"
meta_df2$Technique <- "IMC"

meta_df2$tumor_percentage_grouped035ext_hex_average_method2_CBM <- meta_df2$tumor_percentage_grouped0.35_CBM
meta_df2$Status_simp_2 <- meta_df2$status

# Vector of patient IDs to remove (IMC IDs not present in Xenium)
#patients_to_remove <- c("pt_7290", "08899", "8953", "18574")

# Subset the dataframe
#meta_df2_remove_patients_not_in_Xenium <- meta_df2[!meta_df2$patient_id %in% patients_to_remove, ]




# Step 1: Select relevant columns from both datasets
df_IMC <- meta_df2 %>%
  as.data.frame() %>%
  dplyr::select(
    ObjectName = Status_simp_2,
    common_group,
    tumor_group = tumor_percentage_grouped035ext_hex_average_method2_CBM,
    Technique
  )

df_Xenium <- meta_df %>%
  as.data.frame() %>%
  dplyr::select(
    ObjectName = Status_simp_2,
    common_group,
    tumor_group = tumor_percentage_grouped035ext_hex_average_method2_CBM,
    Technique
  )

# Step 2: Combine both datasets
df_combined <- bind_rows(df_IMC, df_Xenium)

# Step 3: Remove NAs or unqualified entries if needed
df_combined <- df_combined %>% filter(!is.na(tumor_group ), tumor_group  != "unqualified")

# Step 4: Calculate proportions of each common_group per sample
df_plot <- df_combined %>%
  group_by(ObjectName, Technique, common_group) %>%
  tally() %>%
  group_by(ObjectName, Technique) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()


# Factor for clusters according to your colors
df_plot$cluster_factor <- factor(df_plot$common_group, levels = names(common_group_colors))
df_plot$ObjectName_fac <- factor(df_plot$ObjectName, levels = names(status_colors)) # adjust as needed

# Define desired order for the combined bars
combined_order <- as.vector(sapply(levels(df_plot$ObjectName_fac), function(x) {
  paste0(x, "_", c("IMC", "Xenium"))  # IMC first, Xenium second
}))

# Create combined column and make it a factor with proper order
df_plot$ObjectName_Technique <- factor(
  paste0(df_plot$ObjectName_fac, "_", df_plot$Technique),
  levels = combined_order
)


#Plot
plot <- ggplot(df_plot, aes(x = ObjectName_Technique, y = proportion, fill = cluster_factor)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = common_group_colors) +
  #facet_wrap(~Technique, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 15)
  )

plot


pdf(file= file.path(output_dir,"Sup_Figure_5","IMC_and_Xenium_combined_status_barplot_relapse.pdf"), 
    width = 8, height = 12, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Step 4: Calculate proportions of each common_group per sample
df_plot <- df_combined %>%
  group_by(tumor_group , Technique, common_group) %>%
  tally() %>%
  group_by(tumor_group , Technique) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()




# Factor for clusters according to your colors
df_plot$cluster_factor <- factor(df_plot$common_group, levels = names(common_group_colors))
df_plot$ObjectName_fac <- factor(df_plot$tumor_group , levels = names(region_colors)) # adjust as needed

# Define desired order for the combined bars
combined_order <- as.vector(sapply(levels(df_plot$ObjectName_fac), function(x) {
  paste0(x, "_", c("IMC", "Xenium"))  # IMC first, Xenium second
}))

# Create combined column and make it a factor with proper order
df_plot$ObjectName_Technique <- factor(
  paste0(df_plot$ObjectName_fac, "_", df_plot$Technique),
  levels = combined_order
)

#Plot
plot <- ggplot(df_plot, aes(x = ObjectName_Technique, y = proportion, fill = cluster_factor)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = common_group_colors) +
  #facet_wrap(~Technique, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 15)
  )

plot

pdf(file= file.path(output_dir,"Sup_Figure_5","IMC_and_Xenium_combined_tumor_grouped_barplot_horizontal.pdf"), 
    width = 8, height = 12, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


meta_df <- Xenium_Object_Seurat @meta.data

meta_df2 <- colData(spe_IMC)






# For IMC metadata (meta_df2)
meta_df2$common_group <- dplyr::case_when(
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Ki67_likely_erythroid") ~ "Erythroid",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Plasma_Cell") ~ "Plasma Cell",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD105_Structural_Cell", "Stroma") ~ "Stroma",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("MKC_Endothelial") ~ "MKC or Endothelial",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Macrophage") ~ "Macrophage",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Neutrophil") ~ "Neutrophil",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Neutrophil_Progenitor") ~ "Neutrophil Progenitor",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD14_Monocyte") ~ "Monocyte",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD8_T_Cell") ~ "CD8 T Cell",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("CD4_T_Cell") ~ "CD4 T Cell",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("HLA-DR_Cell") ~ "APC",
  meta_df2$cluster_celltype_integrated_202511_good_MKC_CD14_seperate %in% c("Fibrotic_Tissue", "Undefined_Cell") ~ "Other",
  TRUE ~ "Other"
)

# For the other metadata (meta_df)
meta_df$common_group <- dplyr::case_when(
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Erythroid") ~ "Erythroid",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell") ~ "Plasma Cell",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Osteoblast" ,"LEPR MSC","THY1 MSC","APOD MSC","Osteo-Fibroblastic MSC","Fibro MSC","CXCL14 MSC","VSMC", "Adipocyte") ~ "Stroma",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("MKC Lineage","SEC","EC","AEC") ~ "MKC or Endothelial",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Macrophage") ~ "Macrophage",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Mature Neutrophil","LTF Immature Neutrophil 3","MMP9 Immature Neutrophil 4","OLFM4 Immature Neutrophil 2","CEACAM6 Immature Neutrophil 1") ~ "Neutrophil",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Neutrophil Progenitor") ~ "Neutrophil Progenitor",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Classical Monocyte","Non Classical Monocyte") ~ "Monocyte",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("CD8 T Cell","NK or Cytotoxic T Cell","Activated or Exhausted T Cell") ~ "CD8 T Cell",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Regulatory T Cell" ,"Naive or CM CD4 T Cell") ~ "CD4 T Cell",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Plasmacytoid DC","Conventional DC","Plasmacytoid DC") ~ "APC",
  meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Interleukin Producing Cell","Basophil Eosinophil","Mast Cell") ~ "Other",
  TRUE ~ "Other"
)


common_group_colors <- c(
  "Other"                 = "darkgray",       # from Fibrotic_Tissue / Undefined_Cell
  "Erythroid"             = "pink",           # from Ki67_likely_erythroid
  "MKC or Endothelial"    = "#90ee90",        # from MKC_Endothelial
  "Stroma"                = "#1C750C",         # from CD105_Structural_Cell
  "Neutrophil Progenitor" = "#00FFFF",        # from Neutrophil_Progenitor
  "Neutrophil"            = "#00BFFF",        # from Neutrophil
  "Macrophage"            = "#005BFF",        # from Macrophage
  "Monocyte"              = "#c51b8a",        # from CD14_Monocyte
  "APC"                  = "purple",    
  "CD8 T Cell"            = "#FA7921",        # from CD8_T_Cell
  "CD4 T Cell"            = "#FFE66D",        # from CD4_T_Cell
  "Plasma Cell"           = "#FF6666"       # from Plasma_Cell
)



meta_df$Technique <- "Xenium"
meta_df2$Technique <- "IMC"

meta_df2$tumor_percentage_grouped035ext_hex_average_method2_CBM <- meta_df2$tumor_percentage_grouped0.35_CBM
meta_df2$Status_simp_2 <- meta_df2$status

# Vector of patient IDs to remove (IMC IDs not present in Xenium)
#patients_to_remove <- c("pt_7290", "08899", "8953", "18574")

# Subset the dataframe
#meta_df2_remove_patients_not_in_Xenium <- meta_df2[!meta_df2$patient_id %in% patients_to_remove, ]




# Step 1: Select relevant columns from both datasets
df_IMC <- meta_df2 %>%
  as.data.frame() %>%
  dplyr::select(
    ObjectName = Status_simp_2,
    common_group,
    tumor_group = tumor_percentage_grouped035ext_hex_average_method2_CBM,
    Technique
  )

df_Xenium <- meta_df %>%
  as.data.frame() %>%
  dplyr::select(
    ObjectName = Status_simp_2,
    common_group,
    tumor_group = tumor_percentage_grouped035ext_hex_average_method2_CBM,
    Technique
  )

# Step 2: Combine both datasets
df_combined <- bind_rows(df_IMC, df_Xenium)

# Step 3: Remove NAs or unqualified entries if needed
df_combined <- df_combined %>% filter(!is.na(tumor_group ), tumor_group  != "unqualified")

# Step 4: Calculate proportions of each common_group per sample
df_plot <- df_combined %>%
  group_by(ObjectName, Technique, common_group) %>%
  tally() %>%
  group_by(ObjectName, Technique) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()


# Factor for clusters according to your colors
df_plot$cluster_factor <- factor(df_plot$common_group, levels = names(common_group_colors))
df_plot$ObjectName_fac <- factor(df_plot$ObjectName, levels = names(status_colors)) # adjust as needed

# Define desired order for the combined bars
combined_order <- as.vector(sapply(levels(df_plot$ObjectName_fac), function(x) {
  paste0(x, "_", c("IMC", "Xenium"))  # IMC first, Xenium second
}))

# Create combined column and make it a factor with proper order
df_plot$ObjectName_Technique <- factor(
  paste0(df_plot$ObjectName_fac, "_", df_plot$Technique),
  levels = combined_order
)


#Plot
plot <- ggplot(df_plot, aes(x = ObjectName_Technique, y = proportion, fill = cluster_factor)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = common_group_colors) +
  #facet_wrap(~Technique, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 15)
  )

plot


pdf(file= file.path(output_dir,"Sup_Figure_5","IMC_and_Xenium_combined_status_barplot.pdf"), 
    width = 8, height = 12, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Step 4: Calculate proportions of each common_group per sample
df_plot <- df_combined %>%
  group_by(tumor_group , Technique, common_group) %>%
  tally() %>%
  group_by(tumor_group , Technique) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()




# Factor for clusters according to your colors
df_plot$cluster_factor <- factor(df_plot$common_group, levels = names(common_group_colors))
df_plot$ObjectName_fac <- factor(df_plot$tumor_group , levels = names(region_colors)) # adjust as needed

# Define desired order for the combined bars
combined_order <- as.vector(sapply(levels(df_plot$ObjectName_fac), function(x) {
  paste0(x, "_", c("IMC", "Xenium"))  # IMC first, Xenium second
}))

# Create combined column and make it a factor with proper order
df_plot$ObjectName_Technique <- factor(
  paste0(df_plot$ObjectName_fac, "_", df_plot$Technique),
  levels = rev(combined_order)
)

#Plot
plot <- ggplot(df_plot, aes(x = ObjectName_Technique, y = proportion, fill = cluster_factor)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = common_group_colors) +
  coord_flip() +
  #facet_wrap(~Technique, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 15)
  )

plot

pdf(file= file.path(output_dir,"Sup_Figure_5","IMC_and_Xenium_combined_tumor_grouped_barplot_horizontal.pdf"), 
    width = 17, height = 6, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()



# Sup Figure 4F IMC data barplot combined same patients ---------------------------


library(dplyr)


# Convert meta_df2 to data.frame
meta_df2 <- as.data.frame(meta_df2)

# Extract numeric patient IDs
meta_df2$patient_number <- str_extract(meta_df2$patient_id, "\\d+")
meta_df$patient_number <- str_extract(meta_df$ObjectName, "\\d+")

# Identify shared patients
shared_patients <- intersect(meta_df2$patient_number, meta_df$patient_number)
non_shared_patients <- setdiff(meta_df2$patient_number, meta_df$patient_number)

# Add shared_patient column
meta_df2$shared_patient <- ifelse(meta_df2$patient_number %in% shared_patients, "Yes", "No")
meta_df$shared_patient <- ifelse(meta_df$patient_number %in% shared_patients, "Yes", "No")

# Subset both datasets to only shared patients
meta_df2_shared <- meta_df2[meta_df2$patient_number %in% shared_patients, ]
meta_df_shared  <- meta_df[meta_df$patient_number %in% shared_patients, ]


# Vector of patient IDs to remove (IMC IDs not present in Xenium)
#patients_to_remove <- c("pt_7290", "08899", "8953", "18574")

# Subset the dataframe
#meta_df2_remove_patients_not_in_Xenium <- meta_df2[!meta_df2$patient_id %in% patients_to_remove, ]




# Step 1: Select relevant columns from both datasets
df_IMC <- meta_df2_shared %>%
  as.data.frame() %>%
  dplyr::select(
    ObjectName = Status_simp_2,
    common_group,
    tumor_group = tumor_percentage_grouped035ext_hex_average_method2_CBM,
    Technique
  )

df_Xenium <- meta_df_shared %>%
  as.data.frame() %>%
  dplyr::select(
    ObjectName = Status_simp_2,
    common_group,
    tumor_group = tumor_percentage_grouped035ext_hex_average_method2_CBM,
    Technique
  )

# Step 2: Combine both datasets
df_combined <- bind_rows(df_IMC, df_Xenium)

# Step 3: Remove NAs or unqualified entries if needed
df_combined <- df_combined %>% filter(!is.na(tumor_group ), tumor_group  != "unqualified")

# Step 4: Calculate proportions of each common_group per sample
df_plot <- df_combined %>%
  group_by(ObjectName, Technique, common_group) %>%
  tally() %>%
  group_by(ObjectName, Technique) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()


# Factor for clusters according to your colors
df_plot$cluster_factor <- factor(df_plot$common_group, levels = names(common_group_colors))
df_plot$ObjectName_fac <- factor(df_plot$ObjectName, levels = names(status_colors)) # adjust as needed


# Define desired order for the combined bars
combined_order <- as.vector(sapply(levels(df_plot$ObjectName_fac), function(x) {
  paste0(x, "_", c("IMC", "Xenium"))  # IMC first, Xenium second
}))

# Create combined column and make it a factor with proper order
df_plot$ObjectName_Technique <- factor(
  paste0(df_plot$ObjectName_fac, "_", df_plot$Technique),
  levels = combined_order
)


#Plot
plot <- ggplot(df_plot, aes(x = ObjectName_Technique, y = proportion, fill = cluster_factor)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = common_group_colors) +
  #facet_wrap(~Technique, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 15)
  )

plot


pdf(file= file.path(output_dir,"Sup_Figure_5","IMC_and_Xenium_combined_status_barplot_shared.pdf"), 
    width = 8, height = 12, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


# Step 4: Calculate proportions of each common_group per sample
df_plot <- df_combined %>%
  group_by(tumor_group , Technique, common_group) %>%
  tally() %>%
  group_by(tumor_group , Technique) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()




# Factor for clusters according to your colors
df_plot$cluster_factor <- factor(df_plot$common_group, levels = names(common_group_colors))
df_plot$ObjectName_fac <- factor(df_plot$tumor_group , levels = names(region_colors)) # adjust as needed

# Define desired order for the combined bars
combined_order <- as.vector(sapply(levels(df_plot$ObjectName_fac), function(x) {
  paste0(x, "_", c("IMC", "Xenium"))  # IMC first, Xenium second
}))

# Create combined column and make it a factor with proper order
df_plot$ObjectName_Technique <- factor(
  paste0(df_plot$ObjectName_fac, "_", df_plot$Technique),
  levels = combined_order
)

#Plot
plot <- ggplot(df_plot, aes(x = ObjectName_Technique, y = proportion, fill = cluster_factor)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = common_group_colors) +
  #facet_wrap(~Technique, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 15)
  )

plot

pdf(file= file.path(output_dir,"Sup_Figure_5","IMC_and_Xenium_combined_tumor_grouped_barplot_shared.pdf"), 
    width = 8, height = 12, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()


df_plot$ObjectName_fac <- factor(df_plot$tumor_group , levels = rev(names(region_colors))) # adjust as needed
# Define desired order for the combined bars
combined_order <- as.vector(sapply(levels(df_plot$ObjectName_fac), function(x) {
  paste0(x, "_", c("IMC", "Xenium"))  # IMC first, Xenium second
}))

# Create combined column and make it a factor with proper order
df_plot$ObjectName_Technique <- factor(
  paste0(df_plot$ObjectName_fac, "_", df_plot$Technique),
  levels = combined_order
)

#Plot
plot <- ggplot(df_plot, aes(x = ObjectName_Technique, y = proportion, fill = cluster_factor)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = common_group_colors) +
  coord_flip() +  # flip axes
  #facet_wrap(~Technique, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 15)
  )

plot

pdf(file= file.path(output_dir,"Sup_Figure_5","IMC_and_Xenium_combined_tumor_grouped_barplot_shared_flipeed.pdf"), 
    width = 15, height = 5, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()



# Sup Figure 5 A-D Volcanoplots -----------------------------------------------




combined_DEG <- readRDS("~/DEG_not_pseudobulk.rds")

combined_DEG_sub <- combined_DEG[combined_DEG$celltype == "Stroma",]

# Code to use to plot
volcano_ctcvsbm <- make_volcano_plot(
  data        = combined_DEG_sub,
  log2fc_col  = "avg_log2FC",
  padj_col    = "p_val_adj",
  gene_col    = "gene",  # use rownames as gene names, or set to your gene column name
  up_label    = "Dense",
  down_label  = "Sparse & Normal",
  other_label = "Not significant",
  padj_cut    = 0.005,
  fc_cut      = 0.025,
  color_up    = "#D82638",
  color_down  = "#397FB9", 
  color_other = "#D9D9D9",
  title       = "analysis",
  subtitle    = "specific_analysis",
  label_all_sig = FALSE,
  # label_genes   = labels_filtered,   # only label these if significant
  label_top_n   = 40,
  max_overlaps  = 20,
  # or FALSE + label_top_n = 50
  seed        = 123
)
volcano_ctcvsbm

pdf(file= file.path(dir, file), 
    width = 8, height = 6, onefile = FALSE)
print(volcano_ctcvsbm)
dev.off()


# Sup Figure 6D CN SMM ---------------------------------------------------------------

library(pheatmap)



spe_Xenium_SMM <- mcreadRDS("~/spe_Xenium_SMM.rds" )



for_plot <- prop.table(table(spe_Xenium_SMM$cn_celltypes_12_named_SMM_correct, as.character(spe_Xenium_SMM$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv)), margin = 1)

niche_counts <- table(spe_Xenium_SMM$cn_celltypes_12_named_SMM_correct)
neighborhood_percentage <- niche_counts / sum(niche_counts) * 100
# Step 2: Create a data frame for row annotation (with the percentage of the neighborhood size)
row_annotation <- data.frame(Neighborhood_Percentage = neighborhood_percentage)
# Convert row_annotation to a named vector
neighborhood_percentage_vector <- setNames(row_annotation$Neighborhood_Percentage.Freq, row_annotation$Neighborhood_Percentage.Var1)
row_annotation_df <- data.frame(Neighborhood_Percentage = neighborhood_percentage_vector)

plot <- pheatmap( for_plot,
                  color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
                  scale = "column",
                  annotation_row = row_annotation_df,  # tilt column labels
                  border_color = NA,
                  angle_col = 45  # tilt column labels
)
plot

pdf(file= file.path(output_dir,"Sup_Figure_3","CN_SMM.pdf"), width = 12, height = 6)
print(plot)
dev.off()



# Sup Figure 6E Interactions all---------------------------------------------------------------

## interaction with new annot ifndiv:
out_all <- readRDS("~/interactions_all_per_patient.rds")

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(tibble)
library(viridis)

out_all_no_duplicates <- out_all[!(out_all$group_by %in% c( "02390_2_NDMM", "12928_2_PCL", "52998_2_NDMM","17034_2_NDMM","14719_3_CBM", "8667_3_CBM",  
                                                            "9414_3_CBM", "17336_3_CBM")), ]

#out_all_sparse_noCBM <- out_all_sparse[!(out_all_sparse$group_by %in% c( "14719_2_CBM",  "14719_3_CBM", "8667_2_CBM","8667_3_CBM","9141_2_CBM","9414_2_CBM", "9414_3_CBM" ,  "17336_3_CBM", "17336_2_CBM", "14719_2_CBM"  ,   "14719_3_CBM")), ]


out_all_tib <- as_tibble(out_all_no_duplicates)
# Step 2 — summarize interactions
df <- out_all_tib %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE), .groups = "drop")

# 2. Pivot to matrix for clustering
mat <- df %>%
  pivot_wider(names_from = to_label, values_from = sum_sigval, values_fill = 0) %>%
  column_to_rownames("from_label") %>%
  as.matrix()

# 3. Perform hierarchical clustering
row_clust <- hclust(dist(mat))
col_clust <- hclust(dist(t(mat)))

# 4. Extract ordered factor levels from clustering
row_order <- rownames(mat)[row_clust$order]
col_order <- colnames(mat)[col_clust$order]

# 5. Reorder factor levels in original df
df <- df %>%
  mutate(
    from_label = factor(from_label, levels = row_order),
    to_label = factor(to_label, levels = col_order)
  )

# --- Determine symmetric limits for diverging scale ---
limit_val <- max(abs(df$sum_sigval), na.rm = TRUE)

# --- Publication-quality heatmap ---
plot <- ggplot(df, aes(from_label, to_label, fill = sum_sigval)) +
  geom_tile(color = "grey90", size = 0.1) +
  scale_fill_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("red"),
    midpoint = 0,
    limits = c(-limit_val, limit_val),
    name = "Interaction\nstrength"
  ) +
  labs(
    x = "Sender cell type",
    y = "Receiver cell type",
    title = "Interaction analysis whole dataset"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 0, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11)
  )
plot

pdf(file= file.path(output_dir,"Sup_Figure_6","Interactions_all.pdf"), 
    width = 10, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()




# Sup Figure 6F interactions CBM---------------------------------------------------------------

## interaction with new annot ifndiv:
out_all <- readRDS("~/interactions_all_combined_per_patient.rds")

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(tibble)
library(viridis)

out_all_no_duplicates <- out_all[!(out_all$group_by %in% c( "02390_2_NDMM", "12928_2_PCL", "52998_2_NDMM","17034_2_NDMM","14719_3_CBM", "8667_3_CBM",  
                                                            "9414_3_CBM", "17336_3_CBM")), ]

out_all_CBM <- out_all_no_duplicates[out_all_no_duplicates$group_by %in% c( "14719_2_CBM",  "14719_3_CBM", "8667_2_CBM","8667_3_CBM","9141_2_CBM","9414_2_CBM", "9414_3_CBM" ,  "17336_3_CBM", "17336_2_CBM", "14719_2_CBM"  ,   "14719_3_CBM"), ]


out_all_tib <- as_tibble(out_all_CBM)
# Step 2 — summarize interactions
df <- out_all_tib %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE), .groups = "drop")

# 2. Pivot to matrix for clustering
mat <- df %>%
  pivot_wider(names_from = to_label, values_from = sum_sigval, values_fill = 0) %>%
  column_to_rownames("from_label") %>%
  as.matrix()

# 3. Perform hierarchical clustering
row_clust <- hclust(dist(mat))
col_clust <- hclust(dist(t(mat)))

# 4. Extract ordered factor levels from clustering
row_order <- rownames(mat)[row_clust$order]
col_order <- colnames(mat)[col_clust$order]

# 5. Reorder factor levels in original df
df <- df %>%
  mutate(
    from_label = factor(from_label, levels = row_order),
    to_label = factor(to_label, levels = col_order)
  )

# --- Determine symmetric limits for diverging scale ---
limit_val <- max(abs(df$sum_sigval), na.rm = TRUE)

# --- Publication-quality heatmap ---
plot <- ggplot(df, aes(from_label, to_label, fill = sum_sigval)) +
  geom_tile(color = "grey90", size = 0.1) +
  scale_fill_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("red"),
    midpoint = 0,
    limits = c(-limit_val, limit_val),
    name = "Interaction\nstrength"
  ) +
  labs(
    x = "Sender cell type",
    y = "Receiver cell type",
    title = "Interaction analysis CBM"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 0, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )
plot


pdf(file= file.path(output_dir,"Sup_Figure_4","Interactions_CBM.pdf"), 
    width = 10, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

# Sup Figure 6G  Interactions_normal_noCBM.pdf ---------------------------------------------------------------

## interaction with new annot ifndiv:
out_all <- readRDS("~/interactions_combined_wo_duplicated.rds")

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(tibble)
library(viridis)

out_all_no_duplicates <- out_all[!(out_all$group_by %in% c( "02390_2_NDMM", "12928_2_PCL", "52998_2_NDMM","17034_2_NDMM","14719_3_CBM", "8667_3_CBM",  
                                                            "9414_3_CBM", "17336_3_CBM")), ]

out_all_sparse <- out_all_no_duplicates[out_all_no_duplicates$architecture %in% c( "normal PC percentage"), ]
out_all_sparse_noCBM <- out_all_sparse[!(out_all_sparse$group_by %in% c( "14719_2_CBM",  "14719_3_CBM", "8667_2_CBM","8667_3_CBM","9141_2_CBM","9414_2_CBM", "9414_3_CBM" ,  "17336_3_CBM", "17336_2_CBM", "14719_2_CBM"  ,   "14719_3_CBM")), ]


out_all_tib <- as_tibble(out_all_sparse)
# Step 2 — summarize interactions
df <- out_all_tib %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE), .groups = "drop")

# 2. Pivot to matrix for clustering
mat <- df %>%
  pivot_wider(names_from = to_label, values_from = sum_sigval, values_fill = 0) %>%
  column_to_rownames("from_label") %>%
  as.matrix()

# 3. Perform hierarchical clustering
row_clust <- hclust(dist(mat))
col_clust <- hclust(dist(t(mat)))

# 4. Extract ordered factor levels from clustering
row_order <- rownames(mat)[row_clust$order]
col_order <- colnames(mat)[col_clust$order]

# 5. Reorder factor levels in original df
df <- df %>%
  mutate(
    from_label = factor(from_label, levels = row_order),
    to_label = factor(to_label, levels = col_order)
  )

# --- Determine symmetric limits for diverging scale ---
limit_val <- max(abs(df$sum_sigval), na.rm = TRUE)

# --- Publication-quality heatmap ---
plot <- ggplot(df, aes(from_label, to_label, fill = sum_sigval)) +
  geom_tile(color = "grey90", size = 0.1) +
  scale_fill_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("red"),
    midpoint = 0,
    limits = c(-limit_val, limit_val),
    name = "Interaction\nstrength"
  ) +
  labs(
    x = "Sender cell type",
    y = "Receiver cell type",
    title = "Interaction analysis only normal PC percentage (no CBM)"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 0, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11)
  )
plot

pdf(file= file.path(output_dir,"Sup_Figure_4","Interactions_normal_noCBM.pdf"), 
    width = 10, height = 10, useDingbats = FALSE, onefile = FALSE)
print(plot)
dev.off()

