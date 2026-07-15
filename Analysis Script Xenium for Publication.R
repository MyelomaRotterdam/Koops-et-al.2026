
##### ANALYSIS PIPELINE # [ Spatial profiling identifies a distinct and topographically-defined tumor microenvironment that emerges during multiple myeloma evolution ] 



# Read in packages --------------------------------------------------------

pkgs <- unique(c(
  "qs", "Seurat", "imcRtools", "SpatialExperiment", "ggplot2", "dplyr", "BiocStyle", 
  "BiocManager", "remotes", "cytomapper", "sf", "scuttle", "scran", "scater", "uwot",
  "rtracklayer", "patchwork", "pheatmap", "gridExtra", "Rphenograph", "igraph", "dittoSeq",
  "viridis", "bluster", "BiocParallel", "readxl", "readr", "hexbin", "openxlsx", "ggrastr",
  "tibble", "grid", "speckle", "limma", "edgeR", "tidyr", "cowplot", "scales", "RColorBrewer", 
  "stringr", "harmony", "data.table", "arrow"
  
))

# Load all packages
invisible(lapply(pkgs, library, character.only = TRUE))

# Analysis Spatial Transcriptomics Xenium Data ----------------------------


# ---- read in xenium object from scratch ---- 
data.dir <- "~/output-_SMM_MM//"
data.dir <- "~/output-_relapse_MM/"
data.dir <- "~/output-_HBM/"
data.dir <- "~/output-_MM/"

# create CSV file with transcripts

transcripts <- read_parquet(file.path(data.dir, "transcripts.parquet"))
write.csv(transcripts, gzfile(file.path(data.dir, "transcripts.csv.gz")), row.names = FALSE)


list.files <- list.files("~/20250612__085636__202501628B-Hema-12062025/")
data.dir.start <- "~/20250612__085636__202501628B-Hema-12062025/"


for(file in list.files){
  data.dir <- file.path(data.dir.start, file)
  transcripts <- read_parquet(file.path(data.dir, "transcripts.parquet"))
  write.csv(transcripts, gzfile(file.path(data.dir, "transcripts.csv.gz")), row.names = FALSE)
}


# Load the matrix files manually

# Load transcript data
transcript_data <- fread( "~/_SMM_MM//transcripts.csv.gz")
transcript_data_raw <- read_parquet("~/_relapse_MM/transcripts.parquet")

cell_feature_data <- Read10X_h5("~/_SMM_MM/cell_feature_matrix.h5")

total_transcripts <- nrow(transcript_data_raw)
high_quality_transcripts <- sum(transcript_data_raw$qv >= 20)
percentage_high_quality <- (high_quality_transcripts / total_transcripts) * 100

percentage_high_quality


########### manual way to get folder structure, skip if not applicable
data.dir.start <- "~/Run_3"
path <- "~/Run_3"

path <- normalizePath(path)
# Step 1: get immediate subfolders of path
subdirs_lvl1 <- list.dirs(path = path, recursive = FALSE, full.names = TRUE)
subdirs_lvl1 <- subdirs_lvl1[subdirs_lvl1 != path]
# Step 2: for each of those, get only their immediate subfolders (level 2)
list.files <- unlist(lapply(subdirs_lvl1, function(p) {
  sub <- list.dirs(path = p, recursive = FALSE, full.names = TRUE)
  sub[sub != p]
}))
###


xenium.list2 <- list()
names <- c("HBM","MM","SMM_MM", "Relapse_MM")

# Iterate over files
for (i in seq_along(list.files)) {
  # Current name and file
  name <- names[i]
  file <- list.files[i]
  #data.dir <- file.path(data.dir.start, file)
  data <- ReadXenium(
    data.dir = file,
    type = c("centroids", "segmentations"),
  )
  
  assay <- "Xenium"
  segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
  )
  
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay
  )
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  
  data$matrix[["Unassigned Codeword"]] <- data$matrix[["Negative Control Codeword"]]
  data$matrix[["Unassigned Codeword"]][data$matrix[["Unassigned Codeword"]] != 0] <- 0
  xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  
  fov <- "test"
  xenium.obj[[fov]] <- coords
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
  xenium.list2[[name]] <- xenium.obj
}

saveRDS(xenium.list, "~/RDS_files/xenium_list_run3_HBM_andallMM.rds" )




# ---- Divide cells in tissues, use Xenium Explorer to create CSV files with cells per tissue (if you have multiple tissues per slide) ---- 



list.files <- list.files("~/cell_stats_from_explorer/")
# Assuming list.files is your vector of filenames:
files <- list.files
# Folder path where the files are located
folder_path <- "~/cell_stats_from_explorer/"

# 1. Read all cell lists with tissue IDs extracted from filenames
cell_lists <- lapply(files, function(f) {
  cells <- read.csv(file.path(folder_path, f), header = FALSE)[, 1]
  cells[-c(1, 2, 3)]  # Remove first 3 rows, as you did before
})

# 2. Extract tissue IDs from filenames (e.g. "02390" from "02390_cells_stats.csv")
tissue_ids <- sub("_cells_stats.csv", "", files)

# 3. Create a named list for convenience
names(cell_lists) <- tissue_ids

# Vector of tissue IDs corresponding to each object in the list
# Make sure the order matches your list or use named vector
tissue_ids_per_object <- list(
  MM = c("1", "2","3","4","5"),
  HBM = c("1", "2", "3", "4"),
  SMM_MM = c("1", "2", "3", "4", "5", "6"),
  Relapse_MM = c("1", "2", "3")
)

for (sample_name in names(xenium.list)) {
  seurat_obj <- xenium.list[[sample_name]]
  
  raw_cellnames <- rownames(seurat_obj@meta.data)
  
  # Remove everything up to and including the last underscore
  cleaned_cellnames <- sub("^.*_", "", raw_cellnames)
  
  seurat_obj$CleanCellName <- cleaned_cellnames
  
  xenium.list[[sample_name]] <- seurat_obj
}


# 4. Function to assign tissue info to a Seurat object, given the cell lists and tissue IDs to assign
assign_tissue <- function(seurat_obj, cell_lists, tissue_ids) {
  #seurat_obj$CellNames <- colnames(seurat_obj)  # ensure CellNames column exists
  seurat_obj$Tissue <- NA                       # initialize Tissue column
  
  for (tid in tissue_ids) {
    if (tid %in% names(cell_lists)) {
      matching_cells <- seurat_obj$CleanCellName %in% cell_lists[[tid]]
      seurat_obj$Tissue[matching_cells] <- tid
    }
  }
  seurat_obj
}

# Loop over your list, assign tissue IDs, and update list elements
for (obj_name in names(xenium.list)) {
  tissues <- tissue_ids_per_object[[obj_name]]
  xenium.list[[obj_name]] <- assign_tissue(xenium.list[[obj_name]], cell_lists, tissues)
}


# 1. Rename cells to include sample name (from the list names) to make them unique
for (i in names(xenium.list)) {
  xenium.list[[i]] <- RenameCells(xenium.list[[i]], add.cell.id = i)
}



# Remove 2 genes, so Run 2 and Run 3 can be merged ------------------------


# Genes to exclude
exclude_genes <- c("FLT3LG", "IFI6") # these genes had unspecific binding in Run 2



# Function to remove genes from a Seurat object
remove_genes_from_seurat <- function(seurat_obj, genes_to_remove) {
  for (assay_name in names(seurat_obj@assays)) {
    all_genes <- rownames(seurat_obj[[assay_name]])
    genes_to_keep <- setdiff(all_genes, genes_to_remove)
    seurat_obj[[assay_name]] <- subset(seurat_obj[[assay_name]], features = genes_to_keep)
  }
  return(seurat_obj)
}

# Apply to previous Xenium list
prev_xenium_list_filtered <- lapply(xenium_list_run2, remove_genes_from_seurat, genes_to_remove = exclude_genes)

# Apply to new Xenium list
new_xenium_list_filtered <- lapply(xenium.list, remove_genes_from_seurat, genes_to_remove = exclude_genes)

for (obj_name in names(prev_xenium_list_filtered)) {
  seurat_obj <- prev_xenium_list_filtered[[obj_name]]
  
  # Check if 'nCount_RNA' or 'nCount' is present; adjust accordingly
  if ("nCount_Xenium" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$high_nCount <- seurat_obj$nCount_Xenium > 8
  } else if ("nCount" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$high_nCount <- seurat_obj$nCount_Xenium > 8
  } else {
    stop(paste("No nCount column found in", obj_name))
  }
  
  # Save back to the list
  prev_xenium_list_filtered[[obj_name]] <- seurat_obj
}


# QC ----------------------------------------------------------------------



##### beofore merge, loose spatial data, remove nCount < 8  and run SCT on each tissue separately. This will help normalize batch control and help doing DEG later

for (sample_name in names(xenium.list)) {
  seurat_obj <- xenium.list[[sample_name]]
  
  # Drop spatial data to free memory
  seurat_obj@images <- list()
  if ("Spatial" %in% names(seurat_obj@misc)) {
    seurat_obj@misc$Spatial <- NULL
  }
  
  tissues <- unique(seurat_obj$Tissue)
  tissues <- tissues[!is.na(tissues)]
  
  for (tissue in tissues) {
    seurat_sub <- subset(seurat_obj, subset = Tissue == tissue)
    
    counts <- GetAssayData(seurat_sub, assay = "Xenium", slot = "counts")
    new_seurat <- CreateSeuratObject(counts = counts, meta.data = seurat_sub@meta.data)
    new_seurat <- subset(new_seurat, subset = nCount_Xenium > 8) # excluding FLT3LG and IFI6
    new_seurat <- SCTransform(new_seurat, assay = "RNA", verbose = FALSE)
    
    new_name <- paste(sample_name, tissue, sep = "_")
    seurat_list_split[[new_name]] <- new_seurat
  }
}


# Merge -------------------------------------------------------------------



# Merge all new Seurat objects (now cleaned and SCT-normalized)
merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_list_split)
library(pbapply)
# Initialize the merged object as the first in the list
merged_seurat <- seurat_list_split[[1]]
# Merge the rest one by one with a progress bar
for (i in pbsapply(2:length(seurat_list_split), function(x) x)) {
  merged_seurat <- merge(merged_seurat, seurat_list_split[[i]])
}

library(pbapply)
start_idx <- 11
end_idx <- 17
last7_list <- seurat_list_split[start_idx:end_idx]

# Merge last 7 at once with progress bar
merged_last7 <- Reduce(function(x, y) merge(x, y), last7_list)

# Now merge the result with the already merged first 10
merged_seurat_all <- merge(merged_seurat, merged_last7)

saveRDS(merged_seurat_all , "~/seurat_merge_after_sct.rds" )


# Integration -------------------------------------------------------------


DefaultAssay(merged_seurat_all) <- "SCT"
cat("Running PCA...\n")
merged_seurat_all <- RunPCA(merged_seurat_all, assay = "SCT", verbose = T, features = rownames(merged_seurat_all))

cat("Running Harmony...\n")
integrated_seurat_all <- RunHarmony(
  object = merged_seurat_all,
  group.by.vars = c("Tissue"),  # Batch metadata column
  reduction = "pca",
  assay.use = "SCT",  # Specify SCT assay explicitly
  reduction.save = "harmony"
)
## on integrated object
ElbowPlot(integrated_seurat_all, n = 30)
cat("Running UMAP...\n")
integrated_seurat_all <- RunUMAP(integrated_seurat_all, dims = 1:24, reduction.name = "harmony_umap", reduction = "harmony")
integrated_seurat_all <- FindNeighbors(integrated_seurat_all, reduction = "harmony", dims = 1:24)
integrated_seurat_all <- FindClusters(integrated_seurat_all, resolution = 1, cluster.name = "clusters_res1", algorithm = 4)

cat("Saving final object...\n")
saveRDS(integrated_seurat_all , "~/RDS_files/seurat_integrated.rds" )



# FindAllMarkers ----------------------------------------------------------


Idents(integrated_seurat_all  ) <- integrated_seurat_all  $clusters_res1
markers <- FindAllMarkers(integrated_seurat_all , assay = "SCT", only.pos = TRUE)
write.csv(markers , "markers/markers.csv")


# Annotation --------------------------------------------------------------

# Annotation was performed using in house scRNAseq data and literature. 
# For annotation, download the public available barcodes and link to data



# Subclustering -----------------------------------------------------------

## subclustering was done on myeloid cells, T cells and structural cells
## This is an example, with the structural cells

structural_genes <- read.csv("annotated_markers.csv")
structural_genes <- structural_genes[structural_genes$cluster %in% c( "Stromal cells", "Smooth muscle", "Endothelial cells"),]
structural_genes <- structural_genes$gene

structural_cells <- subset(seurat_reference, subset = annot_mostprevalent %in% c("Stromal cells", "Smooth muscle", "Endothelial cells"))

# Run SCTransform on structural cells
structural_cells <- SCTransform(structural_cells)

# Run PCA using structural genes to reduce noise from genes that are not relevant for structural cells, but are influencing because of spillover
structural_cells <- RunPCA(structural_cells, npcs = 50, features = structural_genes, reduction.name = "pca_structural")

ElbowPlot(structural_cells, n = 50)

structural_cells <- RunHarmony(
  object = structural_cells,
  group.by.vars = "Tissue_ext",  # Batch metadata column
  reduction = "pca_structural",
  assay.use = "SCT",  # Specify SCT assay explicitly
  reduction.save = "harmony_structural"
)

# Run UMAP on harmony reduction
structural_cells <- RunUMAP(structural_cells, dims = 1:22, reduction.name = "harmony_umap_structural", reduction = "harmony_structural")
structural_cells <- FindNeighbors(structural_cells, reduction = "harmony_structural", dims = 1:22)
structural_cells <- FindClusters(structural_cells, resolution = 0.6, cluster.name = "cluster_structural_res0.6", algorithm = 4)

saveRDS(structural_cells, "~/RDS_files/structural_cells_run23.rds")

structural_cells <- PrepSCTFindMarkers(structural_cells)
# Find markers per cluster
Idents(structural_cells) <- "cluster_structural_res0.6"
markers <- FindAllMarkers(structural_cells, assay = "SCT", only.pos = TRUE)
write.csv(markers, "structural_cells_run2and3.csv")

# SPATIAL ANALYSES --------------------------------------------------------


# Transfer Spatial data to seurat object ----------------------------------



# Define offsets for slides (assuming names of xenium.list match these keys)
slide_offsets <- c("HBM" = 0, "MM" = 23000, "SMM_MM" = 46000, "Relapse_MM" = 69000)

coords_list <- list()

for (slide_name in names(xenium_list_run3 )) {
  xen_obj <- xenium_list_run3 [[slide_name]]
  
  coords <- GetTissueCoordinates(xen_obj)
  
  # coords rownames must match cell IDs
  coords <- coords %>% as.data.frame()
  coords$CellNames <- coords$cell
  
  # Add offset
  coords <- coords %>%
    mutate(
      x_corrected = x + slide_offsets[slide_name],
      y_corrected = y + slide_offsets[slide_name]
    )
  
  coords_list[[slide_name]] <- coords %>%
    select(CellNames, x_corrected, y_corrected)
}

all_coords <- bind_rows(coords_list)

# Join to Seurat metadata
integrated_seurat_all <- AddMetaData(
  object = merged_run3_annotated,
  metadata = all_coords %>% column_to_rownames("CellNames")
)

integrated_seurat_all <- transfer_metadata_columns_flex(
  from_obj = all_coords,
  to_obj = integrated_seurat_all,
  columns_to_transfer = c("x_corrected", "y_corrected"),
  link_from = "CellNames",
  link_to = "CellName_run3"
)


slide_offsets <- c("MM" = 92000, "HBM" = 115000, "MM_02" = 138000, "HBM_03" = 161000)

coords_list <- list()

for (slide_name in names(xenium_list_run2 )) {
  xen_obj <- xenium_list_run2 [[slide_name]]
  
  coords <- GetTissueCoordinates(xen_obj)
  
  # coords rownames must match cell IDs
  coords <- coords %>% as.data.frame()
  coords$CellNames <- coords$cell
  
  # Add offset
  coords <- coords %>%
    mutate(
      x_corrected = x + slide_offsets[slide_name],
      y_corrected = y + slide_offsets[slide_name]
    )
  coords$Glass_slide <- slide_name
  coords$Cellname_Run2 <- paste0(coords$Glass_slide, "_", coords$CellNames)
  
  coords_list[[slide_name]] <- coords %>%
    select(Cellname_Run2, x_corrected, y_corrected)
}

all_coords2 <- bind_rows(coords_list)
colnames(all_coords2) <- c("Cellname_Run2", "x_corrected_run2", "y_corrected_run2")



# Hexogonal grid  ----------------------------------



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


### i will skip this and do it on the merge
# --- Step 3: Loop over each Seurat object in the list ---
for (obj_name in names(xenium.list)) {
  seurat_obj <- xenium.list[[obj_name]]
  
  # Extract cell centroid coordinates (adjust if your data structure is different)
  centroids <- seurat_obj@images[["test"]]@boundaries[["centroids"]]@coords
  centroids_df <- as.data.frame(centroids)
  cell_points <- st_as_sf(centroids_df, coords = c("x", "y"), crs = NA)
  
  # --- Step 4: Loop over offsets to create shifted hex grids and assign region labels ---
  for (i in seq_along(offsets)) {
    dx <- offsets[[i]][1]
    dy <- offsets[[i]][2]
    
    # Create hexagonal grid over fixed bounding box
    hex_grid <- st_make_grid(bbox, cellsize = 125, square = FALSE, what = "polygons")
    
    # Shift the grid polygons by offset
    hex_grid_shifted <- st_geometry(hex_grid) + c(dx, dy)
    
    # Convert shifted grid back to sf and assign region names
    hex_grid_sf <- st_sf(geometry = hex_grid_shifted)
    hex_grid_sf$region_name <- paste0("region_", seq_len(nrow(hex_grid_sf)))
    
    # Spatial join: assign each cell point to a hex region polygon
    joined <- st_join(cell_points, hex_grid_sf, left = TRUE)
    
    # Define metadata column name for this grid
    region_colname <- paste0("small_region_hex", i)
    
    # Add the region assignment to Seurat metadata
    seurat_obj[[region_colname]] <- joined$region_name
    
    # Optional: print distribution summary for this sample and grid
    message(paste("Sample:", obj_name, "- Hex grid", i))
    print(table(seurat_obj[[region_colname]]))
  }
  
  # Save updated Seurat object back to the list
  xenium.list[[obj_name]] <- seurat_obj
}


# Assuming you have your bounding box 'bbox' and offsets list defined already

# Extract centroid coordinates from metadata
centroids_df <- integrated_seurat_all @meta.data[, c("x_corrected_together", "y_corrected_together")]
cell_points <- st_as_sf(centroids_df, coords = c("x_corrected_together", "y_corrected_together"), crs = NA)

for (i in seq_along(offsets)) {
  dx <- offsets[[i]][1]
  dy <- offsets[[i]][2]
  
  # Create hexagonal grid over fixed bounding box
  hex_grid <- st_make_grid(bbox, cellsize = 125, square = FALSE, what = "polygons")
  
  # Shift the grid polygons by offset
  hex_grid_shifted <- st_geometry(hex_grid) + c(dx, dy)
  
  # Convert shifted grid back to sf and assign region names
  hex_grid_sf <- st_sf(geometry = hex_grid_shifted)
  hex_grid_sf$region_name <- paste0("region_", seq_len(nrow(hex_grid_sf)))
  
  # Spatial join: assign each cell point to a hex region polygon
  joined <- st_join(cell_points, hex_grid_sf, left = TRUE)
  
  # Define metadata column name for this grid
  region_colname <- paste0("small_region_hex", i)
  
  # Add the region assignment to Seurat metadata
  integrated_seurat_all[[region_colname]] <- joined$region_name
  
  # Optional: print distribution summary for this grid
  message(paste("Hex grid", i))
  print(table(integrated_seurat_all[[region_colname]]))
}

#tablecellsperregions <- table(integrated_seurat_all$small_region_hex1, integrated_seurat_all$annot_broad)
tablecellsperregions <- table(integrated_seurat_all$small_region_hex1, integrated_seurat_all$annot_merged_final_with_macro)
tablecellsperregions <- as.data.frame(tablecellsperregions)
tablecellsperregions <- tablecellsperregions %>%
  group_by(Var1) %>%
  mutate(total_counts = sum(Freq))
tablecellsperregions <- tablecellsperregions %>%
  mutate(percentage = Freq/total_counts)

filtered_tumor <- tablecellsperregions[tablecellsperregions$Var2 %in% c("Plasma cells", "Cycling B/Plasma cells"), ]

integrated_seurat_all$tumor_freq <- filtered_tumor$percentage[match(integrated_seurat_all$small_region_hex1, filtered_tumor$Var1)]

hist(filtered_tumor$percentage,
     main = "Distribution of Tumor Frequencies",
     xlab = "Frequency",
     ylab = "Count",
     col = "skyblue",  # Color of bars
     border = "white",
     breaks = 20# Border color of bars
)   
regions_amountofcells <- table(integrated_seurat_all$small_region_hex1)


hist(regions_amountofcells,
     main = "Cell_counts",
     xlab = "Cells",
     ylab = "Region_Count",
     col = "skyblue",  # Color of bars
     border = "white",
     breaks = 50# Border color of bars
)


regions_amountofcells <- as.data.frame(regions_amountofcells)
regionswithnotenoughcells <- regions_amountofcells[regions_amountofcells$Freq < 25,] #30 for 250 by 250, 15 for 125 by 125
regionswithnotenoughcells <- as.character(regionswithnotenoughcells$Var1)

# Select cells with region names matching the desired ones
Idents(integrated_seurat_all) <- "small_region_hex1"
selected_cells <- WhichCells(integrated_seurat_all, ident = regionswithnotenoughcells)

integrated_seurat_all$regionselected <- "qualified"

integrated_seurat_all$regionselected[selected_cells] <- "disqualified"



integrated_seurat_all$tumor_percentage_grouped0.35 <- ifelse(integrated_seurat_all$regionselected == "disqualified", "unqualified",
                                                             ifelse(integrated_seurat_all$tumor_freq > 0.35, "dense",
                                                                    ifelse(integrated_seurat_all$tumor_freq > 0.2, "dispersed","scarce")))



# Tumor regions for all 6 grids -------------------------------------------


hex_columns <- paste0("small_region_hex", 1:6)
meta <- integrated_seurat_all@meta.data

#this works for only plasma cells
for (hex_col in hex_columns) {
  message("Processing ", hex_col)
  
  hex_vec <- as.character(meta[[hex_col]])
  annot_vec <- as.character(meta$annot_merged_final_with_macro)
  
  stopifnot(length(hex_vec) == length(annot_vec))
  
  tab <- table(hex_vec, annot_vec)
  tab_df <- as.data.frame(tab) %>%
    group_by(hex_vec) %>%
    mutate(total_counts = sum(Freq),
           percentage = Freq / total_counts)
  
  filtered <- tab_df[tab_df$annot_vec %in% c("Plasma cells"), ] 
  
  perc_name <- paste0("tumor_freq_ext", hex_col)
  meta[[perc_name]] <- filtered$percentage[match(hex_vec, filtered$hex_vec)]
  
  region_counts <- table(hex_vec)
  region_counts_df <- as.data.frame(region_counts)
  disqualified <- as.character(region_counts_df$hex_vec[region_counts_df$Freq < 25])
  
  status_col <- paste0("regionselected_", hex_col)
  meta[[status_col]] <- ifelse(hex_vec %in% disqualified, "disqualified", "qualified")
  
  group_col <- paste0("tumor_percentage_grouped_035_ext", hex_col)
  meta[[group_col]] <- ifelse(meta[[status_col]] == "disqualified", "unqualified",
                              ifelse(meta[[perc_name]] > 0.50, "extremely_dense",
                                     ifelse(meta[[perc_name]] > 0.35, "dense",
                                            ifelse(meta[[perc_name]] > 0.2, "dispersed",
                                                   ifelse(meta[[perc_name]] > 0.05, "sparse", "normal PC percentage")))))
}

# Write back
integrated_seurat_all@meta.data <- meta


hex_columns <- paste0("small_region_hex", 1:6)
meta <- integrated_seurat_all@meta.data
# this works for plasma cells and cycling and b cells
for (hex_col in hex_columns) {
  message("Processing ", hex_col)
  
  hex_vec <- as.character(meta[[hex_col]])
  annot_vec <- as.character(meta$annot_merged_final_with_macro)
  
  stopifnot(length(hex_vec) == length(annot_vec))
  
  # Frequency table and conversion to percentages
  tab_df <- as.data.frame(table(hex = hex_vec, celltype = annot_vec)) %>%
    group_by(hex) %>%
    mutate(total_counts = sum(Freq),
           percentage = Freq / total_counts)
  
  # Filter for tumor-related cell types and sum their percentages per hex
  tumor_freq_df <- tab_df %>%
    filter(celltype %in% c("Plasma cells",  "Cycling B/Plasma cells")) %>%
    group_by(hex) %>%
    summarise(tumor_freq = sum(percentage), .groups = "drop")
  
  perc_name <- paste0("tumor_freq_ext", hex_col)
  meta[[perc_name]] <- tumor_freq_df$tumor_freq[match(hex_vec, tumor_freq_df$hex)]
  
  # Region qualification
  region_counts <- table(hex_vec)
  region_counts_df <- as.data.frame(region_counts)
  disqualified <- as.character(region_counts_df$hex_vec[region_counts_df$Freq < 25])
  
  status_col <- paste0("regionselected_", hex_col)
  meta[[status_col]] <- ifelse(hex_vec %in% disqualified, "disqualified", "qualified")
  
  group_col <- paste0("tumor_percentage_grouped_035_ext", hex_col)
  meta[[group_col]] <- ifelse(meta[[status_col]] == "disqualified", "unqualified",
                              ifelse(meta[[perc_name]] > 0.50, "extremely_dense",
                                     ifelse(meta[[perc_name]] > 0.35, "dense",
                                            ifelse(meta[[perc_name]] > 0.2, "dispersed",
                                                   ifelse(meta[[perc_name]] > 0.05, "sparse", "normal PC percentage")))))
}

# Write back to Seurat object
integrated_seurat_all@meta.data <- meta

tumor_cells <-  c("Plasma cells", "Cycling B/Plasma cells")
# Compute percentage of tumor cells per tissue
meta <- integrated_seurat_all@meta.data
# Calculate tumor frequency per tissue/sample
tumor_freq_df <- meta %>%
  group_by(Tissue_ext) %>%
  summarise(tumor_freq_total = mean(annot_merged_final_with_macro_noslash %in% tumor_cells)) 
# Join back to metadata by Tissue_ext
meta <- meta %>%
  left_join(tumor_freq_df, by = "Tissue_ext")
# Add to Seurat object's metadata
integrated_seurat_all$tumor_freq_total <- meta$tumor_freq_total
rm(meta)

# Combine all 6 grid classifications into an average call
group_cols <- paste0("tumor_percentage_grouped_035_ext", hex_columns)
df <- integrated_seurat_all@meta.data
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
df$most_frequent <- apply(df[, group_cols], 1, Mode)
integrated_seurat_all$tumor_percentage_grouped035ext_hex_average <- df$most_frequent
# average tumor freq 
integrated_seurat_all$tumor_freq_ext_avg <- rowMeans(integrated_seurat_all@meta.data[, paste0("tumor_freq_extsmall_region_hex", 1:6)], na.rm = TRUE)
# Vector of qualification status column names
qual_cols <- paste0("regionselected_", hex_columns)
# Score: 1 point per "qualified"
# Compute score in the metadata
integrated_seurat_all$qualified_score <- apply(
  integrated_seurat_all@meta.data[qual_cols],
  1,
  function(x) sum(x == "qualified")
)


integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2_with_extreme <-
  ifelse(integrated_seurat_all$qualified_score <  4, "unqualified",
         ifelse(integrated_seurat_all$tumor_freq_ext_avg > 0.50, "extremely_dense",
                ifelse(integrated_seurat_all$tumor_freq_ext_avg > 0.35, "dense",
                       ifelse(integrated_seurat_all$tumor_freq_ext_avg > 0.2, "dispersed",
                              ifelse(integrated_seurat_all$tumor_freq_ext_avg > 0.05, "sparse", "normal PC percentage")))))

integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2 <- integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2_with_extreme
integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2[integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2_with_extreme == "extremely_dense"] = "dense"
integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2_CBM <- integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2
integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2_CBM[integrated_seurat_all$Status_simp == "CBM"] <- "CBM"


integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_withextreme <- integrated_seurat_all$tumor_percentage_grouped035ext_hex_average
integrated_seurat_all$tumor_percentage_grouped035ext_hex_average[integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_withextreme == "extremely_dense"] = "dense"

integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_CBM <- integrated_seurat_all$tumor_percentage_grouped035ext_hex_average
integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_CBM[integrated_seurat_all$Status_simp == "CBM"] <- "CBM"





# Transform to SpatialExperiment for further analysis --------------------------------------------------------


# Extract SCT counts
counts <- GetAssayData(seurat_obj, assay = "SCT", slot = "counts")

# Prepare metadata and gene info
coldata <- seurat_obj@meta.data
rowdata <- data.frame(gene_id = rownames(counts))

# Extract spatial coordinates from metadata
spatial_coords <- as.matrix(coldata[, c("x", "y")])
rownames(spatial_coords) <- rownames(coldata)

# Create the SpatialExperiment object
spe <- SpatialExperiment(
  assays = list(counts = counts),
  colData = coldata,
  rowData = rowdata,
  spatialCoords = spatial_coords
)




# Percent PC interacting with all celltypes -------------------------------


# Extract edges
g <- spe_Xenium@int_colData$colPairs$delaunay_interaction_graph
hits <- g@hits
cell_ids <- colnames(spe_Xenium)
from_cells <- cell_ids[hits@from]
to_cells   <- cell_ids[hits@to]
edges_df <- data.frame(from = from_cells, to = to_cells, stringsAsFactors = FALSE)

# Metadata
meta_df <- colData(spe_Xenium)

# Unique patients
patients <- unique(meta_df$ObjectName_anonymous_grouped)

# Unique cell types (excluding Plasma Cells)
celltypes <- setdiff(unique(meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv), "Plasma Cell")

# Initialize list to store results
results_list <- list()

for (p in patients) {
  
  plasma_ids <- colnames(spe_Xenium)[meta_df$ObjectName_anonymous_grouped == p & meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv %in% c("Plasma Cell")]
  
  if (length(plasma_ids) == 0) next
  
  for (ct in celltypes) {
    
    target_ids <- colnames(spe_Xenium)[meta_df$ObjectName_anonymous_grouped == p & meta_df$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv == ct]
    
    if (length(target_ids) == 0) next
    
    interacting_plasma <- unique(c(
      edges_df$from[edges_df$from %in% plasma_ids & edges_df$to %in% target_ids],
      edges_df$to[edges_df$to %in% plasma_ids & edges_df$from %in% target_ids]
    ))
    
    percent <- length(interacting_plasma) / length(plasma_ids) * 100
    
    results_list[[length(results_list) + 1]] <- data.frame(
      patient = p,
      celltype = ct,
      percent = percent
    )
  }
}


results_list

# Combine all results
results_df_annot_precise <- bind_rows(results_list)

saveRDS(results_df_annot_precise, "results_df_annot_precise_ifndiv.rds")
saveRDS(results_df, "results_df.rds")

results_df <- readRDS("results_df.rds")


results_df <- results_df %>%
  mutate(Disease = sapply(strsplit(as.character(patient), "_"), function(x) x[3]))

order_tumor_freq <- integrated_seurat_all@meta.data %>%
  group_by(Tissue_ext_spe) %>%
  summarise(tumor_freq_total = mean(tumor_freq_total)) %>%
  arrange(desc(tumor_freq_total))
order_tumor_freq <- order_tumor_freq$Tissue_ext_spe


results_df$patient <- factor(results_df$patient, levels = order_tumor_freq)


#### faster code!!
# Metadata
meta_df <- colData(spe_Xenium)


# Metadata
meta_dt <- data.table(
  cell_id = colnames(spe_Xenium),
  patient = meta_df$ObjectName_anonymous_grouped,
  patient = meta_df$Status,
  celltype = meta_df$annot_broad_intbetween_run3_correct_ifndiv,
  architecture = meta_df$tumor_percentage_grouped035ext_hex_average_method2
)


dir.create("interaction_results_architecture", showWarnings = FALSE)

patients <- unique(meta_dt$patient)
architectures <- unique(meta_dt$architecture)

### plasma cell interacting with other cells
for (p in patients) {
  for (arch in architectures) {
    
    outfile <- file.path("interaction_results_architecture", paste0("interactions_broader", p, "_", arch, ".rds"))
    if (file.exists(outfile)) {
      message("Skipping ", p, " - ", arch)
      next
    }
    
    message("Processing ", p, " - ", arch, " ...")
    
    # Subset metadata for this patient x architecture
    meta_sub <- meta_dt[patient == p & architecture == arch]
    if (nrow(meta_sub) == 0) next
    
    # plasma cells in this architecture
    plasma_ids <- meta_sub[celltype == "B/Plasma Cell", cell_id] #Plasma Cell
    if (length(plasma_ids) == 0) next
    
    # target cells per cell type
    target_dt <- meta_sub[celltype != "B/Plasma Cell", .(target_ids = list(cell_id)), by = celltype] #Plasma Cell
    if (nrow(target_dt) == 0) next
    
    # --- Filter edges: only edges where 'from' is a Plasma cell ---
    local_edges <- edges_df[edges_df$from %in% plasma_ids & edges_df$to %in% meta_sub$cell_id, ]
    
    # --- Safety checks ---
    if (!all(local_edges$from %in% plasma_ids)) stop("Some 'from' cells are not Plasma_MM")
    mirrored_edges <- edges_df[edges_df$to %in% plasma_ids & !edges_df$from %in% plasma_ids, ]
    if (nrow(mirrored_edges) > 0) warning("There are mirrored edges that should not be counted")
    
    # Filter edges: only Plasma → Plasma
    local_edges_pc_pc <- local_edges[local_edges$to %in% plasma_ids, ]
    
    # Build neighbor list (only Plasma → other)
    neighbors <- split(local_edges$to, local_edges$from)
    neighbors <- lapply(neighbors, unique)
    
    plasma_in_graph <- plasma_ids[plasma_ids %in% names(neighbors)]
    if (length(plasma_in_graph) == 0) {
      target_dt[, percent := 0]
      target_dt[, `:=`(patient = p, architecture = arch)]
      saveRDS(target_dt[, .(patient, architecture, celltype, percent)], outfile)
      message("No plasma edges for ", p, " - ", arch)
      next
    }
    
    # Compute percent of interacting plasma cells per cell type
    target_dt[, percent := sapply(target_ids, function(target_vec) {
      interacting <- vapply(plasma_in_graph, function(pc) {
        any(neighbors[[pc]] %in% target_vec)
      }, logical(1))
      mean(interacting) * 100
    })]
    
    # Count total number of plasma neighbors per plasma cell
    pc_total_pc_neighbors <- data.table(
      plasma_id = names(neighbors_pc_pc),
      n_neighbors_pc = vapply(neighbors_pc_pc, length, integer(1))
    )
    
    
    target_dt[, patient := p]
    target_dt[, architecture := arch]
    
    # Save result
    saveRDS(target_dt[, .(patient, architecture, celltype, percent)], outfile)
    message("Done ", p, " - ", arch, " (", nrow(target_dt), " cell types)")
  }
}

files <- list.files("~/RDS_files/Interactions/interaction_absolute_results_architecture/Broader_annot_neutro_combined/", full.names = TRUE)
results_df_annot_broad <- rbindlist(lapply(files, readRDS), fill = TRUE)

# Filter out unqualified architectures and order cell types by max percent (optional)
plot_data <- results_df_annot_broad %>%
  filter(architecture != "unqualified") %>%
  group_by(celltype) %>%
  mutate(max_percent = max(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(celltype = factor(celltype, levels = unique(celltype[order(-max_percent)])))

# Plot: cell type on x-axis, bars colored by architecture, facet per patient
ggplot(plot_data, aes(x = celltype, y = percent, fill = architecture)) +
  geom_col(position = "dodge") +
  facet_wrap(~ patient, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Plasma Cell Interactions per Cell Type and Architecture",
       x = "Cell Type",
       y = "Percentage of Plasma Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")



# 1️⃣ Add Disease Stage column
results_df_annot_broad <- results_df_annot_broad %>%
  mutate(DiseaseStage = sapply(strsplit(patient, "_"), `[`, 1))

# 2️⃣ Filter out unqualified architectures
plot_data <- results_df_annot_broad %>%
  filter(architecture != "unqualified")

saveRDS(plot_data,"~/RDS_files/Interactions/interaction_absolute_results_architecture/Broader_annot_neutro_combined/combined.rds")

plot_data <- readRDS("~/RDS_files/Interactions/interaction_absolute_results_architecture/Broader_annot_neutro_combined/combined.rds")

# 3️⃣ Get all cell types
celltypes <- unique(plot_data$celltype)

# Compute mean percent per DiseaseStage / Architecture / Celltype
summary_dt <- plot_data %>%
  group_by(celltype, architecture, DiseaseStage) %>%
  summarise(mean_percent = mean(percent, na.rm = TRUE), .groups = "drop")

# Compute difference between CBM and NDMM per celltype / architecture
diff_dt <- summary_dt %>%
  pivot_wider(names_from = DiseaseStage, values_from = mean_percent) %>%
  mutate(diff_CBM_NDMM = NDMM - CBM)  # or abs(NDMM-CBM) if you want magnitude

# 4️⃣ Create PDF
pdf("PlasmaCell_Interactions_Per_Celltype_broad_test.pdf", width = 10, height = 7)

for (ct in celltypes) {
  data_ct <- plot_data %>% filter(celltype == ct)
  
  p <- ggplot(data_ct, aes(x = DiseaseStage, y = percent)) +
    geom_boxplot(aes(fill = DiseaseStage), outlier.shape = NA) +
    geom_jitter(aes(color = DiseaseStage), width = 0.2, size = 2) +
    facet_wrap(~ architecture, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Plasma Cell Interactions -", ct),
         x = "Disease Stage",
         y = "Percentage of Plasma Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  print(p)
}

dev.off()

# Make a list of ordering for each architecture
arch_ordering <- diff_dt %>%
  group_by(architecture) %>%
  arrange(desc(diff_CBM_NDMM)) %>%
  summarise(celltype_order = list(celltype)) %>%
  deframe()  # named list: architecture -> ordered celltypes


pdf("PlasmaCell_Interactions_Ordered_Per_Architecture_broadxx.pdf", width = 18, height = 7)

architectures <- unique(plot_data$architecture)

for (arch in architectures) {
  data_arch <- plot_data %>% filter(architecture == arch)
  
  # Order celltypes for this architecture
  celltype_levels <- arch_ordering[[arch]]
  data_arch <- data_arch %>%
    mutate(celltype = factor(celltype, levels = celltype_levels))
  
  p <- ggplot(data_arch, aes(x = celltype, y = percent)) +
    geom_boxplot(aes(fill = DiseaseStage), outlier.shape = NA) +
    geom_jitter(aes(color = DiseaseStage), width = 0.2, size = 0) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Plasma Cell Interactions -", arch),
         x = "Cell Type",
         y = "Percentage of Plasma Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}

dev.off()

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
pdf("PlasmaCell_Interactions_Ordered_by_Avg_CBMNDMM_Stages_broadxx.pdf", width = 18, height = 7)

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
}

dev.off()

pdf("PlasmaCell_Interactions_Celltypes_x_Architectures_grouped_perDiseaseStage.pdf",
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



# 1️⃣ Add DiseaseStage column if not already done
plot_data <- results_df_annot_broad %>%
  mutate(DiseaseStage = sapply(strsplit(patient, "_"), `[`, 1)) %>%
  filter(architecture != "unqualified")  # remove unqualified

# 2️⃣ Compute average percent per architecture x disease stage x celltype
avg_dt <- plot_data %>%
  group_by(architecture, DiseaseStage, celltype) %>%
  summarise(mean_percent = mean(percent, na.rm = TRUE), .groups = "drop")

# 3️⃣ Create ordering for each architecture x disease stage
ordering_list <- avg_dt %>%
  group_by(architecture, DiseaseStage) %>%
  arrange(desc(mean_percent)) %>%
  summarise(celltype_order = list(celltype), .groups = "drop") %>%
  as.data.table()

# 4️⃣ Plot
pdf("PlasmaCell_Interactions_Per_Celltype_Arch_Disease_broadxxx.pdf", width = 16, height = 10)

architectures <- unique(plot_data$architecture)
diseasestages <- unique(plot_data$DiseaseStage)

for (arch in architectures) {
  for (ds in diseasestages) {
    
    data_sub <- plot_data %>% filter(architecture == arch, DiseaseStage == ds)
    
    if(nrow(data_sub) == 0) next  # <-- skip if no cells
    
    cell_order <- ordering_list[architecture == arch & DiseaseStage == ds, celltype_order][[1]]
    data_sub <- data_sub %>% mutate(celltype = factor(celltype, levels = cell_order))
    
    p <- ggplot(data_sub, aes(x = celltype, y = percent)) +
      geom_boxplot(aes(fill = celltype), outlier.shape = NA) +
      geom_jitter(aes(color = celltype), width = 0.2, size = 2) +
      theme_minimal(base_size = 14) +
      labs(title = paste0("Plasma Cell Interactions - ", arch, " - ", ds),
           x = "Cell Type",
           y = "Percentage of Plasma Cells") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    print(p)
  }
}
dev.off()



# Plot: stacked bar per patient showing % of plasma cells interacting with each cell type
ggplot(results_df_annot_precise, aes(x = patient, y = percent, fill = celltype)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  labs(title = "Percentage of Plasma Cells Interacting with Each Cell Type per Patient",
       x = "Patient", y = "Percentage of Plasma Cells", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

results_df$patient <- factor(results_df$patient, levels = patient_order)

results_df <- results_df %>%
  mutate(Disease = sapply(strsplit(as.character(patient), "_"), function(x) x[3]))

# Step 2: Calculate the percentage per disease
results_df <- results_df %>%
  group_by(Disease) %>%
  mutate(Total = sum(percent),  # Sum of values per disease
         MD_Percent = percent / Total * 100)  # Percent calculation per disease

# Step 3: Calculate the average percentage per disease
results_df <- results_df %>%
  group_by(Disease, celltype) %>%
  summarise(Avg_Percent_Diseasegroup = mean(percent, na.rm = TRUE))

# Check
head(results_df)
# Plot: one panel per cell type
ggplot(results_df, aes(x = Disease, y = Avg_Percent_Diseasegroup, fill = Disease)) +
  geom_col() +
  facet_wrap(~ celltype_fac, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Plasma Cell Interactions per Cell Type",
       y = "Percentage of Plasma Cells", x = NULL) +
  theme(axis.text.x = element_blank())


facet_order <- results_df %>%
  group_by(celltype) %>%
  summarise(max_val = mean(Avg_Percent_Diseasegroup, na.rm = TRUE)) %>%
  arrange(desc(max_val)) %>%
  pull(celltype)

results_df$celltype_fac <- factor(results_df$celltype, levels = facet_order)


saveRDS(results_df, "results_df.rds")
results_df <- readRDS("results_df.rds")

disease_list <- split(results_df, results_df$Disease)


# Get unique diseases
disease_vec <- unique(results_df$Disease)

# List to store plots (optional)
plots <- list()

# Loop over diseases
for (d in disease_vec) {
  
  # Subset for this Disease
  df_disease <- results_df %>% filter(Disease == d)
  
  # Order celltypes by max value in this Disease
  celltype_order <- df_disease %>%
    group_by(celltype_fac) %>%
    summarise(max_val = max(Avg_Percent_Diseasegroup, na.rm = TRUE)) %>%
    arrange(desc(max_val)) %>%
    pull(celltype_fac)
  
  # Apply ordering
  df_disease$celltype_fac <- factor(df_disease$celltype_fac, levels = celltype_order)
  
  # Create plot
  p <- ggplot(df_disease, aes(x = celltype_fac, y = Avg_Percent_Diseasegroup, fill = celltype_fac)) +
    geom_col() +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0("Plasma Cell Interactions — ", d),
      x = "Cell Type",
      y = "Average % of Plasma Cells"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Store plot in list (optional)
  plots[[d]] <- p
  
  # Show plot
  print(p)
  
  # Optional: save directly
  ggsave(
    filename = paste0("PC_interactions_", d, ".tiff"),
    plot = p,
    width = 10,
    height = 6
  )
}





# plasma cell interacting with plasma cells -------------------------------


patients <- unique(meta_dt$patient)
architectures <- unique(meta_dt$architecture)

for (p in patients) {
  for (arch in architectures) {
    
    outfile <- file.path(
      "RDS_files/Interactions/interaction_absolute_results_architecture/interactions_PC_PC/Per_status_per_archi/",
      paste0("pc_pc_neighbors_", p, "_", arch, ".rds")
    )
    
    if (file.exists(outfile)) {
      message("Skipping ", p, " - ", arch)
      next
    }
    
    message("Processing ", p, " - ", arch, " ...")
    
    meta_sub <- meta_dt[patient == p & architecture == arch]
    if (nrow(meta_sub) == 0) next
    
    plasma_ids <- meta_sub[celltype == "B/Plasma Cell", cell_id]
    if (length(plasma_ids) == 0) next
    
    # Filter edges: only plasma → plasma
    local_edges_pc_pc <- edges_df[edges_df$from %in% plasma_ids & edges_df$to %in% plasma_ids, ]
    
    # Build neighbor list
    neighbors_pc_pc <- split(local_edges_pc_pc$to, local_edges_pc_pc$from)
    neighbors_pc_pc <- lapply(neighbors_pc_pc, unique)
    
    # Count neighbors for plasma cells
    n_neighbors <- rep(0, length(plasma_ids))
    names(n_neighbors) <- plasma_ids
    if (length(neighbors_pc_pc) > 0) {
      n_neighbors[names(neighbors_pc_pc)] <- vapply(neighbors_pc_pc, length, integer(1))
    }
    
    # Summary table
    pc_summary <- as.data.table(table(n_neighbors))
    setnames(pc_summary, c("n_neighbors", "n_cells"))
    pc_summary[, n_neighbors := as.integer(as.character(n_neighbors))]
    
    # Add patient/architecture info
    pc_summary[, patient := p]
    pc_summary[, architecture := arch]
    
    # Save
    saveRDS(pc_summary, outfile)
    message("Done ", p, " - ", arch, " (", nrow(pc_summary), " neighbor bins)")
  }
}


# Folder with all PC-PC RDS files
folder <- "RDS_files/Interactions/interaction_absolute_results_architecture/interactions_PC_PC/Per_status_per_archi/"
files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)

all_summary <- rbindlist(lapply(files, function(f) {
  dt <- readRDS(f)
  
  # Compute total plasma cells in this sample
  total_pcs <- sum(dt$n_cells)
  
  # Compute percentage for each n_neighbors
  dt[, percent := n_cells / total_pcs * 100]
  
  # Compute cumulative percentages
  pct_at_least1 <- sum(dt[n_neighbors >= 1, n_cells]) / total_pcs * 100
  pct_at_least2 <- sum(dt[n_neighbors >= 2, n_cells]) / total_pcs * 100
  
  # Add patient and architecture info
  dt[, `:=`(
    patient = unique(dt$patient),
    architecture = unique(dt$architecture),
    pct_at_least1 = pct_at_least1,
    pct_at_least2 = pct_at_least2
  )]
  
  dt
}))

totals <- all_summary[, .(total_cells = sum(n_cells)), 
                      by = .(patient, architecture)]

valid <- totals[total_cells >= 80]

# Filter pct_wide


# Optionally reshape wide: each n_neighbors becomes a column with the percentage
pct_wide <- dcast(
  all_summary, 
  patient + architecture + pct_at_least1 + pct_at_least2 ~ n_neighbors, 
  value.var = "percent", fill = 0
)

# Rename the columns for clarity
setnames(pct_wide, old = as.character(0:max(all_summary$n_neighbors)), 
         new = paste0("pct_", 0:max(all_summary$n_neighbors)))

# View result
head(pct_wide)
pct_wide <- pct_wide[pct_wide$architecture != "unqualified",]
pct_wide[, patient_arch := paste(patient, architecture, sep = "_")]
pct_wide[, disease_stage := sub("_.*", "", patient)]

pct_wide_filtered <- pct_wide[valid, on = .(patient, architecture)]


summary_dt <- pct_wide_filtered[, .(
  mean_pct = mean(pct_at_least1, na.rm = TRUE),
  sd_pct   = sd(pct_at_least1, na.rm = TRUE),
  n        = .N
), by = .(disease_stage, architecture)]


summary_dt <- pct_wide_filtered[, .(
  mean_pct = mean(pct_at_least2, na.rm = TRUE),
  sd_pct   = sd(pct_at_least2, na.rm = TRUE),
  n        = .N
), by = .(disease_stage, architecture)]




ggplot(summary_dt,
       aes(x = disease_stage, y = mean_pct, fill = architecture)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  theme_bw() +
  labs(x = "Disease Stage",
       y = "Mean pct_at_least1 ± SD")

ggplot(summary_dt,
       aes(x = disease_stage, y = mean_pct)) +
  geom_col(fill = "grey70") +
  geom_errorbar(
    aes(ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  facet_wrap(~ architecture, nrow = 1) +
  theme_bw() +
  labs(x = "Disease Stage",
       y = "Mean pct_at_least1 ± SD")

summary_dt <- pct_wide[, .(
  mean_pct = mean(pct_at_least1, na.rm = TRUE),
  sd_pct   = sd(pct_at_least1, na.rm = TRUE),
  n        = .N
), by = .(architecture)]

ggplot(summary_dt,
       aes(x = architecture, y = mean_pct, fill = architecture)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  theme_bw() +
  labs(x = "Disease Stage",
       y = "Mean pct_at_least1 ± SD")



folder <- "RDS_files/Interactions/interaction_absolute_results_architecture/interactions_PC_PC/Per_status_per_archi/"
files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)

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
  patient + pct_at_least1 + pct_at_least2 ~ n_neighbors, 
  value.var = "percent", fill = 0
)
setnames(pct_patient_wide, old = as.character(0:max(all_summary_patient$n_neighbors)), 
         new = paste0("pct_", 0:max(all_summary_patient$n_neighbors)))

head(pct_patient_wide)



# Copy table
dt <- copy(pct_patient_wide)

# Extract disease stage from patient name
dt[, disease_stage := sub("_.*", "", patient)]

# Compute mean and SD per disease stage
summary_dt <- dt[, .(
  mean_at_least1 = mean(pct_at_least1),
  sd_at_least1   = sd(pct_at_least1),
  mean_at_least2 = mean(pct_at_least2),
  sd_at_least2   = sd(pct_at_least2)
), by = disease_stage]

summary_dt$disease_stage_fac <- factor(summary_dt$disease_stage, 
                                       levels = c("CBM", "SMM", "NDMM", "PCL", "Relapse"))

# ---- Barplot 1: at least 1 PC ----
p1 <- ggplot(summary_dt, aes(x = disease_stage_fac, y = mean_at_least1, fill = disease_stage_fac)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_at_least1 - sd_at_least1, ymax = mean_at_least1 + sd_at_least1), width = 0.2) +
  labs(title = "Percentage of PCs interacting with at least 1 PC", x = "Disease stage", y = "Mean Percentage ± SD") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# ---- Barplot 2: at least 2 PCs ----
p2 <- ggplot(summary_dt, aes(x = disease_stage_fac, y = mean_at_least2, fill = disease_stage_fac)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_at_least2 - sd_at_least2, ymax = mean_at_least2 + sd_at_least2), width = 0.2) +
  labs(title = "Percentage of PCs interacting with at least 2 PCs", x = "Disease stage", y = "Mean Percentage ± SD") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

p1
p2

summary_long <- melt(
  summary_dt,
  id.vars = "disease_stage_fac",
  measure = list(
    value = c("mean_at_least1", "mean_at_least2"),
    sd    = c("sd_at_least1",   "sd_at_least2")
  ),
  variable.name = "interaction_level"
)
summary_long[, interaction_level := factor(interaction_level,
                                           labels = c("1 PC", "2 PCs")
)]
summary_long[, x_label := paste(disease_stage_fac, interaction_level)]

ggplot(summary_long, aes(x = x_label, y = value, fill = disease_stage_fac)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.2) +
  labs(
    x = "Disease stage and interaction level",
    y = "Mean Percentage ± SD",
    title = "Plasma cell–plasma cell interactions"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
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



p_combined <- ggplot(summary_long,
                     aes(x = disease_stage_fac, y = value, fill = disease_stage_fac)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                width = 0.2) +
  facet_wrap(~ interaction_level, nrow = 1) +
  scale_fill_manual(values = status_colors) +
  labs(
    x = "Disease stage",
    y = "Mean Percentage ± SD",
    title = "Plasma cell–plasma cell interactions"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")




### Needed for the first functions
folder <- "RDS_files/Interactions/Interactions_PC_PC_per_cluster/"
files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)

all_summary_patient <- rbindlist(lapply(files, function(f) {
  dt <- readRDS(f)
  patient_id <- gsub("\\.rds$", "", basename(f))  # remove folder path and ".rds"
  dt[, .(pc_cluster, mean_neighbors, patient = patient_id)]
}))
##

all_summary_patient <- readRDS("RDS_files/Interactions/Interactions_PC_PC_per_cluster/all_patients_pc_pc_neighborsMeanTACIposneg2.rds")


all_summary_patient$patient_short <- sapply(strsplit(all_summary_patient$patient, "_"), function(x) paste(tail(x, 2), collapse = "_"))
all_summary_patient$status <- sapply(strsplit(all_summary_patient$patient, "_"), function(x) x[length(x)-1])




pc_numbers <- table(meta_df$annot_B_lineage_true_topmarker, meta_df$ObjectName_anonymous)
# Convert table to data.frame for easy lookup
pc_numbers_df <- as.data.frame(pc_numbers)
colnames(pc_numbers_df) <- c("pc_cluster", "patient_short", "n_cells")


all_summary_patient_NDMM <- all_summary_patient[all_summary_patient$status %in% c("NDMM","PCL"),]
all_summary_patient_NDMM_filtered <- all_summary_patient_NDMM
# Merge with summary
all_summary_patient_NDMM_filtered <- merge(
  all_summary_patient_NDMM,
  pc_numbers_df,
  by.x = c("pc_cluster", "patient_short"),
  by.y = c("pc_cluster", "patient_short")
)

# Keep only cluster × patient with at least 50 cells
all_summary_patient_NDMM_filtered <- all_summary_patient_NDMM_filtered[all_summary_patient_NDMM_filtered$n_cells >= 30,]

# Order pc_cluster by median mean_neighbors
all_summary_patient_NDMM_filtered <- all_summary_patient_NDMM_filtered %>%
  group_by(pc_cluster) %>%
  mutate(median_neighbors = median(mean_neighbors)) %>%
  ungroup() %>%
  mutate(pc_cluster = factor(pc_cluster, levels = unique(pc_cluster[order(median_neighbors)])))


# Plot
ggplot(all_summary_patient_NDMM_filtered, aes(x = pc_cluster, y = mean_neighbors)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Plasma Cell Cluster",
    y = "Average Number of PC Neighbors",
    title = "Distribution of PC–PC Interactions per Cluster (ordered)"
  )




# Cellular neighborhood analysis -----------------------------------------------------------------

spe_Xenium <- buildSpatialGraph(spe_Xenium, img_id = "ObjectName_anonymous_grouped", type = "delaunay", max_dist = 20)


spe_Xenium <- aggregateNeighbors(spe_Xenium , colPairName = "delaunay_interaction_graph", 
                               aggregate_by = "metadata", count_by = "annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv")




set.seed(220705)



cn_10 <- kmeans(spe_Xenium $aggregatedNeighbors, centers = 10)
spe_Xenium$cn_celltypes_10 <- as.factor(cn_10$cluster)


# shows clusters in CN


for_plot <- prop.table(table(spe_Xenium$cn_celltypes_10, as.character(spe_Xenium$annot)), margin = 1)

niche_counts <- table(spe_Xenium$cn_celltypes_10)
neighborhood_percentage <- niche_counts / sum(niche_counts) * 100
# Step 2: Create a data frame for row annotation (with the percentage of the neighborhood size)
row_annotation <- data.frame(Neighborhood_Percentage = neighborhood_percentage)
# Convert row_annotation to a named vector
neighborhood_percentage_vector <- setNames(row_annotation$Neighborhood_Percentage.Freq, row_annotation$Neighborhood_Percentage.Var1)
row_annotation_df <- data.frame(Neighborhood_Percentage = neighborhood_percentage_vector)

pheatmap( for_plot,
          color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
          scale = "column",
          annotation_row = row_annotation_df,  # tilt column labels
          border_color = NA,
          angle_col = 45  # tilt column labels
)

# Convert for_plot to numeric matrix if not already
for_plot_matrix <- as.matrix(for_plot)
storage.mode(for_plot_matrix) <- "numeric"

# Scale columns (z-score per cell type)
for_plot_scaled <- scale(for_plot_matrix)  # scale by column (default)
for_plot_df_neighborhood <- as.data.frame(for_plot_scaled)
# Rename columns for clarity
colnames(for_plot_df_neighborhood) <- c("Neighborhood", "Celltype", "Enrichment_Score")

# Optionally filter to a single neighborhood
selected_neighborhood <- "Tumor-Associated"  # change this as needed, or NULL for all
for_plot_df_neighborhood <- for_plot_df_neighborhood
if (!is.null(selected_neighborhood)) {
  for_plot_df_neighborhood <- for_plot_df_neighborhood %>% filter(Neighborhood == selected_neighborhood)
}

# Order cell types by decreasing proportion (within the filtered neighborhood)
for_plot_df_neighborhood <- for_plot_df_neighborhood %>%
  arrange(Enrichment_Score) %>%
  mutate(Celltype = factor(Celltype, levels = unique(Celltype)))

# Plot
ggplot(for_plot_df_neighborhood, aes(x = Celltype, y = Enrichment_Score, fill = Enrichment_Score)) +
  geom_col() +
  scale_fill_gradient2(low = "darkblue", mid = "gray", high = "darkred", midpoint = 0) +
  coord_flip() +
  theme_minimal() +
  labs(x = "Cell type", y = "Enrichment score") +
  ggtitle(ifelse(is.null(selected_neighborhood), "All neighborhoods", selected_neighborhood)) +
  theme(axis.text.y = element_text(size = 10))



# Interaction analysis ------------------------------------------------------------



# 🔹 Parallel setup (adjust workers depending on your machine)
BPPARAM <- SnowParam(workers = 5, RNGseed = 221029)
library(scales)
system.time(
  out_all <- testInteractions(spe_Xenium, 
                              group_by = "Tissue_ext",
                              label = "annot", 
                              colPairName = "delaunay_interaction_graph",
                              BPPARAM = SerialParam(RNGseed = 221029))
)

system.time(mean(5))
saveRDS(out_all      ,"Y:/RDS_files/interactions_all_new_annot.rds")

library(BiocParallel)

# 🔹 Output directory
output_dir <- "Y:/RDS_files"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)



# 🔹 Unique patients and architectures
patients <- unique(spe_Xenium$ObjectName_anonymous_grouped)
architectures <- unique(spe_Xenium$tumor_percentage_grouped035ext_hex_average_method2)

# 🔹 Loop per patient x architecture
for (p in patients) {
  for (arch in architectures) {
    
    # Output file for this patient × architecture
    outfile <- file.path(output_dir, paste0("interactions_", p, "_", arch, ".rds"))
    if (file.exists(outfile)) {
      message("Skipping ", p, " - ", arch)
      next
    }
    message("Processing ", p, " - ", arch, " ...")
    # Subset SPE for patient × architecture
    spe_subset <- spe_Xenium[, spe_Xenium$ObjectName_anonymous_grouped == p & spe_Xenium$tumor_percentage_grouped035ext_hex_average_method2 == arch]
    if (ncol(spe_subset) == 0) {
      message("  No cells for ", p, " - ", arch, ", skipping.")
      next
    }
    # Run interaction test
    runtime <- system.time({
      interaction_result <- testInteractions(
        spe_subset,
        group_by = "Tissue_ext",
        label = "annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv",
        colPairName = "delaunay_interaction_graph",
        BPPARAM = BPPARAM
      )
    })
    # Save result immediately
    saveRDS(interaction_result, outfile)
    message("✅ Finished ", p, " - ", arch, " in ", runtime[3], " seconds.\n")
  }
}

cat("All patient × architecture results saved!\n")


# Step 1 — convert and rename
out_all_tib <- out_all_no_duplicates %>%
  as_tibble() %>%
  mutate(
    from_label = ifelse(from_label %in% names(rename_vector), 
                        rename_vector[from_label], from_label),
    to_label   = ifelse(to_label %in% names(rename_vector), 
                        rename_vector[to_label], to_label)
  )
# or
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


# 6. Plot the heatmap
ggplot(df, aes(from_label, to_label, fill = sum_sigval)) +
  geom_tile() +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




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

# Plot stacked barplot
ggplot(df_pc, aes(x = interaction, y = sum_sigval, fill = sum_sigval)) +
  geom_col() +
  coord_flip() +  
  scale_fill_gradient2(low = "darkblue", mid = "gray", high = "darkred", midpoint = 0) +
  labs(x = "Interaction", y = "Sum sigval", 
       title = "Plasma Cell Interactions (ranked by sum_sigval)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))




# Combine Spatial Interaction heatmaps of only plasma cells ------------------------------------



# Step 1: Summarize interactions per dataset
plasma_interactions <- lapply(combined, function(df) {
  df %>%
    group_by(group_by, from_label, to_label) %>%   # group_by = patient/sample
    summarize(sum_sigval = sum(sigval, na.rm = TRUE), .groups = "drop") %>%
    filter(from_label == "Plasma Cell")
})



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

# Plot
ggplot(combined_df, aes(x = to_label, y = from_label, fill = avg_sigval)) +
  geom_tile(color = "white") +  # optional: white borders between tiles
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y", switch = "y") +  # datasets on left
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = median(combined_df$avg_sigval, na.rm = TRUE)) +
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

combined_df_archi <- combined_df %>%
  filter(dataset %in% c("normal PC percentage", "sparse", "dispersed", "dense", "all"))

combined_df_archi_no_CBM_with_CBM <- combined_df %>%
  filter(dataset %in% c("CBM", "normal PC percentage", "sparse", "dispersed", "dense", "all"))


combined_df_state <- combined_df %>%
  filter(dataset %in% c("CBM", "SMM", "Relapse", "NDMM", "PCL", "all"))

ggplot(combined_df_archi_no_CBM_with_CBM, aes(x = to_label, y = from_label, fill = avg_sigval)) +
  geom_tile(color = "white") +  # optional: white borders between tiles
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y", switch = "y") +  # datasets on left
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = median(combined_df$avg_sigval, na.rm = TRUE)) +
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



combined_plot_df <- combined_df %>%
  group_by(dataset) %>%
  mutate(n_rows = n()) %>%  # number of rows per dataset
  mutate(row_in_dataset = row_number()) %>%
  ungroup() %>%
  # compute cumulative offset
  mutate(dataset_offset = lag(cumsum(n_rows), default = 0)) %>%
  mutate(y_plot = row_in_dataset + dataset_offset)

# Optional: dataset label positions
dataset_labels <- combined_plot_df %>%
  group_by(dataset) %>%
  summarize(y_mid = mean(y_plot))

# --- Plot ---
ggplot(combined_plot_df, aes(x = to_label, y = y_plot, fill = avg_sigval)) +
  geom_tile(color = "white") +  # white lines between columns
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"),
                       midpoint = median(combined_plot_df$avg_sigval, na.rm = TRUE)) +
  scale_y_reverse() +
  geom_text(data = dataset_labels, aes(x = -0.5, y = y_mid, label = dataset),
            inherit.aes = FALSE, hjust = 1, fontface = "bold", size = 5) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 50),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  labs(
    fill = "Interaction Score",
    title = "Interactions of Plasma Cells with Other Cell Types"
  )


# Combine 3 spatial methods for a niche score -----------------------------

interaction_df <- df_pc[!(df_pc$to_label %in% c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell")), ]
interaction_df <- interaction_df %>%
  arrange(desc(sum_sigval)) %>%
  mutate(order_rank = row_number())  # new column with ranking

for_plot_df_neighborhood_ranked <- for_plot_df_neighborhood[!(for_plot_df_neighborhood$Celltype %in% c("Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell")), ]
for_plot_df_neighborhood_ranked <- for_plot_df_neighborhood_ranked %>%
  arrange(desc(Enrichment_Score)) %>%
  mutate(order_rank = row_number())  # new column with ranking

significance_df_plot_densevsnormal <- significance_df_plot
significance_df_plot_densevsnormal <- significance_df_plot_densevsnormal %>%
  arrange(desc(logp_value)) %>%
  mutate(order_rank = row_number())  # new column with ranking


# Example: df1, df2, df3 with same rownames
# Each has a column "rank"

combined <- interaction_df %>%
  select(celltype = to_label, rank1 = order_rank) %>%
  inner_join(for_plot_df_neighborhood_ranked %>% 
               select(celltype = Celltype, rank2 = order_rank),
             by = "celltype") %>%
  inner_join(significance_df_plot_densevsnormal %>%
               select(celltype = cluster_celltype, rank3 = order_rank),
             by = "celltype") %>%
  mutate(rank_sum = rank1 + rank2 + rank3) %>%
  arrange(rank_sum)

# Plot stacked barplot
combined <- combined %>%
  arrange(desc(rank_sum)) %>%
  mutate(celltype = factor(celltype, levels = celltype))

ggplot(combined, aes(x = celltype, y = rank_sum, fill = rank_sum)) +
  geom_col() +
  coord_flip() +  
  scale_fill_gradient2(low ="darkred", mid = "gray", high =  "darkblue", midpoint = 40) +
  labs(x = "Interaction", y = "Sum sigval", 
       title = "Plasma Cell Interactions (ranked by 3 spatial methods)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Invert rank_sum
combined <- combined %>%
  mutate(rank_sum_inverted = max(rank_sum) - rank_sum + min(rank_sum))  # ensures positive values

# Plot with inverted rank_sum
ggplot(combined, aes(x = celltype, y = rank_sum_inverted, fill = rank_sum_inverted)) +
  geom_col() +
  coord_flip() +  
  scale_fill_gradient2(low ="darkblue", mid = "gray", high =  "darkred", midpoint = median(combined$rank_sum_inverted)) +
  labs(x = "Interaction", y = "Inverted Sum sigval", 
       title = "Plasma Cell Interactions (ranked by 3 spatial methods, inverted)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))




# Cells in neighborhood ---------------------------------------------------


# Combine relevant info into a dataframe
df <- data.frame(
  neighborhood = spe_Xenium$cn_celltypes_12_named,
  celltype = spe_Xenium$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv,
  status = spe_Xenium$Status
)


# Filter only plasma cells
df_pc <- df %>% filter(celltype == "Plasma Cell")

# Count plasma cells per status × CN
pc_counts <- df_pc %>%
  group_by(status, neighborhood) %>%
  summarise(n = n(), .groups = "drop")

# Compute per-status percentages (sum to 100% per status)
pc_percentages <- pc_counts %>%
  group_by(status) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

ggplot(pc_percentages, aes(x = neighborhood, y = percent, fill = status)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    x = "Neighborhood (CN)",
    y = "% of Plasma Cells",
    title = "Distribution of Plasma Cells across Neighborhoods",
    fill = "Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )


# clostest distance to tumor cell -----------------------------------------


################# closest distance to tumor cell

library(RANN)  # For fast nearest neighbor search

# Sample Data: x, y coordinates & cell type
coords <- spatialCoords(spe_Xenium)
coords <- as.data.frame(coords)
spe_Xenium$X <- coords$Pos_X
spe_Xenium$Y <- coords$Pos_Y

# Extract metadata with coordinates
meta_df_raw <- colData(spe_Xenium) %>%
  as.data.frame() %>%
  select(X, Y, annot_merged_final_with_macro_noslash_dis_lowquali_correct)  # Adjust column names as needed
# Separate tumor and non-tumor cells
tumor_cells <- meta_df_raw %>% 
  filter(annot_merged_final_with_macro_noslash_dis_lowquali_correct %in% c("Plasma cells",  "Cycling B/Plasma cells","Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell")) %>% 
  select(X, Y)
non_tumor_cells <- meta_df_raw %>% 
  filter(!(annot_merged_final_with_macro_noslash_dis_lowquali_correct %in% c("Plasma cells", "Cycling B/Plasma cells","Plasma Cell", "B or Plasma Cell", "Cycling B or Plasma Cell"))) %>% 
  select(X, Y)
# Compute minimum distance for each non-tumor cell
if (nrow(tumor_cells) > 0 & nrow(non_tumor_cells) > 0) {
  nn_result <- nn2(data = tumor_cells, query = non_tumor_cells, k = 1)
  meta_df_raw$min_distance <- 0
  meta_df_raw$min_distance[!(meta_df_raw$annot_merged_final_with_macro_noslash_dis_lowquali_correct %in% c("Plasma cells",  "Cycling B/Plasma cells"))] <- nn_result$nn.dists
}
# Compute minimum distance for each tumor cell to other tumor cells
if (nrow(tumor_cells) > 1) {  # Ensure there are at least two tumor cells
  nn_result_tumor <- nn2(data = tumor_cells, query = tumor_cells, k = 2)  # k=2 because the nearest neighbor to itself is excluded
  tumor_cells$min_distance_to_other_tumor <- nn_result_tumor$nn.dists[,2]  # Get the second closest neighbor (ignore self)
  
  # Add the calculated minimum distances back to the meta_df for tumor cells
  meta_df_raw$min_distance_to_other_tumor <- NA
  meta_df_raw$min_distance_to_other_tumor[meta_df_raw$annot_merged_final_with_macro_noslash_dis_lowquali_correct %in% c("Plasma cells", "Cycling B/Plasma cells")] <- tumor_cells$min_distance_to_other_tumor
} 
# Store the result back in the SPE object
colData(spe_Xenium)$min_distance <- meta_df_raw$min_distance

colData(spe_Xenium)$min_distance_to_other_tumor <- meta_df_raw$min_distance_to_other_tumor
colData(spe_Xenium)$min_distance_to_other_tumor_com <- ifelse(colData(spe_Xenium)$min_distance == 0,colData(spe_Xenium)$min_distance_to_other_tumor,colData(spe_Xenium)$min_distance )

colData(spe_HBM)$min_distance_to_other_tumor <- meta_df_raw$min_distance_to_other_tumor
colData(spe_HBM)$min_distance_to_other_tumor_com <- ifelse(colData(spe_HBM)$min_distance == 0,colData(spe_HBM)$min_distance_to_other_tumor,colData(spe_HBM)$min_distance )




# Extract metadata with coordinates
meta_df <- colData(spe_Xenium) %>%
  as.data.frame() %>%
  select(X, Y, annot_even_less_broad_correct, min_distance, annot_merged_final_with_macro_noslash_dis_lowquali_correct,min_distance_to_other_tumor_com, Tissue_ext, cn_celltypes_12_named,tumor_percentage_grouped035ext_hex_average_method2, Status_simp_2 )  ## Adjust column names as needed

meta_df <- meta_df[meta_df$tumor_percentage_grouped035ext_hex_average_method2 %in% c("sparse", "normal PC percentage"),]
meta_df <- meta_df[meta_df$Status_simp_2 != c("CBM"),]
meta_df <- meta_df[meta_df$Status_simp_2 %in% c("PCL"),]
meta_df<- meta_df[meta_df$Status_simp_2 %in% c("NDMM"),]
meta_df<- meta_df[meta_df$Status_simp_2 %in% c("CBM"),]
meta_df<- meta_df[meta_df$Status_simp_2 %in% c("SMM"),]
meta_df<- meta_df[meta_df$Status_simp_2 %in% c("Relapse"),]

# Reorder clusters by median min_distance
meta_df <- meta_df %>%
  mutate(annot_merged_final_with_macro_noslash_dis_lowquali_correct = reorder(annot_merged_final_with_macro_noslash_dis_lowquali_correct,min_distance_to_other_tumor_com, FUN = median, na.rm = TRUE))

# Step 1: Calculate median min_distance per cluster
cluster_order <- meta_df %>%
  group_by(annot_merged_final_with_macro_noslash_dis_lowquali_correct, annot_even_less_broad_correct) %>%  # Keep annot_broad info
  summarise(median_distance = median(min_distance, na.rm = TRUE), .groups = "drop") %>%
  arrange(annot_even_less_broad_correct, median_distance)  # First alphabetically, then by median distance

# Step 2: Apply ordering to the original data
meta_df <- meta_df %>%
  mutate(annot_merged_final_with_macro_noslash_dis_lowquali_correct = factor(annot_merged_final_with_macro_noslash_dis_lowquali_correct, levels = cluster_order$annot_merged_final_with_macro_noslash_dis_lowquali_correct ))


# Create an ordered boxplot
ggplot(meta_df, aes(x = annot_merged_final_with_macro_noslash_dis_lowquali_correct, y = min_distance, fill = annot_even_less_broad_correct)) +
  geom_boxplot(outlier.shape = NA) +
  #stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") +  # Red dot for median
  theme_minimal() +
  labs(title = "Minimum Distance to Tumor Niche",
       x = "Cluster Names (Ordered by Median Distance)",
       y = "Minimum Distance (µm)") +
  coord_cartesian(ylim = c(0, 200))  +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")  # Rotate labels for better readability

meta_df <- meta_df %>%
  mutate(annot_merged_final_with_macro_noslash_dis_lowquali_correct = factor(annot_merged_final_with_macro_noslash_dis_lowquali_correct, levels = sort(unique(as.character(annot_merged_final_with_macro_noslash_dis_lowquali_correct)))))


# Step 1: Compute the median min_distance for each annot_broad group (order annot_broad by median min_distance)
annot_order <- meta_df %>%
  group_by(annot_even_less_broad_correct) %>%
  summarise(median_distance = median(min_distance_to_other_tumor_com, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_distance)  # Order annot_broad by median min_distance

# Step 2: Compute the median min_distance for each ClusterNames and order them within annot_broad
cluster_order <- meta_df %>%
  group_by(annot_even_less_broad_correct , annot_merged_final_with_macro_noslash_dis_lowquali_correct) %>%
  summarise(median_distance_cluster = median(min_distance_to_other_tumor_com, na.rm = TRUE), .groups = "drop") %>%
  arrange(annot_even_less_broad_correct , median_distance_cluster) %>%  # Order within annot_broad by median min_distance of clusters
  distinct(annot_merged_final_with_macro_noslash_dis_lowquali_correct, annot_even_less_broad_correct)  # Get unique ClusterNames for correct order

# Step 3: Apply the ordering to the original data
meta_df <- meta_df %>%
  mutate(
    annot_even_less_broad_correct = factor(annot_even_less_broad_correct, levels = annot_order$annot_even_less_broad_correct),
    annot_merged_final_with_macro_noslash_dis_lowquali_correct = factor(
      annot_merged_final_with_macro_noslash_dis_lowquali_correct,
      levels = unique(cluster_order$annot_merged_final_with_macro_noslash_dis_lowquali_correct)  # Remove duplicates here
    )
  )


ggplot(meta_df, aes(x = annot_merged_final_with_macro_noslash_dis_lowquali_correct, y = min_distance_to_other_tumor_com, fill = annot_even_less_broad_correct )) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() +
  labs(title = "Minimum Distance to Tumor Cell",
       x = "Cluster Names",
       y = "Minimum Distance (µm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate labels for better readability
  coord_cartesian(ylim = c(0, 90))  + 
  #scale_y_log10() + # Limit y-axis
  scale_fill_brewer(palette = "Set3")  # Limit y-axis





# Divide Tissues in regions for pseudobulk-----------------------------------------------


# Assign "dense_dis" for "dense" or "dispersed"
integrated_seurat_all$dense_vs_sparse <- NA  # initialize

integrated_seurat_all$dense_vs_sparse[
  integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2_CBM %in% c("dense", "dispersed")
] <- "dense_dis"

# Assign "sparse_normal" for "normal PC percentage" or "sparse"
integrated_seurat_all$dense_vs_sparse[
  integrated_seurat_all$tumor_percentage_grouped035ext_hex_average_method2_CBM %in% c("normal PC percentage", "sparse")
] <- "sparse_normal"

meta_df <- integrated_seurat_all@meta.data

coords <- meta_df %>% filter(ObjectName_anonymous == "NDMM_5B") %>% select(x_corrected_together, y_corrected_together)
set.seed(123)
km <- kmeans(coords, centers = 4)
meta_df$region_for_pseudobulk[meta_df$ObjectName_anonymous == "NDMM_5B"] <- km$cluster
meta_df_NDMM_2 <- meta_df[meta_df$ObjectName_anonymous == "NDMM_5B",]

meta_df_NDMM_2 <- meta_df_NDMM_2 %>%
  mutate(region_pseudo_dense_or_sparse = ifelse(
    is.na(dense_vs_sparse),
    NA,
    paste0(region_for_pseudobulk, "_", dense_vs_sparse)
  ))


ggplot(meta_df_NDMM_2, aes(x = x_corrected_together, y = y_corrected_together, color = as.factor(region_pseudo_dense_or_sparse))) +
  geom_point(size = 1) +
  scale_color_brewer(palette = "Set1", name = "Region") +
  coord_fixed() +  # keeps aspect ratio square
  theme_minimal(base_size = 14) +
  labs(title = "NDMM_2: Spatial pseudobulk regions") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

### Now a loop for every patient


# Make sure the column for region labels exists
meta_df$region_for_pseudobulk <- NA
meta_df$region_pseudo_dense_or_sparse <- NA

# Loop over all patients
for (pid in unique(meta_df$ObjectName_anonymous)) {
  
  # Subset coordinates for this patient
  coords <- meta_df %>%
    filter(ObjectName_anonymous == pid) %>%
    select(x_corrected_together, y_corrected_together)
  
  # Run k-means to get 4 regions per patient
  set.seed(123)
  km <- kmeans(coords, centers = 4)
  
  # Assign regions back to meta_df
  meta_df$region_for_pseudobulk[meta_df$ObjectName_anonymous == pid] <- km$cluster
  
  # Create combined region + dense/sparse + patient label
  meta_df <- meta_df %>%
    mutate(region_pseudo_dense_or_sparse = ifelse(
      ObjectName_anonymous == pid,
      ifelse(is.na(dense_vs_sparse), 
             NA, 
             paste0(pid, "_region", km$cluster, "_", dense_vs_sparse)),
      region_pseudo_dense_or_sparse
    ))
}

meta_df$region_for_pseudobulk_patient <- paste0(meta_df$ObjectName_anonymous, "_", meta_df$region_for_pseudobulk)

integrated_seurat_all@meta.data <- meta_df

# Check result
table(meta_df$region_pseudo_dense_or_sparse)


min_cells <- 150
meta_df$region_for_pseudobulk_merged <- meta_df$region_pseudo_dense_or_sparse

for(pid in unique(meta_df$ObjectName_anonymous)) {
  df <- meta_df %>% filter(ObjectName_anonymous == pid)
  
  for(status in unique(df$dense_vs_sparse)) {
    df_status <- df %>% filter(dense_vs_sparse == status)
    
    # skip if no cells for this status
    if(nrow(df_status) == 0) next
    
    # compute sizes
    region_sizes <- df_status %>%
      group_by(region_for_pseudobulk_merged) %>%
      summarise(n_cells = n(),
                centroid_x = mean(x_corrected_together),
                centroid_y = mean(y_corrected_together),
                .groups = "drop")
    
    small_regions <- region_sizes %>% filter(n_cells < min_cells)
    
    if(nrow(small_regions) == 0) next  # all good
    
    # merge small regions
    for(sr in small_regions$region_for_pseudobulk_merged) {
      # centroid of small region
      sx <- region_sizes$centroid_x[region_sizes$region_for_pseudobulk_merged == sr]
      sy <- region_sizes$centroid_y[region_sizes$region_for_pseudobulk_merged == sr]
      
      # candidate regions to merge with (larger regions only)
      large_regions <- region_sizes %>% filter(n_cells >= min_cells & region_for_pseudobulk_merged != sr)
      
      # if no large region, merge with nearest other small region
      if(nrow(large_regions) == 0) {
        candidates <- region_sizes %>% filter(region_for_pseudobulk_merged != sr)
      } else {
        candidates <- large_regions
      }
      
      candidates <- candidates %>%
        mutate(dist = sqrt((centroid_x - sx)^2 + (centroid_y - sy)^2))
      
      nearest_region <- candidates$region_for_pseudobulk_merged[which.min(candidates$dist)]
      
      # merge
      df$region_for_pseudobulk_merged[df$region_for_pseudobulk_merged == sr] <- nearest_region
    }
  }
  
  # save back
  meta_df$region_for_pseudobulk_merged[meta_df$ObjectName_anonymous == pid] <- df$region_for_pseudobulk_merged
}

# Recreate the combined label
meta_df <- meta_df %>%
  mutate(region_pseudo_dense_or_sparse_merged = ifelse(
    is.na(dense_vs_sparse),
    NA,
    paste0(region_for_pseudobulk_merged)
  ))

# Check results
table(meta_df$region_pseudo_dense_or_sparse_merged)

# Add the merged region ID
integrated_seurat_all$region_for_pseudobulk_merged <- meta_df$region_for_pseudobulk_merged

# Add the combined region + dense/sparse label
integrated_seurat_all$region_pseudo_dense_or_sparse_merged <- meta_df$region_pseudo_dense_or_sparse_merged

table(integrated_seurat_all$region_pseudo_dense_or_sparse_merged)

meta_df_NDMM_2 <- meta_df[meta_df$ObjectName_anonymous == "SMM_1",]

ggplot(meta_df_NDMM_2, aes(x = x_corrected_together, y = y_corrected_together, color = as.factor(region_pseudo_dense_or_sparse_merged))) +
  geom_point(size = 1) +
  scale_color_brewer(palette = "Set1", name = "Region") +
  coord_fixed() +  # keeps aspect ratio square
  theme_minimal(base_size = 14) +
  labs(title = "NDMM_2: Spatial pseudobulk regions") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

saveRDS(meta_df, "meta_df.rds")
meta_df <- readRDS("meta_df.rds")



hetero_patients <- c("PCL_2", "NDMM_3", "NDMM_4", "NDMM_5", "NDMM_1", "NDMM_2", "PCL_1")
integrated_seurat_all_hetero <- subset(
  integrated_seurat_all,
  subset = ObjectName_anonymous_grouped %in% hetero_patients
)

# Suppose your Seurat object is `integrated_seurat_all`
# Aggregate expression per pseudobulk region
region_expr <- AverageExpression(
  integrated_seurat_all_hetero, 
  group.by = "region_pseudo_dense_or_sparse_merged",
  assays = "RNA", slot = "counts"
)$RNA

region_expr <- as.matrix(region_expr)

# PCA
region_pca <- prcomp(t(region_expr), scale. = TRUE)

# Plot first two PCs
pc_df <- data.frame(region = colnames(region_expr),
                    PC1 = region_pca$x[,1],
                    PC2 = region_pca$x[,2])

pc_df_clean <- pc_df %>% drop_na(PC1, PC2)

ggplot(pc_df_clean, aes(x = PC1, y = PC2, label = region, color = region)) +
  geom_point(size = 3) +
  #geom_text_repel() +
  theme_minimal()

ggplot(pc_df_clean, aes(PC1, PC2, label = region, color = region)) +
  geom_point(size = 3) +
  geom_text_repel(
    max.overlaps = Inf,
    force = 2,
    box.padding = 0.5,
    point.padding = 0.3,
    clip = "off"
  ) +
  theme_minimal() +
  theme(legend.position = "none")




# run for multiple Ks -----------------------------------------------------



assign_pseudobulk_regions <- function(meta_df, K = 4, seed = 123, min_cells = 150) {
  meta_df$region_for_pseudobulk <- NA
  meta_df$region_pseudo_dense_or_sparse <- NA
  meta_df$region_for_pseudobulk_merged <- NA
  for (pid in unique(meta_df$ObjectName_anonymous)) {
    # Subset to current patient AND cells with valid coordinates
    
    df_pid <- meta_df %>%
      filter(ObjectName_anonymous == pid & 
               !is.na(x_corrected_together) & 
               !is.na(y_corrected_together)&
               !is.na(dense_vs_sparse)) %>%
      mutate(cellname = rownames(.))  # store original rownames
    
    
    coords <- df_pid %>%
      select(x_corrected_together, y_corrected_together)
    
    set.seed(seed)
    km <- kmeans(coords, centers = K)
    
    
    df_pid$region_for_pseudobulk[df_pid$ObjectName_anonymous == pid] <- km$cluster
    
    # Only assign regions to cells with dense_vs_sparse not NA
    valid_cells <- !is.na(df_pid$dense_vs_sparse) & df_pid$ObjectName_anonymous == pid
    df_pid$region_pseudo_dense_or_sparse[valid_cells] <- 
      paste0(pid, "_region", km$cluster[valid_cells], "_", df_pid$dense_vs_sparse[valid_cells])
    
    # --- Merge small regions ---
    df <- df_pid[df_pid$ObjectName_anonymous == pid, ]
    df <- df[!is.na(df$dense_vs_sparse), ]
    df$region_for_pseudobulk_merged <- as.character(df$region_pseudo_dense_or_sparse)
    
    repeat_merge <- TRUE
    while(repeat_merge) {
      region_sizes <- df %>%
        group_by(region_for_pseudobulk_merged) %>%
        summarise(
          n_cells = n(),
          centroid_x = mean(x_corrected_together),
          centroid_y = mean(y_corrected_together),
          .groups = "drop"
        )
      
      small_regions <- region_sizes %>% filter(n_cells < min_cells)
      if(nrow(small_regions) == 0) { repeat_merge <- FALSE; break }
      
      for(sr in small_regions$region_for_pseudobulk_merged) {
        sx <- mean(df$x_corrected_together[df$region_for_pseudobulk_merged == sr])
        sy <- mean(df$y_corrected_together[df$region_for_pseudobulk_merged == sr])
        
        candidates <- region_sizes %>% filter(region_for_pseudobulk_merged != sr)
        candidates <- candidates %>% mutate(dist = sqrt((centroid_x - sx)^2 + (centroid_y - sy)^2))
        
        nearest_region <- as.character(candidates$region_for_pseudobulk_merged[which.min(candidates$dist)])
        
        df$region_for_pseudobulk_merged[df$region_for_pseudobulk_merged == sr] <- nearest_region
      }
    }
    
    # Then map back to meta_df using these cellnames
    # Make sure df has rownames as cellnames
    df$cellname <- rownames(df)
    # Map merged regions back to meta_df by matching cellnames
    # ensure column exists
    meta_df$region_for_pseudobulk_merged[match(df$cellname, rownames(meta_df))] <- df$region_for_pseudobulk_merged
  }
  
  meta_df$region_pseudo_dense_or_sparse_merged <- meta_df$region_for_pseudobulk_merged
  return(meta_df)
}

K_values <- c(6)
K_values <- c(3, 4, 5, 6)
seeds_for_k4 <- c(111, 222, 333)

results_list <- list()

meta_df_hetero <- meta_df[meta_df$ObjectName_anonymous_grouped %in% hetero_patients,]

# Ks other than 4
for(K in K_values[K_values != 4]) {
  cat("Running K =", K, "\n")
  df_out <- assign_pseudobulk_regions(meta_df_hetero, K = K, seed = 123)
  results_list[[paste0("K", K, "_seed123")]] <- df_out
}

# K = 4 with 3 seeds
for(s in seeds_for_k4) {
  cat("Running K = 4, seed =", s, "\n")
  df_out <- assign_pseudobulk_regions(meta_df_hetero, K = 4, seed = s)
  results_list[[paste0("K4_seed", s)]] <- df_out
}


meta_df_NDMM_2 <- results_list2$K6_seed123[results_list2$K6_seed123$ObjectName_anonymous == "NDMM_2",]
meta_df_NDMM_2 <- meta_df[meta_df$ObjectName_anonymous == "NDMM_4",]



ggplot(meta_df_NDMM_2, aes(x = x_corrected_together, y = y_corrected_together, color = as.factor(region_for_pseudobulk_merged))) +
  geom_point(size = 1) +
  #scale_color_brewer(palette = "Set1", name = "Region") +
  coord_fixed() +  # keeps aspect ratio square
  theme_minimal(base_size = 14) +
  labs(title = "NDMM_2: Spatial pseudobulk regions") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())#,
#legend.position = "none")

for(name in names(results_list)) {
  df <- results_list[[name]]
  
  # Make sure rownames exist
  if(is.null(rownames(df))) stop(paste0("No rownames in results for ", name))
  
  # Name of the new column in Seurat metadata
  new_col <- paste0("pseudo_regions_", name)
  
  # Initialize column in Seurat metadata if it doesn't exist
  integrated_seurat_all@meta.data[[new_col]] <- NA
  
  # Match by cellnames (rownames)
  matched_idx <- match(rownames(df), rownames(integrated_seurat_all@meta.data))
  
  # Keep only valid matches
  valid <- !is.na(matched_idx)
  
  # Assign the merged region info only to valid matches
  integrated_seurat_all@meta.data[[new_col]][matched_idx[valid]] <- df$region_for_pseudobulk_merged[valid]
}




# Pseudobulk analysis -----------------------------------------------------

# Set cluster identities
Idents(integrated_seurat_all) <- integrated_seurat_all@meta.data$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv

# Prepare results storage
results_all_genes <- list()
results_marker_genes <- list()

# Loop over clusters
for (clust in levels(Idents(integrated_seurat_all))) {
  
  cat("Processing cluster:", clust, "\n")
  
  # Subset cells for cluster
  cells <- WhichCells(integrated_seurat_all, idents = clust)
  seurat_cluster <- subset(integrated_seurat_all, cells = cells)
  
  # Get counts, sample IDs, and group labels
  counts_mat <- GetAssayData(seurat_cluster, assay = "RNA", slot = "counts")
  sample_ids <- seurat_cluster$ObjectName_anonymous_grouped
  sample_group_levels <- seurat_cluster$region_pseudo_dense_or_sparse_merged
  group <- seurat_cluster$dense_vs_sparse
  
  # Remove cells with NA group
  keep <- !is.na(group)
  counts_mat <- counts_mat[, keep, drop = FALSE]
  sample_ids <- sample_ids[keep]
  group <- group[keep]
  
  # --- Pseudobulk per sample × group ---
  #sample_group_levels <- unique(paste(sample_ids, group, sep = "_"))
  
  
  # --- Pseudobulk per region ---
  region_levels <- unique(sample_group_levels)
  pseudobulk <- t(sapply(region_levels, function(reg) {
    idx <- which(sample_group_levels == reg)
    rowSums(counts_mat[, idx, drop = FALSE])
  }))
  rownames(pseudobulk) <- region_levels
  
  
  # Remove pseudobulk samples with zero counts
  keep_samples <- rowSums(pseudobulk) > 0
  pseudobulk <- pseudobulk[keep_samples, , drop = FALSE]
  
  # DGEList
  dge <- DGEList(counts = t(pseudobulk))
  dge <- calcNormFactors(dge)
  
  # 1. Make sure rownames of pseudobulk are samples
  rownames(pseudobulk)   # should be: "NDMM_3_sparse_normal", "NDMM_3_dense_dis", etc.
  
  # 2. Construct metadata
  meta <- data.frame(sample_group = rownames(pseudobulk), stringsAsFactors = FALSE)
  meta <- meta[keep_samples, , drop = FALSE]
  meta$sample <- sapply(strsplit(meta$sample_group, "_"), function(x) paste(x[1:2], collapse = "_"))
  meta$group  <- sapply(strsplit(meta$sample_group, "_"), function(x) paste(x[4:length(x)], collapse="_"))
  
  
  # Make factors
  meta$sample <- factor(meta$sample)
  meta$group  <- factor(meta$group, levels = c("dense_dis", "sparse_normal"))
  
  # Ensure rownames of meta match columns of dge
  rownames(meta) <- meta$sample_group
  all(rownames(meta) == colnames(dge))  # must be TRUE
  
  # 3. Construct design matrix for paired analysis
  design <- model.matrix(~0 + group, data = meta)
  colnames(design) <- gsub("group", "", colnames(design))  # optional, cleaner names
  
  all(colnames(dge) == rownames(meta))
  
  # 1. Voom transformation
  v <- voom(dge, design)  # logCPM + mean-variance weights
  
  # 2. Estimate correlation for paired samples
  # 'block' indicates which samples are paired (e.g., same patient)
  corfit <- duplicateCorrelation(v, design, block = meta$sample)
  
  # 3. Fit linear model accounting for pairing
  fit <- lmFit(v, design, block = meta$sample, correlation = corfit$consensus)
  #fit <- lmFit(v, design)  # no block, no correlation
  
  # 4. Define contrast: Dense vs Sparse
  # Make sure the column names match exactly your design matrix
  colnames(design)
  # Example: "dense_dis" and "sparse_normal"
  contrast <- makeContrasts(DenseVsSparse = dense_dis - sparse_normal, levels = design)
  
  # 5. Fit the contrast
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # 6. Extract DE results
  res_all <- topTable(fit2, coef = "DenseVsSparse", number = Inf, sort.by = "p")
  
}

saveRDS(results_all_genes, "results_all_genes_pseudobulkregionsall_broad.rds")
#saveRDS(results_marker_genes, "results_marker_genes_pseudobulkregionsall_broad.rds")

cat("Pseudobulk DE analysis completed!\n")


# Pseudobulk wo dispersed -------------------------------------------------

heterogenous_samples <- c("PCL_1", "PCL_2", "NDMM_5", "NDMM_4", "NDMM_3")

Xenium_Object_Seurat_sub <- subset(
  Xenium_Object_Seurat,
  subset = ObjectName_anonymous_grouped %in% heterogenous_samples
)

pseudo_cols <- c("pseudo_regions_K4_seed111_wo_dispersed")  # replace with your actual columns

# Remove cells with NA in dense_vs_sparse
Xenium_Object_Seurat_sub <- subset(
  Xenium_Object_Seurat,
  cells = colnames(Xenium_Object_Seurat)[
    !is.na(Xenium_Object_Seurat@meta.data$pseudo_regions_K4_seed111_wo_dispersed)
  ]
)

# Count cells per Sample x Group
cell_counts <- table(Xenium_Object_Seurat_sub@meta.data[[meta_col]])

# Keep only groups with >= 500 cells
keep_groups <- names(cell_counts[cell_counts >= 200])

Xenium_Object_Seurat_sub <- subset(
  Xenium_Object_Seurat_sub,
  cells = rownames(Xenium_Object_Seurat_sub@meta.data)[
    Xenium_Object_Seurat_sub@meta.data$pseudo_regions_K4_seed111_wo_dispersed %in% keep_groups
  ]
)

# Set cluster identities (same as before)
Idents(Xenium_Object_Seurat_sub) <- Xenium_Object_Seurat_sub@meta.data$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv
#Idents(Xenium_Object_Seurat_sub) <- Xenium_Object_Seurat_sub@meta.data$annot_broad_run3_correct_ifndiv

# Loop over clusters
for (clust in levels(Idents(Xenium_Object_Seurat_sub))) {
  
  cat("  Processing cluster:", clust, "\n")
  
  # Subset cells for cluster
  cells <- WhichCells(Xenium_Object_Seurat_sub, idents = clust)
  seurat_cluster <- subset(Xenium_Object_Seurat_sub, cells = cells)
  
  # Get counts
  counts_mat <- GetAssayData(seurat_cluster, assay = "RNA", layer = "counts")
  sample_ids <- seurat_cluster$ObjectName_anonymous_grouped
  
  # Use the current metadata column
  sample_group_levels <- seurat_cluster@meta.data[[meta_col]]
  
  # Remove cells with NA group
  keep <- !is.na(sample_group_levels)
  counts_mat <- counts_mat[, keep, drop = FALSE]
  sample_ids <- sample_ids[keep]
  sample_group_levels <- sample_group_levels[keep]
  
  # --- Pseudobulk per region ---
  region_levels <- unique(sample_group_levels)
  pseudobulk <- t(sapply(region_levels, function(reg) {
    idx <- which(sample_group_levels == reg)
    rowSums(counts_mat[, idx, drop = FALSE])
  }))
  rownames(pseudobulk) <- region_levels
  
  # Remove pseudobulk samples with zero counts
  keep_samples <- rowSums(pseudobulk) > 0
  pseudobulk <- pseudobulk[keep_samples, , drop = FALSE]
  
  # DGEList
  dge <- DGEList(counts = t(pseudobulk))
  dge <- calcNormFactors(dge)
  
  # Construct metadata for voom
  meta <- data.frame(sample_group = rownames(pseudobulk), stringsAsFactors = FALSE)
  meta <- meta[keep_samples, , drop = FALSE]
  meta$sample <- sapply(strsplit(meta$sample_group, "_"), function(x) paste(x[1:2], collapse = "_"))
  meta$group  <- sapply(strsplit(meta$sample_group, "_"), function(x) paste(x[4:length(x)], collapse="_"))
  meta$sample <- factor(meta$sample)
  meta$group  <- factor(meta$group, levels = c("dense", "sparse_normal"))
  rownames(meta) <- meta$sample_group
  
  # Design matrix
  design <- model.matrix(~0 + group, data = meta)
  colnames(design) <- gsub("group", "", colnames(design))
  
  # Voom + duplicateCorrelation
  v <- voom(dge, design)
  corfit <- duplicateCorrelation(v, design, block = meta$sample)
  fit <- lmFit(v, design, block = meta$sample, correlation = corfit$consensus)
  
  # Contrast: Dense vs Sparse
  contrast <- makeContrasts(DenseVsSparse = dense - sparse_normal, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Extract DE results
  res_all <- topTable(fit2, coef = "DenseVsSparse", number = Inf, sort.by = "p")
  results_all_genes[[clust]] <- res_all
  
  # Marker genes
  #marker_genes <- marker_table %>% filter(cluster == clust) %>% pull(gene) %>% intersect(rownames(res_all))
  #res_markers <- res_all[marker_genes, , drop = FALSE]
  #results_marker_genes[[clust]] <- res_markers
}

# Save results per metadata column
saveRDS(results_all_genes, paste0("results_all_genes_specific_wo_dispersed", meta_col, ".rds"))
#saveRDS(results_marker_genes, paste0("results_marker_genes_broad", meta_col, ".rds"))




# Volcano Plots -----------------------------------------------------------



results_DEG_pseudo_hetero_K4_wo_dispersed <- readRDS("~//results_all_genes_broad_wo_dispersedpseudo_regions_K4_seed111.rds")

markers_df_pseudo <- map_df(
  names(results_DEG_pseudo_hetero_K4_wo_dispersed), 
  ~ results_DEG_pseudo_hetero_K4_wo_dispersed[[.x]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    mutate(celltype = .x),
  .id = NULL
)

markers_df_pseudo_per_status <- do.call(rbind, lapply(names(results_DEG_pseudo_between_status_notbroad), function(celltype) {
  
  inner_list <- results_DEG_pseudo_between_status_notbroad[[celltype]]
  
  do.call(rbind, lapply(names(inner_list), function(comparison) {
    df <- inner_list[[comparison]]
    df <- as.data.frame(df)
    df$gene <- rownames(df)
    df$celltype <- celltype
    df$comparison <- comparison
    df
  }))
}))

# Define thresholds
sig_threshold <- 0.05
logfc_threshold <- 0.25

# Create PDF for all volcano plots
pdf("DEG_Volcano_Plots/Volcano_plots_dense_vs_sparse_pseudobulkregion_hetero_filtered_K4_seed111_noDispersed.pdf", width = 10, height = 7)

for (ct in unique(markers_df_pseudo$celltype)) {
  
  cat("Plotting:", ct, "\n")
  
  df_plot <- markers_df_pseudo %>%
    filter(celltype == ct)
  
  pc_genes <- marker_table_broad %>%
    filter(cluster == "B/Plasma Cell", avg_log2FC > 1.2) %>%
    pull(gene)
  
  # Filter genes based on marker_table for this cluster
  cluster_genes <- marker_table_broad %>%
    filter(cluster == ct) %>%
    filter(avg_log2FC> 0.4) %>%
    pull(gene) %>%
    setdiff(pc_genes)
  
  df_plot <- df_plot %>%
    filter(gene %in% cluster_genes)
  
  # Skip if no genes left
  if (nrow(df_plot) == 0) next
  # Add significance column
  df_plot <- df_plot %>%
    mutate(significant = (adj.P.Val    < sig_threshold) & (abs(logFC ) > logfc_threshold))
  # Replace p_val_adj == 0 with a tiny value
  df_plot <- df_plot %>%
    mutate(p_val_adj_plot = ifelse(adj.P.Val       == 0, 1e-300, adj.P.Val))
  
  # Select top up- and down-regulated genes to label
  top_up <- df_plot %>%
    filter(logFC          > 0) %>%
    slice_max(order_by = -log10(adj.P.Val ), n = 30)
  
  top_down <- df_plot %>%
    filter(logFC           < 0) %>%
    slice_max(order_by = -log10(adj.P.Val ), n = 30)
  
  df_labels <- bind_rows(top_up, top_down)
  
  df_labels <- bind_rows(top_up, top_down) %>%
    distinct(gene, .keep_all = TRUE) %>%
    filter(
      significant == TRUE
    )
  
  
  # Volcano plot
  p <- ggplot(df_plot, aes(x = logFC           , y = -log10(p_val_adj_plot))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "black") +
    geom_text_repel(
      data = df_labels,
      aes(x = logFC        , y = -log10(p_val_adj_plot ), label = gene),  # explicitly set x, y, label
      inherit.aes = FALSE,
      size = 3,
      max.overlaps = 20
    )+
    labs(
      title = paste("dense_vs_sparse:", ct),
      x = "log2 Fold Change",
      y = "-log10(adj. p-value)",
      color = "Significant"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  print(p)
}

dev.off()



# IMC data analysis -------------------------------------------------------


library(future)
plan("multisession", workers = 15) #not too much when working with very big objects
options(future.globals.maxSize = 100 * 1024 * 1024 * 1024)  # Set a different value (e.g., 1000 MiB)

setOldClass("ExpData")

options(bitmapType = "cairo")
source('/scripts/readWriteRDS_multiThreaded.R')

library(BiocStyle)
library(imcRtools)
library(BiocManager)
library(remotes)
library(imcRtools)
library(cytomapper) 
library(dplyr)

library(SpatialExperiment)
library(scuttle)
library(scran)
library(scater)
library(uwot)
library(rtracklayer)
library(patchwork)
library(pheatmap)
library(gridExtra)
library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)
library(bluster)
library(BiocParallel)
library(ggplot2)
library(readxl)
library(dplyr)
## read in custom files
library(readr)

cur_features <- read_csv("~/Hyperion_data/Project/analysis/cpout/cell.csv")

dim(cur_features)
head(colnames(cur_features))

counts <- cur_features[,grepl("Intensity_MeanIntensity_", 
                              colnames(cur_features))]
counts2 <- cur_features[, c("AreaShape_Area", "AreaShape_Eccentricity", "AreaShape_MeanRadius")]

library(dplyr)
counts <- select(counts, -starts_with("Intensity_MeanIntensity_ProbabStack"))

counts <- cbind(counts, counts2)

counts[is.nan(counts)] <- 0

meta <- cur_features[,c("ImageNumber", "ObjectNumber", "AreaShape_Area",
                        "AreaShape_Eccentricity", "AreaShape_MeanRadius")]

coords <- cur_features[,c("Location_Center_X", "Location_Center_Y")]


spe <- SpatialExperiment(assays = list(counts = t(counts)),
                         colData = meta, 
                         sample_id = as.character(meta$ImageNumber),
                         spatialCoords = as.matrix(coords))
cur_pairs <- read_csv("~/Hyperion_data/Project/analysis/cpout/Object relationships.csv")

cur_from <- paste(cur_pairs$`First Image Number`, cur_pairs$`First Object Number`)
cur_to <- paste(cur_pairs$`Second Image Number`, cur_pairs$`Second Object Number`)

edgelist <- SelfHits(from = match(cur_from, 
                                  paste(spe$ImageNumber, spe$ObjectNumber)),
                     to = match(cur_to, 
                                paste(spe$ImageNumber, spe$ObjectNumber)),
                     nnode = ncol(spe))

colPair(spe, "neighborhood") <- edgelist

rowData(spe) <- read_csv("~/Hyperion_data/Project/analysis/cpout/panel_rightorder.csv")

colnames(spe@int_colData@listData$spatialCoords) <- c("Pos_X", "Pos_Y")



counts(spe)[1:37, 1:10]
head(colData(spe))
head(spatialCoords(spe))
colPair(spe, "neighborhood")
counts(spe)
head(rowData(spe))
tail(rowData(spe))
rownames(spe) <- rowData(spe)$Target
spe
colnames(spe) <- paste0(spe$sample_id,"cell", spe$ObjectNumber) # create unique identifiers per cell
colnames(spe)
unique(spe$sample_id)


# Suppose spe_list is a list of spe objects per patient
spe_list <- lapply(spe_list, function(spe) {
  
  # 1️⃣ Clip hot pixels (MAD)
  exprs_clipped <- counts(spe)
  for(ch in rownames(exprs_clipped)) {
    x <- exprs_clipped[ch, ]
    med <- median(x)
    mad_val <- mad(x)
    hot <- x > med + 20 * mad_val
    exprs_clipped[ch, hot] <- med + 20 * mad_val
  }
  
  # Store clipped counts
  assay(spe, "counts_clipped") <- exprs_clipped
  
  return(spe)
})


# 1️⃣ Min-max normalization
normalize_minmax <- function(x) {
  if(max(x) == min(x)) return(rep(0, length(x)))
  (x - min(x)) / (max(x) - min(x))
}



# 2️⃣ asinh transformation
assay(spe, "exprs_asinh") <- asinh(assay(spe, "counts_clipped") / 1)
assay(spe, "exprs_norm") <- t(apply(assay(spe, "exprs_asinh"), 1, normalize_minmax))



spe <- runPCA(spe,
                     exprs_values = "exprs_asinh",  
                     subset_row = rowData(spe)$use_channel, 
                     ncomponents = 28,
                     BSPARAM = ExactParam())

set.seed(230616)
out <- RunHarmony(spe, group.by.vars = "patient_id")

# Check that order of cells is the same
stopifnot(all.equal(colnames(spe), colnames(out)))

reducedDim(spe, "harmony") <- reducedDim(out, "HARMONY")

mat <- reducedDim(spe, "harmony")
#out <- Rphenograph(as.matrix(mat), k = 45)
out <- Rphenograph(mat, k = 45) #45  ## how lower K the smaller the clusters

clusters <- factor(membership(out[[2]]))

spe$clusters <- clusters

set.seed(220228)
spe <- runUMAP(spe, dimred = "harmony", name = "UMAP_harmony") 



# IMC image visualization -------------------------------------------------

masks <- readRDS('/masks.rds')

# Sample images
set.seed(220517)
cur_id <- c("NDMM_4_6") 
cur_masks <- masks[names(masks) %in% cur_id]


#without segmentation
plotCells(cur_masks,
          object = spe  ,  
          cell_id = "Cell_number", img_id = "Patient_ImageNumber_correct", 
          colour_by = "annotation_IMC") 

## with visible segmentation
plotCells(cur_masks, image_title=NULL, scale_bar = NULL, #legend = NULL,
          object = spe, 
          cell_id = "ObjectNumber", img_id = "Patient_ImageNumber_correct", thick = F, outline_by = "Patient_name", 
          colour_by = "annotation_IMC", 
          colour = list(annotation_IMC = metadata(annotation_IMC)$color_vectors$celltype, 
                        Patient_name = metadata(spe)$color_vectors$Patient_name ))

# Patch detection and T cell infiltration ---------------------------------


spe_IMC_combined_with202606 <- buildSpatialGraph(spe_IMC_combined_with202606, img_id = "ROI_or_Image", type = "expansion", threshold = 20)


spe_IMC_combined_with202606  <- patchDetection(spe_IMC_combined_with202606 , 
                                               patch_cells = spe_IMC_combined_with202606$celltype_annotation_transfer_with202606_nounderscore == "Plasma Cell",
                                               img_id = "ROI_or_Image",
                                               expand_by = 1,
                                               min_patch_size = 40,
                                               colPairName = "expansion_interaction_graph",
                                               BPPARAM = MulticoreParam())


patch_df <- colData(spe_IMC_combined_with202606) %>%
  as_tibble() %>%
  filter(!is.na(patch_id)) %>%
  group_by(ObjectName_anonymous_grouped, patch_id) %>%
  summarise(
    patch_size = n(),
    Tcell_count = sum(celltype_annotation_transfer_with202606_nounderscore %in%
                        c("CD4 T Cell", "CD8 T Cell")),
    .groups = "drop"
  ) %>%
  mutate(
    Tcell_present = Tcell_count > 0,
    log_Tcells = log1p(Tcell_count)  # avoids log(0) issue
  )


plot_df <- patch_df %>%
  filter(patch_size >= 50) %>%
  mutate(size_bin = cut(
    patch_size,
    breaks = c(25, 50, 100, 250, Inf),
    include.lowest = TRUE
  )) %>%
  group_by(size_bin, Tcell_present) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(size_bin) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

p1 <- ggplot(plot_df, aes(x = size_bin, y = frac, fill = Tcell_present)) +
  geom_col() +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  theme_classic() +
  labs(
    x = "Patch size bin",
    y = "Fraction of patches"
  )
p1

status_colors <- c(  #fir annot even less broad
  
  "CBM" = "#46F0F0",
  "IMC_CBM" = "#46F0F0",
  "MGUS" = "#984EA3" ,
  "SMM" = "#4575B4",              # Macrophages
  "NDMM" = "#FF7F00", 
  "MM" = "#FF7F00", 
  "IMC_NDMM" = "#46F0F0",
  "PCL" = "#D73027",
  "pPCL" = "#D73027",
  "Relapse" = "#984EA3" #New: yellow
)

patch_df <- colData(spe_IMC_combined_with202606) %>%
  as_tibble() %>%
  filter(!is.na(patch_id)) %>%
  group_by(ObjectName_anonymous_grouped, patch_id, status) %>%
  summarise(
    patch_size = n(),
    .groups = "drop"
  ) %>%
  filter(patch_size >= 45)

order_vec <- c(
  "SMM_1",
  "SMM_2","SMM_3",
  "SMM_4",
  "NDMM_1","NDMM_2",
  "NDMM_3","NDMM_4",
  "NDMM_5","NDMM_6",
  "NDMM_7","NDMM_8",
  "NDMM_9",
  "NDMM_10","NDMM_11",
  "NDMM_12",
  "NDMM_13","NDMM_14",
  "PCL_1",
  "PCL_2","PCL_3",
  "PCL_4",
  "Relapse_1",
  "Relapse_2", "Relapse_3",
  "Relapse_4"
)

patch_df <- patch_df %>%
  mutate(
    ObjectName_anonymous_grouped = factor(ObjectName_anonymous_grouped, levels = order_vec)
  )

ggplot(patch_df, aes(
  x = ObjectName_anonymous_grouped,
  y = patch_size,
  color = status
)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  scale_y_log10() +
  scale_color_manual(values = status_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = NULL,
    y = "Patch size (# cells, log10 scale)",
    color = "Status"
  )


p1_b <- ggplot(patch_df, aes(
  x = ObjectName_anonymous_grouped,
  y = patch_size,
  fill = status
)) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.colour = "black",
    alpha = 0.8
  ) +
  scale_y_log10() +
  scale_fill_manual(values = status_colors) +
  scale_color_manual(values = status_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Patch size (# cells, log10 scale)",
    fill = "Status",
    color = "Status"
  )


p1_b

library(colorspace)

dot_colors <- darken(status_colors, amount = 0.35)

p1_b <- ggplot(patch_df, aes(
  x = ObjectName_anonymous_grouped,
  y = patch_size
)) +
  geom_boxplot(
    aes(fill = status),
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = status),
    width = 0.15,
    alpha = 0.6,
    size = 1.5
  ) +
  scale_y_log10() +
  scale_fill_manual(values = status_colors) +
  scale_color_manual(values = dot_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Patch size (# cells, log10 scale)",
    fill = "status",
    color = "status"
  )


library(dplyr)
library(ggplot2)

## Plot 1
df <- colData(spe_IMC_combined) %>% 
  as_tibble() %>%
  filter(!is.na(patch_id)) %>%
  group_by(patch_id, ROI_or_Image) %>%
  summarise(
    Tcell_count = sum(celltype_annotation_transfer_with202606_nounderscore %in%
                        c("CD8 T Cell", "CD4 T Cell")),
    patch_size = n(),
    .groups = "drop"
  ) %>%
  filter(patch_size >= 25)

fit <- lm(log10(Tcell_count + 1) ~ log10(patch_size), data = df)
slope <- coef(fit)[2]
label_text <- paste0("Slope = ", round(slope, 3))

p2 <- ggplot(df, aes(x = patch_size, y = Tcell_count + 1)) +
  geom_point(alpha = 0.6, size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.4) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Patch size",
    y = "T cell count (+1)"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_text,
    hjust = 1.05,
    vjust = 1.4,
    size = 7
  ) +
  theme_classic(base_size = 22) +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.3, "cm")
  )

p2



## Plot 2
df <- colData(spe_IMC_combined) %>% 
  as_tibble() %>%
  filter(!is.na(patch_id)) %>%
  group_by(patch_id, ROI_or_Image, patient_id) %>%
  summarise(
    Tcell_count = sum(celltype_annotation_transfer_with202606_nounderscore %in%
                        c("CD8 T Cell", "CD4 T Cell")),
    patch_size = n(),
    .groups = "drop"
  ) %>%
  filter(patch_size >= 25) %>%
  mutate(
    Tcell_present = Tcell_count > 0,
    Tcell_ratio = Tcell_count / patch_size
  )

p3 <- ggplot(df, aes(x = patch_size, y = Tcell_ratio)) +
  geom_point(
    alpha = 0.5,
    size = 3.5,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 0.15)) +
  labs(
    x = "Patch size",
    y = "T cell ratio"
  ) +
  theme_classic(base_size = 22) +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.3, "cm")
  )



# SEC vs AEC volcanoplot --------------------------------------------------



Idents(integrated_seurat_all      ) <- integrated_seurat_all    @meta.data$annot_merged_final_with_macro_noslash_dis_lowquali_correct_correct_IFNsubdiv_for_s_curves
markers <- FindMarkers(integrated_seurat_all  , assay = "SCT", ident.1 = "SEC", ident.2 = "AEC" )
#markers_umap_clustered_filt <- markers_umap_clustered[!markers_umap_clustered$gene %in% tumor_genes,]
write.csv(markers , "markers_and_otherCSVfiles/SEC_vs_AEC.csv")

# Define thresholds
sig_threshold <- 0.05
logfc_threshold <- 0.25


aec_vs_sec <- read.csv("~/markers_and_otherCSVfiles/SEC_vs_AEC.csv")


p  <- make_volcano_plot(
  data        = aec_vs_sec,
  log2fc_col  = "avg_log2FC",
  padj_col    = "p_val_adj",
  gene_col    = "X",  # use rownames as gene names, or set to your gene column name
  up_label    = "SEC",
  down_label  = "AEC",
  other_label = "Not significant",
  padj_cut    = 0.05,
  fc_cut      = 0.25,
  color_up    = "#D82638",
  color_down  = "#397FB9", 
  color_other = "#D9D9D9",
  title       = "SEC vs AEC",
  #subtitle    = ct,
  label_all_sig = FALSE,
  # label_genes   = labels_filtered,   # only label these if significant
  label_top_n   = 20,
  max_overlaps  = 50,
  # or FALSE + label_top_n = 50
  seed        = 123
)
print(p)


# Create PDF for all volcano plots
pdf("~/Volcano_plots_AEC_vs_SEC_findmarkers.pdf", width = 8, height = 6)
print(p)
dev.off()



# Nodule growth correlation with outcome ----------------------------------




library(survival)
library(survminer)
library(dplyr)



HVN87_merged_combined <- 
  HVN87_merged_combined %>%
  dplyr::rename(
    ISS   = iss,
    pfsi_1 = pfsi
  )



HVN87_merged_combined$architecture_group_3x <- ifelse(
  HVN87_merged_combined$Nodule_SIZE_FINAL_with_3x_scored %in% c(4, 5),
  "Nodular / sheet-like",
  "Diffuse / mini-nodular"
)



# Prepare metadata

HVN87_merged_no_METADATA <- HVN87_merged_combined[!is.na(HVN87_merged_combined$pfs),]


HVN87_merged_no_METADATA_no_lowquality <- HVN87_merged_no_METADATA[!is.na(HVN87_merged_no_METADATA$Nodule_SIZE_FINAL_average),]



HVN87_merged_no_METADATA_no_lowquality$architecture_group_3x <- ifelse(
  HVN87_merged_no_METADATA_no_lowquality$Nodule_SIZE_FINAL_with_3x_scored %in% c(4, 5),
  "Nodular / sheet-like",
  "Diffuse / mini-nodular"
)





library(dplyr)

meta <- HVN87_merged_no_METADATA_no_lowquality


meta$BMPC_10 <- meta$bmpc_path / 10


meta$reviss <- factor(meta$reviss)



meta$architecture_group_average <- factor(meta$architecture_group_3x)



# Cox model — PFS

cox_pfs <- coxph(
  Surv(pfs, pfsi_1) ~
    architecture_group_average +
    #BMPC_10 +
    reviss, # +

  data = meta
)

summary(cox_pfs)

# Forest plot
p <- ggforest(
  cox_pfs,
  data = meta,
  main = "Multivariate Cox Regression PFS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p




# Cox model — OS

cox_os <- coxph(
  Surv(os, osi) ~
    architecture_group_average +
    #BMPC_10 +
    reviss, #+

  
  data = meta
)

summary(cox_os)


# Forest plot
p <- ggforest(
  cox_os,
  data = meta,
  main = "Multivariate Cox Regression OS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p





# Cox model — PFS

cox_pfs <- coxph(
  Surv(pfs, pfsi_1) ~
    architecture_group_average,
  
  data = meta
)

summary(cox_pfs)

# Forest plot
p <- ggforest(
  cox_pfs,
  data = meta,
  main = "Multivariate Cox Regression PFS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p




# Cox model — OS

cox_os <- coxph(
  Surv(os, osi) ~
    architecture_group_average ,
  
  data = meta
)

summary(cox_os)


# Forest plot
p <- ggforest(
  cox_os,
  data = meta,
  main = "Multivariate Cox Regression OS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p



# Cox model — PFS

cox_pfs <- coxph(
  Surv(pfs, pfsi_1) ~
    BMPC_10,
  
  data = meta
)

summary(cox_pfs)

# Forest plot
p <- ggforest(
  cox_pfs,
  data = meta,
  main = "Multivariate Cox Regression PFS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p


# Cox model — OS

cox_os <- coxph(
  Surv(os, osi) ~
    BMPC_10 ,
  
  data = meta
)

summary(cox_os)


# Forest plot
p <- ggforest(
  cox_os,
  data = meta,
  main = "Multivariate Cox Regression OS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p


# Cox model — PFS

cox_pfs <- coxph(
  Surv(pfs, pfsi_1) ~
    #architecture_group_average +
    BMPC_10 +
    reviss, # +
 
  data = meta
)

summary(cox_pfs)

# Forest plot
p <- ggforest(
  cox_pfs,
  data = meta,
  main = "Multivariate Cox Regression PFS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p



# Cox model — OS

cox_os <- coxph(
  Surv(os, osi) ~
    #architecture_group_average +
    BMPC_10 +
    reviss, #+
  
  
  data = meta
)

summary(cox_os)


# Forest plot
p <- ggforest(
  cox_os,
  data = meta,
  main = "Multivariate Cox Regression OS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p






# Cox model — PFS

cox_pfs <- coxph(
  Surv(pfs, pfsi_1) ~
    architecture_group_average +
    BMPC_10 +
    reviss, # +

  
  data = meta
)

summary(cox_pfs)

# Forest plot
p <- ggforest(
  cox_pfs,
  data = meta,
  main = "Multivariate Cox Regression PFS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p


# Cox model — OS

cox_os <- coxph(
  Surv(os, osi) ~
    architecture_group_average +
    BMPC_10 +
    reviss, #+

  
  data = meta
)

summary(cox_os)


# Forest plot
p <- ggforest(
  cox_os,
  data = meta,
  main = "Multivariate Cox Regression OS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1
)
p




library(survival)
library(survminer)

# KM model
km_os_arch <- survfit(
  Surv(os, osi) ~ architecture_group_average,
  data = meta
)

# Plot
p_km_os_arch <- ggsurvplot(
  km_os_arch,
  data = meta,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time",
  ylab = "Overall survival probability",
  legend.title = "Architecture group",
  palette = c("#9ECAE1", "#08519C"),
  ggtheme = theme_bw(base_size = 14),
  risk.table.height = 0.25
)

# Show
p_km_os_arch



library(survival)
library(survminer)

# KM model — PFS
km_pfs_arch <- survfit(
  Surv(pfs, pfsi_1) ~ architecture_group_average,
  data = meta
)

# Plot with at-risk table
p_km_pfs_arch <- ggsurvplot(
  km_pfs_arch,
  data = meta,
  risk.table = TRUE,
  risk.table.y.text = TRUE,
  risk.table.y.text.col = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time",
  ylab = "Progression-free survival probability",
  legend.title = "Architecture group",
  palette = c("#9ECAE1", "#08519C"),
  ggtheme = theme_bw(base_size = 14),
  risk.table.height = 0.30
)

# Show plot + risk table
p_km_pfs_arch


