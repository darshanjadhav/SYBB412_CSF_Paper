# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----   Single cell transcriptomics reveals cerebrospinal fluid immune   -----
# -----  dysregulation during healthy brain aging and cognitive impairment -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 10-11-2021
# Written by: Emma Tapp, Natalie Piehl
# Summary: Remove GEX background contamination with SoupX 
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Hmisc")
  library("Seurat")
  library("Matrix")
  library("SoupX")
  library("DropletUtils")
})

# Define inputs
gex_dir <- "/home/dxj294/Project/new"
soupx_dir <- "/home/mxg979/Project/SYBB412/SoupX"

# Source helper functions
source("00_helper_functions.R")

# Create output directory
dir.create(soupx_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Run SoupX

# Create lists of directories to load as SoupX objects and Seurat objects
sample_dirs <- list.dirs(gex_dir, recursive = FALSE)

# Initialize gene list to estimate contamination fraction
# Selecting for monocyte/dendritic markers which are highly specific and abundant
M.genes <- c("CD14", "CD68", "MS4A7", "FCGR3A")

# Initialize contamination fraction dataframe 
contamination_frac <- setNames(data.frame(matrix(ncol = 2, nrow = 0)),
                               c("sample", "contamination_frac"))

# For each sample...
for (dir in sample_dirs) {
  # Isolate sample name
  sample <- unlist(strsplit(dir, "/")) %>%
    tail(1)
  
  # Move to next iteration if sample F5
  if (sample == "F5") {
    next
  }
  
  if (!(startsWith(sample,"CSF_MCI_AD_E") || startsWith(sample,"CSF_MCI_AD_G") || startsWith(sample,"CSF_MCI_AD_H"))) {
    next
  }
  
  # Print what sample is being processed
  print(paste0("Processing sample ", sample))
  
  # Load in SoupX and Seurat data
  soupx <- load10X(paste0(dir, ""))
  seurat <- Read10X(paste0(dir, "/filtered_feature_bc_matrix")) %>%
    CreateSeuratObject()
  seurat
  
  # Normalize and run PCA on Seurat object
  seurat <- SCTransform(seurat, variable.features.n = 2000, verbose = TRUE) %>%
    RunPCA()
  
  # Generate clusters
  seurat <- RunTSNE(seurat, dims = 1:10) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters(resolution = 0.3)
  
  # Add cluster info to SoupX object
  soupx <- setClusters(soupx, setNames(seurat$seurat_clusters, colnames(seurat)))
  
  # Estimate contamination fraction
  useToEst <- estimateNonExpressingCells(soupx, nonExpressedGeneList = list(M.genes = M.genes))
  soupx <- calculateContaminationFraction(soupx,
                                          list(M.genes = M.genes),
                                          useToEst = useToEst)
  
  # Save contamination fraction
  contamination_frac[nrow(contamination_frac) + 1, ] <- c(sample, soupx$metaData$rho[1])
  
  # Create adjusted counts
  adj_counts <- adjustCounts(soupx)
  
  # Save Corrected Counts
  dir.create(paste0(soupx_dir, "/", sample), showWarnings = FALSE)
  DropletUtils:::write10xCounts(paste0(soupx_dir, "/", sample), adj_counts, overwrite = TRUE)
}

# Print contamination fractions
contamination_frac
