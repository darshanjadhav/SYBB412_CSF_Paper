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
# Date: 06-14-2022
# Written by: Natalie Piehl
# Summary: DE on diagnosis in advanced age only
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggpubr")
  library("ggrepel")
  library("ggthemes")
  library("scales")
  library("grid")
  library("UpSetR")
  library("doMC")
})

# Initialize input parameters
seurat_object <- "/home/mxg979/Project/SYBB412/seurat_object_05s_tcrclean_new"
output_dir <- "/home/mxg979/Project/SYBB412/results/Fig2J/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Source helper functions
source("/home/mxg979/Project/SYBB412/csf_aging-main_local/code/0_preprocessing/00_helper_functions.R")

# Set thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Set core number for parallel model fitting
registerDoMC(cores = 6)

# Load Seurat object
load(seurat_object)

#-------------------------------------------------------------------------------
# Run DE on diagnosis in advanced age and generate volcano plot (Fig 2J)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)


assay5 <- as(s[["RNA"]], Class = "Assay5")
seurat5 <- CreateSeuratObject(assay5, meta.data = s@meta.data)
s <- JoinLayers(seurat5)

# Subset for advanced only
s <- subset(s, age_bin == "advanced")
print(table(s[["age_bin"]]))

# Set Ident to Diagnosis
s <- SetIdent(s, value = "Diagnosis")
print(table(s[["Diagnosis"]]))

# Define Diagnosis Groups
AD_samples <- c("E1", "B4", "B6", "G6", "B7", "C5")
MCI_samples <- c("A6", "B1", "H1", "E4", "E7", "E5", "D2", "D3")

# Add diagnosis to metadata
s@meta.data$Diagnosis2<- rep("HControl", nrow(s@meta.data))
s@meta.data$Diagnosis2[
  which(s@meta.data$ID %in% c(MCI_samples))] <- "MCI"
s@meta.data$Diagnosis2[
  which(s@meta.data$ID %in% c(AD_samples))] <- "AD"

run_de <- function(cell_type) {
  print(cell_type)
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Find DEGs b/w MCI/AD and HC
  degs <-FindMarkers(object = s,
                     ident.1 = "CI",
                     ident.2 = "HC",
                     latent.vars = c("sex"),
                     test.use = "MAST",
                     logfc.threshold = -Inf,
                     min.pct = 0.1,
                     assay = "RNA"
  )
  
  # Remove ribosomal, mitochondrial, and HLA genes
  degs <- degs[-grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = rownames(degs)),]
  
  # Run Benjamini-Hochberg adjustment
  degs$BH <- p.adjust(degs$p_val, method = "BH")
  print(head(degs))
  
  # Write out results
  write.csv(degs, paste0(output_dir, cell_type_label, "_advancedonly_MCIAD_vs_HC_degs.csv"))
  
  # Create volcano plot
  volcano_plot(degs, title = paste0("MCI/AD vs HC (Advanced only)\nin ", cell_type),
               file = paste0(output_dir, cell_type_label, "_advancedonly_MCIAD_vs_HC_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
}

# Run DE on all celltypes
cell_types <- "all_cell"
mclapply(cell_types, run_de, mc.cores = 6)
