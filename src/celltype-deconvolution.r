# since this was run on a cluster, provide it with the condition/sample you are working on 
args <- commandArgs(trailingOnly=TRUE)
SampleCondition <- args[1]

require(Seurat)
require(dplyr)
require(ggplot2)
require(Giotto)
source(file = "scripts/Spatial_Pipeline/Spatial_wrappers.R")

# Read in the samplesheet
samplesheet <- read.csv("good_samplesheet.csv")
samplesheet <- samplesheet[samplesheet$cond_abv == SampleCondition, ]
message("Deconvolution of ", nrow(samplesheet), " ", SampleCondition, " objects is about to begin.")

# Identify scRNA-Seq reference dataset
sc_reference_set <- IdentifySCRNASeqDataset(condition = "merged")

for(SampleName in SampleNames){
  message("Working on ", SampleName)
  load(file = paste0("data/R_objects/spatial/giotto/",SampleName,"_processed_giotto.Rds"))

  # Running Spatial DWLS, Giotto's deconvolution alg 
  deconv_giotto_obj <- SpatialDWLSDeconvolution(giotto_obj, reference_obj = sc_reference_set)
  save(deconv_giotto_obj, file = paste0("data/R_objects/spatial/giotto/", SampleName,"_deconv_giotto.Rds"))
  message("Processing on ", SampleName, " complete. Moving on.")
}
message("Deconvolution of ", nrow(samplesheet), " now complete. Quitting R.")
