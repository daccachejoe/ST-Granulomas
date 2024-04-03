library(Seurat)
library(dplyr)
library(ggplot2)
library(Giotto)
source(file = "scripts/Spatial_Pipeline/Spatial_wrappers.R") # wrapper functions read in

samplesheet <- read.csv("samplesheet.csv")

# load in seurat objects processed in preprocessing.R
# these are qc-filtered, clustered, and annotated seurat objects
object.list <- list()
for(name in samplesheet$SampleName){
  load(file = paste0("data/R_objects/spatial/seurat/", name, "_preprocessed_Seurat.Rds"))
  object.list[[name]] <- object
}

# create giotto object list
giotto.list <- lapply(object.list, function(obj){
  giotto_obj <- PrepareGiottoObject(spatial_obj = obj, obj_name = obj@project.name)
  return(giotto_obj)
})

# run spatial network analysis
giotto.list <- lapply(giotto.list, SpatialNeighborhoodGiotto)

# Detecting spatial correlated genes
giotto.list[-grep("CNTRL", names(giotto.list))] <- lapply(giotto.list[-grep("CNTRL", names(giotto.list))], 
                                                          DetectSpatialCoExpression)

for(name in names(giotto.list)){
  giotto_obj <- giotto.list[[name]]
  save(giotto_obj, file = paste0("data/R_objects/spatial/giotto/", name, "_processed_giotto.Rds"))
}
save(giotto.list, file = "data/R_objects/spatial/giotto/preprocessed_giotto_list.Rds")
