# R script for SpatialPCA analysis of spatial transcriptomic project
library(SpatialPCA)
library(Seurat)
library(future)
library(ggplot2)
library(dplyr)
library(slingshot)

# load in data
object.list <- readRDS("object.list.for.figures.RDS")

# rename meta data column to be less characters
change.var.vec <- c("1. Macrophage", "2. Dermis-INFL","3. Dermis-ADJ","4. Necrobiosis","5. Dermis-CONT","6. Papillary Dermis")
names(change.var.vec) <- c("1M", "2DL", "3DU", "4N", "6C", "5E")
color_vector <- c("#CC0000","#FFCC33", "#CC99FF","#33FF00", "#999999", "#0000CC", "black")
names(color_vector) <- c(names(change.var.vec),"segment") 

object.list <- 
  lapply(object.list, function(obj.list){
    seu <- obj.list[["seurat_obj"]]
    seu$plot.var <- names(change.var.vec)[match(seu$histological_region, change.var.vec)]
    obj.list[["seurat_obj"]] <- seu
    return(obj.list)
  })

conditions.to.test <- c("GA", "NL", "SAR", "NXG")
objects.to.test <-
  names(object.list)[unlist(lapply(conditions.to.test, function(cond) {
    grep(cond, names(object.list))
  }))]



s.pca.list <-
  lapply(objects.to.test,
         function(obj.name) {
           seu <- object.list[[obj.name]][["seurat_obj"]]
           
           rawcount <- seu[["Spatial"]]@counts
           location <- seu[["slice1"]]@coordinates
           location <-
             as.matrix(location[, c("imagerow", "imagecol")])
           
           ST <-
             CreateSpatialPCAObject(
               counts = rawcount,
               location = location,
               project = obj.name,
               gene.type = "spatial",
               sparkversion = "spark",
               gene.number = 3000,
               customGenelist = NULL,
               min.loctions = 20,
               min.features = 20
             )
           
           ST <-
             SpatialPCA_buildKernel(ST, kerneltype = "gaussian", bandwidthtype = "SJ")
           ST <-
             SpatialPCA_EstimateLoading(ST, fast = FALSE, SpatialPCnum = 20)
           ST <- SpatialPCA_SpatialPCs(ST, fast = FALSE)
           
           return(ST)
         })
names(s.pca.list) <- objects.to.test

pt.object.list <-
  lapply(objects.to.test,
         function(obj.name) {
           seu <- object.list[[obj.name]][["seurat_obj"]]
           rawcount <- seu[["Spatial"]]@counts
           location <- seu[["slice1"]]@coordinates
           location <-
             as.matrix(location[, c("imagerow", "imagecol")])
           
           ST <- s.pca.list[[obj.name]]
           
           clusterlabel <-
             seu[["plot.var"]][rownames(ST@location),]
           # trajectory on the whole tissue slice
           sim <- SingleCellExperiment(assays = rawcount)
           reducedDims(sim) <- SimpleList(DRM = t(ST@SpatialPCs))
           colData(sim)$Walktrap <- factor(clusterlabel)
           # in this data we set macrophage region as start cluster
           sim <-
             slingshot(
               sim,
               clusterLabels = 'Walktrap',
               reducedDim = 'DRM',
               start.clus = "1M"
             )
           
           # # focus on macrophage rich region #NOT RUN
           # mac_ind <- which(clusterlabel %in% c("1M"))
           # sim_mac <-
           #   SingleCellExperiment(assays = rawcount[, mac_ind])
           # reducedDims(sim_mac) <-
           #   SimpleList(DRM = t(ST@SpatialPCs[, mac_ind]))
           # colData(sim_mac)$Walktrap <-
           #   factor(clusterlabel[mac_ind])
           # sim_mac <-
           #   slingshot(
           #     sim_mac,
           #     clusterLabels = 'Walktrap',
           #     reducedDim = 'DRM',
           #     start.clus = "1M"
           #   )
           
           
           if(sum(grepl("slingPseudotime_2" ,colnames(colData(sim)))) > 0){
             sim$PT.vec <- 
               ifelse(
                 is.na(sim$slingPseudotime_1), 
                 sim$slingPseudotime_2, 
                 ifelse(
                   is.na(sim$slingPseudotime_2),
                   sim$slingPseudotime_1,
                   ifelse(sim$slingPseudotime_2 > sim$slingPseudotime_1, 
                          sim$slingPseudotime_2,
                          sim$slingPseudotime_1)
                 )
               )
           } else {
             sim$PT.vec <- sim$slingPseudotime_1
           }
           
           # pseudotime_traj1_mac <-
           #   sim_mac@colData@listData$slingPseudotime_1
           # clusterlabels_mac <- clusterlabel[mac_ind]
           
           pseudotime_traj1 <-
             sim$PT.vec
           # pseudotime_traj1[mac_ind] <- NA
           # pseudotime_traj1[mac_ind] <- pseudotime_traj1_mac
           
           return(pseudotime_traj1)
         })
names(pt.object.list) <- objects.to.test

lapply(c(1:10), function(gridnum){
  pt.plot.list <-
    lapply(objects.to.test,
           function(obj.name) {
             pseudotime_traj1 <- pt.object.list[[obj.name]]
             seu <- object.list[[obj.name]][["seurat_obj"]]
             location <- seu[["slice1"]]@coordinates
             location <-
               as.matrix(location[, c("imagerow", "imagecol")])
             ST <- s.pca.list[[obj.name]]
  
             p_traj1 <- plot_trajectory(
               pseudotime_traj1,
               location,
               seu[["plot.var"]][rownames(ST@location),],
               gridnum,
               color_in,
               pointsize = 3.5,
               arrowlength = 0.3,
               arrowsize = 1.3,
               textsize = 15
             )
             return(p_traj1)
           })
  
  pdf(file = paste0("PT-plots-no-ind-mac-gd-",gridnum,".pdf"), height = 8, width = 12)
  lapply(pt.plot.list, function(p_traj1) {
    print(ggpubr::ggarrange(p_traj1[[4]], p_traj1[[1]],
                      ncol = 2, nrow = 1))
  })
  dev.off()
})

# applying it to the seurat objects
seu.list <-
  lapply(objects.to.test,
         function(obj.name){
           ST <- s.pca.list[[obj.name]]
           md.to.add <- as.data.frame(pt.object.list[[obj.name]], row.names = rownames(ST@location))
           colnames(md.to.add) <- "PT"
           seu <- object.list[[obj.name]][["seurat_obj"]]
           seu <- AddMetaData(seu, metadata = md.to.add)
           return(seu)
         })
names(seu.list) <- objects.to.test


p1 <-
  bind_rows(
  lapply(
    seu.list,
    function(seu){
      return(seu@meta.data)
    })) %>%
  mutate(slice.num = sapply(strsplit(orig.ident, split = "_"), getElement, 2)) %>%
  group_by(cond_abv, slice.num) %>%
  mutate(PT = PT/max(PT)) %>%
  filter(plot.var != "5E") %>%
  ggplot(aes(x = PT,y = plot.var, fill = plot.var)) +
  ggridges::geom_density_ridges(alpha = 1) +
  facet_grid(slice.num~cond_abv,
             scales = "free") +
  theme_classic() +
  scale_fill_manual(values = color_vector) 
p1.bulk <- 
  bind_rows(
    lapply(
      seu.list,
      function(seu){
        return(seu@meta.data)
      })) %>%
  group_by(cond_abv) %>%
  mutate(PT = PT/max(PT)) %>%
  filter(plot.var != "5E") %>%
  ggplot(aes(x = PT, y = plot.var, fill = plot.var)) +
  ggridges::geom_density_ridges(alpha = 1) +
  facet_grid(~cond_abv) +
  theme_classic() +
  scale_fill_manual(values = color_vector)

pdf("PT-denstiy-plots.pdf", height = 6, width = 8)
ggpubr::ggarrange(p1, p1.bulk, nrow = 2, ncol = 1, heights = c(0.7, 0.3), align = "hv")
dev.off()

p2 <-
  bind_rows(
  lapply(
    seu.list,
    function(seu){
      md.to.export <- seu[[c("PT","orig.ident", "cond_abv")]]
      md.to.export <- cbind(md.to.export, as.data.frame(t(seu[["dwls"]]@counts)))
      return(md.to.export)
    })) %>%
  tidyr::pivot_longer(cols = c("Myeloid", "Tcell", "Fibroblast")) %>%
  mutate(slice.num = sapply(strsplit(orig.ident, split = "_"), getElement, 2)) %>%
  dplyr::select(PT,slice.num, cond_abv, name, value) %>%
  group_by(cond_abv, slice.num) %>%
  mutate(PT = PT/max(PT)) %>%
  mutate(value = ifelse(value < 0, 0, value)) %>%
  ggplot(aes(y = value, x = PT, fill = name, color = name)) +
  geom_smooth(alpha = 0.5) +
  facet_grid(cond_abv~slice.num,
             scales = "free") +
    theme_classic() +
  scale_x_continuous(labels = scales::percent)
p3 <-
  bind_rows(
  lapply(
    seu.list,
    function(seu){
      md.to.export <- seu[[c("PT", "orig.ident", "cond_abv")]]
      md.to.export <- cbind(md.to.export, as.data.frame(t(seu[["SCT"]]@data[c("CD68", "CD3E", "COL1A1"), ])))
      return(md.to.export)
    })) %>%
  tidyr::pivot_longer(cols = c("CD68", "CD3E", "COL1A1")) %>%
  mutate(slice.num = sapply(strsplit(orig.ident, split = "_"), getElement, 2)) %>%
  dplyr::select(PT,slice.num, cond_abv, name, value) %>%
  group_by(cond_abv, slice.num) %>%
  mutate(PT = PT/max(PT)) %>%
  ggplot(aes(y = value, x = PT, fill = name, color = name)) +
  geom_smooth(alpha = 0.5) +
  facet_grid(cond_abv~slice.num,
             scales = "free") +
  theme_classic() +
  scale_x_continuous(labels = scales::percent)

pdf("PT-gene-plots.pdf", height = 4, width = 6)
p2
p3
dev.off()

