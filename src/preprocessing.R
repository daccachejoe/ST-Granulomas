library(Seurat)
library(dplyr)
library(ggplot2)
library(Giotto)
source(file = "scripts/Spatial_Pipeline/Spatial_wrappers.R")

# read in samplesheet
samplesheet <- read.csv("samplesheet.csv")

# if actually running this script, would plug in T/F values depending on what you want to do
# for now, just setting up the script to run nothing
# some steps require user input so they are all set to FALSE for code sharing purposes
Create_Objects <- FALSE
Combined <- FALSE
Clustering <- FALSE
Clusters_Defined <- FALSE
Annotated <- FALSE
save_objects <- FALSE
post_subset <- FALSE
TestMarkers <- FALSE

# this chunk creates seurat objects from visium data, users filters out spots, and plots some qc metrics
if(Create_Objects){
  # Create our Seurat objects
  object.list <- list()
  for(i in 1:nrow(samplesheet)){
    SampleName <- samplesheet$SampleName[i]
    
    # LoadVisiumData is a UDF in Spatial_wrappers.R
    object.list[[SampleName]] <- LoadVisiumData(dir.path = samplesheet$path[i])
    
    object.list[[SampleName]]@project.name <- SampleName
    object.list[[SampleName]]$orig.ident <- SampleName
    for(col_to_add in md_cols_to_add){
      object.list[[SampleName]][[col_to_add]] <- samplesheet[i, col_to_add]
    }
  }
  
  # Export Alignment QC metrics
  qc_dfs <- lapply(samplesheet$path, LoadAlignmentQC)
  names(qc_dfs) <- samplesheet$SampleName
  df_combined <- bind_rows(qc_dfs, .id = "SampleName")
  write.csv(df_combined, "output/Alignment_qc_DF.csv")
  
  # Filtering out known problematic spots by dragging them out
  # we removed any spots not adjacent to the bulk of the tissue
  object.list <- lapply(object.list, DragToFilter)
  
  # Plot common qc metrics
  # start by combining all the metadata
  combined_md <- bind_rows(
    lapply(object.list, 
           function(x){
             return(x@meta.data)
           }
    )
  )
  
  # recreating the 10X quality control after the R-based spot filtering steps
  export_df <- 
    combined_md %>% 
    group_by(orig.ident) %>%
    select(nCount_Spatial, nFeature_Spatial) %>%
    summarise_all("mean")
  write.csv(export_df, "output/post_R_filter_qc_metrics.csv")
  
  # sequencing depth via nGenes and nReads
  pdf("output/preliminary_qc_plots.pdf", height = 8, width = 20)
  combined_md %>%
    ggplot(aes(x = orig.ident, y = nFeature_Spatial, fill = orig.ident))+
    geom_boxplot() +
    theme_classic() +
    ggtitle("Unique genes per spot") +
    BoldTitle() +
    theme(plot.title = element_text(hjust = 0.5))
  combined_md %>%
    ggplot(aes(x = orig.ident, y = nCount_Spatial, fill = orig.ident))+
    geom_boxplot() +
    theme_classic()+
    ggtitle("total reads per spot") +
    BoldTitle() +
    theme(plot.title = element_text(hjust = 0.5))
  combined_md %>%
    ggplot(aes(x = orig.ident, y = log(nFeature_Spatial)/log(nCount_Spatial), fill = orig.ident))+
    geom_boxplot() +
    theme_classic()+
    ggtitle("Log ratio of unique genes to total reads per spot") +
    BoldTitle() +
    theme(plot.title = element_text(hjust = 0.5))
  combined_md %>%
    ggplot(aes(x = nCount_Spatial, y = nFeature_Spatial, color = orig.ident))+
    geom_point() +
    theme_classic()+
    ggtitle("Relationship between reads and unique reads per spot") +
    BoldTitle() +
    theme(plot.title = element_text(hjust = 0.5))
  dev.off()
}

# if we are not combining samples
# this chunk runs the clustering pipeline on each sample
# to 1) determine the optimal resolution  
# 2) provide plots and a dataframe to annotate clusters  
# 3) read in the annotations and re-cluster if necessary
# 4) generate plots and markers
if(!Combined){
    # this is the clustering step to select optimal resolution
  if(Clustering){ 
    # Normalization: default SCTransform
    object.list <- lapply(object.list, NormalizeMyData)
  
    # Dimensionality Reduction : PCA, Clustering, Visualization (tSNE or UMAP)
    object.list <- lapply(object.list,function(obj){ 
      obj <- RunDimensionalityReduction(obj, resolution = c(0.5, 0.8, 1, 1.5))
    })
    
    # plotting all resolutions to select best ones
    plot.list <- lapply(object.list, function(obj){
      p1 <- SpatialDimPlot(obj, group.by = paste0("SCT_snn_res.", c(0.5, 0.8, 1, 1.5)), pt.size.factor = 3) + ggtitle(obj@project.name)
      p2 <- clustree(obj, prefix = "SCT_snn_res.", layout = "sugiyama")
      p3 <- DimPlot(obj, group.by = paste0("SCT_snn_res.", c(0.5, 0.8, 1, 1.5)))
      return(list(p1, p2, p3))
    })
  
    pdf("output/selecting_resolutions.pdf", height = 12, width = 12)
    plot.list <- unlist(plot.list, recursive = FALSE)
    plot.list
    dev.off()
    
    # selecting a resolution by writing into this csv file for the next step
    # can also edit csv file to remove the slide by selecting resoltuion "rm"
    export_df <- data.frame(sample = samplesheet$SampleName)
    write.csv(export_df, "output/selecting_resolutions.csv")
  
    } else if(Clusters_Defined){ # this plugs in when the user has selected a resolution for each slice 
  res_to_use <- read.csv("output/selecting_resolutions.csv") # read in the now user-edited csv file with resolutions per slide
  res_to_use <- split(res_to_use$resolution, f = res_to_use$sample, drop = T) # convert to list
  
  # assigning the determined appropriate resolution to each object
  # or removing objects flagged for removal with an "rm"
  object.list <- 
    lapply(object.list, function(obj){
      Idents(obj) <- paste0("SCT_snn_res.", res_to_use[[obj@project.name]])
      if(res_to_use[[obj@project.name]] == "rm"){
        warning("Sample ", obj@project.name, " flagged to remove.")
        return(NULL)
      }
      return(obj)
  })
  object.list <- object.list[-which(unlist(lapply(object.list, is.null)) == TRUE)]

  # To annotate the clusters, we first export the plots colored as clusters: tissue pictures must be included in order to not have scale or resolution reduced
  plot.list <- lapply(names(object.list), function(sample){
    message(sample)
    png_image <-  png::readPNG(source = paste0(samplesheet$path[match(sample, samplesheet$SampleName)], "/spatial/tissue_hires_image.png")) # high res image better
    
    plot.md <- as.data.frame(object.list[[sample]]@images$slice1@coordinates[ ,c("imagerow", "imagecol")])
    if(!(all(rownames(plot.md) == rownames(object.list[[sample]]@meta.data)))){
      stop("barcodes not the same")
    }

    plot.md$cluster <- object.list[[sample]]@meta.data[ , paste0("SCT_snn_res.", res_to_use[[sample]])]
    p1 <-
      plot.md %>%
      ggplot(aes(x = imagecol, y = -imagerow, fill = cluster)) +
      geom_point(shape = 21, 
                size = 3, 
                color = "black") +
      theme_classic() +
      ggtitle(object.list[[sample]]@project.name) +
      theme(axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text = element_blank(),
            axis.ticks = element_blank()) + 
      BoldTitle()
    p2 <-
      ggplot(plot.md, aes(x = imagerow, y = imagecol)) +    # Modify image file
      geom_point(col = "white") +
      xlim(- 5, 5) +
      ylim(- 5, 5) +
      theme(axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      annotation_custom(grid::rasterGrob(png_image, width = 1, height = 1),
                        xmin = - 4.5, xmax = Inf,
                        ymin = - Inf, ymax = Inf)
    return(list(p1, p2))
  })
    
  pdf("output/seurat_cluster_dimplots.pdf", height = 24, width = 12)
  ggpubr::ggarrange(plotlist = unlist(plot.list, recursive = FALSE), 
                    nrow = 4, ncol = 2, widths = c(.45, .55))
  dev.off()
  
  plot.list <- lapply(object.list, function(obj){
    p <- SpatialDimPlot(obj, pt.size.factor = 2.5,
                        alpha = 0.5, image.alpha = 1) + 
      ggtitle(obj@project.name) +
      theme(plot.title = element_text(hjust = 0.5))
    return(p)
  })
  pdf("output/seurat_cluster_dimplots_w_images.pdf", height = 12, width = 12)
  ggpubr::ggarrange(plotlist = plot.list, 
                    nrow = 2, ncol = 2)
  dev.off()
  
  plot.list <- lapply(object.list, function(obj){
    p <- VlnPlot(obj, features = c("nCount_Spatial", "nFeature_Spatial"), stack = T, 
                 pt.size = 0.1) + 
      NoLegend() +
      ggtitle(obj@project.name)
    return(p)
  })
  pdf("output/sequencing_violins.pdf", width = 10, height = 6)
  plot.list
  dev.off()
  
  # this "form" is exported so that user can define annotations for each cluster in each sample
  cluster_data_frames <- lapply(object.list, function(seu){
    df <- data.frame(cluster = sort(unique(Idents(seu))),
                    hf = rep("", length((unique(Idents(seu))))))
  })
  clusters_data_frame <- bind_rows(cluster_data_frames, .id = "sample")
  write.csv(clusters_data_frame, file = "output/cluster_iden_df_form.csv",row.names = FALSE)
  
  save(object.list, file = "data/R_objects/spatial/seurat/preprocessed_seurat_list.Rds")  # save the objects

  } else if(Annotated){ # now that the objects are annotated, the dataframes are read in for annotation and subsequent re-clustering if any clusters are removed
  # load in annotations
  if(post_subset){ # if the user has removed any clusters, we read in the post-subsetted dataframe
    clusters_data_frame <- read.csv("output/cluster_iden_df_form_post_subset.csv")
  } else {
    clusters_data_frame <- read.csv("output/cluster_iden_df_form.csv")
  }
    # add our annotations, stored in "histological_region" and "histological_region_indepth" to the objects
  object.list <- lapply(object.list, function(seu){
    md <- seu@meta.data
    ref_df <- clusters_data_frame[clusters_data_frame$sample == seu@project.name, ]
    hf_to_add <- ref_df$hf[match(md[ ,paste0("SCT_snn_res.", res_to_use[[seu@project.name]])], ref_df$cluster)]
    seu <- AddMetaData(seu, hf_to_add, col.name = "histological_region")
    hf_to_add <- ref_df$hf.2[match(md[ ,paste0("SCT_snn_res.", res_to_use[[seu@project.name]])], ref_df$cluster)]
    seu <- AddMetaData(seu, hf_to_add, col.name = "histological_region_indepth")
    return(seu)
  })

    # identify any objects that have clusters annotated with "rm" flagging them for removal
  objects_to_recluster <-
    clusters_data_frame %>%
    filter(hf == "rm") %>%
    pull(sample) %>%
    unique()
  
  # remove any problematic spots post-clustering
  object.list <- lapply(object.list, function(seu){
    seu <- subset(seu, histological_region %in%  
                    c("rm",
                      "Gland",
                      "Eccrine Gland",
                      "Muscle",
                      "Epidermis",
                      "Fat"), invert = TRUE)
    print(SpatialDimPlot(seu, group.by = "histological_region", pt.size.factor = 3))

    # if there are any spots that require manual filtering, we do so here
    question1 <- readline(paste0("Are there any spots that require manual filtering in ", seu@project.name, "? [y/n] "))
    if (regexpr(question1, 'y', ignore.case = TRUE) == 1){
      seu <- DragToFilter(seu, color_def = "histological_region")
      if(!(seu@project.name %in% objects_to_recluster)){
        objects_to_recluster <- c(objects_to_recluster, seu@project.name)
      }
    }
    return(seu)
  })
  
  # In most cases we removed clusters of places we werent interested in
  # So we comment this out bc it was histological reasons
    # re-clustering objects that have been subsetted
  if(length(objects_to_recluster) > 0){
    object.list[objects_to_recluster] <- lapply(object.list[objects_to_recluster], NormalizeMyData)
    object.list[objects_to_recluster] <- lapply(object.list[objects_to_recluster],function(obj){ 
      obj <- RunDimensionalityReduction(obj, resolution = c(0.5, 0.8, 1, 1.5))
    })
  
    clusters_data_frame <- clusters_data_frame[-which(clusters_data_frame$sample %in% objects_to_recluster), ]
    new_clusters_df <- lapply(object.list[objects_to_recluster], function(seu){
    df <- data.frame(cluster = sort(unique(seu$SCT_snn_res.0.8)),
                    hf = rep("", length(unique(seu$SCT_snn_res.0.8))))
    })
    new_clusters_df <- bind_rows(new_clusters_df, .id = "sample")
    clusters_data_frame <- rbind(clusters_data_frame, new_clusters_df)
    write.csv(clusters_data_frame, file = "output/cluster_iden_df_form_post_subset.csv",row.names = FALSE)
    
    plot.list <- lapply(object.list[objects_to_recluster], function(obj){
      p <- SpatialDimPlot(obj) + 
        ggtitle(obj@project.name) +
        theme(plot.title = element_text(hjust = 0.5))
      return(p)
    })
    pdf("output/seurat_cluster_dimplots_w_images_subsetted_objects.pdf", height = 12, width = 12)
    ggpubr::ggarrange(plotlist = plot.list, 
                      nrow = 2, ncol = 2)
    dev.off()
  }
  
  # if the user specified that plots are to be made
  if(ToPlot){
    plot.list <- 
      lapply(object.list, function(obj){
        p1 <- SpatialDimPlot(obj, 
                             group.by = "histological_region",
                             pt.size.factor = 2,
                             image.alpha = 1, 
                             alpha = 0.4,
                             label.size = 2.5,
                             label = T) + 
          ggtitle(obj@project.name)
        p2 <- DimPlot(obj, group.by = "histological_region", label = T, repel = T) + NoLegend()
        return(list(p1, p2))
    })
    
    pdf("output/seurat_annotated_dimplots.pdf", height = 12, width = 12)
    plot.list <- unlist(plot.list,recursive = FALSE)
    ggpubr::ggarrange(plotlist = plot.list, nrow = 2, ncol = 2)
    dev.off()
  
  
  }
  
  # if markers are to be generated
  if(TestMarkers){
    marker.list <- lapply(object.list, function(obj){
      Idents(obj) <- "histological_region"
      M <- FindAllMarkers(obj, only.pos = TRUE)
      return(M)
    })
    WriteXLS::WriteXLS(marker.list, ExcelFileName = "output/histological_region_markers.xlsx", SheetNames = names(marker.list))
  }
  if(save_objects){
    save(object.list, file = "data/R_objects/spatial/seurat/preprocessed_seurat_list.Rds")

    for(name in names(object.list)){
    object <- object.list[[name]]
    save(object, file = paste0("data/R_objects/spatial/seurat/", name, "_preprocessed_Seurat.Rds"))
    }
  }

  }
}

# To be run AFTER the individual approach that filters out spots  
# follows a similar pipeline to the individual approach, but merges like-condition samples together
if(Combined){
  if(Clustering){
    # Merging like-condition samples together
    sample.merged.list <- list(CNTRL = object.list[grep("CNTRL", names(object.list))],
                               GA = object.list[grep("GA", names(object.list))],
                               NL = object.list[grep("NL", names(object.list))],
                               SAR = object.list[grep("SAR", names(object.list))],
                               NXG = object.list[grep("NXG", names(object.list))]])
    sample.merged.list <- lapply(sample.merged.list, function(my_list){
      if(length(my_list) > 1){
        var.feats <- unique(unlist(lapply(my_list, function(obj){
          var_feats <- VariableFeatures(obj)[1:1500] # take the top 1500 variable genes for merging
        })))
        merged.object <- merge(my_list[[1]], my_list[c(2:length(my_list))])
        DefaultAssay(merged.object) <- "SCT"
        VariableFeatures(merged.object) <- var.feats # assign the union of all samples' variable genes as the merged object variable genes
        return(merged.object)
      } else {
        return(my_list[[1]])
      }
    })
    
    # Running Dimensionality Reduction Pipeline on each merged object
    sample.merged.list <- lapply(sample.merged.list, RunDimensionalityReduction)
    
    # Add new cluster defintions to original objects
    object.list <- lapply(names(object.list), function(obj_name){
      obj <- object.list[[obj_name]]
      merged_reference <- sample.merged.list[[obj$cond_abv[1]]]
      md.to.add <- merged_reference@meta.data[ , c("orig.ident", paste0("SCT_snn_res.", c(0.2, 0.5, 0.8, 1, 1.5)))]
      md.to.add <- subset(md.to.add, orig.ident == obj$orig.ident[1])
      rownames(md.to.add) <- sapply(strsplit(rownames(md.to.add), split = "_"), getElement, 1)
      obj <- AddMetaData(obj, metadata = md.to.add)
      return(obj)
    })
    names(object.list) <- unlist(lapply(object.list, function(obj){return(obj@project.name)}))
    
    # plotting all resolutions to select best ones
    plot.list <- lapply(object.list, function(obj){
      plot.list.mini <- lapply(paste0("SCT_snn_res.", 
                                      c(0.2, 0.5, 0.8, 1, 1.5)), 
                               function(res){
                                        p <- SpatialDimPlot(obj,
                                                            group.by = res,
                                                            pt.size.factor = 2.5)
                                        return(p)
                                        })
      names(plot.list.mini) <- c("0.2", "0.5", "0.8", "1", "1.5")
      return(plot.list.mini)
    })
    
    export.plot.list <- list()
    for(cond in names(sample.merged.list)){
      export.plot.list[[cond]] <- list(
        plot.list[grep(cond, names(plot.list))]
      )
      export.plot.list[[cond]] <- unlist(export.plot.list[[cond]], recursive = FALSE)
      export.plot.list[[cond]] <- unlist(export.plot.list[[cond]], recursive = FALSE)
      export.plot.list[[cond]] <- lapply(paste0(c("0.2", "0.5", "0.8", "1", "1.5"),"$"), function(res){
        p.list <-  export.plot.list[[cond]][grep(res, names(export.plot.list[[cond]]))]
        return(p.list)
      })
    }
    
    residual.plot.list <- lapply(sample.merged.list, function(obj){
      p1 <- clustree(obj, prefix = "SCT_snn_res.", layout = "sugiyama")
      p2 <- DimPlot(obj, group.by = paste0("SCT_snn_res.", c(0.2, 0.5, 0.8, 1, 1.5)))
      return(list(p1, p2))
    })
    
    pdf("output/selecting_resolutions_merged.pdf", height = 12, width = 12)
    for(my.cond in names(export.plot.list)){
      
        print(ggpubr::ggarrange(plotlist = lapply(export.plot.list[[my.cond]],
                                            function(p.list){
                                              ggpubr::ggarrange(plotlist = p.list, 
                                                                nrow = 1, 
                                                                ncol = length(p.list), 
                                                                common.legend = TRUE, 
                                                                legend = "right")
                                              }),
                          nrow = 2))
      print(residual.plot.list[[my.cond]][[1]])
      print(residual.plot.list[[my.cond]][[2]])
    }
    dev.off()
    
    # selecting a resolution
    export_df <- data.frame(sample = names(sample.merged.list),
                            resolution = rep("", length(sample.merged.list)))
    write.csv(export_df, "output/selecting_resolutions_merged.csv", row.names = FALSE)
  
    
    } else if(Clusters_Defined){
    res_to_use <- read.csv("output/selecting_resolutions_merged.csv")
    res_to_use <- split(res_to_use$resolution, f = res_to_use$sample, drop = T)
    sample.merged.list <- 
      lapply(sample.merged.list, function(obj){
        Idents(obj) <- paste0("SCT_snn_res.", res_to_use[[obj$cond_abv[1]]])
        if(res_to_use[[obj$cond_abv[1]]] == "rm"){
          warning("Sample ", obj$cond_abv[1], " flagged to remove.")
          return(NULL)
        }
        return(obj)
      })
    
    object.list <- 
      lapply(object.list, function(obj){
        Idents(obj) <- paste0("SCT_snn_res.", res_to_use[[obj$cond_abv[1]]])
        if(res_to_use[[obj$cond_abv[1]]] == "rm"){
          warning("Sample ", obj$cond_abv[1], " flagged to remove.")
          return(NULL)
        }
        return(obj)
      })
    
    # To annotate the clusters, we first export the plots colored as clusters: tissue pictures must be included in order to not have scale or resolution reduced
    color_vector <- scales::hue_pal()(max(unlist(lapply(sample.merged.list, function(obj){length(unique(Idents(obj)))}))))
    names(color_vector) <- 0:(length(color_vector)-1)
    
    plot.list <- lapply(names(object.list), function(sample){
      png_image <-  png::readPNG(source = paste0(samplesheet$path[match(sample, samplesheet$SampleName)], "/spatial/tissue_hires_image.png"))
      plot.md <- as.data.frame(object.list[[sample]]@images$slice1@coordinates[ ,c("imagerow", "imagecol")])
      if(!(all(rownames(plot.md) == rownames(object.list[[sample]]@meta.data)))){
        stop("barcodes not the same")
      }
      
      plot.md$cluster <- object.list[[sample]]@meta.data[ , paste0("SCT_snn_res.", res_to_use[[obj$cond_abv[1]]])]
      p1 <-
        plot.md %>%
        ggplot(aes(x = imagecol, y = -imagerow, fill = cluster)) +
        geom_point(shape = 21, 
                   size = 3, 
                   color = "black") +
        theme_classic() +
        ggtitle(object.list[[sample]]@project.name) +
        theme(axis.title = element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.text = element_blank(),
              axis.ticks = element_blank()) + 
        BoldTitle() +
        scale_fill_manual(values = color_vector)
      p2 <-
        ggplot(plot.md, aes(x = imagerow, y = imagecol)) +    # Modify image file
        geom_point(col = "white") +
        xlim(- 5, 5) +
        ylim(- 5, 5) +
        theme(axis.line = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()) +
        annotation_custom(grid::rasterGrob(png_image, width = 1, height = 1),
                          xmin = - 4.5, xmax = Inf,
                          ymin = - Inf, ymax = Inf)
      return(list(p1, p2))
    })
    
    pdf("output/seurat_cluster_dimplots_merged.pdf", height = 24, width = 12)
    ggpubr::ggarrange(plotlist = unlist(plot.list, recursive = FALSE), 
                      nrow = 4, ncol = 2, widths = c(.45, .55))
    dev.off()
    
    plot.list <- lapply(names(sample.merged.list), function(my.cond){
      obj.list <- object.list[grep(my.cond, names(object.list))]
      p.list <- lapply(obj.list, function(obj){
        p <- 
          SpatialDimPlot(obj, 
                         pt.size.factor = 2.5,
                          alpha = 0.6, 
                          image.alpha = 1,
                         label = T) + 
        ggtitle(obj@project.name) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = color_vector)
      return(p)
      })
      return(p.list)
    })
    
    pdf("output/seurat_cluster_dimplots_w_images_merged.pdf", height = 12, width = 12)
    for(i in 1:length(plot.list)){
      print(
        ggpubr::ggarrange(plotlist = plot.list[[i]],
                         nrow = ifelse(length(plot.list[[i]]) > 2, 2, 1), 
                         ncol = ifelse(length(plot.list[[i]]) > 2, 3, length(plot.list[[i]])), 
                         common.legend = TRUE, 
                         legend = "top")
      )
    }
    dev.off()
    
    plot.list <- lapply(sample.merged.list, function(obj){
      p <- VlnPlot(obj, features = c("nCount_Spatial", "nFeature_Spatial"), stack = T, 
                   pt.size = 0.1) + 
        NoLegend() +
        ggtitle(obj$cond_abv[1])
      return(p)
    })
    pdf("output/sequencing_violins_merged.pdf", width = 10, height = 6)
    plot.list
    dev.off()
    
    cluster_data_frames <- lapply(sample.merged.list, function(seu){
      df <- data.frame(cluster = sort(unique(Idents(seu))),
                       hf = rep("", length((unique(Idents(seu))))))
    })
    clusters_data_frame <- bind_rows(cluster_data_frames, .id = "sample")
    write.csv(clusters_data_frame, file = "output/cluster_iden_df_form_merged.csv",row.names = FALSE)
    
    cluster.marker.list <- lapply(sample.merged.list, function(obj){
      obj <- PrepSCTFindMarkers(obj)
      Idents(obj) <- paste0("SCT_snn_res.", res_to_use[[obj$cond_abv[1]]])
      M <- FindAllMarkers(obj, only.pos = TRUE)
      return(M)
    })
    cluster.marker.list <- cluster.marker.list[!unlist(lapply(cluster.marker.list, length) == 0)]
    openxlsx::write.xlsx(cluster.marker.list, "output/cluster_markers_merged.xlsx")
    
    save(sample.merged.list, file = "data/R_objects/spatial/seurat/preprocessed_merged_seurat_list.Rds")
  
    } else if(Annotated){
      # load in annotations
      clusters_data_frame <- read.csv("output/cluster_iden_df_form_merged.csv")
     
      sample.merged.list <- lapply(sample.merged.list, function(seu){
        md <- seu@meta.data
        ref_df <- clusters_data_frame[clusters_data_frame$sample == seu$cond_abv[1], ]
        hf_to_add <- ref_df$hf[match(md[ ,paste0("SCT_snn_res.", res_to_use[[seu$cond_abv[1]]])], ref_df$cluster)]
        seu <- AddMetaData(seu, hf_to_add, col.name = "histological_region")
        # hf_to_add <- ref_df$hf.2[match(md[ ,paste0("SCT_snn_res.", res_to_use[[seu$cond_abv[1]]])], ref_df$cluster)]
        # seu <- AddMetaData(seu, hf_to_add, col.name = "histological_region_indepth")
        return(seu)
      })
      object.list <- lapply(object.list, function(seu){
        md <- seu@meta.data
        ref_df <- clusters_data_frame[clusters_data_frame$sample == seu$cond_abv[1], ]
        hf_to_add <- ref_df$hf[match(md[ ,paste0("SCT_snn_res.", res_to_use[[seu$cond_abv[1]]])], ref_df$cluster)]
        seu <- AddMetaData(seu, hf_to_add, col.name = "histological_region")
        # hf_to_add <- ref_df$hf.2[match(md[ ,paste0("SCT_snn_res.", res_to_use[[seu$cond_abv[1]]])], ref_df$cluster)]
        # seu <- AddMetaData(seu, hf_to_add, col.name = "histological_region_indepth")
        return(seu)
      })
      
      objects_to_recluster <-
        clusters_data_frame %>%
        filter(hf == "rm") %>%
        pull(sample) %>%
        unique()
      
      # remove any problematic spots post-clustering
      sample.merged.list <- lapply(sample.merged.list, function(seu){
        seu <- subset(seu, histological_region %in%  
                        c("rm",
                          "Gland",
                          "Eccrine Gland",
                          "Muscle",
                          "Epidermis",
                          "Fat"), invert = TRUE)
        return(seu)
      })
      object.list <- lapply(object.list, function(seu){
        seu <- subset(seu, histological_region %in%  
                        c("rm",
                          "Gland",
                          "Eccrine Gland",
                          "Muscle",
                          "Epidermis",
                          "Fat"), invert = TRUE)
        return(seu)
      })
      
      clusters_data_frame <- subset(clusters_data_frame, hf != "rm")
      color_vector <- scales::hue_pal()(length(unique(clusters_data_frame$hf)))
      names(color_vector) <- unique(clusters_data_frame$hf)
      
      plot.list <- lapply(names(sample.merged.list), function(my.cond){
        obj.list <- object.list[grep(my.cond, names(object.list))]
        p.list <- lapply(obj.list, function(obj){
          p <- 
            SpatialDimPlot(obj, 
                           group.by = "histological_region",
                           pt.size.factor = 2.5,
                           alpha = 0.6, 
                           image.alpha = 1,
                           label.size = 3,
                           label = T, repel = T) + 
            ggtitle(obj@project.name) +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_manual(values = color_vector)
          return(p)
        })
        return(p.list)
      })
      # plot.list.2 <- lapply(names(sample.merged.list), function(my.cond){
      #   obj.list <- object.list[grep(my.cond, names(object.list))]
      #   p.list <- lapply(obj.list, function(obj){
      #     p <- 
      #       SpatialDimPlot(obj, 
      #                      group.by = "histological_region_indepth",
      #                      pt.size.factor = 2.5,
      #                      alpha = 0.6, 
      #                      image.alpha = 1,
      #                      label.size = 3,
      #                      label = T, repel = T) + 
      #       ggtitle(obj@project.name) +
      #       theme(plot.title = element_text(hjust = 0.5)) +
      #       scale_fill_manual(values = color_vector)
      #     return(p)
      #   })
      #   return(p.list)
      # })
      
      pdf("output/seurat_annotated_dimplots_merged.pdf", height = 12, width = 12)
      for(i in 1:length(plot.list)){
        print(
          ggpubr::ggarrange(plotlist = plot.list[[i]],
                            nrow = ifelse(length(plot.list[[i]]) > 2, 2, 1), 
                            ncol = ifelse(length(plot.list[[i]]) > 4, 3, 2), 
                            common.legend = TRUE, 
                            legend = "top")
        )
      }
      # for(i in 1:length(plot.list.2)){
      #   print(
      #     ggpubr::ggarrange(plotlist = plot.list.2[[i]],
      #                       nrow = ifelse(length(plot.list.2[[i]]) > 2, 2, 1), 
      #                       ncol = ifelse(length(plot.list.2[[i]]) > 4, 3, 2), 
      #                       common.legend = TRUE, 
      #                       legend = "top")
      #   )
      # }
      dev.off()
      
      plot.list <- lapply(sample.merged.list, function(seu){
        p1 <- DimPlot(seu, group.by = "histological_region", label = T, repel = T) + NoLegend() +scale_color_manual(values = color_vector)
        # p2 <- DimPlot(seu, group.by = "histological_region_indepth", label = T, repel = T) + NoLegend()+scale_color_manual(values = color_vector)
        p3 <- DimPlot(seu, group.by = "orig.ident", label = F)
        # return(list(p1, p2, p3))
        return(list(p1, p3))
      })
      # pdf("output/merged_seurat_umap_plots.pdf", height = 8, width = 16)
      pdf("output/merged_seurat_umap_plots.pdf", height = 8, width = 12)
      for(i in 1:length(plot.list)){
        print(
          # ggpubr::ggarrange(plotlist = plot.list[[i]], nrow = 1, ncol = 3)
          ggpubr::ggarrange(plotlist = plot.list[[i]], nrow = 1, ncol = 2)
        )
      }
      dev.off()
      
      if(TestMarkers){
        marker.list <- lapply(sample.merged.list, function(obj){
          obj <- PrepSCTFindMarkers(obj)
          Idents(obj) <- "histological_region"
          M <- FindAllMarkers(obj, only.pos = TRUE)
          return(M)
        })
        marker.list <- marker.list[!unlist(lapply(marker.list, length) == 0)]
        openxlsx::write.xlsx(marker.list, "output/histological_region_markers_merged.xlsx")
        marker.df <- bind_rows(marker.list, .id = "sample")
        
        # indepth.marker.list <- lapply(sample.merged.list, function(obj){
        #   obj <- PrepSCTFindMarkers(obj)
        #   Idents(obj) <- "histological_region_indepth"
        #   M <- FindAllMarkers(obj, only.pos = TRUE)
        #   return(M)
        # })
        # indepth.marker.list <- indepth.marker.list[!unlist(lapply(indepth.marker.list, length) == 0)]
        # openxlsx::write.xlsx(indepth.marker.list, "output/cluster_markers_merged.xlsx")
        # indepth.marker.df <- bind_rows(indepth.marker.list, .id = "sample")
        
      } else {
        Markers <- list(marker.list = "output/histological_region_markers_merged.xlsx",
                        indepth.marker.list = "output/cluster_markers_merged.xlsx")
        Markers <- lapply(Markers, function(filename){
          sheets <- readxl::excel_sheets(filename)
          x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
          names(x) <- sheets
          x <- bind_rows(x, .id = "sample")
          x$condition <- sapply(strsplit(x$sample, split = "_"), getElement, 1)
          return(x)
        })
        marker.df <- Markers[["marker.list"]]
        indepth.marker.df <- Markers[["indepth.marker.list"]]
        rm(Markers)
      }
      
      plot.list <- lapply(sample.merged.list, function(seu){
        my_cond <- seu$cond_abv[1]
        if(!(my_cond %in% marker.df$sample)){return(NULL)}
        genes.to.plot <- 
          marker.df %>%
          filter(sample == my_cond) %>%
          mutate(PI = -log(p_val_adj)*abs(avg_log2FC)) %>%
          group_by(cluster) %>%
          top_n(n = 20, wt = PI) %>%
          top_n(n = 10, wt = avg_log2FC) %>%
          mutate(cluster = as.character(cluster)) %>%
          arrange(cluster) %>%
          pull(gene)
        seu <- ScaleData(seu, assay = "SCT", features = genes.to.plot)
        seu$plot_var <- paste(seu$histological_region, seu$orig.ident)
        
        p <-
          DoHeatmap(seu, features = genes.to.plot, assay = "SCT",group.by ="plot_var",draw.lines = F,raster = F,
                    size = 4, angle = 90) +
          ggtitle(seu$cond_abv[1])
        return(p)
      })
      plot.list <- plot.list[-which(unlist(lapply(plot.list, is.null)) == TRUE)]

      
      # plot.list.2 <- lapply(sample.merged.list, function(seu){
      #   my_cond <- seu$cond_abv[1]
      #   if(!(my_cond %in% indepth.marker.df$sample)){return(NULL)}
      #   genes.to.plot <- 
      #     indepth.marker.df %>%
      #     filter(sample == my_cond) %>%
      #     mutate(PI = -log(p_val_adj)*abs(avg_log2FC)) %>%
      #     group_by(cluster) %>%
      #     top_n(n = 20, wt = PI) %>%
      #     top_n(n = 10, wt = avg_log2FC) %>%
      #     mutate(cluster = as.character(cluster)) %>%
      #     arrange(cluster) %>%
      #     pull(gene)
      #   seu <- ScaleData(seu, assay = "SCT", features = genes.to.plot)
      #   seu$plot_var <- paste(seu$histological_region_indepth, seu$orig.ident)
        
      #   p <-
      #     DoHeatmap(seu, features = genes.to.plot, assay = "SCT",group.by = "plot_var",draw.lines = F,raster = F, 
      #               size = 4, angle = 90) +
      #     ggtitle(seu$cond_abv[1]) +
      #     NoLegend()
      #   return(p)
      # })
      # plot.list.2 <- plot.list.2[-which(unlist(lapply(plot.list.2, is.null)) == TRUE)]
      
      pdf("output/cluster_heatmaps.pdf", height = 12, width = 12)
      ggpubr::ggarrange(plotlist = plot.list, nrow = 1, ncol = 1)
      # ggpubr::ggarrange(plotlist = plot.list.2, nrow = 1, ncol = 1)
      dev.off()
      
      if(save_objects){
        save(sample.merged.list, file = "data/R_objects/spatial/seurat/preprocessed_merged_seurat_list.Rds")
        
        for(name in names(sample.merged.list)){
          object <- sample.merged.list[[name]]
          save(object, file = paste0("data/R_objects/spatial/seurat/", name, "_merged_preprocessed_Seurat.Rds"))
        }
        for(name in names(object.list)){
          object <- object.list[[name]]
          save(object, file = paste0("data/R_objects/spatial/seurat/", name, "_preprocessed_Seurat.Rds"))
        }
      }
      
    
    } else if(TestDomains){
        # spatially_variable_genes <- lapply(sample.merged.list, function(seu){
        #   regions_of_interest <- c("Fibrotic",
        #                            "Grnauloma",
        #                            "Peri-Granular Space",
        #                            "Inter-Granuloma Space")
        #   if(sum(seu$histological_region %in% regions_of_interest) == 0){return(NULL)}
        #   DefaultAssay(seu) <- "SCT"
        #   seu <- ScaleData(seu, features = VariableFeatures(seu))
        #   # seu <- subset(seu, histological_region %in% regions_of_interest)
        #   seu <- FindSpatiallyVariableFeatures(seu,
        #                                        assay = "SCT",
        #                                        selection.method = "markvariogram",
        #                                        features = VariableFeatures(seu)[1:1000])
        #   var.feats <- SpatiallyVariableFeatures(seu, selection.method = "markvariogram",decreasing = TRUE)
        #   return(var.feats)
        # })
        # spatially_variable_genes <- spatially_variable_genes[
        #   - which(unlist(lapply(spatially_variable_genes, is.null)))]
        # openxlsx::write.xlsx(spatially_variable_genes, file = "output/spatially_variable_genes.xlsx")
        
      filename <- "output/spatially_variable_genes.xlsx"
      sheets <- readxl::excel_sheets(filename)
      spatially_variable_genes <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, col_names = FALSE))
      spatially_variable_genes <- lapply(spatially_variable_genes, function(X) X <- X$...1)
      names(spatially_variable_genes) <- sheets
        
        
        plot.list <- lapply(names(spatially_variable_genes), function(my_cond){
          seu <- sample.merged.list[[my_cond]]
          names(seu@images) <- unique(seu$orig.ident)
          var.feats <- spatially_variable_genes[[my_cond]]
          var.feats <- var.feats[- grep("^KRT", var.feats)]
          genes.to.plot <- var.feats[1:5]
          p1 <- SpatialFeaturePlot(seu, 
                                  features = genes.to.plot, 
                                  pt.size.factor = 2,
                                  alpha = c(0.1, 1))
          genes.to.plot <- var.feats[6:10]
          p2 <- SpatialFeaturePlot(seu, 
                                   features = genes.to.plot, 
                                   pt.size.factor = 2,
                                   alpha = c(0.1, 1))
          chemo_cyto_kines <- var.feats[c(grep("^CCL", var.feats), grep("^CXCL", var.feats))]
          p3 <- SpatialFeaturePlot(seu, 
                                   features = chemo_cyto_kines[1:(ceiling(length(chemo_cyto_kines)/2))], 
                                   pt.size.factor = 2,
                                   alpha = c(0.1, 1))
          p4 <- SpatialFeaturePlot(seu, 
                                   features = chemo_cyto_kines[(ceiling(length(chemo_cyto_kines)/2) +1):length(chemo_cyto_kines)], 
                                   pt.size.factor = 2,
                                   alpha = c(0.1, 1))
          return(list(p1, p2, p3, p4))
        })
        pdf("output/spatially_variable_genes.pdf", height = 15, width = 12)
        plot.list
        dev.off()
        
      }
  
}


message("Pre-processing complete.")
