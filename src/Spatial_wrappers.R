# Wrapper functions for the spatial transcriptomic study
# The goal of these wrapper functions is to streamline the analysis of 
# visium data given a simple .csv samplesheet inputted

# Damsky Lab
# Written by daccache joe
# date of initial write : april 13 2022
# last updated : april 2 2024

# dependant packages
require(Seurat)
require(Giotto)
require(ggplot2)
require(dplyr)
require(patchwork)
require(data.table)
require(clustree)

# Load in protein-coding genes for analysis
annot <- readr::read_delim("data/homo_sapien_gene_info.txt", delim = "\t", escape_double = F, trim_ws = T, progress = F,)
annot <- annot[annot$type_of_gene == "protein-coding", ]
pc.genes <- unique(annot$Symbol_from_nomenclature_authority)
syn.genes <- unlist(strsplit(annot$Synonyms, split = "|", fixed = TRUE))
syn.genes <- syn.genes[-which(syn.genes == "-")]
full.genes <- c(pc.genes, syn.genes)
rm(annot, syn.genes, pc.genes)


# load Visium data function
LoadVisiumData <- function(dir.path){
    object <- Load10X_Spatial(data.dir = dir.path) # seurat function
    all.features <- rownames(object)
    
    # generate feature-level meta data 
    object[["Spatial"]]@meta.features$nCells <- rowSums(object[["Spatial"]]@counts > 0) 
    object[["Spatial"]]@meta.features$gene_depth <- rowSums(object[["Spatial"]]@counts)/object[["Spatial"]]@meta.features$nCells # ratio of total transcript counts of a gene to the number of spots it is present in (i.e. 1 transcript/spot = 1)
    object[["Spatial"]]@meta.features$perc_cells <- object[["Spatial"]]@meta.features$nCells/ncol(object[["Spatial"]])

    # include genes present in >= 20 spots
    genes.to.remove <- rownames(object[["Spatial"]]@meta.features %>% filter(nCells < 20))
    genes.to.remove <- c(genes.to.remove, rownames(object[["Spatial"]]@meta.features %>% filter(gene_depth == 1))) # remove genes that are only present as one transcript/spot
    genes.to.keep <- rownames(object[["Spatial"]])[-which(rownames(object[["Spatial"]]) %in% genes.to.remove)]
    genes.to.keep <- genes.to.keep[genes.to.keep %in% full.genes]
    
    object <- subset(object,subset =nFeature_Spatial > 250, features = genes.to.keep)
    return(object)
}


# load alignment qc function
LoadAlignmentQC <- function(dir.path){
    metrics_path <- paste0(dir.path, "/metrics_summary.csv")
    df <- read.csv(metrics_path)
    return(df)
}

# visually filter spots function
# user-interactive function that enables the user to drag and drop spots to remove them from the analysis
DragToFilter <- function(object, project = NULL, image_name = "slice1", invert = FALSE, color_def = NULL, pt.size.factor = 2.5){
    library('bslib')
    if(is.null(project)){
        project <- object@project.name
    }
  message("Filtering ", project)
  # Create 'df' of the Visium 10X coorindates
  df <- object@images[[image_name]]@coordinates
  if(!is.null(color_def)){
    df$color_def <- object@meta.data[ , color_def]
  }
  # set CONTINUE to TRUE to begin spot selection loop
  CONTINUE <- TRUE
  while (CONTINUE == TRUE) {
    # Run 'CellSelector' and user highlights spots to remove
    CellsToRemove <-
      CellSelector(
        df %>%
          mutate(imagerow = -imagerow) %>%
          ggplot(aes(x = imagecol, y = imagerow, color = color_def)) +
          geom_point(size = 3) +
          ggtitle(sample) +
          theme_classic()
    )
    
  # If the user highlights a non-zero number of spots to remove, the spots are taken out of the image
    if(length(CellsToRemove) != 0){
      df <- df[! rownames(df) %in% CellsToRemove, ]
    }
    
    # the user is asked if they are done removing spots from that object
    question1 <- readline(paste0("Have you completed selecting spots for ", project, "? [y/n] "))
    if (regexpr(question1, 'n', ignore.case = TRUE) == 1) {
      # if 'no' i.e. not done removing spots, they go back to 'CellSelector'
      CONTINUE <- TRUE
    } else{
      # if 'yes', i.e. done removing spots, they move on subsetting the object
      CONTINUE <- FALSE
    } # end question
  } # end 'while' loop
  
  # extract barcodes of spots that made the cut(s)
  CellsToKeep <- rownames(df)
  # subset the Seurat obejct accordingly
  object <- subset(object, cells = CellsToKeep, invert = invert)
  dev.off()
  # Return the seurat object
  return(object)
}

# normalization function, default method = "SCT"
NormalizeMyData <- function(object, method = "SCT", assay = "Spatial"){
    if(method == "SCT"){
        object <- SCTransform(object, 
                            assay = "Spatial", 
                            verbose = FALSE)
    } else if (method == "LogNormalize"){
        object <- NormalizeData(object, verbose = FALSE)
        object <- FindVariableFeatures(object, verbose = FALSE)
        object <- ScaleData(object, verbose = FALSE)
    } else {
        stop("Please enter 'SCT' or 'LogNormalize' as arguement 'method'")
    }
    return(object)
}

# dimensionality reduction function
RunDimensionalityReduction <- function(object, dims = 30, 
                                       assay = "SCT", resolution = c(0.2, 0.5, 0.8, 1, 1.5), 
                                       reduction = "umap", rm.max.genes=TRUE){
  if(rm.max.genes){
    # want to remove genes present in > 95% of spots
    if(sum(VariableFeatures(object) %in% 
        rownames(object[["Spatial"]]@meta.features)[
          object[["Spatial"]]@meta.features$perc_cells > 0.95]) > 0){
      VariableFeatures(object) <- VariableFeatures(object)[- which(VariableFeatures(object) %in% 
                                                                     rownames(object[["Spatial"]]@meta.features)[
                                                                       object[["Spatial"]]@meta.features$perc_cells > 0.95]
      )]
    }
  }
  if(length(VariableFeatures(object)) == 0){
    stop("Lost all var features in ", object@project.name)
  }
    object <- RunPCA(object, dims = dims, assay = assay, verbose = FALSE)
    object <- FindNeighbors(object, reduction = "pca", dims = 1:dims, verbose = FALSE)
    object <- FindClusters(object,resolution= resolution, verbose = FALSE)
    if(reduction == "tsne"){
        object <- RunTSNE(object, dims = 1:dims, verbose = FALSE)
    } else if (reduction == "umap"){
        object <- RunUMAP(object, dims = 1:dims, verbose = FALSE)
    } else {
        stop("Argument 'reduction' must be either 'tsne' or 'umap'")
    }
    return(object)
}

# identify corresponding sc RNA-Seq reference set function
IdentifySCRNASeqDataset <- function(object=NULL, condition=NULL){
  if(is.null(condition)){
    condition <- unique(object@cell_metadata$cond_abv)
    message("object's condition is identifed as: ", condition)
    if(!(condition %in%  c("CNTRL", "SAR", "GA", "NL", "NXG"))){
      stop("object's condition variable not identified.")
    }
    condition <- ifelse(condition %in% c("GA", "NL"), "GA",
                        ifelse(condition %in% c("SAR"),"SAR", 
                               "CNTRL"))
  } else if(!(condition %in%  c("CNTRL", "SAR", "GA", "NL", "merged"))){
      stop("object's condition variable not identified.")
  }


  if(condition == "merged"){
    load(file = "data/R_objects/single_cell/merged_reference_sc_object.Rds")
    scRNA_object <- merged.sc.object
  } else {
    load("data/R_objects/single_cell/reference_sc_object_list.Rds")
    scRNA_object <- integrated_object_list[[condition]]
    scRNA_object <- SCTransform(scRNA_object)
  }
  return(scRNA_object)
}

# Seurat based cell type prediction function, not run on this data
SeuratDeconvolution <- function(object, reference_object, prediction_vector = "cell_type"){
    anchors <- FindTransferAnchors(reference = reference_object, 
                                    query = object, 
                                    normalization.method = normalization.method, 
                                    verbose = FALSE)
    predictions.assay <- TransferData(anchorset = anchors, 
                                        refdata = reference_object[ , prediction_vector], 
                                        prediction.assay = TRUE,
                                        weight.reduction = object[["pca"]], 
                                        dims = 1:30, 
                                        verbose = FALSE)
    object[["predictions"]] <- predictions.assay
    return(object)
}

# merging Seurat objects (one line but more readable when wrapped)
MergeMySeuratObjects <- function(object_list){
    seurat_merged <- merge(object_list[[1]],
                            object_list[c(2:length(object_list))])
    return(seurat_merged)

}

# Seurat integration - there's more options than what it shown here, this is stock CCA integration and was not used in the study
Integrate_My_Data <- function(object.list){
  object.list <- lapply(object.list, function(x){
    x <- NormalizeData(x, verbose = TRUE)
    x <- FindVariableFeatures(x, verbose = TRUE)
  })
  int.features <- SelectIntegrationFeatures(object.list = object.list)
  int.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = int.features)
  # this command creates an 'integrated' data assay
  integrated.object <- IntegrateData(anchorset = int.anchors)
  return(integrated.object)
}

# preperation of a giotto object from the Seurat object
PrepareGiottoObject <- function(spatial_obj, obj_name){
  message("Working on ", obj_name)
  test_Giotto <- createGiottoObject(raw_exprs = spatial_obj[["Spatial"]]@counts,
                                    spatial_locs = spatial_obj@images$slice1@coordinates[ ,c("imagerow", "imagecol")],
                                    cell_metadata = spatial_obj@meta.data)
  message("Normalizing")
  test_Giotto <- normalizeGiotto(test_Giotto, 
                                scale_genes = T, 
                                scale_cells = T,
                                scalefactor = 6000)
  message("Adding Statistics")
  test_Giotto <- addStatistics(gobject = test_Giotto)
  message("Calculating HVG")
  test_Giotto <- calculateHVG(gobject = test_Giotto,
                              show_plot = FALSE)
  gene_metadata <- fDataDT(test_Giotto)
  featgenes <- gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID

  message("Running PCA")
  test_Giotto <- runPCA(gobject = test_Giotto, 
                        genes_to_use = featgenes, 
                        scale_unit = F, center = T, 
                        method="factominer")
  # message("Running tSNE")
  # test_Giotto <- runtSNE(test_Giotto, 
  #                           dimensions_to_use = 1:10)
  message("Running UMAP")
  test_Giotto <- runUMAP(test_Giotto,
                         dimensions_to_use = 1:30)
  ## sNN network (default)
  message("Creating NN")
  test_Giotto <- createNearestNetwork(gobject = test_Giotto, 
                                        dimensions_to_use = 1:30, 
                                        k = 30)
  ## Leiden clustering
  message("Clustering")
  test_Giotto <- doLeidenCluster(gobject = test_Giotto, 
                                    resolution = 0.5,
                                    n_iterations = 1000)
  message("Processing complete on ", obj_name)
  return(test_Giotto)
}

# Spatial DWLS (Giotto Framework) deconvlution method
SpatialDWLSDeconvolution <- function(giotto_obj, reference_obj, assay = "SCT", prediction_vector = "cell_type"){
    message("Making signature matrix.")
    signature_matrix <- makeSignMatrixDWLSfromMatrix(matrix = reference_obj[[assay]]@data, 
                                                   sign_gene = reference_obj[[assay]]@var.features,
                                                   cell_type_vector = reference_obj@meta.data[ , prediction_vector])
    message("Running Deconvolution algorithm")
    # updating Giotto object to be compatible with Giotto Suite
    giotto_obj <- createGiottoObject(expression = list(raw = giotto_obj@raw_exprs,
                                                        normalized = giotto_obj@norm_expr),
                                      spatial_locs = giotto_obj@spatial_locs,
                                      cell_metadata =  list('cell' = list('rna' = giotto_obj@cell_metadata)),
                                      instructions = giotto_obj@instructions)
    message("Starting DWLS Deconvolution now.")
    giotto_obj <- runDWLSDeconv(gobject = giotto_obj, 
                                spat_unit = "cell",
                                cluster_column = "leiden_clus",
                                sign_matrix = signature_matrix)
    message("Deconvolution complete.")
    return(giotto_obj)
}

# transferring data between Seurat and Giotto objects
TransferSpatialData <- function(seurat_obj, giotto_obj, return = "both"){
    message("Adding meta data to ", seurat_obj@project.name)
    giotto_obj <- subsetGiotto(giotto_obj, cell_ids = colnames(seurat_obj))
    giotto_md <- giotto_obj@cell_metadata
    giotto_md_to_keep <- as.data.frame(giotto_md[, c("cell_ID", "nr_genes", "perc_genes", "total_expr", "leiden_clus")])
    rownames(giotto_md_to_keep) <- giotto_md_to_keep$cell_ID
    giotto_md_to_keep <- giotto_md_to_keep[ , - which(colnames(giotto_md_to_keep) == "cell_ID")]
    seurat_obj <- AddMetaData(seurat_obj, metadata = giotto_md_to_keep)
    
    # commented out bc it's better to keep enrichment, LR, neighborhood etc in giotto
    message("Adding Spatial DWLS predictions to ", seurat_obj@project.name)
    new_mat <- as.matrix(giotto_obj@spatial_enrichment$DWLS[ ,-1])
    rownames(new_mat) <- giotto_obj@spatial_enrichment$DWLS$cell_ID
    new_mat <- t(new_mat)
    seurat_obj[["dwls"]] <- CreateAssayObject(counts = new_mat)
    message("Adding Spatial Gene Clusters to ", seurat_obj@project.name)
    
    if(seurat_obj$cond_abv[1] != "CNTRL"){
      new_mat <- as.matrix(giotto_obj@spatial_enrichment$cluster_metagene[ , -"cell_ID"])
      rownames(new_mat) <- giotto_obj@spatial_enrichment$DWLS$cell_ID
      new_mat <- t(new_mat)
      seurat_obj[["spat_cor_genes"]] <- CreateAssayObject(counts = new_mat)
    }
    
    seurat_obj@reductions$giotto <- CreateDimReducObject(embeddings =
                                                        as.matrix(giotto_obj@dimension_reduction$cells$umap$umap$coordinates),
                                                                    assay = "Spatial",
                                                                    key = "giotto_")
    
   
    message("Adding seurat meta data to giotto")
    giotto_colnames <- colnames(giotto_obj@cell_metadata)
    seurat_colnames <- colnames(seurat_obj@meta.data)
    md_cols_to_add <- seurat_colnames[!(seurat_colnames %in% giotto_colnames)]
    md_cols_to_add <- c(md_cols_to_add, "histological_region", "histological_region_indepth")
    giotto_obj@cell_metadata[ , md_cols_to_add] <- seurat_obj@meta.data[ , md_cols_to_add]

  if(return == "giotto"){
      return(giotto_obj)
  } else if(return == "seurat"){
      return(seurat_obj)
  } else if(return == "both"){
    object.list <- list(seurat_obj = seurat_obj, giotto_obj = giotto_obj)
    return(object.list)
  } else { 
      stop("Argument 'return' should be either 'giotto', 'seurat', or 'both'.")
  }
}

# Spatial Neighborhood analysis function
SpatialNeighborhoodGiotto <- function(giotto_obj, stepsize = 30, minimum_k = 1){
  message("Working on ", giotto_obj@cell_metadata$orig.ident[1])
    # Delaunay network method with minimum k = 1
    if(minimum_k == 1){
        spatial_network_name <- "Delaunay_network"
    } else {
        spatial_network_name <- "spatial_network"
    }
    giotto_obj <- createSpatialGrid(gobject = giotto_obj,
                                      sdimx_stepsize = stepsize,
                                      sdimy_stepsize = stepsize,
                                      minimum_padding = 0)
    giotto_obj <- createSpatialNetwork(gobject = giotto_obj,
                                          minimum_k = minimum_k)

    return(giotto_obj)
}

# Detect spatial gene co-ecpression
DetectSpatialCoExpression <- function(giotto_obj, spatial_network_name = "Delaunay_network"){
  kmeans_spatialgenes <- binSpect(giotto_obj, bin_method = "kmeans")
  ext_spatial_genes <- kmeans_spatialgenes[1:1000]$genes

  spat_cor_netw_DT <- detectSpatialCorGenes(giotto_obj,
                                          method = 'network',
                                          spatial_network_name = spatial_network_name,
                                          subset_genes = ext_spatial_genes)
  # spat_cor_netw_DT$cor_DT$to_keep <- rowSums(is.na(spat_cor_netw_DT$cor_DT))
  # unique(spat_cor_netw_DT$cor_DT$variable[is.na(spat_cor_netw_DT$cor_DT$spat_cor)])
  # spat_cor_netw_DT$cor_DT <- subset(spat_cor_netw_DT$cor_DT, to_keep == 0)
  spat_cor_netw_DT <- clusterSpatialCorGenes(spat_cor_netw_DT,
                                          name = 'spat_netw_clus',
                                          k = 8)
  netw_ranks <- rankSpatialCorGroups(giotto_obj, 
                                     spatCorObject = spat_cor_netw_DT,
                                     use_clus_name = "spat_netw_clus")
  
  cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, 
                                         use_clus_name = 'spat_netw_clus',
                                         show_top_genes = 1)
  cluster_genes <- cluster_genes_DT$clus
  names(cluster_genes) = cluster_genes_DT$gene_ID
  
  # create spatial metagenes and visualize
  giotto_obj <- createMetagenes(giotto_obj,
                                     gene_clusters = cluster_genes, 
                                     name = 'cluster_metagene')
  
  giotto_obj@gene_metadata$spat_cor_cluster <- cluster_genes[match(giotto_obj@gene_metadata$gene_ID, names(cluster_genes))]
  giotto_obj@gene_metadata$spat_cor_cluster[is.na(giotto_obj@gene_metadata$spat_cor_cluster)] <- "None"
  return(giotto_obj)
}


# StackedVlnPlot function taken from 
# https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 90), 
          axis.text.y = element_text(size = rel(1), angle = 90), 
          plot.margin = plot.margin )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          Stacked_Title = "",
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90, hjust = 0.8, vjust = 0.5), 
          axis.ticks.x = element_line())
  plot_list[[1]] <- plot_list[[1]] + ggtitle(Stacked_Title) + BoldTitle()
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

# Here we want the middle range of cells to use in the future as reference sets
# Normally we would do this on a per-sample basis but 
Detect_Outliers <- function(df){
  df$logRatio <- log(df$nFeature_RNA)/log(df$nCount_RNA)
  quantiles <- quantile(df$logRatio)  
  upper_outlier <- quantiles[["75%"]]
  lower_outlier <- quantiles[["25%"]]
  
  df$keep_for_feats <- ifelse(df$nFeature_RNA < 500 | df$nFeature_RNA > 5000, FALSE, TRUE)
  df$is.Between <- between(df$logRatio, lower_outlier, upper_outlier)
  df$is.Good <- df$keep_for_feats & df$is.Between
  message("There are ", sum(!df$is.Good), " outliers detected, leaving ", sum(df$is.Good), " good cells.")
  return(df)
}

# volcano plot using ggplot2
VolcanoPlot <- function(df, 
              label = FALSE,
              label.quantile = 0.9,
              unique.labels = FALSE,
              facet = FALSE, 
              only.pos = FALSE,
              plot.title = NULL,
              genes_to_label = NULL,
              color_low="#666666",
              color_high = "#FF0000",
              color_def = NULL,
              discrete_colors = NULL,
              p.val.label.cutoff = 0.05){
  if(unique.labels & !label){
    label <- TRUE
  }
  
  if(!("gene" %in% colnames(df))){
    df$gene <- rownames(df)
  }
  # data prep
  if(only.pos){
    color_def <- color_def[df$avg_log2FC > 0]
    df <- df[df$avg_log2FC > 0, ]
  }
  
  if(label){
    df <- 
      df %>%
      mutate(label = ifelse(abs(avg_log2FC) > quantile(abs(avg_log2FC), label.quantile), gene, NA),
             label = ifelse(p_val_adj < p.val.label.cutoff, label, NA),
             PI = abs(avg_log2FC)*(-log(p_val_adj)))
    if(facet){
      df <- 
        df %>%
        mutate(group = ifelse(avg_log2FC < 0, -1, 1)) %>%
        group_by(cluster, group) %>%
        mutate(rank = rank(-PI),
               label = ifelse(group == 1 & rank <= 30, gene, NA))
               # label = ifelse(group == -1 & rank <= 10, gene, label))
    }
    

    
    if(!is.null(genes_to_label)){
      df <- 
        df %>%
        mutate(label = ifelse(gene %in% genes_to_label, gene, NA))
      if(unique.labels){
        df <- 
          df %>%
          mutate(label = ifelse(avg_log2FC > 0, label, NA),
                 label = ifelse(p_val_adj < 0.05, label, NA))
      } 

      # df <- 
      #   df %>%
      #   mutate(label = ifelse(gene %in% genes_to_label, gene, NA),
      #          label = ifelse(avg_log2FC > 0, label, NA),
      #          label = ifelse(rank <= 50, label, NA))
    } else {
      if(unique.labels){
        genes_to_label <-
          df %>%
          group_by(cluster) %>%
          mutate(to_keep = ifelse(abs(avg_log2FC) > quantile(abs(avg_log2FC), label.quantile), TRUE, FALSE)) %>%
          filter(to_keep) %>%
          ungroup() %>%
          group_by(gene) %>%
          summarise(counts = n()) %>%
          filter(counts == 1) %>%
          pull(gene)
        df <- 
          df %>%
          mutate(label = ifelse(gene %in% genes_to_label, gene, NA))
        df <- 
          df %>%
          mutate(label = ifelse(avg_log2FC > 0, label, NA),
                 label = ifelse(p_val_adj < 0.05, label, NA))
      }
    }
  }
  
  p <- 
    df %>%
    ggplot(aes(x = avg_log2FC, y = -log(p_val_adj), color = color_def)) +
    geom_point() +
    theme_classic() 
  
  if(facet){
    p <- p + facet_grid(~ cluster, scales = "free_x")
  }
  
  if(!is.null(plot.title)){
    p <- p + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5)) + BoldTitle()
  }
  if(label){
    p <- p + ggrepel::geom_label_repel(aes(label = label, 
                                           color = color_def), 
                                       box.padding = 0, 
                                       label.padding = 0.1, 
                                       na.rm = TRUE, 
                                       max.overlaps = 25)
    
  }
  
  if(!is.null(color_def)){
    if(!is.null(discrete_colors)){
      p <- p + scale_color_manual(values = discrete_colors)
    } else {
      my_colors <- scale_color_gradient(low = color_low, high = color_high)
      p <- p  +
        my_colors
    }

  }
  
  return(p)
  
}

# find unique genes from list of cluster markers, I'm sure theres an easier way to do this
FindUniqueGenes <- function(marker_list){
  label_genes <-
    bind_rows(marker_list, .id = "id") %>%
    group_by(gene) %>%
    summarise(counts = n()) %>%
    filter(counts == 1) %>%
    pull(gene)
  return(label_genes)
}
