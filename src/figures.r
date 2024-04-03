library(Seurat)
library(dplyr)
library(ggplot2)
library(Giotto)
source(file = "scripts/Spatial_Pipeline/Spatial_wrappers.R")

# Read in the samplesheet
samplesheet <- read.csv("samplesheet.csv")

# Loading in the processed Seurat objects
seurat.list <- lapply(Samples, function(name){
  load(file = paste0("data/R_objects/spatial/seurat/", name, "_preprocessed_Seurat.Rds"))
  return(object)
})
names(seurat.list) <- Samples

# Transferring giotto and seurat data from each object to each other
# load in each giotto object - both deconvoluted and not
giotto_obj_list <- lapply(Samples, function(name){
  load(file = paste0("data/R_objects/spatial/giotto/",name,"_deconv_giotto.Rds"))
  load(file = paste0("data/R_objects/spatial/giotto/",name,"_processed_giotto.Rds"))
  giotto_obj@spatial_enrichment$DWLS <- deconv_giotto_obj@spatial_enrichment$cell$rna$DWLS
  return(giotto_obj)
})
names(giotto_obj_list) <- Samples

# transferring data between objects and generating comprehensive object.list
object.list <- list()
for(name in Samples){
  object.list[[name]] <- TransferSpatialData(seurat_obj = seurat.list[[name]],
                                             giotto_obj = giotto_obj_list[[name]])
}

# save object list
saveRDS(object.list, "object.list.for.figures.RDS")

# loading the merged objects by condition (from preprocessing but in the merged method)
# we premirality have this for the variable features stored in the condition-specific object
load(file = "data/R_objects/spatial/seurat/preprocessed_merged_seurat_list.Rds")
sample.merged.list <- lapply(sample.merged.list, function(seu){
  my_cond <- seu$cond_abv[1]
  samples <- names(seurat.list)[grep(my_cond, names(seurat.list))]
  if(length(samples) > 1){
    var.feats <- lapply(seurat.list[samples],
                        function(mini_obj){
                          mini.var.feats <- VariableFeatures(mini_obj)[1:1500]
                          return(mini.var.feats)})
  } else {
    var.feats <- VariableFeatures(seurat.list[[samples]])
  }

  var.feats <- unique(unlist(var.feats))
  VariableFeatures(seu) <- var.feats
  return(seu)
})

# creating reference data frame with pt size factors for each object
point_size_df <- data.frame(sample = names(object.list),
                            pt.size = numeric(length = length(object.list)))
point_size_df$pt.size <- unlist(lapply(object.list, function(my_list){
  obj <- my_list[["seurat_obj"]]
  row_max <- max(obj@images$slice1@coordinates[ , c("row")])
  col_max <- max(obj@images$slice1@coordinates[ , c("col")])
  if(row_max > col_max){
    my_range <- row_max - min(obj@images$slice1@coordinates[ , c("col")])
  } else if(row_max < col_max){
    my_range <- col_max - min(obj@images$slice1@coordinates[ , c("row")])
  } else {
    row_min <- min(obj@images$slice1@coordinates[ , c("row")])
    col_min <- min(obj@images$slice1@coordinates[ , c("col")])
    if(row_min < col_min){
      my_range <- row_max - min(obj@images$slice1@coordinates[ , c("row")])
    } else {
      col_min <- min(obj@images$slice1@coordinates[ , c("col")])
    }
  }
  ideal_ratio <-  0.0189067
  ideal_pt_size <- ideal_ratio/Radius(obj@images$slice1)
  return(ideal_pt_size)
}))

# Figure 1B: Spatial Dim Plots of representative samples and UMAPs of full conditions
rep_samples <- c("NL_3", "GA_5", "SAR_1", "NXG_1")
plot.list <- lapply(rep_samples, function(obj_name){
  
  p <- SpatialDimPlot(object.list[[obj_name]][["seurat_obj"]],
                      group.by = "plot.var",
                      pt.size.factor = point_size_df$pt.size[point_size_df$sample == obj_name],
                      image.alpha = 0,
                      # stroke = NA,
                      label = FALSE) +
    scale_fill_manual(values = color_vector) +
    theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(4, 'cm'),
          legend.title = element_blank())
  return(p)
})
pdf("output/annotated_slide_images_representative.pdf", height = 8, width = 12)
ggpubr::ggarrange(plotlist = plot.list, nrow = 2, ncol = 2, common.legend = TRUE)
dev.off()

# plotting the UMAPs of each condition
plot.list <- lapply(sample.merged.list, function(seu){
  p1 <- 
    DimPlot(seu,
            reduction = "umap",
            group.by = "plot.var",
            label = T) + 
    scale_color_manual(values = color_vector) +
    ggtitle(seu$cond_abv[1]) +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 12))
  p2 <- 
    DimPlot(seu,
            reduction = "umap",
            group.by = "orig.ident") + 
    ggtitle(seu$cond_abv[1]) +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          legend.position = "top")
  my_list <- list(p1, p2)
  names(my_list) <- c("CT", "batch")
  return(my_list)
})
plot.list <- unlist(plot.list, recursive = FALSE)

pdf("output/annotated_umap_plots.pdf", height = 7, width = 28)
ggpubr::ggarrange(plotlist = 
                    plot.list[grep("CT", names(plot.list))],
                  nrow = 1, common.legend = TRUE)
ggpubr::ggarrange(plotlist = 
                    plot.list[grep("batch", names(plot.list))],
                  nrow = 1
dev.off()


# Figure 1F
# spatial feature plots of important genes
# first list is of rep samples then the rest are ordered normally
genes_to_plot <- c( "CD68", "CD3E", "COL1A1")
f1.samples <- c("GA_1", "NL_3", "SAR_1", "NXG_1")

plot.list.by.gene <- lapply(genes_to_plot, function(gene){
  max_cutoff <- max(unlist(lapply(object.list[f1.samples], function(my_list){
    my_max <- max(my_list[["seurat_obj"]][["SCT"]]@data[gene, ])
    return(my_max)
  })))
  
  plot.list.mini <- lapply(object.list[f1.samples], function(my_list){
    obj <- my_list[["seurat_obj"]]
    p <- SpatialFeaturePlot(obj,combine = F,
                                         features = gene,
                                         alpha = 1,
                            image.alpha = 0,
                            stroke = NA,
                                         pt.size.factor = point_size_df$pt.size[point_size_df$sample == obj@project.name]) 
    my_palette <- c("#330099",
                    "#00FF00",
                    "#FFFF00",
                    "#FF9933",
                    "#990000")
    p <- p[[1]]  +
      scale_fill_gradientn(
        limits = c(0, max_cutoff),
        colors = alpha(my_palette,
                       alpha = c(1, 1, 1, 1, 1)))
  return(p)
  })
  names(plot.list.mini) <- unlist(lapply(object.list[f1.samples], function(my_list){
    my_cond <- my_list[["seurat_obj"]]$cond_abv[1]
    return(my_cond)
  }))
  return(plot.list.mini)
})
names(plot.list.by.gene) <- genes_to_plot

pdf("output/cannonical_genes_spatial_plots.pdf", height = 4, width = 12)
ggpubr::ggarrange(plotlist = plot.list.by.gene[[1]],
                  align = "hv",
                  ncol = 4,
                  nrow = 1, 
                  common.legend = TRUE, 
                  legend = "top")
ggpubr::ggarrange(plotlist = plot.list.by.gene[[2]],
                  align = "hv",
                  ncol = 4,
                  nrow = 1, 
                  common.legend = TRUE, 
                  legend = "top")
ggpubr::ggarrange(plotlist = plot.list.by.gene[[3]],
                  align = "hv",
                  ncol = 4,
                  nrow = 1, 
                  common.legend = TRUE, 
                  legend = "top")
dev.off()

# Figure 1D: dell type deconvolution heatmap 
df.list <- lapply(names(sample.merged.list), function(my_cond){
  message(my_cond)
  my_samples <- Samples[grep(my_cond, Samples)]
  if(my_cond %in% c("CNTRL", "XANTH")){
    return(NULL)
  }
  
  obj_list <- object.list[my_samples]
  obj_list <- unlist(obj_list, recursive = FALSE)
  obj_list <- obj_list[grep("seurat", names(obj_list))]
  if(length(my_samples) > 1){
    obj <- merge(obj_list[[1]], obj_list[c(2:length(obj_list))])
  } else {
    obj <- obj_list[[1]]
  }
  
  obj$plot_var <- paste(obj$histological_region, obj$cond_abv)
  plot_mat <- AverageExpression(obj, 
                                assays = "dwls", 
                                group.by = c("plot_var"))
  
  plot_mat <- as.data.frame(t(plot_mat$dwls))
  df_list <- lapply(plot_mat, function(vec){vec <- scales::rescale(vec, to = c(-1, 1))})
  new_df <- as.data.frame(bind_rows(df_list))
  new_df <- new_df[ , -which(colnames(new_df) %in% c("Mast", "Melano"))] 
  new_rownames <- rownames(plot_mat)
  row_to_remove <- grep("^6. Papillary", new_rownames)
  new_df <- new_df[-row_to_remove, ]
  rownames(new_df) <- new_rownames[-row_to_remove]
  return(new_df)
})
df.list <- df.list[!unlist(lapply(df.list, is.null))]
plot_df.4 <- bind_rows(df.list)
plot_df.4 <- plot_df.4[sort(rownames(plot_df.4)), ]
p4 <- 
  pheatmap::pheatmap(t(plot_df.4),
                     cluster_cols = F, 
                     angle_col = 90,
                     scale = "none", 
                     color = colorRampPalette(c("#0000CC","#FFFFFF", "#CC0000"))(100),
                     border_color = NA, 
                     treeheight_row = 0,
                     fontsize = 12,
                     silent = T)


pdf("output/merged_deconv_heatmap_for_pub.pdf", height = 6, width = 4.5)
print(ggpubr::ggarrange(p4[[4]]))
dev.off()


# Figure 2A
# merging into one big object
merged.seu <- merge(object.list[[1]][["seurat_obj"]],
                    lapply(object.list[c(2:length(object.list))], 
                           function(obj){
                             return(obj[["seurat_obj"]])
                           }))
merged.seu$cond_abv <- factor(merged.seu$cond_abv,
                              levels = c("CNTRL", "NXG", "NL", "SAR", "GA"))
cond.color.vector <- c('#999999', '#996600', '#990000', '#330099', '#339900')
names(cond.color.vector) <- c("CNTRL", "NXG", "NL", "SAR", "GA")

# heatmapping genes
heatmap_genes <- c(
  "CD3E",
  "CD4",
  "CD8A",
  "FOXP3",
  "CTLA4",
  "PDCD1",
  "KLRG1",
  "LAG3",
  "HAVCR2",
  "TIGIT",
  "TOX",
  "CD28",
  "ICOS",
  "LEF1",
  "TCF7",
  "TBX21",
  "IFNG",
  "IL15",
  "IL12A",
  "IL12B",
  "CSF2",
  "OSM",
  "IL6",
  "CXCL9",
  "CXCL10",
  "CXCL11",
  "IL23A",
  "IL17A",
  "IL17F",
  "IL4",
  "IL5",
  "IL13",
  "MS4A1",
  "CD19",
  "CD79A",
  "CD68",
  "ITGAE",
  "CD14",
  "FCGR3A",
  "FCGR3B",
  "CCR2",
  "IL18",
  "IL1B",
  "CD80",
  "CD86",
  "NOS2",
  "IRF3",
  "IRF5",
  "STAT1",
  "NFKBIA",
  "NFKB1",
  "CD274",
  "PDCD1LG2",
  "IL18BP",
  "CX3CR1",
  "MERTK",
  "MARCO",
  "CD163",
  "ARG1",
  "STAT3",
  "STAT6",
  "PPARG",
  "IRF4",
  "IL10")

genes.to.plot <- heatmap_genes[(heatmap_genes %in% rownames(merged.seu))]

plot_mat <- AverageExpression(merged.seu, 
                              features = genes.to.plot,
                              assays = "Spatial", 
                              slot = "data",
                              group.by = c("orig.ident"))
plot_mat <- as.data.frame(plot_mat$Spatial)
pheatmap::pheatmap(t(plot_mat),
                    cluster_cols = F,
                    cluster_rows = F, 
                    angle_col = 90,
                    fontsize_col = 8,
                    scale = "column", 
                    color = colorRampPalette(c("#0000CC","#FFFFFF", "#CC0000"))(100),
                    border_color = NA, 
                    treeheight_col = 0,
                    fontsize = 12, main = "Column Scaled",
                    silent = T)[[4]]

# Figure 2B: UMAP of Mac spots only
# clustering macrophage spots alone
mac.spots <- SCTransform(mac.spots, vars.to.regress = "orig.ident")
DefaultAssay(mac.spots) <- "SCT"
mac.spots <- RunPCA(mac.spots, npcs = 30)
mac.spots <- FindNeighbors(mac.spots)
mac.spots <- RunUMAP(mac.spots, dims = 1:30)

pdf("umap.plot.mac.spots.only.pdf", height = 6, width = 6)
DimPlot(mac.spots, group.by = "cond_abv", cols = cond.color.vector)
DimPlot(mac.spots, group.by = "orig.ident", cols = c("#006600", "#339900", "#33CC00", "#00FF33",
                                                     "#990000","#993333", "#CC3333", "#CC3300","#FF3333",
                                                     "#996600", 
                                                     "#330099", "#0000FF"))
dev.off()

# Figure 2C: Volcano plot of mac to mac comparisons
Idents(mac.spots) <- "cond_abv"
mac.spots <- PrepSCTFindMarkers(mac.spots)
cond.markers <- FindAllMarkers(mac.spots)

volcano.color.vector <- c(cond.color.vector, "#CCCCCC")
names(volcano.color.vector) <- c(names(cond.color.vector), "ns")
volcano.color.vector <- volcano.color.vector[-1]
multi.way.volcanos <- 
  lapply(list(cond.markers, mac.cond.markers), 
         function(M){
           p <- 
             M %>%
             mutate(cluster = as.character(cluster),
                    xvar = ifelse(cluster %in% c("NXG", "NL"), -1 * avg_log2FC, avg_log2FC),
                    yvar = ifelse(-log(p_val_adj) == Inf, 600, -log(p_val_adj)),
                    yvar = ifelse(cluster %in% c("NXG", "SAR"), -1 * yvar, yvar),
                    color = ifelse(p_val_adj < 0.05, cluster, "ns"),
                    gene.info = abs(avg_log2FC) * -log(p_val_adj)) %>%
             group_by(cluster) %>%
             mutate(gene.rank = rank(-gene.info),
                    label = ifelse(gene.rank <= 30, gene, NA),
                    label = ifelse(gene %in% c("IRF1", "ME1", "PLD3", "CXCR4", "CCL13", "CCL17"), 
                                   gene, 
                                   label)) %>%
             filter(avg_log2FC > 0) %>%
             ggplot(aes(x = xvar, y = yvar, color = color)) +
             geom_point(size = 0.25) + 
             ggrepel::geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 100) +
             scale_color_manual(values = volcano.color.vector) +
             ylim(c(-700,700)) +
             xlab("ABS(avg_log2FC)") +
             ylab("ABS(-log(p_val_adj))") +
             scale_x_continuous(labels = abs) +
             scale_y_continuous(labels = abs) +
             geom_hline(yintercept = 0, lty = 2, color = "gray") +
             geom_vline(xintercept = 0, lty = 2, color = "gray") +
             theme_classic()
           return(p)
         })
multi.way.volcanos[[1]] <- multi.way.volcanos[[1]] + 
  ggtitle("All spots") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
multi.way.volcanos[[2]] <- multi.way.volcanos[[2]] + 
  ggtitle("1M spots only") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
pdf("four.way.volcano.plot.pdf", height = 8, width = 8)
multi.way.volcanos
dev.off()

# Figure 2D, GSEA of mac vs mac markers presented in the volcano plots above
gene_list <-
  mac_vs_mac_markers %>%
  filter(p_val_adj < 0.05,
         avg_log2FC > 0)
gene_list <- split(gene_list, f = gene_list$cluster)
gene_list <- lapply(gene_list, function(df){return(df$gene)})

enriched.list <- lapply(gene_list, function(genes){
  res <-
    enrichr(genes, dbs)
  res <-
    bind_rows(res, .id = "GSEA_type")
  return(res)
})

df <-
  enriched.list %>%
  bind_rows(.id = "sample") %>%
  mutate(term_length = sapply(strsplit(Overlap, split = "/"), getElement, 2))
# remove confounding and redundant pathways from plot
pathways_to_remove <- c("Leishmaniasis",
                        "Inflammatory Response",
                        "Rheumatoid arthritis",
                        "Type I diabetes mellitus",
                        "Staphylococcus aureus infection",
                        "Allograft rejection",
                        "Asthma",
                        "Tuberculosis",
                        "Malaria",
                        "Graft-versus-host disease",
                        "Th17 cell differentiation",
                        "Viral myocarditis",
                        "Autoimmune thyroid disease",
                        "Cell adhesion molecules",
                        "Intestinal immune network for IgA production",
                        "Inflammatory bowel disease",
                        "Hematopoietic cell lineage",
                        "UV Response Dn",
                        "Complement and coagulation cascades",
                        "Epstein-Barr virus infection",
                        "Chagas disease",
                        "Human T-cell leukemia virus 1 infection",
                        "PI3K-Akt signaling pathway",
                        "Apical Junction",
                        "Human papillomavirus infection",
                        "Myogenesis")
pathways_to_remove[!(pathways_to_remove %in% unique(unlist(lapply(enriched.list, function(df){return(df$Term)}))))]

# some data clean up
plot_df <-
  df %>%
  group_by(sample) %>%
  filter(Adjusted.P.value < 0.001) %>%
  top_n(n = 20, wt = Combined.Score) %>%
  filter(!(Term %in% pathways_to_remove))
df_to_add <- df[df$Term %in% plot_df$Term, ]
plot_df <-
  plot_df %>%
  bind_rows(.,
            df_to_add) %>%
  distinct() %>%
  group_by(Term) %>%
  mutate(max_score = max(Combined.Score),
         nhits = n()) %>%
  arrange(
    # desc(nhits),
    desc(max_score), .by_group = FALSE)
plot_df$Term <- factor(plot_df$Term, levels = rev(unique(plot_df$Term)))

pdf("output/GSEA_barchart_mac_vs_mac.pdf", height = 3, width = 6)
plot_df %>%
  ggplot(aes(x = -log(Adjusted.P.value), y = Term, fill = sample))+
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme_classic() +
  xlab("-log(p_val_adj)") +
  facet_grid(~sample) +
  scale_fill_manual(values = c("#FF9999", "#CC99CC",
                               "#3399FF","#00CC99",
                               "#FFCC33")) +
  theme(
    axis.text.y = element_text(color = "black", size = 6.5),
    axis.text.x = element_text(color = "black", size = 6),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(color = "black", size = 7),
    axis.ticks = element_blank())  +
  scale_x_continuous(expand = c(0,0)) +
  NoLegend()
plot_df %>%
  mutate(Term = factor(Term,
         levels = rev(levels(plot_df$Term)))) %>%
  ggplot(aes(y = Combined.Score, x = Term, fill = sample))+
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme_classic() +
  xlab("Combined Score") +
  scale_fill_manual(values = c("#FF9999", "#CC99CC",
                               "#3399FF","#00CC99",
                               "#FFCC33")) +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black",
                               angle = 90,
                               hjust = 1,
                               vjust = 1),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top")  +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(drop = FALSE)
dev.off()

# Figure 2D, 3D-G, 4D-G, 5D-I, 6D-G, 7A&B: Violin and Spatail Feature plots of key genes
my_genes <- c("INSERT GENE VECTOR HERE, MULTIPLE FIGURES GENERATED THIS WAY")
VlnPlot(full_merged,
        features = my_genes,
        group.by = "histological_region",
        split.by = "cond_abv",
        pt.size = 0,
        stack = T)

# or if mac.spots only
VlnPlot(mac.spots,
        features = my_genes,
        group.by = "cond_abv",
        pt.size = 0,
        cols = c("#669900", "#FFCC33","#FF9933", "#990099", "#99CCFF")
        stack = T)

# for spatial feature plots
plot.list <- lapply(object.list, function(my_list){
  obj <- my_list[["seurat_obj"]]
  
  plot.list.mini <- SpatialFeaturePlot(obj,
                                      combine = F,
                                      features = my_genes,
                                      alpha = c(0.1, 1),
                                      pt.size.factor = point_size_df$pt.size[point_size_df$sample == obj@project.name])
  names(plot.list.mini) <- my_genes[my_genes %in% rownames(obj)]
  if(!all(my_genes %in% rownames(obj))){
    missing_genes <- my_genes[!my_genes %in% rownames(obj)]
    
    missing_plots <- SpatialFeaturePlot(obj,
                                        combine = F,
                                        features = rep("FBP1",length(missing_genes)), # insert a blank gene to plot as placeholder
                                        alpha = c(0),
                                        pt.size.factor = point_size_df$pt.size[point_size_df$sample == obj@project.name])
    missing_plots <- lapply(missing_plots, function(p){p <- p + NoLegend()})
    names(missing_plots) <- missing_genes
    plot.list.mini <- c(plot.list.mini, missing_plots)
  }
  return(plot.list.mini)
})

for(my_gene in my_genes){
  print(
    ggpubr::ggarrange(plotlist = lapply(plot.list, function(mini_list){return(mini_list[[my_gene]])}),
                      align = "hv",
                      ncol = length(plot.list))
  )
}

# cell type deconvolution spatial feature plots
# deconvolution results as spatial feature plots
pdf("output/deconvolution-results-spat-feat-plots.pdf", height = 12, width = 10)
lapply(names(object.list)[-grep("CNTRL", names(object.list))],
      function(obj_name){
        seu <- object.list[[obj_name]][["seurat_obj"]]
        DefaultAssay(seu) <- "dwls"
        
        Lymphs <- c("Tcell", "NK", "Tcell-Cycling", "Bcell")
        seu$Lymphocytes <- colSums(seu[["dwls"]]@data[Lymphs, ])
        
        Fibroblasts <- c("Myofibroblast", "Fibroblast")
        seu$Fibroblasts <- colSums(seu[["dwls"]]@data[Fibroblasts, ])
        
        
        feats.to.plot <- c("Myeloid", "Fibroblasts", "Lymphocytes", "Keratinocyte")
        plot.list <- 
          SpatialFeaturePlot(seu,
                             features = feats.to.plot,
                             image.alpha = 0,
                             combine = FALSE,
                             alpha = 1,
                             stroke = NA,
                             pt.size.factor = point_size_df$pt.size[point_size_df$sample == seu@project.name])
        names(plot.list) <- feats.to.plot
        plot.list <-
          lapply(names(plot.list), function(p_name){
            p <- plot.list[[p_name]]
            p <- p  +
              scale_fill_gradientn(
                limits = c(0, 1),
                colors = alpha(my_palette,
                               alpha = c(1, 1, 1, 1, 1)
                )) +
              ggtitle(p_name) +
              theme(plot.title = element_text(size = 12,
                                              face = "bold",
                                              color = "black",
                                              hjust = 0.5,
                                              vjust = 0),
                    legend.title = element_blank(),
                    plot.margin=unit(c(0,-10,2,-10),"pt"),
                    legend.position = "top")
          return(p)
            })
        ggpubr::ggarrange(plotlist = plot.list, nrow = 2, ncol = 2, common.legend = T)
})
dev.off()