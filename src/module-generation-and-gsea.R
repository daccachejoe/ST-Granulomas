# this script was used to generate condition-specifc gene modules and then run GSEA on those gene modules
# this script spans panels C in Figures 3-6 

runMultipleSlicesGiotto <- function(condition.to.do){
seu <- sample.merged.list[[condition.to.do]]

seu$heatmap.var <- paste(seu$plot.var,
                         seu$orig.ident,
                         sep = "--")


# make a location df that doesnt create any overlapping regions and does not scale up distances between points within samples
location_df <- bind_rows(
  lapply(seu@images, 
         function(my.list){
           return(
             my.list@coordinates[,c("imagerow", "imagecol")]
           )
         })
)

location_df$sample <- as.integer(sapply(strsplit(rownames(location_df), split = "_"), getElement, 2))
location_df$imagerow <- location_df$imagerow+(location_df$sample*2000)
location_df %>% ggplot(aes(x = imagerow, y = imagecol, color = sample)) + geom_point()

# create the giotto object
test_Giotto <- createGiottoObject(raw_exprs = seu[["Spatial"]]@counts,
                                  instructions = object.list[[8]][["giotto_obj"]]@instructions,
                                  spatial_locs =location_df[, c("imagerow", "imagecol")],
                                  cell_metadata = seu@meta.data)

# normalize, find HVG, and add gene level stats to giotto objects
test_Giotto <- normalizeGiotto(test_Giotto, 
                               scale_genes = T, 
                               scale_cells = T,
                               scalefactor = 6000)
test_Giotto <- addStatistics(gobject = test_Giotto)
test_Giotto <- calculateHVG(gobject = test_Giotto,
                            show_plot = FALSE)


# create the giotto network
spatial_network_name <- "Delaunay_network"
stepsize = 30
minimum_k = 1

test_Giotto <- createSpatialGrid(gobject = test_Giotto,
                                sdimx_stepsize = stepsize,
                                sdimy_stepsize = stepsize,
                                minimum_padding = 0)
test_Giotto <- createSpatialNetwork(gobject = test_Giotto,
                                   minimum_k = minimum_k)
# plot the grid to make sure no links across samples
spatPlot(gobject = test_Giotto, show_network = T,
         network_color = 'blue', spatial_network_name = 'Delaunay_network')

# whatever this means
kmeans_spatialgenes <- binSpect(test_Giotto, bin_method = "kmeans")
ext_spatial_genes <- kmeans_spatialgenes[1:1000]$genes

spat_cor_netw_DT <- detectSpatialCorGenes(test_Giotto,
                                          method = 'network',
                                          spatial_network_name = spatial_network_name,
                                          subset_genes = ext_spatial_genes)
# spat_cor_netw_DT$cor_DT$to_keep <- rowSums(is.na(spat_cor_netw_DT$cor_DT))
# unique(spat_cor_netw_DT$cor_DT$variable[is.na(spat_cor_netw_DT$cor_DT$spat_cor)])
# spat_cor_netw_DT$cor_DT <- subset(spat_cor_netw_DT$cor_DT, to_keep == 0)
spat_cor_netw_DT <- clusterSpatialCorGenes(spat_cor_netw_DT,
                                           name = 'spat_netw_clus',
                                           k = 8)
netw_ranks <- rankSpatialCorGroups(test_Giotto, 
                                   spatCorObject = spat_cor_netw_DT,
                                   use_clus_name = "spat_netw_clus")

cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, 
                                       use_clus_name = 'spat_netw_clus',
                                       show_top_genes = 1)
cluster_genes <- cluster_genes_DT$clus
names(cluster_genes) = cluster_genes_DT$gene_ID

# create spatial metagenes and visualize
test_Giotto <- createMetagenes(test_Giotto,
                              gene_clusters = cluster_genes, 
                              name = 'cluster_metagene')

test_Giotto@gene_metadata$spat_cor_cluster <- cluster_genes[match(test_Giotto@gene_metadata$gene_ID, names(cluster_genes))]
test_Giotto@gene_metadata$spat_cor_cluster[is.na(test_Giotto@gene_metadata$spat_cor_cluster)] <- "None"


# add Giotto information to seurat object across objects
new_mat <- as.matrix(test_Giotto@spatial_enrichment$cluster_metagene[ , -"cell_ID"])
rownames(new_mat) <- test_Giotto@spatial_enrichment$cluster_metagene$cell_ID
colnames(new_mat) <- paste("Module", colnames(new_mat))
new_mat <- t(new_mat)
seu[["spat_cor_genes"]] <- CreateAssayObject(counts = new_mat)

# prep the heatmap matrix
DefaultAssay(seu) <- "spat_cor_genes"
heatmap.matrix <- AverageExpression(seu,
                                    assays = "spat_cor_genes",
                                    group.by = "heatmap.var")
heatmap.matrix <- heatmap.matrix$spat_cor_genes
heatmap.matrix <- t(scale(t(heatmap.matrix)))

# heatmap annotation information
orig.ident.color.vector <- c(RColorBrewer::brewer.pal(name = "Pastel2", n = length(unique(seu$orig.ident)))[1:length(unique(seu$orig.ident))])
names(orig.ident.color.vector) <- unique(seu$orig.ident)

col_annot_df <- sapply(strsplit(colnames(heatmap.matrix), split = "-"), getElement, 1)
col_annot <- HeatmapAnnotation(celltype = sapply(strsplit(colnames(heatmap.matrix), split = "-"), getElement, 1),
                               orig.ident = sapply(strsplit(colnames(heatmap.matrix), split = "-"), getElement, 3),
                               col = list("celltype" = color_vector,
                                          "orig.ident" = orig.ident.color.vector
                               ))

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))

# make the heatmap
pdf(paste0("output/",condition.to.do,"-heatmap-combined.objects.pdf"), height = 8, width = 5)
print(
  ComplexHeatmap::Heatmap(heatmap.matrix,
                        # split = length(unique(seu$plot.var)),
                        column_split = col_annot_df,
                        # col = col_fun,
                        # use_raster = F,
                        cluster_columns = FALSE,
                        top_annotation = col_annot,
                        show_column_names = FALSE)
)
dev.off()

plot.list <-                     
  SpatialFeaturePlot(seu, 
                    features = rownames(seu[["spat_cor_genes"]])[4], 
                    combine = FALSE,
                    image.alpha = 0,
                    stroke = NA)
names(plot.list) <- unique(seu$orig.ident)
cond.plot.list <- 
  lapply(rownames(new_mat),
         function(feat.name){
           plot.list <- 
             lapply(unique(seu$orig.ident),
                    function(p.name){
                      mini.seu <- subset(seu, orig.ident == p.name)
                      slice.to.keep <- unlist(lapply( mini.seu@images, function(coord){nrow(coord@coordinates)}))
                      slice.to.keep <- which(slice.to.keep > 0)
                      mini.seu@images <- mini.seu@images[slice.to.keep]
                      max_cutoff <- max(seu[["spat_cor_genes"]]@data[feat.name, rownames(seu@meta.data)[seu$orig.ident == p.name]])
                      p <- 
                        SpatialFeaturePlot(mini.seu, 
                                           features = feat.name, 
                                           image.alpha = 0,
                                           pt.size.factor = point_size_df$pt.size[point_size_df$sample == p.name],
                                           stroke = NA) +
                        scale_fill_gradientn(
                          limits = c(0, max_cutoff),
                          colors = alpha(my_palette,
                                         alpha = c(1, 1, 1, 1, 1)
                          ))
                      return(p)
                    })
         })


pdf(paste0("output/",condition.to.do,"-module-plots-combined.objects.pdf"), height = 8, width = 12)
lapply(cond.plot.list, function(p.list){
  print(
    ggpubr::ggarrange(plotlist = p.list,
                    nrow = 2,
                    ncol = 3)
  )
})
dev.off()

# run enrichment on those genes
library(enrichR)
gene.df.list <- 
  test_Giotto@gene_metadata %>%
  filter(spat_cor_cluster != "None") %>%
  group_by(spat_cor_cluster) %>%
  select(gene_ID) %>%
  group_split()

dbs <- c(
  # "GO_Biological_Process_2021", 
  #        "GO_Molecular_Function_2021", 
  "KEGG_2021_Human",
  "MSigDB_Hallmark_2020"
  # "TF_Perturbations_Followed_by_Expression",
  # "Jensen_DISEASES",
  # "Disease_Signatures_from_GEO_up_2014",
  # "Reactome_2016"
)
enriched.list <- lapply(gene.df.list, function(my.df){
  genes <- as.data.frame(my.df)
  genes <- genes$gene_ID
  res <- enrichr(genes, dbs)
  Sys.sleep(1)
  res <-
    bind_rows(res, .id = "GSEA_type")
  return(res)
})
names(enriched.list) <- paste("Module", c(1:8), sep = "_")
enriched.list <- enriched.list[unlist(lapply(enriched.list, nrow)) > 0]
pdf(file = paste0("output/",condition.to.do,"-combined-modules-enrichment-plots.pdf"), height = 12, width = 7)
print(
  enriched.list %>%
  bind_rows(., .id = "Module") %>%
  mutate(term_length = sapply(strsplit(Overlap, split = "/"), getElement, 2),
         term_overlap = sapply(strsplit(Overlap, split = "/"), getElement, 1),
         module = sapply(strsplit(Module, split = "_"), getElement, 2)) %>%
  filter(Adjusted.P.value < 0.05,
         term_length > 15,
         term_overlap > 1
         ) %>%
  group_by(Module) %>%
  top_n(n = 15, wt = -log(Adjusted.P.value)) %>%
  group_by(Term) %>%
  mutate(nhits = length(unique(Module)),
         group = paste(Module, collapse = ";")) %>%
  ungroup() %>%
  arrange(nhits, group, Combined.Score) %>%
  mutate(Term = factor(Term, levels = unique(Term)),
         color.var = ifelse(-log(Adjusted.P.value) >= 50, 50, -log(Adjusted.P.value))) %>%
  ggplot(aes(x = module, 
             y = Term,
             color = color.var, 
             size = log(Combined.Score))) +
  geom_point(alpha = 1) +
  scale_color_gradient(low = "#99CCFF", high = "#FF0000") +
  scale_size_continuous() +
  # ggtitle(toupper(cond)) +♣
  theme_classic() +
  BoldTitle() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        # axis.text.x = element_text(angle = 45, hjust =1, vjust = 1, color = "black"),
        axis.title.y = element_blank(),
        axis.ticks = element_blank())
)
dev.off()


return(list(test_Giotto, seu, enriched.list))
}



# output.list <- lapply(c("GA", "NL", "SAR"), function(cond){
output.list <- lapply(c("SAR"), function(cond){
  my.list <- runMultipleSlicesGiotto(cond)
  return(my.list)
})
