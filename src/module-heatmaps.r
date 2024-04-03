# Script used for generating module heatmaps per condition (panel A in figures 3-6)
# script for generating module heatmaps per condition
library(ComplexHeatmap)
conditions <- c("GA", "NL", "SAR", "NXG")

heatmap.list <-
  lapply(conditions,
         function(my_cond){
           message(my_cond)
           objects.to.plot <-
             names(object.list)[grep(my_cond, names(object.list))]
           lapply(object.list[objects.to.plot], 
                  function(obj.list){
                    seu <- obj.list[["seurat_obj"]]
                    plot.df <- seu[["spat_cor_genes"]]@data
                    seu <- ScaleData(seu, assay = "spat_cor_genes")
                    
                    # plot.df <- t(scale(t(as.matrix(plot.df))))
                    plot.df <- seu[["spat_cor_genes"]]@scale.data
                    rownames(plot.df) <- paste("Module", rownames(plot.df))
                    barcodes <- rownames(seu@meta.data)[order(seu$plot.var)]
                    
                    col_annot_df <- seu@meta.data[barcodes, "plot.var"]
                    
                    col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
                    col_fun(seq(-2, 2))
                    
                    ht <-Heatmap(plot.df[,barcodes],
                            cluster_columns = F,
                            col = col_fun,
                            column_split = col_annot_df,
                            show_column_names = F,
                            show_row_names = T)
                    draw(ht, 
                         column_title = (seu@project.name),
                         column_title_gp = gpar(fontsize = 20, fontface = "bold"))
                  })
         })

pdf("output/module-heatmaps-complex-.pdf", height = 6, width = 5)
heatmap.list
dev.off()





