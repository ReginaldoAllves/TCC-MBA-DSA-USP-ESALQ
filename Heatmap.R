rm(list=ls())

## Libraries

library(tidyverse)
library("ComplexHeatmap")

## Reading files

ucec.table = read.table("./data/utero.data.txt",
                        header = T, sep = "\t", quote = "", fill = T)

## Selecting only molecular data

ucec.molecular = dplyr::select(ucec.table, all_of(c(1, 12:27)))

## Transforming data

ucec.molecular$msi = ifelse(ucec.molecular$msi == 0, "mss", "msi")

## Heatmap substitutions

ucec.sub = data.frame(dplyr::select(ucec.molecular, 1,11:17), row.names = 1)
matrix.sub = as.matrix(ucec.sub[-7])

matrix.sub = t(matrix.sub / 100)


# Annotation

col = list(msi.status = c("mss" = "green", "msi" = "black"))
annotation <- HeatmapAnnotation(
              msi.status = ucec.sub$msi,
              col = col
)

Heatmap(matrix.sub, 
        name = "Contribution (%)", #title of legend
        row_names_gp = gpar(fontsize = 12), # Text size for row names
        cluster_rows = F,
        show_column_names = F,
        bottom_annotation = annotation,
        column_dend_height = unit(4, "cm"),
        clustering_method_columns = "complete",
        column_km = 3, border = TRUE
)