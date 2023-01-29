# Conduct gene set over-representation analysis for gene clusters from time course heatmap
# This script is inspired by a script by Dr. Dennis Goldfarb from Washington University (https://github.com/GoldfarbLab/H522_paper_figures/blob/master/Figure5.R).



# ---
# title: "Enrichment Heatmap"
# author: "Yiqing Wang"
# date: "1/19/2023"
# input: filtered CPM table and gene cluster table from time course heatmap analysis
# output: heatmap demonstrating gene sets over-represented in each cluster, in PNG format
# ---



library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(tidyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)



### SETUP ###

# The base path should contain an enrichment_heatmap/ folder and a time_course_heatmap/ folder, which was created in the time course heatmap script.

basepath <- "path/of/choice"

path <- file.path(basepath, "er_heatmap")
tc_path <- file.path(basepath, "tc_heatmap", "outputs")

dir.create(path = file.path(path, "enrichment_results"))
dir.create(path = file.path(path, "heatmap"))
dir.create(path = file.path(path, "selected_term"))


prefix <- "prefix_of_choice"



### DATA PREPARATION ###

# filtered cpm table from time course cpm heatmap for gene universe
cpm <- read.delim(file.path(tc_path, paste0(prefix, "_cpm")),
                  row.names = 1)

outputname <- prefix


gene_cluster <- read.delim(file.path(tc_path, paste0(prefix, "_clusterInfo")))
colnames(gene_cluster) <- c("gene_name", "cluster")

# total number of clusters
num_cluster <- max(gene_cluster$cluster)

# clusters selected to be plotted in the enrichment heatmap, will modify the variable later
selected_clusters <- seq(num_cluster)



geneNames.universe <- row.names(cpm)


stopifnot(gene_cluster$gene_name %in% geneNames.universe)




## Retrieving info from Molecular Signatures database ##

msig_h_bp_c6 <- msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "H" | gs_cat == "C6" | (gs_cat == "C5" & gs_subcat == "GO:BP")) %>%
  select(gs_name, human_gene_symbol)

save(msig_h_bp_c6, file = file.path(path, "enrichment_results", "msig_h_bp_c6.RData"))




### OVER-REPRESENTATION ANALYSIS ###

enrich <- function(cluster_table, universe, database, selected_clusters) {
  enrichment.results <- NULL

  ## Over-representation analysis for each cluster, filtering q-value <= 0.1 ##
  for (cluster.id in selected_clusters) {

    geneNames.cluster <- cluster_table[cluster_table$cluster == cluster.id, ]$gene_name

    em <- enricher(gene = geneNames.cluster,
                   pAdjustMethod = "fdr",
                   universe = universe,
                   TERM2GENE = database,
                   minGSSize = 10,
                   qvalueCutoff = 0.1)


    if (!is.null(em)) {
      result <- em@result
      result$cluster <- cluster.id
      result <- filter(result, qvalue <= 0.1)
      if (nrow(result) > 0) {
        enrichment.results <- rbind(enrichment.results, result)
      }
    }
  }

  return(enrichment.results)
}


enrichment.results <- enrich(cluster_table = gene_cluster,
                             universe = geneNames.universe,
                             database = msig_h_bp_c6,
                             selected_clusters = selected_clusters)



enrichment.results$cluster <- as.integer(enrichment.results$cluster)

# clusters that contain terms with q-value <= 0.1
enriched_cluster <- unique(enrichment.results$cluster)

enrichment.results <- rename(enrichment.results, genes = geneID)


## Getting number of sig terms enriched for each cluster
sig.term.in.cluster <- data.frame(table(enrichment.results$cluster))
colnames(sig.term.in.cluster) <- c("cluster", "n_sig_terms")
write.csv(sig.term.in.cluster,
          file.path(path, "enrichment_results", paste0(outputname, "_sig_term_in_cluster.csv")),
          row.names = F)


enrichment.results$category <- ""
enrichment.results$category[str_detect(enrichment.results$Description, pattern = "^HALLMARK")] <- "Hallmark"
enrichment.results$category[str_detect(enrichment.results$Description, pattern = "^GOBP")] <- "Biological Process"
enrichment.results$category[!str_detect(enrichment.results$Description, pattern = "(^HALLMARK)|(^GOBP)")] <- "Oncogenic"


# function to enrich individual cluster with no q-value cutoff (for diagnostic purposes)
enrich_cluster <- function(cluster_table, universe, database, cluster_no) {
  result <- NULL
  geneNames.cluster <- cluster_table[cluster_table$cluster == cluster_no, ]$gene_name
  em <- enricher(gene = geneNames.cluster,
                 pAdjustMethod = "fdr",
                 universe = universe,
                 TERM2GENE = database,
                 minGSSize = 10,
                 qvalueCutoff = 0.1)
  if (!is.null(em)) {
    result <- em@result
    result$cluster <- cluster_no
  }
  return(result)
}

# as an example, enrich cluster 7 with enrich_cluster function
c7_msig <- enrich_cluster(cluster_table = gene_cluster, universe = geneNames.universe,
                          database = msig_h_bp_c6, cluster_no = 7)




## Getting number of genes present in each cluster ##
num.in.cluster <- data.frame(table(gene_cluster$cluster))
colnames(num.in.cluster) <- c("cluster", "n")
write.csv(num.in.cluster,
          file.path(path, "enrichment_results", paste0(outputname, "_num_in_cluster.csv")),
          row.names = F)


## Getting the ratio of genes in a cluster that belong to a term against the total number of genes in the cluster ##
enrichment.results$ratio <- enrichment.results$Count / num.in.cluster[match(enrichment.results$cluster, num.in.cluster$cluster), ]$n
write.csv(enrichment.results,
          file.path(path, "enrichment_results", paste0(outputname, "_results.csv")),
          row.names = F)


## "Flattening" the results by showing results for each cluster in parallel, merging terms shared by different clusters ##
enrichment.results_wide <- pivot_wider(enrichment.results,
                                       id_cols = c(Description, category),
                                       names_from = cluster,
                                       values_from = c(pvalue, p.adjust, qvalue, Count,
                                                       GeneRatio, BgRatio, ratio, genes))
enrichment.results_wide <- arrange(enrichment.results_wide, category)
write.csv(enrichment.results_wide,
          file.path(path, "enrichment_results", paste0(outputname, "_results_wide.csv")),
          row.names = F)
save(enrichment.results_wide,
     file = file.path(path, "enrichment_results", paste0(outputname, "_results_wide.RData")))




## Manual selection of interesting terms ##

# The terms shown on the eventual heatmap are manually selected from the terms in enrichment_results_wide.
# The indices of the selected terms should be specified in a "selected_term_index.txt" file and saved in the "selected_term" folder.
# Each index should take up one line in the file.


selected_term_index <- read.table(file = file.path(path, "selected_term", "selected_term_index.txt"),
                                  header = F)[, 1]

enrichment.results.pruned <- enrichment.results_wide[selected_term_index, ]

enrichment.results.pruned <- arrange(enrichment.results.pruned, category)


# clusters that have selected terms
# Note that it is possible that not all clusters with significantly enriched terms have a term that is selected to be in the final heatmap.
# In other words, enriched_cluster and selected_clusters may not be the same.
selected_clusters <- c(1, 3, 4, 5, 8, 9, 10, 11)


## Changing the style of the Description column ##
# convert to lower case; replace "_" with " ".

selected_terms <- tolower(enrichment.results.pruned$Description) %>%
  gsub(pattern = "_", replacement = " ")

write.table(data.frame(selected_terms),
            file = file.path(path, "selected_term", paste0(outputname, "_selected_terms.txt")),
            col.names = F, row.names = F, quote = F)

# Here, some manual style editing of the selected terms is required.
# For example, change the first letter of each word to upper case, and change all abbreviations to upper case.
# Save the edited selected term file to "[outputname]_selected_terms_edited.txt".

## After manual editing of descriptions ##
selected_terms_edited <- read.table(file = file.path(path, "selected_term", paste0(outputname, "_selected_terms_edited.txt")),
                                    header = F,
                                    sep = "\t")[, 1]

enrichment.results.pruned$Description <- selected_terms_edited



save(enrichment.results.pruned,
     file = file.path(path, "enrichment_results", paste0(outputname, "_pruned.RData")))
load(file = file.path(path, "enrichment_results", paste0(outputname, "_pruned.RData")))

write.csv(enrichment.results.pruned,
          file = file.path(path, "enrichment_results", paste0(outputname, "_pruned.csv")),
          row.names = F)



## Selecting q-values for plotting ##
em.q.values <- select(enrichment.results.pruned, grep(pattern = "qvalue",
                                                      colnames(enrichment.results.pruned),
                                                      value = T))

# getting rid of q-value columns whose cluster does not contain selected terms (columns that are not all NA)
em.q.values <- em.q.values[, colSums(is.na(em.q.values)) < nrow(em.q.values)]

em.q.values <- as.data.frame(em.q.values)
row.names(em.q.values) <- as.character(enrichment.results.pruned$Description)
colnames(em.q.values) <- selected_clusters


# taking log10 of the q-values
em.q.values <- -log10(as.matrix(em.q.values))


write.csv(em.q.values,
          file = file.path(path, "enrichment_results", paste0(outputname, "_qvalues.csv")))

save(em.q.values,
     file = file.path(path, "enrichment_results", paste0(outputname, "_qvalues.RData")))
load(file = file.path(path, "enrichment_results", paste0(outputname, "_qvalues.RData")))




## Selecting ratio values for plotting ##
ratio.values <- select(enrichment.results.pruned, grep(pattern = "ratio",
                                                       colnames(enrichment.results.pruned),
                                                       value = T))

ratio.values <- ratio.values[, colSums(is.na(ratio.values)) < nrow(ratio.values)]

ratio.values <- as.data.frame(ratio.values)
row.names(ratio.values) <- as.character(enrichment.results.pruned$Description)
colnames(ratio.values) <- selected_clusters


write.csv(ratio.values,
          file = file.path(path, "enrichment_results", paste0(outputname, "_ratio.csv")))
save(ratio.values,
     file = file.path(path, "enrichment_results", paste0(outputname, "_ratio.RData")))
load(file = file.path(path, "enrichment_results", paste0(outputname, "_ratio.RData")))






### ENRICHMENT HEATMAP ###

## col_fun function for heatmap to plot different shades of colors of dots according to q-values ##
max.q <- max(em.q.values, na.rm = T)
min.q <- min(em.q.values, na.rm = T)

# below is to fix the problem of legend not plotting the whole range of values from min to max
col_fun <- colorRamp2(breaks = c((min.q %/% 20) * 20, (max.q %/% 20 + 1) * 20), colors = c("black", "red"))


## Column annotation blocks with cluster information ##
colPalette <- c("#AEC7E87F", "#98DF8A7F", "#FF98967F",
                "#C49C947F", "#C5B0D57F", "#FFBB787F", "#D627287F",
                "#FF7F0E7F", "#9467BD7F", "#8C564B7F", "#E377C27F",
                "#F7B6D27F")


col.annot <- HeatmapAnnotation(cluster = anno_simple(selected_clusters,
                                                     height = unit(3, "mm"),
                                                     col = structure(colPalette[selected_clusters],
                                                                     names = selected_clusters)
                                                     ),
                               show_annotation_name = F,
                               show_legend = F)



## Heatmap ##

# Graphic parameters can be customized for the atheistics of the plot.

h <- Heatmap(em.q.values,

             # labels
             column_title = "Enriched Gene Sets",
             column_title_gp = gpar(fontsize = 13),
             column_names_gp = gpar(fontsize = 11),
             column_names_rot = 0,
             row_title_gp = gpar(fontsize = 11),
             row_names_gp = gpar(fontsize = 10),

             # annotation
             bottom_annotation = col.annot,

             # legends (custom legends are introduced later)
             show_heatmap_legend = F,

             show_row_names = T,
             show_column_names = T,

             # clustering
             cluster_rows = F,
             cluster_columns = F,

             split = enrichment.results.pruned$category,

             # circles
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height,
                         gp = gpar(fill = "white", col = "#EEEEEE")) # white rectangles in the background
               grid.circle(x = x, y = y,
                           r = sqrt(ratio.values[i, j]) * 4, # <<< change size of circles according to ratio distribution, so that all circles fit in the grids
                           default.units = "mm",

                           # coloring the circles according to the q-values
                           gp = gpar(fill = col_fun(em.q.values[i, j]), col = NA))
             }
)



## Q-value legend ##

q_legend <- Legend(col_fun = col_fun,
                   title = "-log10(q-value)",
                   title_gp = gpar(fontsize = 11),
                   labels_gp = gpar(fontsize = 10),
                   grid_height = unit(4, "mm"),
                   direction = "horizontal")



## Ratio value legend ##

# ratio legend breaks: customize according to the actual distribution of ratio values
ratio_breaks <- c(0.5, 0.4, 0.3, 0.2, 0.1)

ratio_legend <- Legend(labels = ratio_breaks,
                       grid_height = unit(5.5, units = "mm"),
                       grid_width = unit(5.5, units = "mm"),
                       title = "Percent\nof Cluster",
                       title_gp = gpar(fontsize = 11),
                       labels_gp = gpar(fontsize = 10),
                       type = "points",
                       pch = 1, # circle

                       # circle size of ratio legend breaks
                       # IMPORTANT: change the multiplier number below so that the legend circle sizes match the sizes of circles in the heatmap grids.
                       # 2.5 times seems to be good.
                       size = unit(sqrt(ratio_breaks) * 10, units = "mm"),

                       background = "white"
)



## Graphic Devices Settings ##
# Can customize the numbers for aesthetics

png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1900, res = 300)
draw(h, padding = unit(c(2, 12, 8, 25), "mm")) # padding: empty space around heatmap, for legends, row names...
draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
draw(ratio_legend, x = unit(0.91, "npc"), y = unit(0.88, "npc"))
dev.off()
