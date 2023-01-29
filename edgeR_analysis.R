# Differential Gene Expression (DGE) Analysis using edgeR
# The DGE analysis result table produced is used in subsequent time-course heatmap and TE analyses.



# ---
# title: "DGE Analysis"
# author: "Kasyap Tenneti, NaKyung Lee, Yiqing Wang"
# date: "6/23/2022"
# input: count files of two conditions being compared
# output: MDS plot, BCV plot, DGE analysis result table, volcano plot, CPM heatmap, collapsed heatmap, logFC heatmap
# ---



library(edgeR)
library(RColorBrewer)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)



## Data Preparation ##

directory <- "path/to/count/files"

# count file names, separated by comma
# note that the count files need to be ordered by replicates, then by conditions, with the reference condition first
files <- "cond1_rep1,cond2_rep1,cond1_rep2,cond2_rep2,cond1_rep3,cond2_rep3"
files <- unlist(strsplit(files, split = ","))

# replicate names of the count files
donor_list <- "rep1,rep1,rep2,rep2,rep3,rep3"
donor_list <- unlist(strsplit(donor_list, split = ","))

# conditions of the count files
conditions <- "cond1,cond2,cond1,cond2,cond1,cond2"
conditions <- unlist(strsplit(conditions, split = ","))

prefix <- "prefix_of_choice"

# path to a gene-only GTF whose order of genes is the same as in count files
gtf_path <- "path/to/gtf"

output <- "output/path"

dir.create(output)


print(files)
print(donor_list)
print(conditions)


setwd(directory)

sample_design <- data.frame(files, donor_list, conditions)
write.table(sample_design, file = file.path(output, paste0(prefix, "_sample_design")),
			col.names = NA, quote = F, sep = "\t")

labels <- paste(conditions, donor_list, sep = "_") # <<< may change the labels here
print(labels)
DG <- readDGE(files, header = FALSE, labels = labels)

write.table(DG$counts, file = file.path(output, paste0(prefix, "_original_counts")),
			col.names = NA, quote = F, sep = "\t")



## Gene ID Conversion ##

# Gene ID conversion using gtf used in mapping. The total number and order of genes in the gtf should be exactly the same as in the individual count files.
# Thus, can probably directly replace gene ids in DG with list of gene names in the gtf.
if (grepl(pattern = "ENSG", row.names(DG)[1])) {
	gtf <- read.table(gtf_path, sep = " ")

	gene_id <- as.character(gtf[, 2])
	gene_id <- substr(gene_id, start = 1, stop = nchar(gene_id) - 1)
	stopifnot(row.names(DG) == gene_id)

	gene_names <- as.character(gtf[, 6])
	gene_names <- substr(gene_names, start = 1, stop = nchar(gene_names) - 1) # getting rid of the ";" at the end
	row.names(DG) <- gene_names
}



## Filtering ##

keep <- rowSums(cpm(DG) > 1) >= length(unique(donor_list)) # filtering by expression level
print(paste0("before filtering: ", nrow(DG)), quote = F)
DG <- DG[keep, ,keep.lib.sizes = FALSE]
print(paste0("after filtering: ", nrow(DG)), quote = F)


## Deduplicating ##

duplicate_genes <- row.names(DG)[duplicated(row.names(DG))]
DG <- DG[!(row.names(DG) %in% duplicate_genes), ,keep.lib.sizes = F]
print(paste0(length(duplicate_genes) + length(unique(duplicate_genes)),
			 " duplicated genes removed."), quote = F)



## Normalization ##

DG <- calcNormFactors(DG)
DGgroups <- conditions #conditions is a comma separated list of conditions i.e. c("IFN+", "IFN+", "IFN-", "IFN-")
DG <- DGEList(counts = DG, group = factor(DGgroups))
DG <- calcNormFactors(DG)
Donor <- factor(donor_list) #donor_list is list of the donor of each file;
Condition <- factor(conditions) #conditions is same as conditions for DGgroups;


## MDS Plot ##

if(length(labels) > 9) {
  colors <- c(brewer.pal(12,"Paired"), brewer.pal(8,"Dark2"))[seq_along(labels)] # max 20 colors
} else if(length(labels) > 2) {
	colors <- brewer.pal(length(labels),"Set1")
} else {
  colors <- c("deepskyblue", "red")
}

names(colors) <- labels


mdsfilename <- paste0(prefix, "_MDSplot.pdf")
pdf(file.path(output, mdsfilename))

par(mar = c(5, 4, 4, 6)) # leave more margin to the right for the legend
plotMDS(DG, col = colors, pch = 18, cex = 2)
legend(x = 'topright', title="Experiment", legend=labels,
	   col=colors, cex=.75, pch=18, xpd = T, bty="n", inset=c(-0.16, 0))

dev.off()


# MDS dim2 vs dim3 #

mdsfilename <- paste0(prefix, "_MDSplot_2-3.pdf")
pdf(file.path(output, mdsfilename))

par(mar = c(5, 4, 4, 6))
plotMDS(DG, dim.plot = c(2, 3), col = colors, pch = 18, cex = 2)
legend(x = 'topright', title="Experiment", legend=labels,
	   col=colors, cex=.75, pch=18, xpd = T, bty="n", inset=c(-0.16, 0))

dev.off()


## Dispersion Estimation, BCV Plot, and Exact Test ##

design <- model.matrix(~ Donor + Condition)
DG <- estimateDisp(DG, design, robust = TRUE)
print(paste0("common dispersion: ", DG$common.dispersion), quote = F)

dis <- data.frame(DG$AveLogCPM, DG$tagwise.dispersion, DG$trended.dispersion)
row.names(dis) <- row.names(DG$counts)
write.table(dis, file = file.path(output, paste0(prefix, "_dispersion")),
			col.names = NA, sep = "\t", quote = F)


bcv_file_name <- paste0(prefix, "_BCVplot.pdf")
pdf(file.path(output, bcv_file_name)) #name of file to write bcv graph to;
plotBCV(DG)
dev.off()


de <- exactTest(DG, pair = unique(conditions))
de.genes <- row.names(topTags(de, n = nrow(de$table))$table)
diffex_file <- paste0(prefix, "_diffexgenes")

# DGE analysis result table
diffex_file_full <- paste0(prefix, "_diffexgenes_full")
write.table(de.genes, file = file.path(output, diffex_file),
			quote = FALSE, row.names = FALSE) #file to write de.genes to
write.table(topTags(de, n = nrow(de$table))$table, file = file.path(output, diffex_file_full),
			quote = FALSE, col.names = NA, sep = '\t')


## Smear Plot ##

smearplotname <- paste0(prefix, "_SmearPlot.pdf")
pdf(file.path(output, smearplotname))
plotSmear(DG, de.tags = de.genes) #Generates a smear plot for differentially expressed genes, upregulated and downregulated;
dev.off()

print('>>>>>>>>>>>>')


## Volcano Plot ##

detags <- topTags(de, n = n) #Choose n to be how many genes you want to show on the graph (not working currently)

volcanodata <- cbind(detags$table$logFC, -log10(detags$table$PValue), detags$table$FDR)
row.names(volcanodata) <- row.names(detags$table)
colnames(volcanodata) <- c("logFC", "PValue", "FDR")
write.table(volcanodata, file = file.path(output, paste0(prefix, "_volcanodata")),
			quote = F, col.names = NA, sep = "\t")


volcanodata_insig <- volcanodata[volcanodata[,2] < -log10(.05),]
volcanodata_sig <- volcanodata[volcanodata[,2] >= -log10(.05),]
volcano_logFC <- volcanodata_sig[abs(volcanodata_sig[,1]) >= 2,]
volcanoplotname <- paste0(prefix, "_VolcanoPlot.png")


windowsize <- 15
print('>>>>>>>>>>>>')
print(max(volcanodata[is.finite(volcanodata[ ,2]), 2]))
upperlim <- max(volcanodata[is.finite(volcanodata[ ,2]), 2]) + 1

png(file.path(output, volcanoplotname), height = 1000, width = 1000)
par(mar = c(5, 5, 5, 5))
plot(volcanodata_insig, col = 'darkgray', pch = 19, cex = 1.5,
	 ylim = c(0,upperlim), xlim = c(-windowsize, windowsize),
	 cex.lab = 2, cex.axis = 1.5, font.lab = 2, font = 2)
points(volcanodata_sig, col = 'firebrick', cex = 1.5, pch = 19)
points(volcano_logFC, col = 'dodgerblue', cex = 1.5, pch = 19)

sig_FDR <- volcanodata[volcanodata[,3] < .05 & volcanodata[,2] < upperlim, ]
sig_FDR <- sig_FDR[order(abs(sig_FDR[,1]), decreasing = T), ]
sig_FDR <- sig_FDR[sig_FDR[,1] < windowsize, ]
outlier <- volcanodata[volcanodata[,2] > 10, ] # <<< Adjust this number depending on the number of DE genes with higher -log10(P Value)

print(dim(outlier))

num_gene <- min(15, length(rownames(sig_FDR)))
num_out <- min(20, length(rownames(outlier)))

text(sig_FDR[1:num_gene,1], sig_FDR[1:num_gene,2] ,
	 labels=row.names(sig_FDR[1:num_gene,]),
	 cex=1.2, pos=2, font=2)
text(outlier[1:num_out,1], outlier[1:num_out,2] ,
	 labels=row.names(outlier[1:num_out,]),
	 cex=1.2, pos=2, font=2)
abline(h=-log10(.05), col="blue", lty=2, lwd=3)

dev.off()


## CPM Heatmap ##

y <- cpm(DG, log = TRUE, prior.count = 1)
cpmname <- paste0(prefix, "_cpm")
write.table(y, file.path(output, cpmname), col.names = NA, quote = F, sep = "\t")

volcdata_sub <- volcanodata[rownames(volcanodata) %in% rownames(y),]

# filter genes with FDR <0.05, and order according to logFC
max_rows <- rownames(volcdata_sub[order(abs(volcdata_sub[volcdata_sub[,3] < .05, ][,1]), decreasing=TRUE), ])
head(max_rows)

selY <- y[max_rows, ]
selY <- selY[1:min(30, length(rownames(selY))),] # selects top 30 genes with the most logFC

htmaptablename <- paste0(prefix, "_heatmap_table")
write.table(selY, file.path(output, htmaptablename), col.names = NA, quote=FALSE, sep = "\t")

heatmapname <- paste0(prefix, "_HeatMap.png")
png(file.path(output, heatmapname), height=650, width=500)
Heatmap(selY,
		column_order=colnames(selY),
		heatmap_height=unit(20,'cm'),
		heatmap_width=unit(12,'cm'),
		heatmap_legend_param = list(title = "log2(CPM)"))
dev.off()


## Collapsed Heatmap ##

cols <- unique(conditions)
collapsed <- matrix(nrow = nrow(selY), ncol = length(cols))
if("mock" %in% cols) {
	cols <- c("mock", cols[cols != "mock"])
}
colnames(collapsed) <- cols
cols <- colnames(collapsed)
rownames(collapsed) <- rownames(selY)

for(i in seq(1, length(cols))) {
	collapsed[,i] <- rowMeans(selY[,grepl(cols[i], colnames(selY))])
}

collapsename <- paste0(prefix, "_collapsed_heatmap.png")
png(file.path(output, collapsename), height=650, width=350)
Heatmap(collapsed,
		column_order=colnames(collapsed),
		heatmap_height=unit(20,'cm'),
		heatmap_width=unit(8,'cm'),
		heatmap_legend_param = list(title = "log2(CPM)"))
dev.off()


## LogFC Heatmap ##

donors <- unique(donor_list)
logFC_df <- matrix(nrow = nrow(selY), ncol = length(donors), dimnames = list(rownames(selY), donors))
for (i in seq(length(donors))) {
	logFC_df[, i] <- selY[, 2 * i] - selY[, 2 * i - 1] # logFC calculation depends the order of variables (columns) in selY, which depends on filelist order.
}

limit <- max(abs(logFC_df))
col_fun <- colorRamp2(breaks = c(-limit, 0, limit), colors = c("blue", "white", "red")) # to make sure logFC=0 is white

logFC_name <- paste0(prefix, "_logFC_HeatMap.png")
png(file.path(output, logFC_name), height = 650, width = 350)
Heatmap(logFC_df,
		col = col_fun,
		column_order = colnames(logFC_df),
		heatmap_height = unit(20,'cm'),
		heatmap_width = unit(8,'cm'),
		heatmap_legend_param = list(title = "log2(FC)"))
dev.off()
