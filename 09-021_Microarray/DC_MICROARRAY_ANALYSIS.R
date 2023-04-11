library(oligo)
library(Biobase)
library(oligoClasses)
library(ArrayExpress)
library(pd.hugene.2.0.st)
#library(arrayQualityMetrics)
library(ggplot2)
library(limma)
library(stringr)
library(dplyr)
library(pheatmap)
library(colorRamps)
library(colorspace)

normalized_module <- read.table("DCExpression_Normalized_Matrix.txt", row.names = 1, header=TRUE)
dim(normalized_module) # 34659    99
# this is just an additinal section in case the names in expression set don't match the names in pheno set
################ I took out the column names to make sure they match the names in pheno data table
subset <- normalized_module[1,]
subset
write.table(subset, file="subset.txt", sep="\t", quote=F, col.names=NA,)
# after writing out the table I combined it with pheno table names
head(normalized_module)
################ 

pheno_data_set <- read.table("DCExpression_PhenoData.txt", row.names = 1, header=TRUE)
summary(pheno_data_set) 
dim(pheno_data_set) #99 59

# convert pheno_data_set_1 to Annotated data frame format
pheno_data_set_1 <- AnnotatedDataFrame(pheno_data_set)
pheno_data_set_1

# convert normalized_module "data frame" in matrix_data "matrix"

matrix_data <- as.matrix(normalized_module)
colnames(matrix_data)
## create ExpressionSet object:
eset_normalied_pat <- new("ExpressionSet", exprs = matrix_data, phenoData = pheno_data_set_1)

# check some characteristics to make sure things match 
pData(eset_normalied_pat)
fData(eset_normalied_pat)
eset_normalied_pat$Cell_Type

# we can see the boxplot of normalized data:
oligo::boxplot(eset_normalied_pat, target = "core",
               main = "Boxplot of log2-intensitites for the raw data", cex.axis=0.5, las=2)

########## PCA anaysis of normalized data 
exp_raw_J <- log2(Biobase::exprs(eset_normalied_pat))

PCA_raw_J <- prcomp(t(exp_raw_J), scale. = FALSE)
percentVar_n <- round(100*PCA_raw_J$sdev^2/sum(PCA_raw_J$sdev^2),1)
sd_ratio_n <- sqrt(percentVar_n[2] / percentVar_n[1])

pheno_data_set

dataGG_n <- data.frame(PC1 = PCA_raw_J$x[,1], PC2 = PCA_raw_J$x[,2],PC3 = PCA_raw_J$x[,3],
                       Celltype = pData(pheno_data_set_1)$Cell_Type,
                       Outcome.1 = pData(pheno_data_set_1)$Outcome_1,
                       Outcome = pData(pheno_data_set_1)$Outcome_2,
                       Clinical_Response = pData(pheno_data_set_1)$Response)
dataGG_n                       
ggplot(dataGG_n, aes(PC1, PC2)) +
  geom_point(aes(shape = Celltype, colour = Clinical_Response), size = 2.5, alpha=.5) + 
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar_n[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar_n[2], "%")) + ylim(-20,20)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  coord_fixed(ratio = sd_ratio_n) +
  scale_shape_manual(values = c(16,15,17)) +
  scale_color_manual(values = c("slateblue4", "darkgreen","red","purple","green"))+ theme_bw()


###########################################################################################
###########################################################################################
###########################################################################################

library(ggrepel)
library("scatterplot3d")  
library(FactoClass)
require(RCurl)
library(Rgb)
######## DEFINE COLORS 
mycol <- rgb(0,100,0, max = 255, alpha = 110, names = "darkgreen")
mycol1 <- rgb(72,61,139, max = 255, alpha = 150, names = "blue50")
mycol2 <- rgb(255,0,0, max = 255, alpha = 130, names = "red")
colors
colors <- c(mycol1,mycol,mycol2)
colors <- colors[as.numeric(dataGG_n$Celltype)]
shapes = c(16,15, 17)
shapes <- shapes[as.numeric(dataGG_n$Celltype)]
colors


library(pca3d)
colnames(dataGG_n)
pca <- prcomp(t(exp_raw_J), scale.= FALSE) 
gr <- factor(dataGG_n[,4])
summary(gr)

pca3d( pca, group=gr, col=colors,
       bg= "white", show.group.labels = F, show.plane = F,
       axes.color= "gray",legend="bottomleft")
pca3d( pca, group=gr, col=colors,
       bg= "white", show.group.labels = F, show.plane = F,
       axes.color= "gray", show.centroids =TRUE)
snapshotPCA3d(file="Lymph-3D_PCA_1.png") # save picture
pca3d( pca, group=gr, col=colors,
       bg= "white",  show.group.labels = F, show.plane = F,
       axes.color= "gray50", show.centroids =TRUE, show.ellipses=TRUE, ellipse.ci=0.75)
snapshotPCA3d("testfile.png")

########################### OUTCOME vs ###########################
########################### OUTCOME vs ###########################


#### Heatmap clustering analysis #### Heatmap clustering analysis #### Heatmap clustering analysis
#### Heatmap clustering analysis #### Heatmap clustering analysis #### Heatmap clustering analysis

#### Heatmap clustering analysis #### Heatmap clustering analysis #### Heatmap clustering analysis
#### Heatmap clustering analysis #### Heatmap clustering analysis #### Heatmap clustering analysis

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
dists <- as.matrix(dist(t(exp_raw_J), method = "manhattan"))
colnames(dists) <- NULL
diag(dists) <- NA
rownames(dists) <-  pData(pheno_data_set_1)$Cell_Type
hmcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(255)
########## cool color 
hmcol <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
pheatmap(dists, col = rev(hmcol), clustering_distance_rows = "manhattan",clustering_distance_cols = "manhattan",fontsize_row = 7)
                                                          

####################################################################################
library(stringr)



$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
str(pData(eset_normalied_pat))


Cell_type <-as.character(pheno_data_set_1$Cell_Type)
Cell_type

Outcome <-as.character(pheno_data_set_1$Outcome_1)
Outcome

Outcome_2 <-as.character(pheno_data_set_1$Outcome_2)
Outcome_2

CD8.Resp <-as.character(pheno_data_set_1$CD8_Resp_Per_Patient)
CD8.Resp

CD4.Resp <-as.character(pheno_data_set_1$CD4_Resp_Per_Patient)
CD4.Resp


annotation_for_heatmap <-
  data.frame(Cell.type = Cell_type, Outcome = Outcome, CD8=CD8.Resp,CD4=CD4.Resp)

row.names(annotation_for_heatmap) <- row.names(pData(eset_normalied_pat))

dists <- as.matrix(dist(t(exp_raw_J), method = "manhattan"))

rownames(dists) <- row.names(pData(eset_normalied_pat))


hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
hmcol <- colorRampPalette(rev(brewer.pal(9, "YlOrBr")))(255)
colnames(dists) <- Nano_SampleInfo$CART_Patient_Code_J
diag(dists) <- NA
Nano_SampleInfo$ARM
ann_colors <- list(Cell.type = c(iDC = "grey60", mDC = "white", TMM2 = "black"),
                   Outcome = c(Bad = "black", Good = "dodgerblue"),
                   CD8 = c(No = "grey20", NA. = "white", Yes = "khaki1"),
                   CD4 = c(No = "grey20", NA. = "white", Yes = "khaki1"))

colous <- colorRampPalette(c("white","green", "black","red"))(50)
colous <- colorRampPalette(c("lightblue","dodgerblue","black","gold4","yellow"))(100)

png("myplot_M_MART1.PD1.png", height=7.21, width=5.93, units="in", res=300)
pheatmap(dists, col = colous,
         annotation_col = annotation_for_heatmap,
         annotation_colors = ann_colors,
         cluster_rows = T,
         border_color= "black",
         show_rownames = F,
         show_colnames = F,
         legend = TRUE,
         treeheight_row = 15,
         treeheight_col = 15,
         fontsize_row =8,
         scale="none",
         cutree_rows = 1,
         cutree_cols = 1,
         clustering_distance_rows="manhattan",
         clustering_distance_cols="manhattan",
         clustering_method="complete",
         legend_breaks = c(min(dists, na.rm = TRUE),
                           max(dists, na.rm = TRUE)),
         legend_labels = (c("small distance", "large distance")),
         main = "DC Clustering heatmap ")
dev.off()

########## Filtering based on intensity
DC_medians <- rowMedians(Biobase::exprs(eset_normalied_pat))

######## Histogram of the median intensities per gene
hist_res <- hist(DC_medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

emp_mu <- hist_res$breaks[which.max(hist_res$density)]
emp_sd <- mad(DC_medians)/2
prop_cental <- 1

lines(sort(DC_medians), prop_cental*dnorm(sort(DC_medians),
                                          mean = emp_mu , sd = emp_sd),
      col = "grey10", lwd = 4)

# We visually set a cutoff line man_threshold to the left of the histogram peak in order not to exclude too many genes.
man_threshold <- 1.5

hist_res <- hist(DC_medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)

# Transcripts that do not have intensities larger than the threshold in at least as many arrays as the 
# smallest experimental group are excluded.
# In order to do so, we first have to get a list with the number of samples (=arrays) (no_of_samples) 
# in the experimental groups:
no_of_samples <-
  table(paste0(pData(eset_normalied_pat)$Cell_Type, "_",
               pData(eset_normalied_pat)$Cell_Type))
no_of_samples

iDC_iDC   mDC_mDC TMM2_TMM2 
33        33        33 

# We then create a table of idx_man_threshold to summarize the results and get an overview over 
# how many genes are filtered out.
samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(Biobase::exprs(eset_normalied_pat), 1.5,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
idx_man_threshold
FALSE  TRUE 
16 34643 

DC_eset<- subset(eset_normalied_pat, idx_man_threshold)
dim(DC_eset)
Features  Samples 
34643       99 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

write.table(data.frame(fData(DC_eset),exprs(DC_eset)),file="Expression_for_GSEA_7_26_2019.txt",row.names=FALSE,sep="\t")

write.table(data.frame(exprs(DC_eset)),file="Matrix_MSigDB_7_26_2019.txt",row.names=T,sep="\t")


write.table(data.frame(fData(DC_eset),exprs(DC_eset)),file="Expression_for_GSEA_7_26_2019.txt",row.names=T,sep="\t")


# create a design matrix based on the sample attributes
design <- model.matrix(~factor(DC_eset$Cell_Type))
# print the design matrix
design
colnames(design)

samples <- DC_eset$Cell_Type
samples <- as.factor(samples)
samples
design <- model.matrix(~0 + samples)
design
colnames(design) <- c("iDC", "mDC", "TMM2")
design
# fit the linear model to the filtered expression set
fit <- lmFit(DC_eset, design)
# set up a contrast matrix to compare tissues v cell line

contr.matrix <- makeContrasts(mDCvsiDC = mDC-iDC, 
                              TMM2vsmDC = TMM2-mDC,
                              TMM2vsiDC = TMM2-iDC,
                              levels = colnames(design))


contr.matrix
Contrasts
Levels mDCvsiDC TMM2vsmDC TMM2vsiDC
iDC        -1         0        -1
mDC         1        -1         0
TMM2        0         1         1

# Now the contrast matrix is combined with the per-probeset linear model fit.
DC_fits <- contrasts.fit(fit, contr.matrix)
DC_ebFit <- eBayes(DC_fits)
topTable(DC_ebFit, number=10, coef=1)
topTable(DC_ebFit, number=10, coef=2)
topTable(DC_ebFit, number=10, coef=3)

nrow(topTable(DC_ebFit, coef=1, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=2, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=3, number=Inf, lfc=1))


Genes_COEF_1 <- (topTable(DC_ebFit, coef=1, lfc=1, number=Inf))
Genes_COEF_1 <- Genes_COEF_1[order(-Genes_COEF_1$logFC), ]
head(Genes_COEF_1)
dim(Genes_COEF_1)
write.table(Genes_COEF_1, file="DE_coef1_GENES_mDCvsiDC.txt", sep="\t", quote=F, col.names=NA)

Genes_COEF_2 <- (topTable(DC_ebFit, coef=2, lfc=1, number=Inf))
Genes_COEF_2 <- Genes_COEF_2[order(-Genes_COEF_2$logFC), ]
head(Genes_COEF_2)
dim(Genes_COEF_2)
write.table(Genes_COEF_2, file="DE_coef2_GENES_TMM2vsmDC.txt", sep="\t", quote=F, col.names=NA)

Genes_COEF_3 <- (topTable(DC_ebFit, coef=3, lfc=1, number=Inf))
Genes_COEF_3 <- Genes_COEF_3[order(-Genes_COEF_3$logFC), ]
head(Genes_COEF_3)
dim(Genes_COEF_3)
write.table(Genes_COEF_3, file="DE_coef3_GENES_TMM2vsiDC.txt", sep="\t", quote=F, col.names=NA)

#################### MA Plot Visualization 
################# COEF = 1 ################# ################# ################# ################# 
################# COEF = 1 ################# ################# ################# ################# 
################# COEF = 1 ################# ################# ################# ################# 

A_MA_Plot <- (topTable(DC_ebFit, coef=1, number=Inf))
A_MA_Plot["Expression"] <- "NotSignificant"
A_MA_Plot[which(A_MA_Plot['adj.P.Val'] < 0.05 & abs(A_MA_Plot['logFC']) < 1 ),"Expression"] <- "Significant"
A_MA_Plot[which(A_MA_Plot['adj.P.Val'] > 0.05 & abs(A_MA_Plot['logFC']) > 1 ),"Expression"] <- "FoldChange"
A_MA_Plot[which(A_MA_Plot['adj.P.Val'] < 0.05 & abs(A_MA_Plot['logFC']) > 1 ),"Expression"] <- "Significant 2 FoldChange"

# and generate plot 
UA <- ggplot(data = A_MA_Plot, aes(x=AveExpr, y=logFC)) + geom_point(aes(color=Expression), size=2, alpha=0.2) +xlim(c(0, 14)) + ylim(c(-6.5, 9))+scale_color_manual(values=c("navy","limegreen","red")) + theme_bw()+geom_hline(yintercept=0, color = "black",lwd = 0.5,)+geom_hline(yintercept=1, color = "black",lwd = 0.3,linetype="dashed") +geom_hline(yintercept=-1, color = "black",lwd = 0.3,linetype="dashed")+ggtitle("mDC vs iDC")
UA + theme(plot.title = element_text(color="black", size=14, face="bold"))

##### VOLACON PLOT #####################
##### VOLACON PLOT #####################
##### VOLACON PLOT #####################

library("ggrepel") #Avoid overlapping labels
mutateddfA <- A_MA_Plot[order(-A_MA_Plot$logFC), ]  #Will have different colors depending on significance
mutateddf_A <- mutateddfA[1:10, ]
mutateddf_A

head(A_MA_Plot)

mutateddfA_d <- A_MA_Plot[order(A_MA_Plot$logFC), ]  #Will have different colors depending on significance
mutateddf_Ad <- mutateddfA_d[1:10, ]
mutateddf_Ad


volc = ggplot(A_MA_Plot, aes(logFC, -log10(adj.P.Val))) +xlim(c(-9, 11)) + ylim(c(-1, 87))+ #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Expression), alpha=0.4) + #add points colored by significance
  scale_color_manual(values=c("navy","limegreen","red")) + 
  ggtitle("mDC vs iDC") #e.g. 'Volcanoplot DESeq2'
A <- volc+geom_text_repel(data=mutateddf_A, aes(label=row.names(mutateddf_A)),size=3.5)+geom_vline(xintercept = -1, col="black", lwd = 0.3, linetype="dashed")+ geom_vline(xintercept = 1, col="black", lwd = 0.3, linetype="dashed") +geom_hline(yintercept=0, col="black", lwd = 0.3)+geom_point(data=mutateddfA[1:10, ], colour="orange", size=2, alpha=0.5)+theme_bw() #adding text for the top 10 genes
A +geom_text_repel(data=mutateddf_Ad, aes(label=row.names(mutateddf_Ad)),size=3.5)+geom_vline(xintercept = -1, col="black", lwd = 0.3, linetype="dashed")+ geom_vline(xintercept = 1, col="black", lwd = 0.3, linetype="dashed") +geom_hline(yintercept=0, col="black", lwd = 0.3)+geom_point(data=mutateddfA_d[1:10, ], colour="orange", size=2, alpha=0.5)+theme_bw() #adding text for the top 10 genes


################# COEF = 2 ################# ################# ################# ################# 
################# COEF = 2 ################# ################# ################# ################# 
################# COEF = 2 ################# ################# ################# ################# 

B_MA_Plot <- (topTable(DC_ebFit, coef=2, number=Inf))
B_MA_Plot["Expression"] <- "NotSignificant"
B_MA_Plot[which(B_MA_Plot['adj.P.Val'] < 0.05 & abs(B_MA_Plot['logFC']) < 1 ),"Expression"] <- "Significant"
B_MA_Plot[which(B_MA_Plot['adj.P.Val'] > 0.05 & abs(B_MA_Plot['logFC']) > 1 ),"Expression"] <- "FoldChange"
B_MA_Plot[which(B_MA_Plot['adj.P.Val'] < 0.05 & abs(B_MA_Plot['logFC']) > 1 ),"Expression"] <- "Significant 2 FoldChange"

# and generate plot 
UB <- ggplot(data = B_MA_Plot, aes(x=AveExpr, y=logFC)) + 
  geom_point(aes(color=Expression), size=2, alpha=0.35) +
  xlim(c(0, 14)) + ylim(c(-2, 5))+
  scale_color_manual(values=c("navy","limegreen","red")) + 
  theme_bw()+geom_hline(yintercept=0, color = "black",lwd = 0.5,)+
  geom_hline(yintercept=1, color = "black",lwd = 0.3,linetype="dashed") +
  geom_hline(yintercept=-1, color = "black",lwd = 0.3,linetype="dashed")+ 
  annotate("text", x = 4.2, y=4.8, label = "TYR")+ggtitle("TMM2 vs mDC")

UB + theme(plot.title = element_text(color="black", size=14, face="bold"))

##### VOLACON PLOT #####################
##### VOLACON PLOT #####################
##### VOLACON PLOT #####################

library("ggrepel") #Avoid overlapping labels
mutateddfB <- B_MA_Plot[order(-B_MA_Plot$logFC), ]  #Will have different colors depending on significance
mutateddf_B <- mutateddfB[1:10, ]
mutateddf_B

mutateddfB_d <- B_MA_Plot[order(B_MA_Plot$logFC), ]  #Will have different colors depending on significance
mutateddf_Bd <- mutateddfB_d[1:10, ]
mutateddf_Bd

volc = ggplot(B_MA_Plot, aes(logFC, -log10(adj.P.Val))) +xlim(c(-3, 5)) + ylim(c(-1, 40))+ #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Expression), alpha=0.7) + #add points colored by significance
  scale_color_manual(values=c("navy","limegreen","red")) + 
  ggtitle("TMM2 vs mDC") #e.g. 'Volcanoplot DESeq2'

B <- volc+geom_text_repel(data=mutateddf_B, aes(label=row.names(mutateddf_B)),size=3.5)+geom_vline(xintercept = -1, col="black", lwd = 0.3, linetype="dashed")+ geom_vline(xintercept = 1, col="black", lwd = 0.3, linetype="dashed") +geom_hline(yintercept=0, col="black", lwd = 0.3)+geom_point(data=mutateddf_B[1:10, ], colour="orange", size=2, alpha=0.5)+theme_bw() #adding text for the top 10 genes
B +geom_text_repel(data=mutateddf_Bd, aes(label=row.names(mutateddf_Bd)),size=3.5)+geom_vline(xintercept = -1, col="black", lwd = 0.3, linetype="dashed")+ geom_vline(xintercept = 1, col="black", lwd = 0.3, linetype="dashed") +geom_hline(yintercept=0, col="black", lwd = 0.3)+geom_point(data=mutateddf_Bd[1:10, ], colour="orange", size=2, alpha=0.5)+theme_bw() #adding text for the top 10 genes

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ labeling points in graph 
B_MA_Plot <- (topTable(DC_ebFit, coef=2, number=Inf))
B_MA_Plot["Expression"] <- "Expression"
B_MA_Plot[which(B_MA_Plot['adj.P.Val'] < 0.05 & abs(B_MA_Plot['logFC']) < 1 ),"Expression"] <- "Significant"
B_MA_Plot[which(B_MA_Plot['adj.P.Val'] > 0.05 & abs(B_MA_Plot['logFC']) > 1 ),"Expression"] <- "FoldChange"
B_MA_Plot[which(B_MA_Plot['adj.P.Val'] < 0.05 & abs(B_MA_Plot['logFC']) > 1 ),"Expression"] <- "Significant 2 FoldChange"

TYR_MGAE_MLA <- subset(B_MA_Plot, rownames(B_MA_Plot) %in% c("TYR","MLANA","MAGEA6"))
TYR_MGAE_MLA

png("MAPlot_TMM2_vs_mDC_MA_labeled.png", height=3.75, width=5.60, units="in", res=600)
ggplot(data = B_MA_Plot, aes(x=AveExpr, y=logFC)) + 
  geom_point(aes(color=Expression), size=2.5, alpha=0.5) +
  xlim(c(0, 14)) + ylim(c(-2, 5))+
  scale_color_manual(values=c("navy","limegreen","red")) + 
  theme_bw()+geom_hline(yintercept=0, color = "black",lwd = 0.5,)+
  geom_hline(yintercept=1, color = "black",lwd = 0.3,linetype="dashed") +
  geom_hline(yintercept=-1, color = "black",lwd = 0.3,linetype="dashed")+ 
  ggtitle("TMM2 vs mDC")+
  geom_text_repel(data=TYR_MGAE_MLA, aes(label=row.names(TYR_MGAE_MLA)),min.segment.length = 1,size=3.5,segment.size = 0.2,nudge_x = -0.5,nudge_y = 0.2)+
  geom_point(data=TYR_MGAE_MLA, colour="orange", size=2.5, alpha=0.75)
dev.off()


volc = ggplot(B_MA_Plot, aes(logFC, -log10(adj.P.Val))) +xlim(c(-3, 5)) + ylim(c(-1, 40))+ #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Expression),size=2.5, alpha=0.5) + #add points colored by significance
  scale_color_manual(values=c("navy","limegreen","red")) + 
  ggtitle("TMM2 vs mDC") #e.g. 'Volcanoplot DESeq2'

png("MVolcano_TMM2_vs_mDC_MA_labeled.png", height=3.75, width=5.60, units="in", res=600)
volc+geom_text_repel(data=TYR_MGAE_MLA, aes(label=row.names(TYR_MGAE_MLA)),size=3.5)+
  geom_vline(xintercept = -1, col="black", lwd = 0.3, linetype="dashed")+ 
  geom_vline(xintercept = 1, col="black", lwd = 0.3, linetype="dashed") +
  geom_hline(yintercept=0, col="black", lwd = 0.3)+
  geom_point(data=TYR_MGAE_MLA, colour="orange", size=2.5, alpha=0.9)+
  theme_bw() #adding text for the top 10 genes
dev.off()





################# COEF = 3 ################# ################# ################# ################# 
################# COEF = 3 ################# ################# ################# ################# 
################# COEF = 3 ################# ################# ################# ################# 

C_MA_Plot <- (topTable(DC_ebFit, coef=1, number=Inf))
C_MA_Plot["Expression"] <- "NotSignificant"
C_MA_Plot[which(A_MA_Plot['adj.P.Val'] < 0.05 & abs(C_MA_Plot['logFC']) < 1 ),"Expression"] <- "Significant"
C_MA_Plot[which(A_MA_Plot['adj.P.Val'] > 0.05 & abs(C_MA_Plot['logFC']) > 1 ),"Expression"] <- "FoldChange"
C_MA_Plot[which(A_MA_Plot['adj.P.Val'] < 0.05 & abs(C_MA_Plot['logFC']) > 1 ),"Expression"] <- "Significant 2 FoldChange"
# and generate plot 
UC <- ggplot(data = C_MA_Plot, aes(x=AveExpr, y=logFC)) + geom_point(aes(color=Expression), size=2, alpha=0.2) +xlim(c(0, 14)) + ylim(c(-6.5, 9))+scale_color_manual(values=c("navy","limegreen","red")) + theme_bw()+geom_hline(yintercept=0, color = "black",lwd = 0.5,)+geom_hline(yintercept=1, color = "black",lwd = 0.3,linetype="dashed") +geom_hline(yintercept=-1, color = "black",lwd = 0.3,linetype="dashed")+ggtitle("iDC vs TMM2")
UC + theme(plot.title = element_text(color="black", size=14, face="bold"))

##### VOLACON PLOT #####################
##### VOLACON PLOT #####################
##### VOLACON PLOT #####################

library("ggrepel") #Avoid overlapping labels
mutateddfC <- C_MA_Plot[order(-C_MA_Plot$logFC), ]  #Will have different colors depending on significance
mutateddf_C <- mutateddfC[1:10, ]
mutateddf_C

mutateddfC_d <- C_MA_Plot[order(C_MA_Plot$logFC), ]  #Will have different colors depending on significance
mutateddf_Cd <- mutateddfC_d[1:10, ]
mutateddf_Cd


volc = ggplot(C_MA_Plot, aes(logFC, -log10(adj.P.Val))) +xlim(c(-9, 11)) + ylim(c(-1, 87))+ #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Expression), alpha=0.4) + #add points colored by significance
  scale_color_manual(values=c("navy","limegreen","red")) + 
  ggtitle("TMM2 vs iDC") #e.g. 'Volcanoplot DESeq2'
C <- volc+geom_text_repel(data=mutateddf_C, aes(label=row.names(mutateddf_C)),size=3.5)+geom_vline(xintercept = -1, col="black", lwd = 0.3, linetype="dashed")+ geom_vline(xintercept = 1, col="black", lwd = 0.3, linetype="dashed") +geom_hline(yintercept=0, col="black", lwd = 0.3)+geom_point(data=mutateddfC[1:10, ], colour="navy", size=2, alpha=0.5)+theme_bw() #adding text for the top 10 genes
C +geom_text_repel(data=mutateddf_Cd, aes(label=row.names(mutateddf_Cd)),size=3.5)+geom_vline(xintercept = -1, col="black", lwd = 0.3, linetype="dashed")+ geom_vline(xintercept = 1, col="black", lwd = 0.3, linetype="dashed") +geom_hline(yintercept=0, col="black", lwd = 0.3)+geom_point(data=mutateddfC_d[1:10, ], colour="navy", size=2, alpha=0.5)+theme_bw() #adding text for the top 10 genes

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# create a design matrix based on the sample attributes
design <- model.matrix(~factor(DC_eset$Cell_Type))
# print the design matrix
design
colnames(design)

samples <- DC_eset$Cell_Type
samples <- as.factor(samples)
samples

OUTCOME <- as.factor(DC_eset$Outcome_1)
OUTCOME

z <- factor(paste(OUTCOME, samples, sep = "_"))
design <- model.matrix(~0+z)
colnames(design) <- gsub("z", "", colnames(design))
design

# fit the linear model to the filtered expression set
fit <- lmFit(DC_eset, design)
# set up a contrast matrix to compare tissues v cell line

contr.matrix <- makeContrasts(Good_mDCvsGood_iDC = Good_mDC-Good_iDC, 
                              Good_TMM2vsGood_mDC = Good_TMM2-Good_mDC,
                              Bad_mDCvsBad_iDC = Bad_mDC-Bad_iDC,
                              Bad_TMM2vsBad_mDC = Bad_TMM2-Bad_mDC,
                              Good_iDCsBad_iDC = Good_iDC-Bad_iDC,
                              Good_mDCvsBad_mDC = Good_mDC-Bad_mDC,
                              Good_TMM2vsBad_TMM2 = Good_TMM2-Bad_TMM2,
                              levels = colnames(design))

contr.matrix
Contrasts
Levels      Good_mDCvsGood_iDC Good_TMM2vsGood_mDC Bad_mDCvsBad_iDC Bad_TMM2vsBad_mDC
Bad_iDC                    0                   0               -1                 0
Bad_mDC                    0                   0                1                -1
Bad_TMM2                   0                   0                0                 1
Good_iDC                  -1                   0                0                 0
Good_mDC                   1                  -1                0                 0
Good_TMM2                  0                   1                0                 0

# Now the contrast matrix is combined with the per-probeset linear model fit.
DC_fits <- contrasts.fit(fit, contr.matrix)
DC_ebFit <- eBayes(DC_fits)
topTable(DC_ebFit, number=10, coef=1)
topTable(DC_ebFit, number=10, coef=2)
topTable(DC_ebFit, number=10, coef=3)
topTable(DC_ebFit, number=10, coef=4)
topTable(DC_ebFit, number=10, lfc=1, coef=5)
topTable(DC_ebFit, number=10, lfc=1, coef=6)
topTable(DC_ebFit, number=10, lfc=1, coef=7)

nrow(topTable(DC_ebFit, coef=1, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=2, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=3, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=4, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=5, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=6, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=7, number=Inf, lfc=1))

"Good_mDCvsGood_iDC"
Genes_COEF_1 <- (topTable(DC_ebFit, coef=1, lfc=1, number=Inf))
Genes_COEF_1 <- Genes_COEF_1[order(-Genes_COEF_1$logFC), ]
head(Genes_COEF_1)
dim(Genes_COEF_1)
write.table(Genes_COEF_1, file="DE_GENES_Good_mDCvsGood_iDC.txt", sep="\t", quote=F, col.names=NA)

"Good_TMM2vsGood_mDC"
Genes_COEF_2 <- (topTable(DC_ebFit, coef=1, lfc=1, number=Inf))
Genes_COEF_1 <- Genes_COEF_1[order(-Genes_COEF_1$logFC), ]
head(Genes_COEF_1)
dim(Genes_COEF_1)
write.table(Genes_COEF_1, file="DE_coef1_GENES_mDCvsiDC.txt", sep="\t", quote=F, col.names=NA)

"Bad_mDCvsBad_iDC"
Genes_COEF_3 <- (topTable(DC_ebFit, coef=3, lfc=1, number=Inf))
Genes_COEF_3 <- Genes_COEF_3[order(-Genes_COEF_3$logFC), ]
head(Genes_COEF_3)
dim(Genes_COEF_3)
write.table(Genes_COEF_3, file="DE_GENES_Bad_mDCvsBad_iDC.txt", sep="\t", quote=F, col.names=NA)

"Bad_TMM2vsBad_mDC"
Genes_COEF_4 <- (topTable(DC_ebFit, coef=1, lfc=1, number=Inf))
Genes_COEF_1 <- Genes_COEF_1[order(-Genes_COEF_1$logFC), ]
head(Genes_COEF_1)
dim(Genes_COEF_1)
write.table(Genes_COEF_1, file="DE_coef1_GENES_mDCvsiDC.txt", sep="\t", quote=F, col.names=NA)


#################################################################################################################################
PREVIOUS THERAPY
####################################################################
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  # create a design matrix based on the sample attributes
design <- model.matrix(~factor(DC_eset$Cell_Type))
# print the design matrix
design
colnames(design)

table(DC_eset$Prev_Therapy)
samples <- DC_eset$Cell_Type
samples <- as.factor(samples)
samples

PREV_THY <- as.factor(DC_eset$Prev_Therapy)
PREV_THY

z <- factor(paste(PREV_THY, samples, sep = "_"))
design <- model.matrix(~0+z)
colnames(design) <- gsub("z", "", colnames(design))
design

# fit the linear model to the filtered expression set
fit <- lmFit(DC_eset, design)
# set up a contrast matrix to compare tissues v cell line

contr.matrix <- makeContrasts(Checkpoint_iDCvsNoneiDC  = Checkpoint_iDC-None_iDC, 
                              Checkpoint_mDCsNone_mDC = Checkpoint_mDC-None_mDC,
                              Checkpoint_TMM2vsNone_TMM2 = Checkpoint_TMM2 -None_TMM2,
                              levels = colnames(design))

contr.matrix
Contrasts
Levels      Good_mDCvsGood_iDC Good_TMM2vsGood_mDC Bad_mDCvsBad_iDC Bad_TMM2vsBad_mDC
Bad_iDC                    0                   0               -1                 0
Bad_mDC                    0                   0                1                -1
Bad_TMM2                   0                   0                0                 1
Good_iDC                  -1                   0                0                 0
Good_mDC                   1                  -1                0                 0
Good_TMM2                  0                   1                0                 0

# Now the contrast matrix is combined with the per-probeset linear model fit.
DC_fits <- contrasts.fit(fit, contr.matrix)
DC_ebFit <- eBayes(DC_fits)
topTable(DC_ebFit, number=10, coef=1)
topTable(DC_ebFit, number=10, coef=2)
topTable(DC_ebFit, number=10, coef=3)

nrow(topTable(DC_ebFit, coef=1, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=2, number=Inf, lfc=1))
nrow(topTable(DC_ebFit, coef=3, number=Inf, lfc=1))

Genes_COEF_1 <- topTable(DC_ebFit, coef=1, number=Inf, lfc=1)
write.table(Genes_COEF_1, file="AdvTMM2_Previous_Therapy.txt", sep="\t", quote=F, col.names=NA)



##### PLOTS ##################### ##### PLOTS ##################### ##### PLOTS #####################
##### PLOTS ##################### ##### PLOTS ##################### ##### PLOTS #####################
### CONTINUOS VARIABLE 
# create a design matrix based on the sample attributes

Subsetting eset object:

head(exprs(DC_eset[ ,(1:33)]))
head(exprs(DC_eset[ ,(34:66)]))
head(exprs(DC_eset[ ,(67:99)]))


pData(DC_eset)
exprs(DC_eset)
dim(DC_eset)
AdvTMM2_eset <- DC_eset[ ,(67:99)]
exprs(AdvTMM2_eset)

dim(pData(DC_eset[ ,(67:99)]))
dim(pData(AdvTMM2_eset))

iDC_eset <-DC_eset[ ,(1:33)]
mDC_eset <-DC_eset[ ,(34:66)]


#############################################################################################################################3
TYR

samples <- AdvTMM2_eset$Cell_Type
samples <- as.factor(samples)
samples

TYROSINASE <- AdvTMM2_eset$TYR
mm <- model.matrix(~samples+TYROSINASE)

fit <- lmFit(AdvTMM2_eset, mm)
tmp <- contrasts.fit(fit, coef = 4) # test "MAGE.A6" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- topTable(tmp, sort.by = "logFC", n = Inf)
head(top.table, 20)

write.table(top.table, file="TMM2_DE_Gene_Corrrelation_TYROSINASE.txt", sep="\t", quote=F, col.names=NA,)

#############################################################################################################################3
MLANA

pData(DC_eset)
exprs(DC_eset)
dim(DC_eset)
AdvTMM2_eset <- DC_eset[ ,(67:99)]
exprs(AdvTMM2_eset)

dim(pData(DC_eset[ ,(67:99)]))
dim(pData(AdvTMM2_eset))

samples <- AdvTMM2_eset$Cell_Type
samples <- as.factor(samples)
samples

MLN <- AdvTMM2_eset$MLANA.
mm <- model.matrix(~samples+MLN)

fit <- lmFit(AdvTMM2_eset, mm)
tmp <- contrasts.fit(fit, coef = 4) # test "MAGE.A6" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- topTable(tmp, sort.by = "logFC", n = Inf)
head(top.table, 20)

write.table(top.table, file="TMM2_DE_Gene_Corrrelation_MLANA.txt", sep="\t", quote=F, col.names=NA,)
  
#############################################################################################################################3
MAGEA6

pData(DC_eset)
exprs(DC_eset)
dim(DC_eset)
AdvTMM2_eset <- DC_eset[ ,(67:99)]
exprs(AdvTMM2_eset)

dim(pData(DC_eset[ ,(67:99)]))
dim(pData(AdvTMM2_eset))

samples <- AdvTMM2_eset$MAGEA6
samples <- as.factor(samples)
samples

MLN <- AdvTMM2_eset$MLANA.
mm <- model.matrix(~samples+MLN)

fit <- lmFit(AdvTMM2_eset, mm)
tmp <- contrasts.fit(fit, coef = 4) # test "MAGE.A6" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- topTable(tmp, sort.by = "logFC", n = Inf)
head(top.table, 20)

write.table(top.table, file="TMM2_DE_Gene_Corrrelation_MLANA.txt", sep="\t", quote=F, col.names=NA,)

#############################################################################################################################3
O Survival 

samples <- AdvTMM2_eset$Cell_Type
samples <- as.factor(samples)
samples

O.Survival <- AdvTMM2_eset$O_Survival
mm <- model.matrix(~samples+O.Survival)

fit <- lmFit(AdvTMM2_eset, mm)
tmp <- contrasts.fit(fit, coef = 4) # test "MAGE.A6" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- topTable(tmp, sort.by = "logFC", n = Inf)
head(top.table, 20)

write.table(top.table, file="TMM2_DE_Gene_Corrrelation_MLANA.txt", sep="\t", quote=F, col.names=NA,)

#############################################################################################################################3
PF Survival 

samples <- AdvTMM2_eset$Cell_Type
samples <- as.factor(samples)
samples

PF.Survival <- AdvTMM2_eset$PF_Survival
mm <- model.matrix(~samples+O.Survival)

fit <- lmFit(AdvTMM2_eset, mm)
tmp <- contrasts.fit(fit, coef = 4) # test "MAGE.A6" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- topTable(tmp, sort.by = "logFC", n = Inf)
head(top.table, 20)

write.table(top.table, file="TMM2_DE_Gene_Corrrelation_MLANA.txt", sep="\t", quote=F, col.names=NA,)
  
#############################################################################################################################3
CD8 T cells: 

samples <- AdvTMM2_eset$Cell_Type
samples <- as.factor(samples)
samples

CD8_Resp <- AdvTMM2_eset$CD8_Resp_Per_Patient
mm <- model.matrix(~samples+CD8_Resp)

fit <- lmFit(AdvTMM2_eset, mm)
tmp <- contrasts.fit(fit, coef = 4) # test "MAGE.A6" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- topTable(tmp, sort.by = "logFC", n = Inf)
head(top.table, 20)

write.table(top.table, file="TMM2_DE_Gene_Corrrelation_CD8_Resp.txt", sep="\t", quote=F, col.names=NA,)

  
#############################################################################################################################3
VEGFA: 
  
samples <- AdvTMM2_eset$Cell_Type
samples <- as.factor(samples)
samples

VEGFA <- AdvTMM2_eset$PDL2
mm <- model.matrix(~samples+VEGFA)

fit <- lmFit(AdvTMM2_eset, mm)
tmp <- contrasts.fit(fit, coef = 4) # test "MAGE.A6" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- topTable(tmp, sort.by = "logFC", n = Inf)
head(top.table, 20)

write.table(top.table, file="TMM2_DE_Gene_Corrrelation_VEGFA.txt", sep="\t", quote=F, col.names=NA,)

#############################################################################################################################
#############################################################################################################################
iDC and mDC cells
Overall Survival:
  
iDC
samples <- mDC_eset$Cell_Type
samples <- as.factor(samples)
samples

PF_Survival <- mDC_eset$PF_Survival
mm <- model.matrix(~samples+PF_Survival)

fit <- lmFit(mDC_eset, mm)
tmp <- contrasts.fit(fit, coef = 4) # test "MAGE.A6" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- topTable(tmp, sort.by = "logFC", n = Inf)
head(top.table, 20)

write.table(top.table, file="TMM2_DE_Gene_Corrrelation_VEGFA.txt", sep="\t", quote=F, col.names=NA,)

iDC_eset <-DC_eset[ ,(1:33)]
mDC_eset <-DC_eset[ ,(34:66)]

#######################################################
note on fit object subsetting:
fit.day <- fit1[,11:13]
colnames(fit.day)
[1] "day2" "day3" "day4"
> head(fit.day$coef)


