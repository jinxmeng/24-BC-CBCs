#### Zhewen Qin, 20240805 ####

###加载所需要的包
library(tidydr)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(msigdbr)

temp <- c('GSE117783_Crypts_processed.csv')
scRNAlist <- list()
for(i in 1:length(temp)){
  counts <- read.csv(temp[i],header = TRUE, row.names = 1)
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells = 3, min.features =200)
}

pbmc <- scRNAlist[[1]]

pbmc[["group"]] <- ifelse(str_detect(pbmc@meta.data$orig.ident,"C05"),"Normal","Irradiated")

#标准化，默认是log2,针对每个细胞的
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
##建议只对高变基因进行标准化(官方推荐方法)
all.genes <-  VariableFeatures(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
##线性降维（PCA）,必须用高变基因集,但也可通过features参数自己指定；
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

ElbowPlot(pbmc, ndims=50, reduction="pca")
##细胞聚类
##非线性降维（UMAP/tSNE)基于PCA空间中的欧氏距离计算nearest neighbor graph,优化任意两个细胞间的距离权重（输入上一步得到的PC维数）。
pbmc <- FindNeighbors(pbmc, dims = 1:15)

##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
pbmc <- FindClusters(pbmc, resolution = 0.8)##0.5可以往上下调整，以增加或减少分群的cluster数

##如需要计算UMAP，适合大型数据

pbmc <- RunTSNE(pbmc, dims = 1:15)

DimPlot(pbmc,reduction = "tsne",label = TRUE,pt.size = 0.5,cols = col.pal)

FeaturePlot(pbmc,reduction = "tsne",features = c("Lgr5","Ppara","Hmgcs2","Hopx"))
FeaturePlot(pbmc,reduction = "tsne",features = c("Dclk1","Trpm5","Chgb","Tac1","Defa24","Defa22","Lyz1",
                                                 'Ptprc',"Tff3","Agr2","Muc2","Ascl2","Fabp1","Alpi","Stmn1","Cd44"))

pbmc <- subset(pbmc,idents = "19",invert=T)
pbmc <- RenameIdents(pbmc,"3" = "Immune cells","10" = "Immune cells","17" = "Immune cells",
                     "18"="Tuft","7" = "CBCs","16" = "CBCs","15" = "Entero-endocrine",
                     "4" = "Enterocytes","5" = "Enterocytes","11" = "Enterocytes","13" = "Enterocytes",
                     "0" = "Enterocytes","1" = "Enterocytes","2" = "Enterocytes","14" = "Enterocytes",
                     "9" = "Unknown","8" = "Goblets","12" = "Goblets","6" = "Paneth")               

col.pal <- c("#19937d","#4ca8bf","#cf4b35","#c98d62","#709356","#58869c","#71b666","#ae9a84","#8bc3b6","#41529b","#bb3e31","#709356","#e2da84","#44657d","#af5f39","#76215f","#71b666","#77181c","#d4d4cc","#edab4c","#a7ab7b","#c98d62","#a8706b","#878a78","#6f5461","#e26b21","#b11f23","#7b3d87","#568d90","#e68b9c","#040000")
DimPlot(pbmc,reduction = "tsne",label = TRUE,pt.size = 0.5,cols =col.pal )


DoHeatmap(pbmc,features = c("Dclk1","Trpm5","Chgb","Tac1","Defa24","Defa22","Lyz1",
                            'Ptprc',"Tff3","Agr2","Muc2","Ascl2","Fabp1","Alpi","Stmn1","Cd44") )
DotPlot(pbmc,features = c("Dclk1","Trpm5","Chgb","Tac1","Defa24","Defa22","Lyz1",
                          'Ptprc',"Tff3","Agr2","Muc2","Fabp1","Alpi","Ascl2","Stmn1","Cd44"),cols = "PiYG",dot.scale =10)
