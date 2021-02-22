#### step2.1-2.6（对应如下的1-6） ####

### 1、下载、探索、整理数据----
## 1.1 下载、探索数据
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465
sessionInfo()
# 读取文件耗时比较长，请耐心等待
a <- read.table("../../rawdata/GSE84465_GBM_All_data.csv.gz")
a[1:4,1:4]
#行名为symbol ID
#列名为sample，看上去像是两个元素的组合。
summary(a[,1:4]) 
boxplot(a[,1:4])
head(rownames(a))
tail(rownames(a),10)
# 可以看到原文的counts矩阵来源于htseq这个计数软件，所以有一些不是基因的行需要剔除：
#  "no_feature"           "ambiguous"            "too_low_aQual"        "not_aligned"          "alignment_not_unique"
tail(a[,1:4],10)

a=a[1:(nrow(a)-5),]

#原始counts数据

#3,589 cells of 4 human primary GBM samples, accession number GSE84465
#2,343 cells from tumor cores and 1,246 cells from peripheral regions
b <- read.table("../../rawdata/SraRunTable.txt",
                sep = ",", header = T)
b[1:4,1:4]
table(b$Patient_ID) # 4 human primary GBM samples
table(b$TISSUE) # tumor cores and peripheral regions
table(b$TISSUE,b$Patient_ID)


## 1.2 整理数据 
# tumor and peripheral 分组信息
head(colnames(a))
head(b$plate_id)
head(b$Well)
#a矩阵行名（sample）并非为GSM编号，而主要是由相应的plate_id与Well组合而成

b.group <- b[,c("plate_id","Well","TISSUE","Patient_ID")]
paste0("X",b.group$plate_id[1],".",b.group$Well[1])
b.group$sample <- paste0("X",b.group$plate_id,".",b.group$Well)
head(b.group)
identical(colnames(a),b.group$sample)

# 筛选tumor cell
index <- which(b.group$TISSUE=="Tumor")
length(index)
group <- b.group[index,] #筛选的是行
head(group)

a.filt <- a[,index] #筛选的是列
dim(a.filt)
identical(colnames(a.filt),group$sample)

sessionInfo()


### 2、构建seurat对象，质控绘图----
# 2.1 构建seurat对象，质控
#In total, 2,343 cells from tumor cores were included in this analysis.
#quality controlstandards: 
#1) genes detected in < 3 cells were excluded; 筛选基因
#2) cells with < 50 total detected genes were excluded; 筛选细胞 
#3) cells with ≥ 5% of mitochondria-expressed genes were excluded. 筛选细胞
sessionInfo()
library("Seurat")
?CreateSeuratObject
sce.meta <- data.frame(Patient_ID=group$Patient_ID,
                       row.names = group$sample)
head(sce.meta)
table(sce.meta$Patient_ID)
# 这个函数 CreateSeuratObject 有多种多样的执行方式
scRNA = CreateSeuratObject(counts=a.filt,
                           meta.data = sce.meta,
                           min.cells = 3, 
                           min.features = 50)
#counts:a matrix-like object with unnormalized data with cells as columns and features as rows 
#meta.data:Additional cell-level metadata to add to the Seurat object
#meta.data: features detected in at least this many cells. 
#min.features:cells where at least this many features are detected.
head(scRNA@meta.data)
#nCount_RNA：the number of cell total counts
#nFeature_RNA：the number of cell's detected gene
summary(scRNA@meta.data)
scRNA@assays$RNA@counts[1:4,1:4]
# 可以看到，之前的counts矩阵存储格式发生了变化：4 x 4 sparse Matrix of class "dgCMatrix"

dim(scRNA)
#  20050  2342 仅过滤掉一个细胞
#接下来根据染色体基因筛选低质量细胞
#Calculate the proportion of transcripts mapping to mitochondrial genes
table(grepl("^MT-",rownames(scRNA)))
#FALSE 



#20050 没有染色体基因
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
head(scRNA@meta.data)
summary(scRNA@meta.data)
#结果显示没有线粒体基因，因此这里过滤也就没有意义，但是代码留在这里
# 万一大家的数据里面有线粒体基因，就可以如此这般进行过滤啦。
pctMT=5 #≥ 5% of mitochondria-expressed genes
scRNA <- subset(scRNA, subset = percent.mt < pctMT)
dim(scRNA)

table(grepl("^ERCC-",rownames(scRNA)))
#FALSE  TRUE 
#19961    86  发现是有ERCC基因
#External RNA Control Consortium，是常见的已知浓度的外源RNA分子spike-in的一种
#指标含义类似线粒体含量，ERCC含量大，则说明total sum变小
scRNA[["percent.ERCC"]] <- PercentageFeatureSet(scRNA, pattern = "^ERCC-")
head(scRNA@meta.data)
summary(scRNA@meta.data)
rownames(scRNA)[grep("^ERCC-",rownames(scRNA))]
#可以看到有不少ERCC基因

sum(scRNA$percent.ERCC< 40)
#较接近原文过滤数量2149，但感觉条件有点宽松了，先做下去看看
#网上看了相关教程，一般ERCC占比不高于10%
sum(scRNA$percent.ERCC< 10)   #就只剩下461个cell，明显低于文献中的数量
pctERCC=40
scRNA <- subset(scRNA, subset = percent.ERCC < pctERCC)
dim(scRNA)
# 20047  2142   原文为19752  2149
dim(a.filt)
#23460  2343 未过滤前


# 2.2 可视化
#图A：观察不同组cell的counts、feature分布
col.num <- length(unique(scRNA@meta.data$Patient_ID))
library(ggplot2)

p1_1.1 <- VlnPlot(scRNA,
                features = c("nFeature_RNA"),
                group.by = "Patient_ID",
                cols =rainbow(col.num)) +
  theme(legend.position = "none") +
  labs(tag = "A")
p1_1.1
p1_1.2 <- VlnPlot(scRNA,
                features = c("nCount_RNA"),
                group.by = "Patient_ID",
                cols =rainbow(col.num)) +
  theme(legend.position = "none") 
p1_1.2
p1_1 <- p1_1.1 | p1_1.2
p1_1
VlnPlot(scRNA,
        features = c("nFeature_RNA","nCount_RNA","percent.ERCC"))
#图B：nCount_RNA与对应的nFeature_RNA关系
p1_2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                       group.by = "Patient_ID",pt.size = 1.3) +
  labs(tag = "B")
p1_2
FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ERCC")
sessionInfo()

### 3、挑选hvg基因，可视化----
#highly Variable gene:简单理解sd大的
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 1500) 
#根据文献原图，挑选变化最大的1500个hvg
top10 <- head(VariableFeatures(scRNA), 10) 
top10
plot1 <- VariableFeaturePlot(scRNA) 
#标记top10 hvg
p1_3 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) +
  theme(legend.position = c(0.1,0.8)) +
  labs(tag = "C")
p1_3

#看看ERCC
ERCC <- rownames(scRNA)[grep("^ERCC-",rownames(scRNA))]
LabelPoints(plot = plot1, points = ERCC, repel = TRUE, 
            size=2.5,colour = "blue") +
  theme(legend.position = c(0.1,0.8)) +
  labs(tag = "C")
#可以直观看到ERCC均不是高变基因，而且部分的ERCC基因表达量确实很高
?LabelPoints
p1_1 | p1_2 | p1_3 #上图

# 这里展开介绍一下 scater 
# https://bioconductor.org/packages/release/bioc/html/scater.html
# Single-cell analysis toolkit工具箱 for expression 
#scater contains tools to help with the analysis of single-cell transcriptomic data, 
#focusing on low-level steps such as quality control, normalization and visualization.
#based on the SingleCellExperiment class (from the SingleCellExperiment package)
#关于sce对象，https://www.jianshu.com/p/9bba0214844b
library(scater)
ct=as.data.frame(scRNA@assays$RNA@counts)
pheno_data=scRNA@meta.data
sce <- SingleCellExperiment(
  assays = list(counts = ct), 
  colData = pheno_data
)
#SingleCellExperiment是SingleCellExperiment包的函数；在加载scater包时会一起加载
sce
?stand_exprs
stand_exprs(sce) <- log2(
  calculateCPM(sce ) + 1)  #只对自己的文库的标准化
assays(sce)
sum(counts(sce)[,1])
head(counts(sce)[,1])
log2(1*10^6/422507+1)
#logcounts(sce)[1:4,1:4]
#exprs(sce)[1:4,1:4]


sce <- logNormCounts(sce) #可以考虑不同细胞的文库差异的标准化
assays(sce)
logcounts(sce)[1:4,1:4]

#https://osca.bioconductor.org/normalization.html#spike-norm
#基于ERCC的标准化方式也有许多优势（不同cell的量理论上是一样的），详见链接
#关于一些常见的FPKM等方式在番外篇会有简单的介绍与学习



#观察上面确定的top10基因在四个样本的分布比较
plotExpression(sce, top10 ,
               x = "Patient_ID",  colour_by = "Patient_ID", 
               exprs_values = "logcounts") 
# 下面的绘图非常耗时：(保存为本地文件查看比较高效，建议为pdf文件)
p1 <- plotHighestExprs(sce, exprs_values = "logcounts")
ggsave("../../out/2.3HighestExprs.pdf", plot = p1, width = 15, height = 18) 
#如果按照ERCC 40%的过滤标准，ERCC表达量也十分大
?plotHighestExprs
# Sometimens few spike-in transcripts may also be present here, 
# though if all of the spike-ins are in the top 50, 
# it suggests that too much spike-in RNA was added

#后续步骤暂时还按照pctERCC=40的过滤标准的结果进行分析

#tips:后面主要还是基于Seurat对象，此处可以删除该变量，节约内存。


### 4、降维，PCA分析，可视化----
#先进行归一化（正态分布）
scRNA <- ScaleData(scRNA, features = (rownames(scRNA)))
#储存到"scale.data"的slot里
GetAssayData(scRNA,slot="scale.data",assay="RNA")[1:8,1:4]
#对比下原来的count矩阵
GetAssayData(scRNA,slot="counts",assay="RNA")[1:8,1:4]
#scRNA@assays$RNA@
#PCA降维，利用之前挑选的hvg，可提高效率
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
#挑选第一，第二主成分对cell可视化
DimPlot(scRNA, reduction = "pca", group.by="Patient_ID")
#发现与原文献中颠倒了
?DimPlot
?RunPCA
#seed.use	:Set a random seed. By default, sets the seed to 42. 
#Setting NULL will not set a seed.
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA),seed.use=3)
#尝试了seed.use的不同取值发现图形只有四种变化（四个拐角），其中以seed.use=3为代表的一类与原文文献一致
DimPlot(scRNA, reduction = "pca", group.by="Patient_ID")
#与文献一致了。个人觉得颠倒与否如果只是随机种子的差别的话，对后续分析应该没影响
p2_1 <- DimPlot(scRNA, reduction = "pca", group.by="Patient_ID")+
  labs(tag = "D")
p2_1

DimPlot(scRNA, reduction = "pca",  split.by = 'Patient_ID')

#挑选主成分，RunPCA默认保留了前50个
scRNA <- JackStraw(scRNA,reduction = "pca", dims=20)
scRNA <- ScoreJackStraw(scRNA,dims = 1:20)
#鉴定前20个主成分的稳定性?
p2_2 <- JackStrawPlot(scRNA,dims = 1:20, reduction = "pca") +
  theme(legend.position="bottom") +
  labs(tag = "E")
p2_2
p2_3 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
p2_3
#结果显示可挑选前20个pc

p2_1| (p2_2 | p2_3) #中图

### 5、聚类，筛选marker基因，可视化----
#5.1 聚类
pc.num=1:20
#基于PCA数据
scRNA <- FindNeighbors(scRNA, dims = pc.num) 
# dims参数，需要指定哪些pc轴用于分析；这里利用上面的分析，选择20
scRNA <- FindClusters(scRNA, resolution = 0.5)
table(scRNA@meta.data$seurat_clusters)
#  0   1   2   3   4   5   6   7   8   9  10  11  12 
#  374 373 329 293 219 127 113  84  84  41  37  37  31 

scRNA = RunTSNE(scRNA, dims = pc.num)
DimPlot(scRNA, reduction = "tsne",label=T)
?RunTSNE
p3_1 <- DimPlot(scRNA, reduction = "tsne",label=T ) +
  labs(tag = "F")
p3_1

colnames(scRNA@meta.data)
DimPlot(scRNA, reduction = "tsne",label=T,split.by = "Patient_ID" ) 
table(scRNA@meta.data$Patient_ID,scRNA@meta.data$seurat_clusters)


#5.2 marker gene
#进行差异分析，一般使用标准化数据
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize")
#结果储存在"data"slot里
GetAssayData(scRNA,slot="data",assay="RNA")[1:8,1:4]
#if test.use is "negbinom", "poisson", or "DESeq2", slot will be set to "counts"


diff.wilcox = FindAllMarkers(scRNA)##默认使用wilcox方法挑选差异基因，大概4-5min
head(diff.wilcox)
dim(diff.wilcox)
library(tidyverse)
all.markers = diff.wilcox %>% select(gene, everything()) %>%
  subset(p_val<0.05 & abs(diff.wilcox$avg_logFC) > 0.5)
#An adjusted P value < 0.05  and | log 2 [fold change (FC)] | > 0.5 
#were considered the 2 cutoff criteria for identifying marker genes.
dim(all.markers)
summary(all.markers)
save(all.markers,file = "../../tmp/markergene.Rdata")
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(scRNA)) 
top10
length(top10)
length(unique(sort(top10)))

p3_2 <- DoHeatmap(scRNA, features = top10, group.by = "seurat_clusters") +
  labs(tag = "G")
p3_2
p3_1 | p3_2 #下图

### 6、拼图，比较----
p <- (p1_1 | p1_2 | p1_3 ) /
  ((p2_1| p2_2 | p2_3) /
     (p3_1 | p3_2))
ggsave("../../out/my_try.pdf", plot = p, width = 15, height = 18) 
save(scRNA,file = "../../tmp/scRNA.Rdata")

