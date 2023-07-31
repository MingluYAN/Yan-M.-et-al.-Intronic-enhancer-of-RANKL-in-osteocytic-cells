# Yan-M.-et-al.-Intronic-enhancer-of-RANKL-in-osteocytic-cells
# Analysis of scRNA-seq data GSE154719 Jialiang S. Wang et al. “Control of osteocyte dendrite formation by Sp7 and its target gene osteocrin”. https://doi.org/10.1038/s41467-021-26571-7
# load data Dmp1-Cre; Sp7+/+; tdTm+
Con_Td <- Read10X(data.dir = "/Users/yanminglu/Desktop/")

# QC
Con_Td <- CreateSeuratObject(counts = Con_Td, project = "Con_Td", min.cells = 3, min.features = 200)
Con_Td <- PercentageFeatureSet(object = Con_Td, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(Con_Td, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Con_Td <- subset(Con_Td, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
# analysis
Con_Td@assays$RNA 
Con_Td <- NormalizeData(Con_Td, normalization.method = "LogNormalize", scale.factor = 10000)
Con_Td <- FindVariableFeatures(Con_Td, selection.method = "vst", nfeatures = 2000)
top15 <- head(VariableFeatures(Con_Td), 15)
plot1 <- VariableFeaturePlot(Con_Td)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
plot1 + plot2
all.genes <- rownames(Con_Td)
Con_Td <- ScaleData(Con_Td, features = all.genes)
Con_Td <- RunPCA(Con_Td, features = VariableFeatures(object = Con_Td))

Con_Td <- FindNeighbors(Con_Td, dims = 1:20)
Con_Td <- FindClusters(Con_Td, resolution = 1)
Con_Td <- RunUMAP(Con_Td, dims = 1:20)
DimPlot(Con_Td, reduction = "umap", label = T)
markers <- FindAllMarkers(Con_Td, min.pct = 0.5)
write.csv(markers, "markers1.csv")
# remove contaminated hematopoietic cells and endothelial cells 
Con_Td <- subset(Con_Td, idents = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","17","19","20","0"))
Con_Td <- RunUMAP(object = Con_Td, dims = 1:20, verbose = F)
Con_Td <- FindNeighbors(Con_Td, dims = 1:20)
Con_Td <- FindClusters(Con_Td, resolution = 0.5)
# remove red blood cells
Con_Td <- subset(Con_Td, idents = c("1","2","3","4","5","6","7","8","0"))

# Figure 2a 
pdf("Fig.2a_control_umap.pdf", useDingbats = F, height = 3, width = 4)
DimPlot(object = Con_Td, label = T, reduction = "umap", label.size = 4) + theme(aspect.ratio = 0.8) + NoLegend()
dev.off()

pdf("Fig.2a_control_dotplot.pdf", useDingbats = F, height = 4.5, width = 5)
DotPlot(Con_Td, features = c("Sost","Ackr3", "Fbln7", "Dmp1","Irx5", "Dkk1"), dot.scale = 10, cols = c("lightgrey","red"), col.min = 0, col.max = 3)
dev.off()

# load data Dmp1-Cre; Sp7flox/flox; tdTm+
Sp7_OcyKO <- Read10X(data.dir = "/Users/yanminglu/Desktop/")

#QC
Sp7_OcyKO <- CreateSeuratObject(counts = Sp7_OcyKO, project = "Sp7_OcyKO", min.cells = 3, min.features = 200) 
Sp7_OcyKO <- PercentageFeatureSet(object = Sp7_OcyKO, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(Sp7_OcyKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Sp7_OcyKO <- subset(Sp7_OcyKO, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# analysis
Sp7_OcyKO <- NormalizeData(Sp7_OcyKO, normalization.method = "LogNormalize", scale.factor = 10000) 
Sp7_OcyKO <- FindVariableFeatures(Sp7_OcyKO, selection.method = "vst", nfeatures = 2000)
top15 <- head(VariableFeatures(Sp7_OcyKO), 15)
plot1 <- VariableFeaturePlot(Sp7_OcyKO)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
plot1 + plot2
all.genes <- rownames(Sp7_OcyKO)
Sp7_OcyKO <- ScaleData(Sp7_OcyKO, features = all.genes)
Sp7_OcyKO <- RunPCA(Sp7_OcyKO, features = VariableFeatures(object = Sp7_OcyKO))
ElbowPlot(Sp7_OcyKO)

Sp7_OcyKO <- FindNeighbors(Sp7_OcyKO, dims = 1:20)
Sp7_OcyKO <- FindClusters(Sp7_OcyKO, resolution = 1)
Sp7_OcyKO <- RunUMAP(Sp7_OcyKO, dims = 1:20)

DimPlot(Sp7_OcyKO, reduction = "umap", pt.size = 1, label = T)
markers2 <- FindAllMarkers(Sp7_OcyKO, min.pct = 0.5)
write.csv(markers2, "markers2.csv")
# remove contaminated hematopoietic cells and red blood cells
Sp7_OcyKO <- subset(Sp7_OcyKO, idents = c("3","4","5","6","7","10","8","0","1","2"))

#Figure 2b 
pdf("Fig.2b_Sp7_OcyKO_umap.pdf", useDingbats = F, height = 3, width = 4)
DimPlot(object = Sp7_OcyKO, label = T, reduction = "umap", label.size = 0) + theme(aspect.ratio = 0.8) + NoLegend()
dev.off()

pdf("Fig.2b_ dotplot.pdf", useDingbats = F, height = 4.5, width = 5)
DotPlot(Sp7_OcyKO, features = c("Sost","Ackr3", "Fbln7", "Dmp1","Irx5", "Dkk1"), dot.scale = 10, cols = c("lightgrey","red"), col.min = 0, col.max = 3)
dev.off()

# Seurat data intergration
Con_Td $dataset <- "Con_Td"
Sp7_OcyKO $dataset <- "Sp7_OcyKO"
InteOcy <- merge(x= Con_Td, y = Sp7_OcyKO)
InteOcy.list <- SplitObject(object = InteOcy, split.by = "dataset")
InteOcy.list <- lapply(X = InteOcy.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = InteOcy.list)
InteOcy.anchors <- FindIntegrationAnchors(object.list = InteOcy.list, anchor.features = features)
InteOcy.integration <- IntegrateData(anchorset = InteOcy.anchors)

InteOcy.integration@assays #26769 cells
DefaultAssay(InteOcy.integration) <- "integrated"

InteOcy.integration <- ScaleData(InteOcy.integration, verbose = FALSE)
InteOcy.integration <- RunPCA(InteOcy.integration, npcs = 20, verbose = FALSE)
InteOcy.integration <- RunUMAP(InteOcy.integration, reduction = "pca", dims = 1:20)
InteOcy.integration <- FindNeighbors(InteOcy.integration, reduction = "pca", dims = 1:20)
InteOcy.integration <- FindClusters(InteOcy.integration, resolution = 0.5)
# Figure 2c
p1 <- DimPlot(InteOcy.integration, reduction = "umap", group.by = "dataset")
p2 <- DimPlot(InteOcy.integration, reduction = "umap", label = TRUE, repel = TRUE)
p1+p2

DefaultAssay(InteOcy.integration) <- "RNA"

pdf("Figure 2c_Sost.pdf", useDingbats = F, height = 3, width = 4)
FeaturePlot(InteOcy.integration,"Sost",pt.size = 0.5, min.cutoff = 0, max.cutoff = 1)
dev.off()

pdf("Figure 2c_Ackr3.pdf", useDingbats = F, height = 3, width = 4)
FeaturePlot(InteOcy.integration,"Ackr3",pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.8)
dev.off()

pdf("Figure 2c_Fbln7.pdf", useDingbats = F, height = 3, width = 4)
FeaturePlot(InteOcy.integration,"Fbln7",pt.size = 0.5, min.cutoff = 0, max.cutoff = 1.5)
dev.off()

pdf("Figure 2c_Dmp1.pdf", useDingbats = F, height = 3, width = 4)
FeaturePlot(InteOcy.integration,"Dmp1",pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

pdf("Figure 2c_Irx5.pdf", useDingbats = F, height = 3, width = 4)
FeaturePlot(InteOcy.integration,"Irx5",pt.size = 0.5, min.cutoff = 0, max.cutoff = 1.5)
dev.off()

pdf("Figure 2c_Dkk1.pdf", useDingbats = F, height = 3, width = 4)
FeaturePlot(InteOcy.integration,"Dkk1",pt.size = 0.5, min.cutoff = 0, max.cutoff = 1.5)
dev.off()

# Comparative analysis of osteocytic cells derived from Con_Td and Sp7_OcyKO
Idents(InteOcy.integration) 
InteOcy.integration$dataset 
InteOcy.integration$cluster.dataset<-paste(Idents(InteOcy.integration), InteOcy.integration$dataset, sep = "-")
InteOcy.integration$cluster.dataset
Idents(InteOcy.integration) <- "cluster.dataset"
Ocy.compare.marker <- FindMarkers(InteOcy.integration, ident.1 = "9-Con_Td", ident.2 = "9-Sp7_OcyKO", min.pct = 0,verbose = FALSE)
write.csv(Ocy.compare.marker, "Ocy.compare.csv") # Metascape Figure 2d,e

# Figure 2f

pdf("Tnfsf11.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Tnfsf11"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,2))
dev.off()

pdf("Cebpb.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Cebpb"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,1.5))
dev.off()

pdf("Cebpa.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Cebpa"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,1))
dev.off()

pdf("Trp53.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Trp53"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,2))

pdf("Sost.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Sost"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,4.5))
dev.off()

pdf("Ackr3.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Ackr3"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,3))
dev.off()

pdf("Fbln7.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Fbln7"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,3))
dev.off()

pdf("Dmp1.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Dmp1"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,5))
dev.off()

pdf("Irx5.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Irx5"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,4))
dev.off()

pdf("Dkk1.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Dkk1"), cols =c("pink", "grey"), pt.size = 0, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,5))
dev.off()

# supplementary figure
my <- c("0","1","7","10","4","3","2","5","6","8","9")
factor(Idents(InteOcy.integration), levels= my)
Idents(InteOcy.integration) <- factor(Idents(InteOcy.integration), levels= my)
# Supplementary figure 1A
pdf("Supplementary Figure 1A.pdf", useDingbats = F, height = 5, width = 7)
DotPlot(InteOcy.integration, features = c("Sost","Fbln7","Dkk1","Dmp1","Pdgfrl","Bglap3","Bglap","Bglap2","Lpl","Kcnk2","Cxcl12","Igfbp5"), dot.scale = 10) + RotatedAxis()
dev.off()

# Supplementary figure 1B
pdf("Supplementary Figure 1B.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features = c("Sp7"), cols =c("pink", "grey"), pt.size = 1, idents = c("9"), split.by = "dataset") + scale_y_continuous(limits = c(0,1.5))
dev.off()

# Supplementary figure 1C

pdf("Supplementary figure 1C_Tnfsf11_c9.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features = c("Tnfsf11"), cols =c("pink", "grey"), pt.size = 1, idents = c("9"), split.by = "dataset") + scale_y_continuous(limits = c(0,2))
dev.off()

pdf("Supplementary figure 1C_Tnfsf11_c4.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features = c("Tnfsf11"), cols =c("pink", "grey"), pt.size = 1, idents = c("4"), split.by = "dataset") + scale_y_continuous(limits = c(0,2))
dev.off()

pdf("Supplementary figure 1C_Tnfrsf11b_c9.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Tnfrsf11b"), cols =c("pink", "grey"), pt.size = 1, idents = "9", split.by = "dataset") + scale_y_continuous(limits = c(0,2))
dev.off()

pdf("Supplementary figure 1C_Tnfrsf11b_c4.pdf", useDingbats = F, height = 3, width = 3)
VlnPlot(InteOcy.integration, features =c("Tnfrsf11b"), cols =c("pink", "grey"), pt.size = 1, idents = "4", split.by = "dataset") + scale_y_continuous(limits = c(0,2))
dev.off()

# Supplementary figure 1D
pdf("Supplementary figure 1D.pdf", useDingbats = F, height = 7, width = 3.5)
DotPlot(InteOcy.integration, features = c("Tnfsf11","Tnfrsf11b"))
dev.off()
