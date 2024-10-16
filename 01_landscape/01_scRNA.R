#01####
# 
library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)
library(cols4all)
library(data.table)
library(stringr)
library(ComplexHeatmap)
library(scales)
library(ggsci)
#
dir.create("01_landscape")
## 
### 
assays <- dir("00_origin_datas/GEO/GSE208653_RAW/")
dir <- paste0("00_origin_datas/GEO/GSE208653_RAW/", assays)
assays
dir
# 
samples_name = assays
samples_name
# 
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}

### 
names(scRNAlist) <- samples_name

# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])


#2.####
mycolor <- c4a('poly.alphabet2',26)#c4a("brewer.paired", 12)
show_col(mycolor)
colors <- mycolor
#p = VlnPlot(scRNA, features=c("nFeature_RNA","nCount_RNA",'percent.mito'), pt.size=0.0, cols=colors)
VlnPlot_nFeature_RNA_before = VlnPlot(scRNA,features=c("nCount_RNA",'nFeature_RNA','percent.mito'),pt.size = 0)
VlnPlot_nFeature_RNA_before
ggsave("01_landscape/VlnPlot_nFeature_RNA_before.pdf", VlnPlot_nFeature_RNA_before, 
       width=12, height=5)


### 
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<8000&percent.mito<10)
VlnPlot_nFeature_RNA_after = VlnPlot(scRNA,features=c("nCount_RNA",'nFeature_RNA','percent.mito'),pt.size = 0)
VlnPlot_nFeature_RNA_after
ggsave("01_landscape/VlnPlot_nFeature_RNA_after.pdf", VlnPlot_nFeature_RNA_after, 
       width=12, height=5)
# NormalizeData, FidnVariableFeatures, ScaleData 
scRNA = SCTransform(scRNA, vars.to.regress = "percent.mito", verbose = F)
#
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
colnames(scRNA@meta.data)
########
library(harmony)
scRNA = RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony=50, lambda=0.5,assay.use = "SCT")

# scRNA ##########
save(scRNA,file = 'scRNA.Rdata')
# load('scRNA.Rdata')

#3.####
ElbowPlot(scRNA,ndims = 50)
scRNA <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims = 1:30,reduction="harmony")
p1 = DimPlot(scRNA, reduction="umap", group.by="orig.ident", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("01_landscape/UMAP_Sample.pdf", p1, height=6, width=6)
###3.1################################ 
#c4a_gui()
#colors <- sample(c4a('poly.sky24',24),size = 20,replace = F)
mydata <- FindClusters(scRNA, resolution=0.1)
p2 <- UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
ggsave("01_landscape/cluster_umap.pdf",p2,height = 10,width = 10)
markers <- FindAllMarkers(mydata,group_by='seurat_clusters',only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "01_landscape/All_cluster_type.txt", col.names=T, row.names=T, quote=F, sep="\t")
table(mydata$seurat_clusters)
# 0     1     2     3     4     5     6     7     8     9    10    11    12    13 
# 13445  5205  3449  3412  3350  2368  1374  1207  1054   720   680   623   537    81  
###3.2 ##################

VlnPlot(mydata, features=c("CPA3"), pt.size=0, cols=colors,group.by ='seurat_clusters')+NoLegend()+theme(axis.title.x=element_blank())
#0: NK/T cells: CD3D	CD3E	NKG7		
#1: Neutrophil cells: FCGR3B	CSF3R	CXCR2
#2: Epithelial cells : KRT8 KRT19 KRT18 CDH1
#3: Epithelial cells : EPCAM,KRT8 KRT19 KRT18 CDH1
#4: Macrophage cells	:CD68	CD163	MRC1
#5: Epithelial cells: KRT8 KRT19 KRT18 CDH1
#6: Plasma cells: CD79A	TNFRSF17	
#7: Fibroblast cells: COL1A1	DCN	COL1A2	ACTA2
#8: Plasma cells: CD79A	TNFRSF17	
#9: Mast cells:TPSAB1	CPA3	TPSB2	MS4A2
#10 Epithelial cells：EPCAM,KRT19 KRT18 CDH1
#11 B cells：MS4A1	CD79A	CD79B	CD40
#12 Endothelial cells ：VWF	CDH5	CLDN5
#13 Epithelial cells：EPCAM,KRT19 KRT18 CDH1
#mydata = subset(mydata, seurat_clusters %in% c(0,1,2,3,4,6,7,8,9))
#mydata@meta.data$seurat_clusters = droplevels(mydata@meta.data$seurat_clusters, exclude=c(5, 8, 9, 11, 12, 13))
### 3.3################
cell_label = c("NK/T cells", "Neutrophil cells", "Epithelial cells", "Epithelial cells", "Macrophage cells",
               "Epithelial cells", "Plasma cells", "Fibroblast cells","Plasma cells",'Mast cells',
               'Epithelial cells','B cells','Endothelial cells','Epithelial cells')
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
# c4a_gui()
colors = c4a('carto.pastel',11)

p2 = UMAPPlot(mydata, pt.size=1, label=T, label.size=9)+NoLegend()+scale_color_manual(values=colors)
library(ggpubr)
p2=ggarrange(AugmentPlot(p2,dpi = 300,width = 13,height = 13), ncol = 1, nrow = 1,
             legend ='right',legend.grob=get_legend(p2),
             font.label = list(size = 15, face = "bold",family ='Times')) # 

ggsave("01_landscape/UMAP_cell_type2.pdf", p2, width=8, height=8)
genes = c('CD3D',	'CD3E','NKG7',
          'FCGR3B',	'CSF3R',	'CXCR2',
          'KRT8', 'KRT19', 'KRT18', 'CDH1',
          'EPCAM',
          'CD68',	'CD163',	'MRC1',
          'CD79A',	'TNFRSF17',
          'COL1A1',	'DCN',	'COL1A2',	'ACTA2',
          'TPSAB1',	'CPA3',	'TPSB2',	'MS4A2',
          'MS4A1',	'CD79B',	#'CD40',
          'VWF',	'CDH5',	'CLDN5')



p3 = DotPlot(mydata, features=genes)+coord_flip()+scale_color_gradientn(colors=c("black", "dodgerblue", "white", "orange", "firebrick1"))+theme_minimal()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("01_landscape/dotplot_gene_marker.pdf",p3,height = 8,width = 10)


###marker####
genes2 <-c('CD3D',	#'CD3E','NKG7',
           'FCGR3B',#	'CSF3R',	'CXCR2',
           'KRT8', #'KRT19', 'KRT18', 'CDH1',
           'EPCAM',
           'CD68',	#'CD163',	'MRC1',
           'CD79A',#	'TNFRSF17',
           'COL1A1',#	'DCN',	'COL1A2',	'ACTA2',
           'TPSAB1',#	'CPA3',	'TPSB2',	'MS4A2',
           'MS4A1',	#'CD79B',	#'CD40',
           'VWF')	#'CDH5',	'CLDN5')

p3_2 <- VlnPlot(mydata, features=genes2, pt.size=0, cols=colors,group.by ='cell_type',ncol = 5)+NoLegend()+theme(axis.title.x=element_blank())
ggsave("01_landscape/volin_marker.pdf", p3_2, width=12, height=8)
##
save(mydata,file = 'mydata.Rdata')
# 4#######

###4.2##########################################################  
mydata[['site']]=ifelse(substring(mydata$orig.ident,1,2)=='CA','Tumor','Normal')
Type_label = c("Normal", "Tumor")
bar$site = factor(bar$site, levels=Type_label)
bar = bar %>% group_by(site) %>% mutate(percent=100*Freq/sum(Freq))

p4 <- ggplot(data=bar, aes(x=cell_type, y=percent, fill=site,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c('#4793AF','#DD5746'))+theme_classic()+
  ggtitle("Percent(%)")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave("01_landscape/barplot_pair_number.pdf", p4, width=8, height=4)


