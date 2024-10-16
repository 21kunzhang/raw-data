# 
library(tidyverse)
library(Seurat)
library(cols4all)
library(glmGamPoi) # SCT
library(cowplot)
library(clustree)
library(AUCell)
library(clusterProfiler) # 
library(ggsignif)
library(GSVA)
#1.####
#
dir.create("02_AUCell")
# 
mydata = readRDS("mydata.rds")
table(mydata@meta.data$cell_type)


Fibroblast_cells = subset(mydata, cell_type == "Fibroblast cells")

table(Fibroblast_cells@meta.data$cell_type)
saveRDS(Fibroblast_cells,"Fibroblast_cells.rds")

#Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(Fibroblast_cells@assays$RNA@data))

# load gene set
hallmark <- read.gmt("02_AUCell/h.all.v2023.2.Hs.symbols.gmt")
head(hallmark)

#
geneSets<-lapply(unique(hallmark$ont),function(x){hallmark$gene[hallmark$ont==x]})
names(geneSets) <- unique(hallmark$ont)
saveRDS(geneSets,file ='02_AUCell/geneSets.RDS' )

#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(geneSets, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)

# 
geneSet <- unique(hallmark$ont)
geneSet
aucs <- getAUC(cells_AUC)[geneSet, ]
aucs = t(aucs) %>%
  as.data.frame() %>%
  write.csv("02_AUCell/aucs.csv")


#2.####
# 
####hallmark####
aucs = read.csv("02_AUCell/aucs.csv", header = T,row.names = 1)
site <- Fibroblast_cells$site
aucs <- aucs[names(site),]
aucs$site <- site
head(aucs)
dim(aucs)
#1207   51
write.csv(aucs,file = "04_immu_metabolism/aucs_site.csv",row.names = T)
#aucs_go = read.csv("02_AUCell/go_aucs.csv", header = T,row.names = 1)

# 
tumor = aucs %>%
  subset(site == "Tumor")
dim(tumor)
#444  51
normal = aucs %>%
  subset(site == "Normal")
dim(normal)
#763  51
# 
sig_signaling <- data.frame() # 
for (i in colnames(tumor)[1:50]){
  result = wilcox.test(tumor[,i], normal[,i])
  sig_signaling[i,'logFC'] = log2(mean(tumor[,i]/mean(normal[,i])))
  sig_signaling[i,'P_Value'] = result$p.value
}

# write.csv(sig_signaling,"02_AUCell/sig_signaling.csv")

sig_signaling = read_delim("02_AUCell/sig_signaling.csv")
head(sig_signaling)
colnames(sig_signaling )[1] <- 'hallmark'
# 
p1 = ggplot(subset(sig_signaling, `P_Value` < 0.01), aes(reorder(hallmark, logFC), logFC, col = `P_Value`, size = abs(logFC))) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#667BC6", high = "#FFB4C2") +
  coord_flip() +
  labs(title = "AUCell score")
p1
ggsave("02_AUCell/sig_signaling.pdf",p1,height = 6,width = 6)
p2 = ggplot(subset(sig_signaling, `P_Value` < 0.01), 
            aes(x = reorder(hallmark, logFC), y = logFC, fill = `P_Value`)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_gradient(low = "#667BC6", high = "#FFB4C2") +
  coord_flip() +
  labs(title = "AUCell score")
p2
ggsave("02_AUCell/sig_signaling_AUCell_score.pdf",p2,height = 6,width = 6)
#3.####
#
signaling_rows <- grepl("signaling", sig_signaling$hallmark, ignore.case = TRUE)

# 
filtered_signaling<- sig_signaling[signaling_rows, ]
View(filtered_signaling)
# HALLMARK_TGF_BETA_SIGNALING
TGF_BETA = aucs[,c("site","HALLMARK_TGF_BETA_SIGNALING")]
dim(TGF_BETA)
TGF_BETA$TGF_BETA_type <- ifelse(TGF_BETA$HALLMARK_TGF_BETA_SIGNALING>median(TGF_BETA$HALLMARK_TGF_BETA_SIGNALING),'High','Low')

Fibroblast_cells[['TGF_BETA_type']] <- TGF_BETA[rownames(Fibroblast_cells@meta.data),'TGF_BETA_type']
Fibroblast_cells[['TGF_score']] <- TGF_BETA[rownames(Fibroblast_cells@meta.data),'HALLMARK_TGF_BETA_SIGNALING']

saveRDS(Fibroblast_cells,file = "02_AUCell/Fibroblast_cells.RDS")
TGF <- aucs[,c("site","HALLMARK_TGF_BETA_SIGNALING" )]
p3 = ggplot(TGF, aes(site, HALLMARK_TGF_BETA_SIGNALING, fill = site)) +
  geom_violin(alpha = 0.5, width = 0.5) +
  #geom_jitter(aes(color = site), width = 0.1) +
  geom_boxplot(fill="white", width=0.1, aes(color = site))+
  scale_fill_manual(values = c("#667BC6","#FFB4C2"))+
  scale_color_manual(values = c("#667BC6","#FFB4C2"))+
  #theme_bw()+
  #theme(axis.text=element_text(face="bold", size=13))+
  geom_signif(comparisons=list(c("Normal", "Tumor")), 
              textsize=5, 
              test=wilcox.test, 
              step_increase=0.2,
              map_signif_level=T) +
  labs(x = "site", y = "AUCell Score", title = "HALLMARK_TGF_BETA_SIGNALING")
p3
ggsave("02_AUCell/TGF_signaling_tumor_normal.pdf",p3,height = 4,width = 4)
#4.####
# 
table(Fibroblast_cells@meta.data$site)
DEG = FindMarkers(Fibroblast_cells, group.by="site", 
                  ident.1="Tumor", ident.2="Normal", 
                  logfc.threshold=0, min.pct=0)
dim(DEG)
head(DEG)
write.csv(DEG,file = '02_AUCell/deg.csv',quote = F,row.names =T)
DEG = read.csv("02_AUCell/deg.csv", header = T, row.names = 1)
deg = subset(DEG, p_val_adj < 0.05 & avg_logFC > 0.25)
dim(deg)
#106   5
# 
 write.csv(deg, "02_AUCell/deg_005_025.csv",row.names = T)

# HALLMARK_TGF_BETA_SIGNALING
inter_gene2 = subset(hallmark, ont == "HALLMARK_TGF_BETA_SIGNALING")$gene
inter_gene2
dput(inter_gene2)
length(inter_gene2)
inter_gene2 <- inter_gene2[c(19,22,32,34,42,54)]#

# 
Idents(Fibroblast_cells) = Fibroblast_cells@meta.data$site
p4 = VlnPlot(Fibroblast_cells, features=inter_gene2,ncol = 3, 
             cols=c("#FFB4C2", "#667BC6"), 
             pt.size=0)+
  NoLegend()
p4
ggsave("02_AUCell/TGF_signaling_gene.pdf",p4,height = 8,width = 6)
