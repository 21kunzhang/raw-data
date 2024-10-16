
library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
library(cols4all)
##
#####
dir.create("03_fibroblast_Pseudotime")

# 
#BiocManager::install("monocle")
# 
library(monocle) #


Fibroblast_cells = readRDS("02_AUCell/Fibroblast_cells.rds")
#Fibroblast_cells <- subset(Fibroblast_cells,site=="Tumor")
table(Fibroblast_cells$site)
# 1.

# 
expr_matrix = as(as.matrix(Fibroblast_cells@assays$RNA@counts), 'sparseMatrix')
# 
p_data = Fibroblast_cells@meta.data
p_data$cell_type = Fibroblast_cells@active.ident # 
f_data = data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
## （cell number）
dim(expr_matrix)
dim(p_data)
dim(f_data)

# 2.####
pd = new('AnnotatedDataFrame', data = p_data)
fd = new('AnnotatedDataFrame', data = f_data)
# 
cds = newCellDataSet(expr_matrix,
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

# 
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

# 
cds = detectGenes(cds, min_expr = 0.1) #
num_cells_expressed
print(head(fData(cds)))
expressed_genes = row.names(subset(fData(cds),
                                   num_cells_expressed >= 10))
length(expressed_genes)
# 
# 

# 3.1####
# diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~site", cores = 1)
# head(diff)
# # 
# deg = subset(diff, qval < 0.01)
# deg = deg[order(deg$qval, decreasing = F),]
# 
# # 
#write.csv(deg, "03_fibroblast_Pseudotime/deg.csv")
# 
# # 
# ordergene = rownames(deg)
# cds = setOrderingFilter(cds, ordergene)
# # [["use_for_ordering"]]，table(cds@featureData@data[["use_for_ordering"]])
# pdf("03_fibroblast_Pseudotime/train.ordergenes.pdf")
# plot_ordering_genes(cds)
# dev.off()
# # 
# cds = reduceDimension(cds, max_components = 2, method = "DDRTree")
# 
# # 
# cds = orderCells(cds) # 
# 
# # 
# ## 
# plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = T)
# plot_cell_trajectory(cds, color_by = "site", size = 1, show_backbone = T)



#3.2FindMarkers Tumor,NORMAL ####
# 
# deg.cluster = FindMarkers(Fibroblast_cells, group.by="Group", ident.1="Tumor", ident.2="Normal")
# head(deg.cluster)
deg.cluster <- read.csv("02_AUCell/deg.csv",header = T,row.names = 1)
express_genes = rownames(subset(deg.cluster, p_val_adj < 0.05, abs(avg_logFC) > 0))
head(express_genes)
cds = setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)
# 
cds = reduceDimension(cds, max_components = 2, method = "DDRTrees")
# 4####
cds = orderCells(cds) # 
saveRDS(cds,file = "03_fibroblast_Pseudotime/cds.RDS")
p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime", size = 2, show_backbone = T)
p1
ggsave("03_fibroblast_Pseudotime/cell_trajectory_Pseudotime.pdf", p1, width = 6, height = 3)
p2 <- plot_cell_trajectory(cds, color_by = "site", size = 2, show_backbone = T)
ggsave("03_fibroblast_Pseudotime/cell_trajectory_site.pdf", p2, width = 6, height = 3)
library(ggpubr)
df <- pData(cds) 
## pData(cds)
View(df)
p2.1 <- ggplot(df, aes(Pseudotime, colour = site, fill=site)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
p2.1 
ggsave("03_fibroblast_Pseudotime/cell_trajectory_site_density.pdf", p2.1, width = 6, height = 3)

colors = c4a("seaborn.bright",9) 
p3 = plot_cell_trajectory(cds, color_by = "State", size = 2, show_backbone = T) +
  scale_color_manual(values = colors) 
  
p3
ggsave("03_fibroblast_Pseudotime/cell_trajectory_state.pdf", p3, width = 6, height = 3)

p3.1 = plot_cell_trajectory(cds, color_by = "TGF_BETA_type", size = 2, show_backbone = T) +
  scale_color_manual(values = colors) 

p3.1
ggsave("03_fibroblast_Pseudotime/cell_trajectory_TGF_BETA_type.pdf", p3.1, width = 6, height = 3)


# 
Time_diff = differentialGeneTest(cds[ordergene,], cores = 1,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
head(Time_diff)

# 5.####
geneSets <- readRDS("02_AUCell/geneSets.RDS") 
hh_genes <- unlist(geneSets['HALLMARK_TGF_BETA_SIGNALING'])
inter_genes = intersect(Time_diff$gene_short_name, hh_genes)
head(inter_genes)
inter_genes



Time_diff_hh = subset(Time_diff, gene_short_name %in% inter_genes)
dim(Time_diff_hh)
Time_genes = Time_diff_hh %>% pull(gene_short_name) %>%as.character()

p4 = plot_pseudotime_heatmap(cds[Time_genes,], 
                             num_clusters = 3, 
                             show_rownames = T, 
                             return_heatmap = T)
p4
ggsave("03_fibroblast_Pseudotime/pseudotime_heatmap_TGF.pdf", p4, width = 6, height = 6)
####fibroblast_proliferation####
dput(Time_genes )
neg_apoptotic <- c( 'DDX3X', 'HSPA5',  'HSPA1B', 'HSPA1A' )#'DNAJA1', 'GREM1', 'ARF4', 'CD74','IFI6', 'DPEP1', 'IGF1','MCL1'
cell_proliferation <- c('NFKBIA', 'CXCL10', 'CXCL11', 'CXCL9', 'JUN', 'JUND', 'STAT1', 'JUNB')
cell_cycle <- c( 'JUND', 'IRF1', 'CCNL1', 'JUNB', 'ACTB')#'HSPA8', 'JUN',
fibroblast_proliferation <- c('CD74',  'FN1')#,'S100A6','IGF1', 'ESR1')

p7 <- plot_genes_in_pseudotime(cds[fibroblast_proliferation,], color_by = "Pseudotime")
ggsave("03_fibroblast_Pseudotime/fibroblast_proliferation_genes_in_pseudotime.pdf", p7, width = 6, height = 6)

p8 <- plot_genes_in_pseudotime(cds[cell_proliferation [c(2,4,7)],], color_by = "Pseudotime")
p8
ggsave("03_fibroblast_Pseudotime/cell_proliferation_genes_in_pseudotime.pdf", p8, width = 6, height = 9)
p9 <- plot_genes_in_pseudotime(cds[cell_cycle[2],], color_by = "Pseudotime")
p9
ggsave("03_fibroblast_Pseudotime/cell_cycle_genes_in_pseudotime.pdf", p9, width = 6, height = 3)

# #######
# p8<-ggplot(data = df, mapping = aes(x = Pseudotime, y = TGF_score)) +
#   geom_point(alpha = 0.7, aes(color=site)) +
#   scale_color_brewer(palette="Paired")  +
#   labs(
#     x = "Pseudotime",
#     y = "TGF_score",color = "site")+ 
#   geom_smooth(method = "lm",se=T,color='#FEA82F') +
#   stat_cor(method = "pearson",digits = 3,size=6)+
#   theme_bw()


####TGF_BETA_type####
p10<- VlnPlot(Fibroblast_cells, features=neg_apoptotic[c(3,4)] , pt.size=0, cols=colors,group.by ='TGF_BETA_type')+NoLegend()+theme(axis.title.x=element_blank())
ggsave("03_fibroblast_Pseudotime/neg_apoptotic_genes_in_TGF_BETA_type.pdf", p10, width = 6, height = 3)
####TGF_BETA_type分组中cell_proliferation####
p11<- VlnPlot(Fibroblast_cells, features=cell_proliferation[c(5,6,8)] , pt.size=0, cols=colors,group.by ='TGF_BETA_type')+NoLegend()+theme(axis.title.x=element_blank())
ggsave("03_fibroblast_Pseudotime/cell_proliferation_genes_in_TGF_BETA_type.pdf", p11, width = 8, height = 3)
####TGF_BETA_type####
p12<- VlnPlot(Fibroblast_cells, features=fibroblast_proliferation , pt.size=0, cols=colors,group.by ='TGF_BETA_type')+NoLegend()+theme(axis.title.x=element_blank())

####TGF_BETA_type####
p13<- VlnPlot(Fibroblast_cells, features=cell_cycle[c(1,3,4)] , pt.size=0, cols=colors,group.by ='TGF_BETA_type')+NoLegend()+theme(axis.title.x=element_blank())
ggsave("03_fibroblast_Pseudotime/cell_cycle_genes_in_TGF_BETA_type.pdf", p13, width = 8, height = 3)





