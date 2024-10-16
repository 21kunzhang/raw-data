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
library(dplyr)


dir.create('04_immu_metabolism')
# HALLMARK_TGF_BETA_SIGNALING
aucs = read.csv("04_immu_metabolism/aucs_site_.csv", header = T,row.names = 1)
TGF_BETA = aucs[,c("site","HALLMARK_TGF_BETA_SIGNALING")]
dim(TGF_BETA)
#1.####
# TGF_BETA$TGF_BETA_type <- ifelse(TGF_BETA$HALLMARK_TGF_BETA_SIGNALING>median(TGF_BETA$HALLMARK_TGF_BETA_SIGNALING),'High','Low')
# 
# Fibroblast_cells[['TGF_BETA_type']] <- TGF_BETA[rownames(Fibroblast_cells@meta.data),'TGF_BETA_type']
DEG.TGF_BETA_type = FindMarkers(Fibroblast_cells, group.by="TGF_BETA_type", 
                              ident.1="High", ident.2="Low", 
                              logfc.threshold=0, min.pct=0)
dim(DEG.TGF_BETA_type)
head(DEG.TGF_BETA_type)
write.csv(DEG.TGF_BETA_type,file = '04_immu_metabolism/DEG.TGF_BETA_type.csv',quote = F,row.names =T)
DEG.TGF_BETA_type = read.csv("04_immu_metabolism/DEG.TGF_BETA_type.csv", header = T, row.names = 1)
deg.TGF_BETA_type = subset(DEG.TGF_BETA_type, p_val_adj < 0.05 & avg_logFC > 0.25)
dim(deg.TGF_BETA_type)
# 
write.csv(deg.TGF_BETA_type, "04_immu_metabolism/deg.TGF_BETA_type_005_025.csv",row.names = T)
signaling.gene.enrichment=mg_clusterProfiler(genes =rownames(DEG.TGF_BETA_type))
fig4a=list()
fig4a[[1]]=enrichplot::dotplot(signaling.gene.enrichment$KEGG)+ggtitle('KEGG')
fig4a[[2]]=enrichplot::dotplot(signaling.gene.enrichment$GO_BP)+ggtitle('Biological Process')
# fig4a[[3]]=enrichplot::dotplot(signaling.gene.enrichment$GO_CC)+ggtitle('Cellular Component')
# fig4a[[4]]=enrichplot::dotplot(signaling.gene.enrichment$GO_MF)+ggtitle('Molecular Function')
fig4=mg_merge_plot(fig4a[[1]],fig4a[[2]],labels = c('A','B'),
                   #mg_merge_plot(fig4a[[3]],fig4a[[4]],ncol=2,nrow=1,labels = LETTERS[3:4],heights = c(1,1),widths = c(1,2)),
                   ncol=2,widths  = c(1,1))

#########
# 
aucs = read.csv("02_AUCell/aucs.csv", header = T, row.names = 1)
colnames(aucs)
colnames(aucs) <- gsub('HALLMARK_',"",colnames(aucs))
dput(colnames(aucs))
singal <- c("ADIPOGENESIS", "ALLOGRAFT_REJECTION", "ANDROGEN_RESPONSE", 
  "ANGIOGENESIS", "APICAL_JUNCTION", "APICAL_SURFACE", "APOPTOSIS", 
  "BILE_ACID_METABOLISM", "CHOLESTEROL_HOMEOSTASIS", "COAGULATION", 
  "COMPLEMENT", "DNA_REPAIR", "E2F_TARGETS", "EPITHELIAL_MESENCHYMAL_TRANSITION", 
  "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE", "FATTY_ACID_METABOLISM", 
  "G2M_CHECKPOINT", "GLYCOLYSIS", "HEDGEHOG_SIGNALING", "HEME_METABOLISM", 
  "HYPOXIA", "IL2_STAT5_SIGNALING", "IL6_JAK_STAT3_SIGNALING", 
  "INFLAMMATORY_RESPONSE", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", 
  "KRAS_SIGNALING_DN", "KRAS_SIGNALING_UP", "._SPINDLE", 
  "MTORC1_SIGNALING", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "MYOGENESIS", 
  "NOTCH_SIGNALING", "OXIDATIVE_PHOSPHORYLATION", "P53_PATHWAY", 
  "PANCREAS_BETA_CELLS", "PEROXISOME", "PI3K_AKT_MTOR_SIGNALING", 
  "PROTEIN_SECRETION", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "SPERMATOGENESIS", 
  "TGF_BETA_SIGNALING", "TNFA_SIGNALING_VIA_NFKB", "UNFOLDED_PROTEIN_RESPONSE", 
  "UV_RESPONSE_DN", "UV_RESPONSE_UP", "WNT_BETA_CATENIN_SIGNALING", 
  "XENOBIOTIC_METABOLISM")
# 
## tgf
Hh_group = c("TGF_BETA_SIGNALING")

## Inflammatory Group

infla_group = toupper(c(
  "TNFA_SIGNALING_VIA_NFKB", "IL6_JAK_STAT3_signaling", "Interferon_alpha_response",
  "Interferon_gamma_response","Inflammatory_response", "IL2_STAT5_signaling")
)

# Proliferation Group
prolif_group = toupper(c(
  "Mitotic_spindle","G2M_checkpoint",
  "PI3K_AKT_mTOR_signaling", "mTORC1_signaling", "E2F_targets",
  "Kras_signaling_up","p53_pathway","Kras_signaling_dn")
)

# Metabolism Group
meta_group = toupper(c(
  "Oxidative_phosphorylation", "Glycolysis", "Adipogenesis",
  "Fatty_acid_metabolism","Bile_acid_metabolism"
))

cor = data.frame()
for(i in infla_group){
  result = cor.test(aucs[,Hh_group], aucs[,c(i)])
  cor[i,"Inflammatory cor"] = result$estimate
  cor[i,"Inflammatory Pvalue"] = result$p.value
}

for(i in prolif_group){
  result = cor.test(aucs[,Hh_group], aucs[,c(i)])
  cor[i,"Proliferation cor"] = result$estimate
  cor[i,"Proliferation Pvalue"] = result$p.value
}

for(i in meta_group){
  result = cor.test(aucs[,Hh_group], aucs[,c(i)])
  cor[i,"Meatbolism cor"] = result$estimate
  cor[i,"Meatbolism Pvalue"] = result$p.value
}
cor

cor$Group = rownames(cor)
cor
saveRDS(cor,file = "04_immu_metabolism/cor.RDS")
p2 = ggplot(subset(cor[infla_group,],`Inflammatory Pvalue`<0.05) ,aes(reorder(Group, `Inflammatory cor`), 
                                   `Inflammatory cor`, 
                                   color = `Inflammatory Pvalue`, 
                                   size = abs(`Inflammatory cor`)))+
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#1679AB", high = "#FFB1B1") +
  labs(y = "Person Correlation", x = "Hallmark", title = "Inflammatory") +
  coord_flip()
p2

p3 = ggplot(subset(cor[prolif_group,], `Proliferation Pvalue` <0.05), aes(reorder(Group, `Proliferation cor`),
                                                                          `Proliferation cor`, 
                                                                          color = `Proliferation Pvalue`, 
                                                                          size = abs(`Proliferation cor`)))+
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#1679AB", high = "#FFB1B1") +
  labs(y = "Person Correlation", x = "Hallmark", title = "Proliferation") +
  coord_flip()
p3


p4 = ggplot(subset(cor[meta_group,], `Meatbolism Pvalue`<0.05), aes(reorder(Group, `Meatbolism cor`), 
                                                                    `Meatbolism cor`, 
                                                                    color = `Meatbolism Pvalue`, 
                                                                    size = abs(`Meatbolism cor`)))+
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#1679AB", high = "#FFB1B1") +
  labs(y = "Person Correlation", x = "Hallmark", title = "Metabolism") +
  coord_flip()
p4
pdf("04_immu_metabolism/cor_tgf_immu_metabolism_cellcycle.pdf",height = 6,width = 18)
plot_grid(p2,p3,p4,labels = c("A","B","C"), ncol = 3)
dev.off()
