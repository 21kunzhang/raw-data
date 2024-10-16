dir.create("05_CellChat")

dim(mydata@meta.data)
mydata[['site']]=ifelse(substring(mydata$orig.ident,1,2)=='CA','Tumor','Normal')
mydata_tumor<-subset(mydata,site=="Tumor")
dim(mydata_tumor@meta.data)
devtools::install_github("sqjin/CellChat")
save(mydata_tumor,file = 'mydata_tumor.Rdata')
Sys.setenv(RETICULATE_PYTHON="/usr/bin/python3")
# 
library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
# install.packages("mindr")
options(stringsAsFactors = FALSE)

###1. 
#
#⚠️：
cellchat_tumor <- createCellChat(object=mydata_tumor,group.by = "cell_type")
#cellchat_tumor <- createCellChat(ma_cm_disease@assays$RNA@data, meta = ma_cm_disease@meta.data, group.by = "cell_type")
cellchat_tumor
summary(cellchat_tumor)
str(cellchat_tumor)
levels(cellchat_tumor@idents)
#cellchat_tumor <- setIdent(cellchat_tumor, ident.use = "cell_type")
groupSize <- as.numeric(table(cellchat_tumor@idents))  
#
groupSize
# [1] 8480 3192 3623 2472  374  444  192  451  135
###2.
CellChatDB <- CellChatDB.human
#ellChatDB <- CellChatDB.mouse
str(CellChatDB) #
#interaction、complex、cofactor geneInfo dataframe
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new() #
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)#
#
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat_tumor@DB <- CellChatDB.use # set the used database in the object
####
cellchat_tumor <- subsetData(cellchat_tumor)
#future::plan("multiprocess", workers = 4)
#
cellchat_tumor <- identifyOverExpressedGenes(cellchat_tumor)
cellchat_tumor <- identifyOverExpressedInteractions(cellchat_tumor) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#cellchat@LR$LRsig
cellchat_tumor <- projectData(cellchat_tumor, PPI.human) 
#
#
cellchat_tumor <- computeCommunProb(cellchat_tumor, raw.use = FALSE, population.size = TRUE) #
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_tumor <- filterCommunication(cellchat_tumor, min.cells = 10)
df.net <- subsetCommunication(cellchat_tumor)
write.csv(df.net, "05_CellChat/mydata_tumor_secreted_net.csv")

cellchat_tumor <- computeCommunProbPathway(cellchat_tumor)
df.netp <- subsetCommunication(cellchat_tumor, slot.name = "netP")
write.csv(df.netp, "05_CellChat/mydata_tumor_secreted_net_pathway.csv")
save(cellchat_tumor ,file = "05_CellChat/cellchat_tumor.Rdata")
###
cellchat_tumor <- aggregateNet(cellchat_tumor)
#
groupSize <- as.numeric(table(cellchat_tumor@idents))
pdf("05_CellChat/cell_type_number.pdf", width = 12, height =6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_tumor@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_tumor@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()

##
mat <- cellchat_tumor@net$count
par(mfrow = c(4,4), xpd=TRUE)
#for (i in 1:nrow(mat)) {
   i = 6
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf("05_CellChat/cell_type_count.pdf",cell_type_count,width = 4,height = 4)
   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
#}
dev.off()
#
mat <- cellchat_tumor@net$weight
#par(mfrow = c(3,3), xpd=T)
#for (i in 1:nrow(mat)) {
  i=6
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf("05_CellChat/cell_type_weight.pdf",cell_type_count,width = 4,height = 4)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
#}
# save as TIL/net_strength_individual.pdf
  dev.off()
######
cellchat_tumor@netP$pathways  #
pathways.show <- c("TGFb")
#（Hierarchy plot）
levels(cellchat_tumor@idents)    # show all celltype
vertex.receiver = c(1,2,3,4,5,7,8,9) # define a numeric vector giving the index of the celltype as targets
#par(mar=c(5.1,4.1,4.1,2.1))
pdf("05_CellChat/netVisual_TGF.pdf",cell_type_count,width = 4,height = 4)
netVisual_aggregate(cellchat_tumor, signaling = pathways.show,  vertex.receiver = vertex.receiver,sources.use = "Fibroblast cells" )
dev.off()
#
pdf("05_CellChat/netAnalysis_contribution_TGF.pdf",cell_type_count,width = 6,height = 6)
netAnalysis_contribution(cellchat_tumor, signaling = pathways.show)
dev.off()
pairLR.TGFb <- extractEnrichedLR(cellchat_tumor, signaling = pathways.show, geneLR.return = FALSE) #
# save as TIL/CXCL_LR_contribution.pdf
p_netVisual_individual <- list()
for (i in  rownames(pairLR.TGFb)) {
  
LR.show <- pairLR.TGFb[i,] 
vertex.receiver =  c(1,2,3,4,5,7,8,9) # a numeric vector
netVisual_individual(cellchat_tumor, signaling = pathways.show, pairLR.use = LR.show, layout = "circle",sources.use = "Fibroblast cells")

}
pdf("05_CellChat/etVisual_individual_TGF3.pdf",width = 6,height = 6)
LR.show <- pairLR.TGFb[3,] 
vertex.receiver =  c(1,2,3,4,5,7,8,9) # a numeric vector
netVisual_individual(cellchat_tumor, signaling = pathways.show, pairLR.use = LR.show, layout = "circle",sources.use = "Fibroblast cells")
dev.off()
pdf("05_CellChat/etVisual_individual_TGF1.pdf",width = 6,height = 6)
LR.show <- pairLR.TGFb[1,] 
vertex.receiver =  c(1,2,3,4,5,7,8,9) # a numeric vector
netVisual_individual(cellchat_tumor, signaling = pathways.show, pairLR.use = LR.show, layout = "circle",sources.use = "Fibroblast cells")
dev.off()
pdf("05_CellChat/etVisual_individual_TGF9.pdf",cell_type_count,width = 6,height = 6)
LR.show <- pairLR.TGFb[9,] 
vertex.receiver =  c(1,2,3,4,5,7,8,9) # a numeric vector
netVisual_individual(cellchat_tumor, signaling = pathways.show, pairLR.use = LR.show, layout = "circle",sources.use = "Fibroblast cells")
dev.off()
pdf("05_CellChat/netVisual_bubble_TGF.pdf",cell_type_count,width = 8,height = 8)
pairLR.use <- extractEnrichedLR(cellchat_tumor, signaling = c("TGFb"))
netVisual_bubble(cellchat_tumor, sources.use = c(6), targets.use = c(1,2,3,4,5,7,8,9), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)
dev.off()
