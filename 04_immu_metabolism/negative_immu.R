###hallmarker
options(stringsAsFactors = F)
#go####
#Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(Fibroblast_cells@assays$RNA@data))

# load gene set
hallmark_go <- read.gmt("02_AUCell/c5.go.bp.v2023.2.Hs.symbols.gmt")
head(hallmark_go)

#
geneSets_go<-lapply(unique(hallmark_go$ont),function(x){hallmark_go$gene[hallmark_go$ont==x]})
names(geneSets_go) <- unique(hallmark_go$ont)
saveRDS(geneSets_go,file ='02_AUCell/go_geneSets.RDS' )
#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC_go <- AUCell_calcAUC(geneSets_go, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)

# 
geneSet_go <- unique(hallmark_go$ont)
geneSet_go
aucs_go <- getAUC(cells_AUC_go)#[geneSet_go, ]
aucs_go = t(aucs_go) %>%
  as.data.frame() %>%
  write.csv("02_AUCell/go_aucs.csv")

###go####
site <- Fibroblast_cells$site
aucs_go <- aucs_go[names(site),]
aucs_go$site <- site
head(aucs_go)
dim(aucs_go)
#1207 7629
write.csv(aucs_go,file = "04_immu_metabolism/aucs_site_go.csv",row.names = T)

#########
# 
aucs_go = read.csv("02_AUCell/go_aucs.csv", header = T, row.names = 1)
colnames(aucs_go)
colnames(aucs_go) <- sub('^[^_]*_', '', colnames(aucs_go))
dput(colnames(aucs_go))
signaling_name <- as.data.frame(colnames(aucs_go))
filtered_signaling_name <- signaling_name[grep("negative", signaling_name$`colnames(aucs_go)`, ignore.case = TRUE), ]
filtered_signaling_name <-filtered_signaling_name[ grep("CELL|IMMUNE", filtered_signaling_name, ignore.case = TRUE)]
write_csv(as.data.frame(filtered_signaling_name),file = '04_immu_metabolism/filtered_signaling_name.csv')
filtered_signaling_name <- read.csv("04_immu_metabolism/filtered_signaling_name_by_hand.csv")
singal <- filtered_signaling_name$filtered_signaling_name
aucs = read.csv("04_immu_metabolism/aucs_site.csv", header = T,row.names = 1)
aucs_go <- cbind(aucs_go,aucs[,c('HALLMARK_TGF_BETA_SIGNALING',"HALLMARK_GLYCOLYSIS")])
aucs_go[,Hh_group]
# 
## tgf
Hh_group = c("HALLMARK_TGF_BETA_SIGNALING")

## Inflammatory Group

infla_group = singal

cor = data.frame()
for(i in infla_group){
  
  result = cor.test(aucs_go[,Hh_group], aucs_go[,c(i)])
  cor[i,"Inflammatory cor"] = result$estimate
  cor[i,"Inflammatory Pvalue"] = result$p.value

}
cor

cor$Group = rownames(cor)
cor
saveRDS(cor,file = "04_immu_metabolism/cor_GO.RDS")

p2 = ggplot(subset(cor[infla_group,],subset = cor$`Inflammatory Pvalue`<0.05&cor$`Inflammatory cor`>0) ,aes(reorder(Group, `Inflammatory cor`), 
                                                                      `Inflammatory cor`, 
                                                                      color = `Inflammatory Pvalue`, 
                                                                      size = abs(`Inflammatory cor`)))+
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#1679AB", high = "#FFB1B1") +
  labs(y = "Person Correlation", x = "Hallmark", title = "Inflammatory") +
  coord_flip()
p2
ggsave("04_immu_metabolism/negative_immu.pdf",p2,width = 10,height = 6)
