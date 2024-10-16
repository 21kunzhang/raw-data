#####
hallmark <- read.gmt("02_AUCell/LIN_TUMOR_ESCAPE_FROM_IMMUNE_ATTACK.v2023.2.Hs.gmt")
head(hallmark)

#
geneSets<-lapply(unique(hallmark$ont),function(x){hallmark$gene[hallmark$ont==x]})
names(geneSets) <- unique(hallmark$ont)
saveRDS(geneSets,file ='02_AUCell/geneSets_escape.RDS' )

#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(geneSets, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)

# 
geneSet <- unique(hallmark$ont)
geneSet
aucs <- getAUC(cells_AUC)[geneSet, ]
aucs = t(aucs) %>%
  as.data.frame() %>%
  write.csv("02_AUCell/aucs_escape.csv")
# options(stringsAsFactors = F)



#########
# 
aucs_escape = read.csv("02_AUCell/aucs_escape.csv", header = T, row.names = 1)
aucs_escape <- data_frame(t(aucs_escape))
colnames(aucs_escape)[1]='TUMOR_ESCAPE_FROM_IMMUNE_ATTACK'
colnames(aucs_escape) <- sub('^[^_]*_', '', colnames(aucs_escape))
dput(colnames(aucs_escape))
signaling_name <- as.data.frame(colnames(aucs_escape))
filtered_signaling_name <- signaling_name
write_csv(as.data.frame(filtered_signaling_name),file = '04_immu_metabolism/filtered_signaling_name_escape.csv')
filtered_signaling_name <- read.csv("04_immu_metabolism/filtered_signaling_name_escape.csv")
singal <- filtered_signaling_name$colnames.aucs_escape.
aucs = read.csv("04_immu_metabolism/aucs_site.csv", header = T,row.names = 1)
aucs_escape <- cbind(aucs_escape,aucs[,c('HALLMARK_TGF_BETA_SIGNALING',"HALLMARK_GLYCOLYSIS")])
aucs_escape[,Hh_group]
# 
## tgf
Hh_group = c("HALLMARK_TGF_BETA_SIGNALING")

## METABOLIC Group

infla_group = singal

cor = data.frame()
for(i in infla_group){
  
  result = cor.test(aucs_escape[,Hh_group], aucs_escape[,c(i)])
  cor[i,"Inflammatory cor"] = result$estimate
  cor[i,"Inflammatory Pvalue"] = result$p.value
  
}
cor

cor$Group = rownames(cor)
cor
saveRDS(cor,file = "04_immu_metabolism/cor_GO_Inflammatory_escape.RDS")
#
cor_GO <- readRDS("/home/pub252/users/liy/20240715_CESC/04_immu_metabolism/cor_GO.RDS")
filtered_signaling_name_immu <- read.csv("04_immu_metabolism/filtered_signaling_name_by_hand.csv")
colnames(cor) <- dput(colnames(cor_GO))
cor_all <- rbind(cor_GO[filtered_signaling_name_immu$filtered_signaling_name,],cor)

p2 = ggplot(subset(cor_all,subset = cor_all$`Inflammatory Pvalue`<0.05&cor_all$`Inflammatory cor`>0) ,aes(reorder(Group, `Inflammatory cor`), 
                                                                                                      `Inflammatory cor`, 
                                                                                                      color = `Inflammatory Pvalue`, 
                                                                                                      size = abs(`Inflammatory cor`)))+
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "#1679AB", high = "#FFB1B1") +
  labs(y = "Person Correlation", x = "Hallmark", title = "Inflammatory") +
  coord_flip()
p2
ggsave("04_immu_metabolism/cor_Inflammatory_escape.pdf",p2,width = 10,height = 6)
