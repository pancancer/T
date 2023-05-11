# Code for manuscript: Multiomics Profiling Uncovers Pan-Cancer T-Cell Dynamics in Immunotherapy Response

library(Seurat)
library(ComplexHeatmap)
library(AUCell)
library(scFunctions)
library(rcartocolor)

x <- readRDS("Pan_Cancer_T_Seurat.RDS")
cd4 <- subset(x, subset = Lineage == "CD4")
cd8 <- subset(x, subset = Lineage == "CD8")

somePDFPath = paste("UMAP_CD4_CD8_Cell_Types.pdf", sep = "")
pdf(file=somePDFPath, width=10, height=8,pointsize=12)
DimPlot(cd4, group.by = "Cell_Type")
DimPlot(cd8, group.by = "Cell_Type")
dev.off()

files <- list.files(path = "../DEGs/", pattern = "RDS", full.names = T)

allave <- NULL
degs <- NULL
for(i in 1:length(files)){
  cname <- NULL
  cname <- gsub(".*(CD.*)_DEGs.*","\\1",files[i])
  x <- NULL
  x <- readRDS(files[i])
  x <- x[which(x$p_val_adj < 0.05 & x$avg_log2FC > 0.25),]
  x <- x[order(x$p_val_adj, decreasing = F),]
  x <- x[order(x$avg_log2FC, decreasing = T),]
  x <- split(x, x$cluster)
  x <- lapply(x, function(y){
    y <- y[1:ifelse(nrow(y) > 10, 10, nrow(y)),]
  })
  x <- do.call(rbind.data.frame, x)
  if(length(grep("CD4", cname, ignore.case = T)) > 0){
    cd4$ID <- paste(cd4$Cell_Type, cd4$Treatment_Group, sep = "_")
    Idents(cd4) <- "ID"
    ave <- NULL
    ave <- AverageExpression(cd4)$RNA
    allave[[i]] <- ave
    names(allave)[i] <- cname
  }else if(length(grep("CD8", cname, ignore.case = T)) > 0){
      cd8$ID <- paste(cd8$Cell_Type, cd8$Treatment_Group, sep = "_")
      Idents(cd8) <- "ID"
      ave <- NULL
      ave <- AverageExpression(cd8)$RNA
      allave[[i]] <- ave
      names(allave)[i] <- cname
  }
  
  ave <- ave[which(row.names(ave) %in% unique(x$gene)),]
  plotx <- NULL
  plotx <- scale(ave)
  plotx <- t(scale(t(plotx)))
  
  ctype <- NULL
  ctype <- unique(gsub("(CD4|CD8).*","\\1",colnames(plotx)))

  ca = HeatmapAnnotation(show_legend = c(T,T),
                         Class = gsub(".*(Pre|Post)","\\1",colnames(plotx)),
                         Cell_Type_Class = gsub("(CD4|CD8).*","\\1",colnames(plotx)),
                         col = list(
                           Class = c(Pre = "#D75B58", Post = "#F5BFD3"),
                           Cell_Type_Class = c(CD4 = color_conditions$cold[i], CD8 = color_conditions$cold[i])))
  
  cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
  # cols = colorRampPalette(color_conditions$gradient)(100)
  
  csplit <- NULL
  csplit <- gsub(".*(Pre|Post)","\\1",colnames(plotx))
  somePDFPath = paste("DEGs_Heatmap_",cname,".pdf", sep = "")
  pdf(file=somePDFPath, width=6, height=4,pointsize=12)
  print(Heatmap(plotx, column_title = "", name = "Expression", column_split = csplit,
                top_annotation = ca, 
                cluster_columns = T,
                col = cols, row_names_gp = gpar(fontsize = 0),
                column_names_gp = gpar(fontsize = 8),
                show_column_dend = TRUE, show_row_dend = TRUE,
                show_row_names = TRUE, show_column_names = TRUE))
  dev.off()
  
}

cts <- c("CD4","CD8")

pdf(file=paste(cdir,"Signature_Genes_CD4_CD8_DotPlot.pdf",sep = ""), width=10, height=14,pointsize=10)
for(j in 1:length(cts)){
    cmarkers <- NULL
    cmarkers <- toupper(markers[[cts[j]]])
    cref <- NULL
    cref <- read.table(paste(cts[j],"_Markers.txt",sep = ""), header = T, sep = "\t")
    x <- NULL
    if(cts[j] == "CD4"){
      x <- cd4
    }else if(cts[j] == "CD8"){
      x <- cd8
    }
    Idents(x) <- "Cell_Type"
    x$Cell_Type <- factor(x$Cell_Type, levels = sort(unique(x$Cell_Type)))
    p <- NULL
    p <- DotPlot(x, features =  cmarkers, group.by = "Cell_Type", scale = T, cluster.idents = T) +
      theme_light(base_size = 13)+
      coord_flip()+RotatedAxis()+ggtitle(paste(cts[j], sep = ""))+
      scale_color_gradientn(colors = rev(color_conditions$RedYellowBlue))+ #ccols
      scale_size(range = c(1,6))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA), strip.text.x = element_blank(), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 16))+
      xlab(paste(cts[j], ":Signature Markers", sep = ""))+ylab("Celltype")
    p$data$pct.exp <- abs(p$data$avg.exp.scaled)
    p+guides(size = guide_legend(title = "Abs(Average Expression)"))
    print(p)
}
dev.off()

rm(x)

prop <- readRDS("Pre_Post_Prop.RDS")
somePDFPath = paste("Pre_Post_Barplot.pdf", sep = "")
pdf(file=somePDFPath, width=14, height=6,pointsize=12)
ggplot(prop, aes(Cancer, Proportion, fill = Cell_Type, label = Cell_Type_Index))+
  geom_bar(position = "stack", stat = "identity")+
  facet_wrap(~IIID, ncol = 4)+
  scale_fill_manual(values = ccols)+theme_linedraw(base_size = 20)+
  RotatedAxis()+theme(legend.position = "bottom", legend.title= element_blank()) +
  guides(color=guide_legend(ncol=2))+xlab("")
dev.off()

aucresults <- NULL
allaucs <- readRDS("PrePost_CD4CD8_Regulon_AUC_Celltypes.RDS")
groups <- sort(unique(allaucs$Group))
prepost <- c("Pre","Post")
somePDFPath = paste("PrePost_CD4CD8_Regulon_AUC_Tree_Plot.pdf", sep = "")
pdf(file=somePDFPath, width=14, height=14,pointsize=12)
  for(j in 1:length(prepost)){
    # print(i)
    print(j)
    plotx <- allaucs[which(allaucs$Treatment_Group == prepost[j]),]
    plotx$SID <- paste(plotx$Treatment_Group, plotx$Cancer, plotx$Cell_Type, sep = "|")
    plotx$Regulon <- gsub("(.*) \\(.*\\)$","\\1",plotx$Regulon, ignore.case = T)
    plotx$Regulon <- gsub("_extended","",plotx$Regulon, ignore.case = T)
    mplotx <- NULL
    mplotx <- plotx
    plotx <- reshape2::dcast(mplotx,SID ~ Regulon, value.var = "AUC")
    row.names(plotx) <- plotx$SID
    plotx <- plotx[,-1]
    plotx[is.na(plotx)] <- 0
    plotx <- scale(plotx)
    plotx <- t(scale(t(plotx)))
    cdist <- cor(t(plotx))
    ctree <- ape::bionj(dist(cdist, method="euclidean"))
    groupcols <- c(CD4 = color_conditions$cold[1], CD8 = color_conditions$cold[2])
    clabel <- data.frame(SID = row.names(plotx),
                         Treatment_Group = gsub("(.*?)\\|(.*?)\\|(.*)","\\1",row.names(plotx), ignore.case = T),
                         Cancer = gsub("(.*?)\\|(.*?)\\|(.*)","\\2",row.names(plotx), ignore.case = T),
                         Cell_Type = gsub("(.*?)\\|(.*?)\\|(.*)","\\3",row.names(plotx), ignore.case = T)
    )
    clabel <- clabel[match(ctree$tip.label, clabel$SID),]
    clabel$node <- 1:nrow(clabel)
    cctree <- NULL
    cctree <- full_join(ctree, clabel, by = 'node')
    p <- ggtree(cctree, layout="daylight", branch.length = 'none')
    clist <- NULL
    clist <- list(ctree = ctree, cctree = cctree, mplotx = mplotx, plotx = plotx, auctree = p)
    aucresults[[length(aucresults)+1]] <- clist
    names(aucresults)[length(aucresults)] <- paste("AUC", prepost[j], sep = "_")
    
    p <- p + geom_tippoint(aes(color=Cancer, shape = Cell_Type), size=5, alpha=0.8)+
      scale_color_manual(values = cancercols)  +
      scale_shape_manual(values = 1:length(mplotx$Cell_Type))+
      geom_tiplab(aes(label = Cell_Type), size=2)+ggtitle(paste("Regulon Activity: ",prepost[j], sep = ""))+
      theme(plot.title = element_text(hjust = 0.5, size = 20))
    print(p)
}

dev.off()

tcrdata <- readRDS("Frequency_CD4CD8_Clonal_Expansion_Test_Statistics.RDS")
plotx <- tcrdata
p1 <- ggplot(plotx[which(plotx$Cutoff == 0),], aes(Cell_Type, Cancer, color = FC_Group, size = Sig_Level))+
  geom_point() + 
  ggtitle("Post(R) vs Post(NR)")+
  theme_classic(base_size = 20)+RotatedAxis()+xlab("")+
  theme(plot.title = element_text(size = 20, face = "bold"))
p2 <- ggplot(plotx[which(plotx$Cutoff == 0),], aes(Cell_Type, Cancer, color = Clonal_log2FC_G2vsG1, size = Sig_Level))+
  geom_point() +
  ggtitle("Post(R) vs Post(NR)")+
  theme_classic(base_size = 20)+RotatedAxis()+xlab("")+
  theme(plot.title = element_text(size = 20, face = "bold"))

pdf(file=paste("Clonal_Homeostasis_Post_RNR_Specific_TBCRs.pdf", sep = ""), width=12, height=6,pointsize=10)
print(p1)
print(p2)
dev.off()

ccompare <- list(c("Post(R)", "Post(NR)"))

chosen_cutoff <- 0
samplepostrnr <- readRDS(paste("../TBCR/Tuning_New_Category/Cutoff_", chosen_cutoff,"_samplepostrnr.RDS", sep = ""))
pdf(file=paste(cdir, "Figure4A2_Boxplot_Expansion0_Clonal_Homeostasis_Post_RNR_Specific_TBCRs.pdf", sep = ""), width=8, height=8,pointsize=12)
# for(i in 1:length(selected)){

  p <- NULL
  for(i in 1:length(selected)){
  p[[i]] <- ggbetweenstats(current[which(current$TBCR == "TCR" & current$Celltype %in% selected[i]),], title = selected[i], Group, Cell_Num, p.adjust.method = "BH", xlab = "Post-Treatment Groups", ylab = "Clonal Frequency")
  }
  
  combine_plots(p,plotgrid.args = list(nrow = 2), annotation.args = list(title = paste("TCR Clonal Expansion",sep = ""))) # caption = "PanCancer"
  
  p <- NULL
  p <- ggplot(current[which(current$TBCR == "TCR" & current$Celltype %in% selected),], aes(x=Group, y=Cell_Num, color=Group)) + facet_wrap(~Celltype, ncol = 2, scales = "free")+ 
    geom_boxplot(alpha = 0) +  ggtitle(paste("TCR Clonal Expansion",sep = ""))+
    geom_point(size = 2, position = position_jitterdodge()) +
    scale_color_manual(values = ccols$general)+  ylab("Frequency")+xlab("Post-Treatment Groups")+
    stat_compare_means(comparisons = ccompare, method="t.test", label="p.format", step.increase = 0.2, hide.ns = T, hjust = -1, vjust = 1.2, tip.length = 0.1, color = "black", size = 4) + # , aes(label = paste0("p = ", ..p.format..))
    # stat_compare_means(aes(group = Group), method="wilcox.test", label="p.format", hjust = 1, vjust = 0.8, tip.length = 0.1, color = "black")+
    theme_classic(base_size = 20)+
    theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  print(p)
dev.off()

library(liana)

files <- NULL
files <- list.files(path = "CellCell/", full.names = T)
for(i in 1:length(files)){
  x <- NULL
  x <- readRDS(files[i])
  x <- liana_aggregate(x)
  cname <- NULL
  cname <- gsub(".*Immune_(.*)_Liana_Output(.*).RDS","\\1\\2",files[i])
  cname <- gsub("Non-responder_Post_RvsNR","Post_NR",cname, ignore.case = T)
  cname <- gsub("Responder_Post_RvsNR","Post_R",cname, ignore.case = T)
  x$source <- gsub(".*_(.*)","\\1",x$source, ignore.case = T)
  x$target <- gsub(".*_(.*)","\\1",x$target, ignore.case = T)
  x$source <- gsub("(.*)_(Pre|Post)$","\\1",x$source, ignore.case = T)
  x$target <- gsub("(.*)_(Pre|Post)$","\\1",x$target, ignore.case = T)
  
  x$Label <- cname
  saveRDS(x, paste("Figure5C_CRC_CellCell_",cname,".RDS", sep = ""))
  
  somePDFPath = paste("CellCell_",cname,"_Top_20_Interactions.pdf", sep = "")
  pdf(file=somePDFPath, width=20, height=8,pointsize=12)
    plotx <- NULL
    plotx <- filter(x, aggregate_rank < 0.01)
    cts <- NULL
    cts <- sort(unique(x$source))
    for(k in 1:length(cts)){
      if(nrow(plotx[which(plotx$source == cts[k]),]) > 0){
        p <- NULL
        p <- liana_dotplot(plotx, source_groups = cts[k],target_groups = unique(c(plotx$source,plotx$target)),ntop=20)+RotatedAxis()+theme(axis.text.x=element_text(colour="black", size = 18, face = "plain"), plot.title = element_text(hjust = 0.5, size = 30))+scale_color_gradientn(colours = color_conditions$gradient)+ggtitle(paste("CRC: ",cname,sep = ""))
        print(p)
      }
    }
  dev.off()
  
}

################################################################################
################################################################################
# Mono2: State DEGs

library(monocle)
library(Seurat)

files <- NULL
files <- list.files(pattern = "_Figure3D_Mono2_Pseudotime_")
cdir <- "Mono2/"
setwd(cdir)

for(j in 1:length(cts)){
data <- NULL
data <- cts[[j]]

for(i in 1:length(files)){
  print(i)
  cname <- NULL
  cname <- gsub("(.*)_Figure3D_Mono2_Pseudotime_(.*).RDS","\\1_\\2",files[i], ignore.case = T)
  print(i)
  m <- NULL
  m <- readRDS(files[i])
  if(length(which(unique(data$Cell_Type) %in% unique(pData(m)$Cell_Type))) > 0){
    x <- NULL
    x <- subset(data, cells = row.names(pData(m)))
    x$State <- pData(m)[match(colnames(x), row.names(pData(m))),"State"]
    Idents(x) <- "State"
    cmarkers <- NULL
    cmarkers <- FindAllMarkers(x, min.pct = 0.25)
    cmarkers <- cmarkers[which(cmarkers$p_val_adj < 0.05),]
    print(cname)
    saveRDS(cmarkers, paste(cdir, "State_DEGs_",cname, ".RDS", sep = ""))
  }
}
}

files <- NULL
files <- list.files(path = cdir, pattern = "State_DEGs_", full.names = T)

library(msigdbr)
library(fgsea)
library(org.Hs.eg.db)
library(limma)

msigref <- msigdbr(species = "Homo sapiens", category = "H")
msigref <- msigref %>% split(x = .$gene_symbol, f = .$gs_name)

goall <- NULL
keggall <- NULL
gseaall <- NULL

for(i in 1:length(files)){
  print(i)
  x <- NULL
  x <- readRDS(files[i])
  cname <- NULL
  cname <- gsub(".*_State_DEGs_(.*).RDS","\\1",files[i], ignore.case = T)
  
  cclust <- NULL
  cclust <- sort(unique(x$cluster))
  table(x$cluster)
  
  for(j in 1:length(cclust)){
    print(j)
    r <- NULL
    r <- x[which(x$cluster == cclust[j]),]
    r <- r[grep("RPL[0-9]+|^RPS|^MT-", r$gene, ignore.case = T, invert = T),] # 
    r <- r[order(r$avg_log2FC, decreasing = T),]
    
    crank <- NULL
    crank <- r$avg_log2FC
    names(crank) <- r$gene
    Rhallmark <- NULL
    Rhallmark <- fgseaMultilevel(msigref, stats = crank, nPermSimple = 1000, eps = 0)
    if(nrow(Rhallmark) > 0){
      gseaall <- rbind(gseaall, data.frame(Group = cname, State = cclust[j], Rhallmark))
    }
    
    r <- r[which(r$avg_log2FC > 0),]
    ctop <- NULL
    ctop <- data.frame(gene=unique(r$gene)[1:ifelse(nrow(r)>200,200,nrow(r))])
    entrez_out <- NULL
    entrez_out <- AnnotationDbi::select(org.Hs.eg.db, keys = ctop$gene, keytype = 'SYMBOL', columns = 'ENTREZID')
    ctop$ENTREZID <- entrez_out[match(ctop$gene, entrez_out$SYMBOL),"ENTREZID"]
    ctop <- ctop[which(!is.na(ctop$ENTREZID)),]
    
    Rkegga <- NULL
    Rkegga <- kegga(unique(ctop$ENTREZID), species="Hs")
    Rkegga <- Rkegga[order(Rkegga$P.DE, decreasing = F),]
    keggall <- rbind(keggall, data.frame(Group = cname, State = cclust[j], Rkegga))
    
    Rgoana <- NULL
    Rgoana <- goana(unique(ctop$ENTREZID), species="Hs")
    Rgoana <- Rgoana[order(Rgoana$P.DE, decreasing = F),]
    goall <- rbind(goall, data.frame(Group = cname, State = cclust[j], Rgoana))
    
  }
  
}

allkegg <- readRDS("Source_Data/Figure3H_KEGG_Combined_BH005.RDS")
allkegg <- allkegg[which(allkegg$Treatment_Group == "Post" & allkegg$Response %in% c("Responder","Non-responder") & allkegg$Cancer == "NonSpecific"),]
allkegg <- allkegg[which(allkegg$BH < 0.05),]
allkegg <- allkegg[which(allkegg$Regulation == "Up"),]
row.names(allkegg) <- NULL

# 16S rDNA
rdna <- readRDS("PICRUSt2_KEGG_level3_Group.P_vs_S_diff.RDS")

str(rdna)
rdna$Response[which(rdna$Group == "P")] <- "Responder"
rdna$Response[which(rdna$Group == "S")] <- "Non-responder"
table(rdna$Response)
rdna <- split(rdna, rdna$description)
rdna <- lapply(rdna, function(x){
  if((x[which(x$Response == "Responder"),"value.sum.mean"] - x[which(x$Response == "Non-responder"),"value.sum.mean"]) > 0){
    x <- x[which(x$Response == "Responder"),]
  }else{
    x <- x[which(x$Response == "Non-responder"),]
  }
  return(x)
})
rdna <- do.call(rbind.data.frame, rdna)
rdna <- rdna[which(rdna$p_value < 0.05),]

rdnasummary <- NULL
allkegg$ID <- paste(allkegg$Cell_Type, allkegg$Response, sep = "_")

cts <- NULL
cts <- unique(allkegg$ID)

for(i in 1:length(cts)){
  x <- NULL
  x <- allkegg[which(allkegg$ID == cts[i]),]
  y <- NULL
  if(unique(x$Response) == "Responder"){
    y <- rdna[which(rdna$Response == "Responder"),]
  }else if(unique(x$Response) == "Non-responder"){
    y <- rdna[which(rdna$Response == "Non-responder"),]
  }
  
  current <- NULL
  current <- data.frame(unique(x[,c("DB","Cell_Type","Response","Treatment_Group")]))
  current$Overlap <- sum(length(which(unique(x$Pathway) %in% unique(y$description))))
  current$Celltype_KEGG <- length(unique(x$Pathway))
  current$rDNA_KEGG <- length(unique(y$description))
  current$Overlap_Prop <- current$Overlap/current$Celltype_KEGG
  current$Overlap_KEGG <- paste(unique(x[which(unique(x$Pathway) %in% unique(y$description)),"Pathway"]),collapse = ",")
  rdnasummary <- rbind(rdnasummary, current)
}

rdnasummary <- rdnasummary[order(rdnasummary$Overlap, decreasing = T),]
saveRDS(rdnasummary, paste("../Metagenomics_Metabolomics/Sjalv/Overlaps_KEGG_rDNA_vs_Figure3H_KEGG_Combined_BH005_RNR_CancerNonSpecific.RDS", sep = ""))
write.table(rdnasummary, "../Metagenomics_Metabolomics/Sjalv/Overlaps_KEGG_rDNA_vs_CD4CD8_RNA_CancerNonSpecific.txt", row.names = F, quote = F, sep = "\t") # sheetName = "Celltype_vs_16S_rDNA", row.names = F

# Metabolomics
ref <- read.table("/Users/luadmpan/Xuexinli_Group\ Dropbox/Lu\ Pan/KI/Studies/LXXLZ/Website/LXX_GSEA/LXX_HT/LXX_Immune/Data/Metagenomics_Metabolomics/Unzipped_Original/metabolism_r\ nr/4.MetaboliteComparison\ 2/P_S/P_S.significant_R_Input.csv", quote = "\"", sep = ";", header = T)
ref$Response <- NULL
ref[which(ref$regulated == "up"),"Response"] <- "Responder"
ref[which(ref$regulated == "down"),"Response"] <- "Non-responder"
table(ref$Response)

metabolites <- NULL
metabolites <- read.table("/Users/luadmpan/Xuexinli_Group\ Dropbox/Lu\ Pan/KI/Studies/LXXLZ/Website/LXX_GSEA/LXX_HT/LXX_Immune/Data/Metagenomics_Metabolomics/Unzipped_Original/metabolism_r\ nr/4.MetaboliteComparison\ 2/P_S/idms2.pathway/P_S.xls", header = T, sep = "\t")

metabolitessummary <- NULL

for(i in 1:length(cts)){
  x <- NULL
  x <- allkegg[which(allkegg$ID == cts[i]),]
  y <- NULL
  if(unique(x$Response) == "Responder"){
    y <- rdna[which(rdna$Response == "Responder"),]
  }else if(unique(x$Response) == "Non-responder"){
    y <- rdna[which(rdna$Response == "Non-responder"),]
  }
  
  current <- NULL
  current <- data.frame(unique(x[,c("DB","Cell_Type","Response","Treatment_Group")]))
  current$Overlap <- sum(length(which(unique(x$Pathway) %in% unique(y$description))))
  current$Celltype_KEGG <- length(unique(x$Pathway))
  current$rDNA_KEGG <- length(unique(y$description))
  current$Overlap_Prop <- current$Overlap/current$Celltype_KEGG
  current$Overlap_KEGG <- paste(unique(x[which(unique(x$Pathway) %in% unique(y$description)),"Pathway"]),collapse = ",")
  rdnasummary <- rbind(rdnasummary, current)
}

rdnasummary <- rdnasummary[order(rdnasummary$Overlap, decreasing = T),]
saveRDS(rdnasummary, paste("Overlaps_KEGG_rDNA_vs_Figure3H_KEGG_Combined_BH005_RNR_CancerNonSpecific.RDS", sep = ""))
