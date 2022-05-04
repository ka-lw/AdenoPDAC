suppressPackageStartupMessages(c(library(ggfortify),
                                 library(DESeq2),
                                 library(tidyverse),
                                 library(ggrepel),
                                 library(ggfortify), 
                                 library(edgeR),
                                 library(readxl),
                                 library(EnhancedVolcano),
                                 library(biomaRt),
                                 library(ggbreak),
                                 library(pheatmap)))

# Input parameters:
# raw_counts			integer dataframe of samples (columns) and genes (rows)
# metaData		    dataframe with sample name and group (columns) of individual mice (rows)
# Ggenes/Kgenes		character vector of genes exported from GO Term / Kegg Pathway analysis
# enrichr_out 		character dataframe of 4 columns summarising Enrichr output; Type (GO/Kegg), Pathway, padj, log(padj)

keep<-filterByExpr(raw_counts, group=metaData$Group) 
raw_counts <- raw_counts[keep,]; rm(keep)

mraw_counts<-as.matrix(raw_counts) #vst from DESEQ package needs either matrix or dds object
vst_counts <- vst(mraw_counts, blind = T)
pcDat <- prcomp(t(vst_counts))
autoplot(pcDat,
         labSize=3.5,
         data = metaData, 
         colour="Group",
         shape=FALSE)+ expand_limits(x = c(0.75,-0.75), y = c(0.5,-0.5))

# remove outlier 23729
raw_counts<-raw_counts[,-c(which(colnames(raw_counts) == "sample_23729" ))]
metaData<-metaData[-c(which(metaData$SampleName == "23729" )),]

dds <- DESeqDataSetFromMatrix(countData=raw_counts, 
                              colData=metaData, 
                              design= ~ Group, tidy = F)
vsd <- vst(dds, blind =T )

# DESeq2
dds <- DESeq(dds)

res <- results(dds) %>%
  data.frame() %>%
  arrange(padj)

ensembl=useMart("ensembl")

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
res$geneid<-rownames(res)
genenames<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = res$geneid, 
                 mart = ensembl)

genenames$geneid<-genenames$ensembl_gene_id
res_named<-merge(genenames, res, by="geneid"); res_named<-res_named[,-2]; colnames(res_named)[2]<-"gene"

res_named<- res_named%>%
  arrange(padj)


######## VISUALISE #########

#### VOLCANO ####
res_named<-res_named[!(res_named$gene == ""),]
pval.set<-0.05

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_table_thres <- res_named %>% 
  mutate(threshold = padj < pval.set & abs(log2FoldChange) >= 0.58)

res_table_thres_ordered <- res_table_thres[order(res_table_thres$padj), ] 

res_table_thres_ordered$genelabels <- ""

#label genes above threshold with lowest adj pvalue
for (i in 1:25){
  res_table_thres_ordered$genelabels[i]<-ifelse(res_table_thres_ordered$threshold[i] == TRUE, res_table_thres_ordered$gene[i], "" )
}

res_table_thres_ordered<-arrange(res_table_thres_ordered ,log2FoldChange)
for (i in 1:15){ #label genes above threshold with highest negative logfoldchange
  res_table_thres_ordered$genelabels[i]<-ifelse(res_table_thres_ordered$threshold[i] == TRUE, res_table_thres_ordered$gene[i], "" )
}

res_table_thres_ordered<-arrange(res_table_thres_ordered,desc(log2FoldChange))
for (i in 1:15){ #label genes above threshold with highest positive logfoldchange
  res_table_thres_ordered$genelabels[i]<-ifelse(res_table_thres_ordered$threshold[i] == TRUE, res_table_thres_ordered$gene[i], "" )
}

g.names<-res_table_thres_ordered$genelabels[which(res_table_thres_ordered$genelabels != "")]; length(g.names)
paste("'",as.character(g.names),"'",collapse=", ",sep="")
res_named$anno<-ifelse(res_named$gene %in% g.names, 1, 0 )
res1<-res_named
res1 <- res1 %>%
  arrange(anno)

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
keyvals <- ifelse(
  res1$log2FoldChange < -0.58 & res1$padj < pval.set, '#78ABD0',
  ifelse(res1$log2FoldChange > 0.58 & res1$padj < pval.set, '#E1A46B',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == '#E1A46B'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'none'
names(keyvals)[keyvals == '#78ABD0'] <- 'down'

x<-bquote(~Log[2]~ 'fold change')
y<-bquote(~-Log[10]~ '(adj' ~italic (P)~')')
pval.set<-0.05
labSize = 5
pointSize = 2
axisLabSize = 15
axistextsize = 15

EnhancedVolcano(res1,
                lab = res1$gene,
                selectLab = g.names,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = x,
                ylab = y,
                pCutoff = pval.set,
                FCcutoff = 0.58,
                cutoffLineCol = "grey48",
                cutoffLineWidth = 0.4,
                cutoffLineType="dashed",
                pointSize = pointSize,
                axisLabSize = axisLabSize,
                colCustom = keyvals,
                labCol = 'black',
                labSize = labSize,
                title = "",
                caption="",
                subtitle = "",
                colAlpha = ifelse(res1$gene %in% g.names, 1,0.35),
                legendPosition = '',
                gridlines.major = TRUE,
                gridlines.minor = FALSE,               
                border = 'partial',
                borderWidth = 0.3,
                borderColour = 'black'
) +
  theme(panel.grid.major = element_line(size=0.5), 
        axis.ticks=element_line(colour = "black", size = 0.2),
        axis.ticks.length=unit(3, "pt"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) + scale_y_continuous(breaks = c(0, 5, 10, 15, 20), limits = c(-0.5,21.1), position="left")+
  scale_y_break(c(15.5, 19.5), expand = F, space = 0.5) +
  scale_x_continuous(breaks = c(-4,-2,0,2,4,6), limits = c(-4.8,6), expand=c(0,0))

#### HEATMAP ####
# Generate normalised counts 
dds <- estimateSizeFactors(dds)
normalised<-as.data.frame(counts(dds, normalized=TRUE))
genenames<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = rownames(normalised), 
                 mart = ensembl)
genenames$geneid<-genenames$ensembl_gene_id; normalised$geneid<-rownames(normalised)
norm_named<-merge(genenames, normalised, by="geneid");norm_named<-norm_named[,-c(1,2)]; colnames(norm_named)[1]<-"gene"
colnames(norm_named)[2:12]<-sub(".*_", "", colnames(norm_named)[2:12]) ; norm_named<-norm_named[norm_named$gene != "NA",]



GOgenes<-as.data.frame(Ggenes) %>% unlist(use.names = F) %>% unique()
Kegggenes<-as.data.frame(Kgenes) %>% unlist(use.names = F) %>% unique()


hmapgenes<-GOgenes #Kegggenes

# Subset normalised counts by the Kegg/GO terms 
m<-match(hmapgenes,norm_named$gene ); mat<-as.data.frame(norm_named[m,]) #this method keeps the same order as the original GO list
rownames(mat)<-mat$gene ; mat<- mat[,-c(1)]

anno <- as.data.frame(metaData[,3]); rownames(anno)<-metaData$SampleName;colnames(anno)<-"Group"; anno<-anno %>% arrange(Group); anno

mat<-mat[,match(rownames(anno), colnames(mat))] 

my_colour =list(
  "Group" = c(Treatment = "azure4", Control = "azure3")
)

pheatmap(mat, 
         cluster_cols = F,
         cluster_rows = F,
         scale="row", 
         annotation_col = anno, 
         fontsize_col = 15,
         fontsize_row = 15,
         fontsize = 15,
         annotation_colors = my_colour,
         annotation_legend = F,
         legend_breaks = c(-2,2), 
         main = "", legend_labels = c("Low Expression", "High Expression"),
         legend = T, 
         cutree_cols = 3) 


#### BARPLOT ####

enrichr_out<-enrichr_out %>% group_by(Type) %>% arrange((logpadj), .by_group = TRUE)

ggplot(enrichr_out, aes(x = Pathway, y = logpadj))+
  geom_col(aes(fill = Type), width = 0.7, show.legend = F) + coord_flip()+
  scale_x_discrete(limits=enrichr_out1$Pathway) +
  scale_y_continuous(breaks = c(0,2,4,6,10), limits=c(0, 10.7), expand = c(0, 0)) +
  ylab(bquote(~-Log[10]~ '(adj' ~italic (P)~')')) +
  xlab('') +
  scale_fill_manual(values=c("grey25", "grey50")) + #"#009f97", "#530058"
  scale_y_break(c(6.3, 9.7), expand = F, space = 0.3) +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15))