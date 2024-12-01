
#???ð?
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05      #pֵ????????
qvalueFilter=0.05      #????????pֵ????????

#??????ɫ
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
setwd("D:\\bioinfor\\sepsis\\glycolysis\\21.3 KEGG")            #???ù???Ŀ¼
rt=read.table("interGenes.txt", header=T, sep="\t", check.names=F)     #??ȡ?????ļ?

#????????ת??Ϊ????id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #ȥ??????idΪNA?Ļ???
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg????????
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#?????????????Ľ???
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#??????ʾͨ·????Ŀ
showNum=7
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#??״ͼ
pdf(file="kegg barplot.pdf", width=9, height=10)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, color=colorSel)
dev.off()

#????ͼ
pdf(file="kegg bubble.pdf", width = 9, height = 10)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, color=colorSel)
dev.off()



