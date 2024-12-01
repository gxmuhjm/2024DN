library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
scRNA <- readRDS("D:\\bioinfor\\sepsis\\DKD_identify.rds")

# ×Ô¶¨ÒåÑÕÉ«ÏòÁ¿£¬È·±£28ÖÖÑÕÉ«¸ß¶Ô±ÈÇÒÊÊºÏSCIENCEÆÚ¿¯·ç¸ñ
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", 
            "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff6347", 
            "#4682b4", "#32cd32", "#ffa500", "#a52a2a", "#daa520", "#8a2be2", "#7fff00", 
            "#dc143c", "#00ced1", "#9400d3", "#ff1493", "#00bfff", "#696969", "#ff4500")

# Ê¹ÓÃ×Ô¶¨ÒåÑÕÉ«»æÖÆUMAPÍ¼
p <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", label = TRUE, 
             cols = colors, 
             label.size = 3.5,  # µ÷Õû±êÇ©×ÖÌå´óÐ¡
             pt.size = 0.6) +   # µ÷ÕûÊý¾Ýµã´óÐ¡
  theme_minimal(base_size = 14) +  # Ê¹ÓÃ¼ò½àµÄÖ÷Ìâ
  theme(legend.position = "right",  # Í¼Àý·ÅÖÃÔÚÓÒ²à
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white", color = NA))

# Èç¹û±êÇ©ÓëµãÖØµþ£¬Ê¹ÓÃgeom_text_repelÀ´±ÜÃâÖØµþ
library(ggrepel)
p <- p + geom_text_repel(aes(label = seurat_clusters), size = 3.5)

# ÏÔÊ¾Í¼ÐÎ
print(p)

library(ggplot2)
library(ggplot2)

# ´´½¨ DimPlot
p <- DimPlot(scRNA, reduction = "umap", group.by = "cell_identify", label = TRUE, 
             cols = c("PT" = "#1f77b4", "LOH" = "#ff7f0e", "DCT" = "#2ca02c", 
                      "PC" = "#d62728", "IC" = "#9467bd", "Podo"="#8c564b",
                      "Fib" = "#e377c2" , "Endo" = "#7f7f7f" , "T cells" = "#bcbd22", "NKT" = "#17becf",
                      "B cells" = "#aec7e8", "Mac" = "#ffbb78", "Mono" = "#98df8a"
                      ),
             label.size = 4, 
             pt.size = 0.5)

# ÐÞ¸ÄµãµÄÍ¸Ã÷¶È
p <- p + scale_alpha(range = c(0, 0))

# ½øÒ»²½µ÷ÕûÍ¼ÐÎµÄÃÀ¹Û¶È
p + theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white", color = NA))


library(ggplot2)
library(dplyr)
library(ggplot2)
library(dplyr)

# ¼ÆËãÃ¿¸ö×é±ðÖÐÃ¿ÖÖÏ¸°ûÀàÐÍµÄÊýÁ¿
cell_counts <- scRNA@meta.data %>%
  group_by(sample_type, cell_identify) %>%
  summarise(count = n()) %>%
  group_by(sample_type) %>%
  mutate(percentage = count / sum(count) * 100)  # ¼ÆËãÃ¿¸öÏ¸°ûÀàÐÍÔÚÃ¿¸ö×éÖÐµÄ±ÈÀý

# »æÖÆ¶ÑµþÌõÐÎÍ¼£¬ÏÔÊ¾Ã¿¸ö×é±ðÖÐÃ¿ÖÖÏ¸°ûÀàÐÍµÄ±ÈÀý
ggplot(cell_counts, aes(x = sample_type, y = percentage, fill = cell_identify)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Cell Type Composition by Sample Type", x = "Sample Type", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values = c("PT" = "#1f77b4", "LOH" = "#ff7f0e", "DCT" = "#2ca02c", 
                               "PC" = "#d62728", "IC" = "#9467bd", "Podo"="#8c564b",
                               "Fib" = "#e377c2" , "Endo" = "#7f7f7f" , "T cells" = "#bcbd22", "NKT" = "#17becf",
                               "B cells" = "#aec7e8", "Mac" = "#ffbb78", "Mono" = "#98df8a"))

# ±£´æÍ¼ÐÎÊ±ÉèÖÃ¸ß¶ÈºÍ¿í¶È
ggsave("Cell_Composition.pdf", plot = last_plot(), width = 5, height = 6, dpi = 300)
ggplot(cell_counts, aes(x = sample_type, y = percentage, fill = cell_identify)) +
  geom_bar(stat = "identity", position = "fill", width = 0.5) +  # ½« width ÉèÖÃÎª0.5
  theme_minimal() +
  labs(title = "Cell Type Composition by Sample Type", x = "Sample Type", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values = c("PT" = "#1f77b4", "LOH" = "#ff7f0e", "DCT" = "#2ca02c", 
                               "PC" = "#d62728", "IC" = "#9467bd", "Podo"="#8c564b",
                               "Fib" = "#e377c2" , "Endo" = "#7f7f7f" , "T cells" = "#bcbd22", "NKT" = "#17becf",
                               "B cells" = "#aec7e8", "Mac" = "#ffbb78", "Mono" = "#98df8a"))



        
features = c("GATM", "PCBD1", "F11", "HRSP12", "G6PC")

FeaturePlot(scRNA, features = features)


features = c("GATM", "PCBD1", "F11", "HRSP12", "G6PC")

features = c("GATM", "PCBD1", "F11", "HRSP12", "G6PC")

FeaturePlot(scRNA, features = features, reduction = "umap", 
            cols = c("lightgrey", "skyblue", "dodgerblue", "firebrick"), # Ê¹ÓÃìÅ¿áµÄÑÕÉ«ÌÝ¶È
            pt.size = 1.5, # Ôö´óµãµÄ´óÐ¡
            order = TRUE) + # ±£Ö¤¸ß±í´ïµÄµãÔÚÉÏ²ãÏÔÊ¾
  theme_void() + # È¥³ý±³¾°
  theme(legend.position = "bottom", # Í¼ÀýÎ»ÖÃÉèÖÃÔÚµ×²¿
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "darkred"), # ÉèÖÃìÅ¿áµÄ±êÌâÑùÊ½
        legend.text = element_text(size = 12, color = "darkblue"), # ÉèÖÃÍ¼ÀýÎÄ±¾ÑùÊ½
        legend.key.size = unit(1.2, "cm"), # µ÷ÕûÍ¼Àý´óÐ¡
        plot.margin = margin(10, 10, 10, 10)) # µ÷ÕûÍ¼ÐÎµÄ±ß¾à
Idents(scRNA) <- scRNA$sample_type
setwd("D:\\bioinfor\\sepsis\\2024.8.31")
Idents(scRNA)= "cell_identify"
DimPlot(scRNA)
degs <- FindAllMarkers(scRNA, logfc.threshold = 0.5,
                       test.use = "roc", 
                       return.thresh = 0.25, 
                       min.pct = 0.3, only.pos = T) 
write.csv(degs, file = "degs1.csv")




###»ùÒòµÄÐ¡ÌáÇÙÍ¼
scRNA <- readRDS("D:\\bioinfor\\sepsis\\DKD_identify.rds")
features = c("GATM", "PCBD1", "F11", "HRSP12", "G6PC")
VlnPlot(scRNA, features = features, pt.size = 0, ncol = 2)+
  scale_x_discrete("")+
  theme(
    axis.text.x.bottom = element_blank()
  )

# ´Ó²îÒì·ÖÎö½á¹ûÖÐÌáÈ¡Äã¸ÐÐËÈ¤µÄ»ùÒòµÄÏÔÖøÐÔÐÅÏ¢
degs_of_interest <- degs[degs$gene %in% features,]
library(Seurat)
library(ggplot2)
library(ggpubr)

# Éú³ÉÐ¡ÌáÇÙÍ¼
p <- VlnPlot(scRNA, features = features, pt.size = 0, ncol = 2) +
  scale_x_discrete("") +
  theme(axis.text.x.bottom = element_blank())

# Ìí¼ÓÏÔÖøÐÔ±ê¼Ç
p + stat_compare_means(aes(group = sample_type), 
                       label = "p.signif", 
                       method = "wilcox.test", 
                       hide.ns = TRUE, 
                       comparisons = list(c("DKD", "HC")), 
                       label.y = c(1.5, 2, 2.5, 3, 3.5))  # ¸ù¾ÝÊý¾Ýµ÷Õû y ÖáÎ»ÖÃ
VlnPlot(object=scRNA,features = features,slot = "data",group.by = 'cell_identify',split.by = 'sample_type',pt.size = 0)+stat_compare_means(method = 't.test')+stat_boxplot()
VlnPlot(object=scRNA,features = features,slot = "data",group.by = 'sample_type',pt.size = 0)+stat_compare_means(method = 't.test')+stat_boxplot()
getwd()
VlnPlot(object=scRNA, features = features, slot = "data", group.by = 'cell_identify', pt.size = 0) +
  stat_compare_means(method = 'wilcox.test', label = "p.format", hide.ns = F)
# ½« cell_identify ÉèÎª»î¶¯±êÊ¶·û
scRNA <- SetIdent(scRNA, value = "cell_identify")

# ×Ô¶¨ÒåÑÕÉ«Ó³Éä
cell_colors <- c("PT" = "#1f77b4", "LOH" = "#ff7f0e", "DCT" = "#2ca02c", 
                 "PC" = "#d62728", "IC" = "#9467bd", "Podo"="#8c564b",
                 "Fib" = "#e377c2" , "Endo" = "#7f7f7f" , "T cells" = "#bcbd22", 
                 "NKT" = "#17becf", "B cells" = "#aec7e8", 
                 "Mac" = "#ffbb78", "Mono" = "#98df8a")

# »æÖÆÐ¡ÌáÇÙÍ¼
VlnPlot(object = scRNA, features = "G6PC", slot = "data", pt.size = 0)  +
  scale_fill_manual(values = cell_colors) +
  theme_minimal()


library(ggplot2)
library(dplyr)

# ¼ÙÉèÄãÓÐÒ»¸ö°üº¬»ùÒò±í´ïºÍpÖµµÄÊý¾Ý¿òdf£¬ÐÎÈç:
 df <- data.frame(
   cell_identify = rep(unique(scRNA$cell_identify), each = length(features)),
   gene = rep(features, times = length(unique(scRNA$cell_identify))),
   expression = runif(length(features) * length(unique(scRNA$cell_identify))),
   p_value = runif(length(features) * length(unique(scRNA$cell_identify)))
 )
 # ´´½¨ÆøÅÝÍ¼
 ggplot(df, aes(x = cell_identify, y = gene)) +
   geom_point(aes(size = expression, color = -log10(p_value))) +  # ÓÃ±í´ïÖµ¿ØÖÆÆøÅÝ´óÐ¡£¬pÖµ¿ØÖÆÑÕÉ«
   scale_size(range = c(1, 10)) +  # µ÷ÕûÆøÅÝ´óÐ¡·¶Î§
   scale_color_viridis_c(option = "E") +  # Ê¹ÓÃ Viridis ÅäÉ«·½°¸
   theme_minimal() +
   xlab("cell_identify") +
   ylab("gene") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # µ÷ÕûxÖá±êÇ©½Ç¶È
         axis.text.y = element_text(size = 10))  # µ÷ÕûyÖá±êÇ©´óÐ¡
# ×Ô¶¨ÒåScience·ç¸ñµÄÅäÉ«
library(ggplot2)
library(ggpubr)

library(ggplot2)
library(ggpubr)

VlnPlot(object = scRNA, features = "PCBD1", slot = "data", group.by = 'sample_type', pt.size = 0) +
  stat_compare_means(method = 'wilcox.test', label = "p.format", hide.ns = FALSE) +
  scale_fill_manual(values = c("DKD(CKD)" = "#E6846D", "Health kidney" = "#8DCDD5"))  # ÐÞ¸ÄÑÕÉ«


library(ggplot2)
library(ggpubr)

####AUC
memory.limit(32000000)
##ç³»ç»ŸæŠ¥é”™æ”¹ä¸ºè‹±æ–‡
Sys.setenv(LANGUAGE = "en")
##ç¦æ­¢è½¬åŒ–ä¸ºå› å­?
options(stringsAsFactors = FALSE)
##æ¸…ç©ºçŽ¯å¢ƒ
rm(list=ls())

setwd("D:\\bioinfor\\sepsis\\DN_old_mito\\28.AUCell")
###åŠ è½½æ‰€éœ€è¦çš„åŒ?
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(ggplot2)

##BiocManager::install("AUCell")
library(AUCell)

##mn BiocManager::install("clusterProfiler")
library(clusterProfiler)

sc.id=sample(colnames(scRNA),10000)
sc2=scRNA[,sc.id]
##install.packages("doParallel")
##install.packages("doRNG")
slotNames(sc2@assays$SCT)
cells_rankings <- AUCell_buildRankings(sc2@assays$SCT@data,  nCores=6, plotStats=TRUE) 

cells_rankings

c2 <- read.gmt("DN.gmt") 
geneSets <- lapply(unique(c2$term), function(x){print(x);c2$gene[c2$term == x]})
names(geneSets) <- unique(c2$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)

geneSet <- "DN"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
sc2$AUC <- aucs
df<- data.frame(sc2@meta.data, sc2@reductions$umap@cell.embeddings)
colnames(df)
class_avg <- df %>%
  group_by(seurat_clusters) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

ggplot(df, aes(umap_1, umap_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="E") +
  ggrepel::geom_label_repel(aes(label = seurat_clusters),
                            data = class_avg,
                            size = 5,
                            label.size = 1,
                            segment.color = NA
  )+   theme(legend.position = "none") + theme_bw() 


library(ggplot2)
library(viridis)
library(ggrepel)

# åˆ›å»ºä¸€ä¸ªæ–°åˆ—ï¼Œå°? "Type" è½¬æ¢ä¸ºå› å­å˜é‡?
df$Type <- factor(df$Type)
library(ggplot2)
library(viridis)
library(ggrepel)



# åˆ›å»ºä¸€ä¸ªæ–°åˆ—ï¼Œå°? "Type" è½¬æ¢ä¸ºå› å­å˜é‡?
df$Type <- factor(df$Type)

# åˆ†åˆ«ç»˜åˆ¶ä¸¤ä¸ªå­å›¾
gg <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = AUC)) + 
  viridis::scale_color_viridis(option="E") +
  ggrepel::geom_label_repel(
    aes(label = ""),  # è®¾ç½®ä¸€ä¸ªç©ºå­—ç¬¦ä¸²ä½œä¸? label æ˜ å°„
    data = class_avg,
    size = 5,
    label.size = 1,
    segment.color = NA
  ) +
  theme(legend.position = "none") + 
  theme_bw() +
  facet_wrap(~ Type, ncol = 1)  # æ ¹æ® "Type" åˆ—åˆ†é¢ï¼Œæ¯åˆ—ä¸€ä¸ªå­å›?

print(gg)

library(ggplot2)
library(viridis)
library(ggrepel)

# åˆ›å»ºä¸€ä¸ªæ–°åˆ—ï¼Œå°? "Type" è½¬æ¢ä¸ºå› å­å˜é‡?
df$Type <- factor(df$Type)

# åˆ†åˆ«ç»˜åˆ¶ä¸¤ä¸ªå­å›¾
gg <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = AUC)) + 
  viridis::scale_color_viridis(option="E", limits = c(0, 0.1)) +  # è®¾ç½® AUC å°ºåº¦èŒƒå›´
  ggrepel::geom_label_repel(
    aes(label = ""),  # è®¾ç½®ä¸€ä¸ªç©ºå­—ç¬¦ä¸²ä½œä¸? label æ˜ å°„
    data = class_avg,
    size = 5,
    label.size = 1,
    segment.color = NA
  ) +
  theme(legend.position = "none") + 
  theme_bw() +
  facet_wrap(~ Type, ncol = 1)  # æ ¹æ® "Type" åˆ—åˆ†é¢ï¼Œæ¯åˆ—ä¸€ä¸ªå­å›?

print(gg)

library(ggplot2)

# È·±£ df ÖÐÓÐ seurat_clusters ºÍ AUC ÁÐ
# »æÖÆÐ¡ÌáÇÙÍ¼
ggplot(df, aes(x = factor(cell_identify), y = AUC)) +
  geom_violin(trim = FALSE, fill = "lightblue") + 
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  theme_bw() +
  xlab("cell_identify") +
  ylab("AUC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggplot2)

library(ggplot2)

library(ggplot2)

library(ggplot2)

# ¶¨ÒåÑÕÉ«Ó³Éä
cell_colors <- c("PT" = "#1f77b4", "LOH" = "#ff7f0e", "DCT" = "#2ca02c", 
                 "PC" = "#d62728", "IC" = "#9467bd", "Podo"="#8c564b",
                 "Fib" = "#e377c2" , "Endo" = "#7f7f7f" , "T cells" = "#bcbd22", 
                 "NKT" = "#17becf", "B cells" = "#aec7e8", 
                 "Mac" = "#ffbb78", "Mono" = "#98df8a")

# »æÖÆÐ¡ÌáÇÙÍ¼£¬²¢Ó¦ÓÃÑÕÉ«Ó³Éä
ggplot(df, aes(x = factor(cell_identify), y = AUC, fill = factor(cell_identify))) +
  geom_violin(trim = FALSE, color = "black", size = 0.3) + 
  scale_fill_manual(values = cell_colors) + 
  theme_classic() +
  xlab("Cell Type") +
  ylab("AUC Value") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.line = element_line(size = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "none"
  )

