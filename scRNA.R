library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
scRNA <- readRDS("D:\\bioinfor\\sepsis\\DKD_identify.rds")

# 自定义颜色向量，确保28种颜色高对比且适合SCIENCE期刊风格
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", 
            "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff6347", 
            "#4682b4", "#32cd32", "#ffa500", "#a52a2a", "#daa520", "#8a2be2", "#7fff00", 
            "#dc143c", "#00ced1", "#9400d3", "#ff1493", "#00bfff", "#696969", "#ff4500")

# 使用自定义颜色绘制UMAP图
p <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", label = TRUE, 
             cols = colors, 
             label.size = 3.5,  # 调整标签字体大小
             pt.size = 0.6) +   # 调整数据点大小
  theme_minimal(base_size = 14) +  # 使用简洁的主题
  theme(legend.position = "right",  # 图例放置在右侧
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white", color = NA))

# 如果标签与点重叠，使用geom_text_repel来避免重叠
library(ggrepel)
p <- p + geom_text_repel(aes(label = seurat_clusters), size = 3.5)

# 显示图形
print(p)

library(ggplot2)
library(ggplot2)

# 创建 DimPlot
p <- DimPlot(scRNA, reduction = "umap", group.by = "cell_identify", label = TRUE, 
             cols = c("PT" = "#1f77b4", "LOH" = "#ff7f0e", "DCT" = "#2ca02c", 
                      "PC" = "#d62728", "IC" = "#9467bd", "Podo"="#8c564b",
                      "Fib" = "#e377c2" , "Endo" = "#7f7f7f" , "T cells" = "#bcbd22", "NKT" = "#17becf",
                      "B cells" = "#aec7e8", "Mac" = "#ffbb78", "Mono" = "#98df8a"
                      ),
             label.size = 4, 
             pt.size = 0.5)

# 修改点的透明度
p <- p + scale_alpha(range = c(0, 0))

# 进一步调整图形的美观度
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

# 计算每个组别中每种细胞类型的数量
cell_counts <- scRNA@meta.data %>%
  group_by(sample_type, cell_identify) %>%
  summarise(count = n()) %>%
  group_by(sample_type) %>%
  mutate(percentage = count / sum(count) * 100)  # 计算每个细胞类型在每个组中的比例

# 绘制堆叠条形图，显示每个组别中每种细胞类型的比例
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

# 保存图形时设置高度和宽度
ggsave("Cell_Composition.pdf", plot = last_plot(), width = 5, height = 6, dpi = 300)
ggplot(cell_counts, aes(x = sample_type, y = percentage, fill = cell_identify)) +
  geom_bar(stat = "identity", position = "fill", width = 0.5) +  # 将 width 设置为0.5
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
            cols = c("lightgrey", "skyblue", "dodgerblue", "firebrick"), # 使用炫酷的颜色梯度
            pt.size = 1.5, # 增大点的大小
            order = TRUE) + # 保证高表达的点在上层显示
  theme_void() + # 去除背景
  theme(legend.position = "bottom", # 图例位置设置在底部
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "darkred"), # 设置炫酷的标题样式
        legend.text = element_text(size = 12, color = "darkblue"), # 设置图例文本样式
        legend.key.size = unit(1.2, "cm"), # 调整图例大小
        plot.margin = margin(10, 10, 10, 10)) # 调整图形的边距
Idents(scRNA) <- scRNA$sample_type
setwd("D:\\bioinfor\\sepsis\\2024.8.31")
Idents(scRNA)= "cell_identify"
DimPlot(scRNA)
degs <- FindAllMarkers(scRNA, logfc.threshold = 0.5,
                       test.use = "roc", 
                       return.thresh = 0.25, 
                       min.pct = 0.3, only.pos = T) 
write.csv(degs, file = "degs1.csv")




###基因的小提琴图
scRNA <- readRDS("D:\\bioinfor\\sepsis\\DKD_identify.rds")
features = c("GATM", "PCBD1", "F11", "HRSP12", "G6PC")
VlnPlot(scRNA, features = features, pt.size = 0, ncol = 2)+
  scale_x_discrete("")+
  theme(
    axis.text.x.bottom = element_blank()
  )

# 从差异分析结果中提取你感兴趣的基因的显著性信息
degs_of_interest <- degs[degs$gene %in% features,]
library(Seurat)
library(ggplot2)
library(ggpubr)

# 生成小提琴图
p <- VlnPlot(scRNA, features = features, pt.size = 0, ncol = 2) +
  scale_x_discrete("") +
  theme(axis.text.x.bottom = element_blank())

# 添加显著性标记
p + stat_compare_means(aes(group = sample_type), 
                       label = "p.signif", 
                       method = "wilcox.test", 
                       hide.ns = TRUE, 
                       comparisons = list(c("DKD", "HC")), 
                       label.y = c(1.5, 2, 2.5, 3, 3.5))  # 根据数据调整 y 轴位置
VlnPlot(object=scRNA,features = features,slot = "data",group.by = 'cell_identify',split.by = 'sample_type',pt.size = 0)+stat_compare_means(method = 't.test')+stat_boxplot()
VlnPlot(object=scRNA,features = features,slot = "data",group.by = 'sample_type',pt.size = 0)+stat_compare_means(method = 't.test')+stat_boxplot()
getwd()
VlnPlot(object=scRNA, features = features, slot = "data", group.by = 'cell_identify', pt.size = 0) +
  stat_compare_means(method = 'wilcox.test', label = "p.format", hide.ns = F)
# 将 cell_identify 设为活动标识符
scRNA <- SetIdent(scRNA, value = "cell_identify")

# 自定义颜色映射
cell_colors <- c("PT" = "#1f77b4", "LOH" = "#ff7f0e", "DCT" = "#2ca02c", 
                 "PC" = "#d62728", "IC" = "#9467bd", "Podo"="#8c564b",
                 "Fib" = "#e377c2" , "Endo" = "#7f7f7f" , "T cells" = "#bcbd22", 
                 "NKT" = "#17becf", "B cells" = "#aec7e8", 
                 "Mac" = "#ffbb78", "Mono" = "#98df8a")

# 绘制小提琴图
VlnPlot(object = scRNA, features = "G6PC", slot = "data", pt.size = 0)  +
  scale_fill_manual(values = cell_colors) +
  theme_minimal()


library(ggplot2)
library(dplyr)

# 假设你有一个包含基因表达和p值的数据框df，形如:
 df <- data.frame(
   cell_identify = rep(unique(scRNA$cell_identify), each = length(features)),
   gene = rep(features, times = length(unique(scRNA$cell_identify))),
   expression = runif(length(features) * length(unique(scRNA$cell_identify))),
   p_value = runif(length(features) * length(unique(scRNA$cell_identify)))
 )
 # 创建气泡图
 ggplot(df, aes(x = cell_identify, y = gene)) +
   geom_point(aes(size = expression, color = -log10(p_value))) +  # 用表达值控制气泡大小，p值控制颜色
   scale_size(range = c(1, 10)) +  # 调整气泡大小范围
   scale_color_viridis_c(option = "E") +  # 使用 Viridis 配色方案
   theme_minimal() +
   xlab("cell_identify") +
   ylab("gene") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # 调整x轴标签角度
         axis.text.y = element_text(size = 10))  # 调整y轴标签大小
# 自定义Science风格的配色
library(ggplot2)
library(ggpubr)

library(ggplot2)
library(ggpubr)

VlnPlot(object = scRNA, features = "PCBD1", slot = "data", group.by = 'sample_type', pt.size = 0) +
  stat_compare_means(method = 'wilcox.test', label = "p.format", hide.ns = FALSE) +
  scale_fill_manual(values = c("DKD(CKD)" = "#E6846D", "Health kidney" = "#8DCDD5"))  # 修改颜色


library(ggplot2)
library(ggpubr)

####AUC
memory.limit(32000000)
##绯荤粺鎶ラ敊鏀逛负鑻辨枃
Sys.setenv(LANGUAGE = "en")
##绂佹杞寲涓哄洜瀛?
options(stringsAsFactors = FALSE)
##娓呯┖鐜
rm(list=ls())

setwd("D:\\bioinfor\\sepsis\\DN_old_mito\\28.AUCell")
###鍔犺浇鎵�闇�瑕佺殑鍖?
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

# 鍒涘缓涓�涓柊鍒楋紝灏? "Type" 杞崲涓哄洜瀛愬彉閲?
df$Type <- factor(df$Type)
library(ggplot2)
library(viridis)
library(ggrepel)



# 鍒涘缓涓�涓柊鍒楋紝灏? "Type" 杞崲涓哄洜瀛愬彉閲?
df$Type <- factor(df$Type)

# 鍒嗗埆缁樺埗涓や釜瀛愬浘
gg <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = AUC)) + 
  viridis::scale_color_viridis(option="E") +
  ggrepel::geom_label_repel(
    aes(label = ""),  # 璁剧疆涓�涓┖瀛楃涓蹭綔涓? label 鏄犲皠
    data = class_avg,
    size = 5,
    label.size = 1,
    segment.color = NA
  ) +
  theme(legend.position = "none") + 
  theme_bw() +
  facet_wrap(~ Type, ncol = 1)  # 鏍规嵁 "Type" 鍒楀垎闈紝姣忓垪涓�涓瓙鍥?

print(gg)

library(ggplot2)
library(viridis)
library(ggrepel)

# 鍒涘缓涓�涓柊鍒楋紝灏? "Type" 杞崲涓哄洜瀛愬彉閲?
df$Type <- factor(df$Type)

# 鍒嗗埆缁樺埗涓や釜瀛愬浘
gg <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = AUC)) + 
  viridis::scale_color_viridis(option="E", limits = c(0, 0.1)) +  # 璁剧疆 AUC 灏哄害鑼冨洿
  ggrepel::geom_label_repel(
    aes(label = ""),  # 璁剧疆涓�涓┖瀛楃涓蹭綔涓? label 鏄犲皠
    data = class_avg,
    size = 5,
    label.size = 1,
    segment.color = NA
  ) +
  theme(legend.position = "none") + 
  theme_bw() +
  facet_wrap(~ Type, ncol = 1)  # 鏍规嵁 "Type" 鍒楀垎闈紝姣忓垪涓�涓瓙鍥?

print(gg)

library(ggplot2)

# 确保 df 中有 seurat_clusters 和 AUC 列
# 绘制小提琴图
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

# 定义颜色映射
cell_colors <- c("PT" = "#1f77b4", "LOH" = "#ff7f0e", "DCT" = "#2ca02c", 
                 "PC" = "#d62728", "IC" = "#9467bd", "Podo"="#8c564b",
                 "Fib" = "#e377c2" , "Endo" = "#7f7f7f" , "T cells" = "#bcbd22", 
                 "NKT" = "#17becf", "B cells" = "#aec7e8", 
                 "Mac" = "#ffbb78", "Mono" = "#98df8a")

# 绘制小提琴图，并应用颜色映射
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

