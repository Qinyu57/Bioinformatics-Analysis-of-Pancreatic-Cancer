options(future.globals.maxSize = 64000 * 1024^2)  # 64GB
options(future.rng.onMisuse = "ignore")
Sys.setenv('R_MAX_VSIZE'=64000000000)  # 64GB向量内存
invisible(gc())
install.packages("hdf5r", type = "binary")
if (!requireNamespace("harmony", quietly = TRUE)) {
  install.packages("harmony")
}
library(harmony)
library(hdf5r)
library(Seurat)
library(ggplot2)
library(ggrepel)
# 创建函数来加载数据并立即删除原始.data对象
load_and_create_seurat <- function(data_dir, project_name, is_h5 = FALSE) {
  if (is_h5) {
    data <- Read10X_h5(data_dir)
  } else {
    data <- Read10X(data.dir = data_dir)
  }
  seurat_obj <- CreateSeuratObject(counts = data, project = project_name)
  rm(data)  # 删除原始数据以释放内存
  gc()      # 强制垃圾回收
  return(seurat_obj)
}



# 加载PDAC_PBMC样本

PDAC_PBMC_1_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_1_QC/", "PDAC_PBMC_1_QC")
PDAC_PBMC_2_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_2_QC/", "PDAC_PBMC_2_QC")
PDAC_PBMC_3_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_3_QC/", "PDAC_PBMC_3_QC")
PDAC_PBMC_4_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_4_QC/", "PDAC_PBMC_4_QC")
PDAC_PBMC_5_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_5_QC/", "PDAC_PBMC_5_QC")
PDAC_PBMC_6_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_6_QC/", "PDAC_PBMC_6_QC")
PDAC_PBMC_7_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_7_QC/", "PDAC_PBMC_7_QC")
PDAC_PBMC_8_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_8_QC/", "PDAC_PBMC_8_QC")
PDAC_PBMC_9_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_9_QC/", "PDAC_PBMC_9_QC")
PDAC_PBMC_10A_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_10A_QC/", "PDAC_PBMC_10A_QC")
PDAC_PBMC_10B_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_10B_QC/", "PDAC_PBMC_10B_QC")
PDAC_PBMC_11_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_11_QC/", "PDAC_PBMC_11_QC")
PDAC_PBMC_12_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_12_QC/", "PDAC_PBMC_12_QC")
PDAC_PBMC_13_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_13_QC/", "PDAC_PBMC_13_QC")
PDAC_PBMC_14_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_14_QC/", "PDAC_PBMC_14_QC")
PDAC_PBMC_15_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_15_QC/", "PDAC_PBMC_15_QC")
PDAC_PBMC_16_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/PDAC_PBMC_16_QC/", "PDAC_PBMC_16_QC")

PDAC_PBMC_1_QC$condition <-"Cancer"
PDAC_PBMC_2_QC$condition <-"Cancer"
PDAC_PBMC_3_QC$condition <-"Cancer"
PDAC_PBMC_4_QC$condition <-"Cancer"
PDAC_PBMC_5_QC$condition <-"Cancer"
PDAC_PBMC_6_QC$condition <-"Cancer"
PDAC_PBMC_7_QC$condition <-"Cancer"
PDAC_PBMC_8_QC$condition <-"Cancer"
PDAC_PBMC_9_QC$condition <-"Cancer"
PDAC_PBMC_10A_QC$condition <-"Cancer"
PDAC_PBMC_10B_QC$condition <-"Cancer"
PDAC_PBMC_11_QC$condition <-"Cancer"
PDAC_PBMC_12_QC$condition <-"Cancer"
PDAC_PBMC_13_QC$condition <-"Cancer"
PDAC_PBMC_14_QC$condition <-"Cancer"
PDAC_PBMC_15_QC$condition <-"Cancer"
PDAC_PBMC_16_QC$condition <-"Cancer"


# 加载Healthy_PBMC样本
Healthy_PBMC_1_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/Healthy_PBMC_1_QC/", "Healthy_PBMC_1_QC")
Healthy_PBMC_2_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/Healthy_PBMC_2_QC/", "Healthy_PBMC_2_QC")
Healthy_PBMC_3_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/Healthy_PBMC_3_QC/", "Healthy_PBMC_3_QC")
Healthy_PBMC_4_QC <- load_and_create_seurat("~/Desktop/基因癌症预测/QC/Healthy_PBMC_4_QC/", "Healthy_PBMC_4_QC")

Healthy_PBMC_1_QC$condition <-"Healthy"
Healthy_PBMC_2_QC$condition <-"Healthy"
Healthy_PBMC_3_QC$condition <-"Healthy"
Healthy_PBMC_4_QC$condition <-"Healthy"


rm(list = ls(pattern = "\\.data$"))  # 删除所有.data对象
invisible(gc())  # 强制垃圾回收

# 合并所有样本


PBMC.merged <- merge(PDAC_PBMC_1_QC, y=c(PDAC_PBMC_2_QC,PDAC_PBMC_3_QC,PDAC_PBMC_4_QC,PDAC_PBMC_5_QC,PDAC_PBMC_6_QC,PDAC_PBMC_7_QC,PDAC_PBMC_8_QC,PDAC_PBMC_9_QC,PDAC_PBMC_10A_QC,PDAC_PBMC_10B_QC,PDAC_PBMC_11_QC,PDAC_PBMC_12_QC,PDAC_PBMC_13_QC,PDAC_PBMC_14_QC,PDAC_PBMC_15_QC,PDAC_PBMC_16_QC,Healthy_PBMC_1_QC,Healthy_PBMC_2_QC,Healthy_PBMC_3_QC,Healthy_PBMC_4_QC),
                     add.cell.ids=c("PDAC_PBMC_1_QC","PDAC_PBMC_2_QC","PDAC_PBMC_3_QC","PDAC_PBMC_4_QC","PDAC_PBMC_5_QC","PDAC_PBMC_6_QC","PDAC_PBMC_7_QC","PDAC_PBMC_8_QC","PDAC_PBMC_9_QC","PDAC_PBMC_10A_QC","PDAC_PBMC_10B_QC","PDAC_PBMC_11_QC","PDAC_PBMC_12_QC","PDAC_PBMC_13_QC","PDAC_PBMC_14_QC","PDAC_PBMC_15_QC","PDAC_PBMC_16_QC","Healthy_PBMC_1_QC","Healthy_PBMC_2_QC","Healthy_PBMC_3_QC","Healthy_PBMC_4_QC"),
                     project = "SampleIntegrated")


# 删除单个样本对象以释放内存
rm(PDAC_PBMC_1_QC, PDAC_PBMC_2_QC, PDAC_PBMC_3_QC, PDAC_PBMC_4_QC, PDAC_PBMC_5_QC, PDAC_PBMC_6_QC, PDAC_PBMC_7_QC, PDAC_PBMC_8_QC, PDAC_PBMC_9_QC, PDAC_PBMC_10A_QC, PDAC_PBMC_10B_QC, PDAC_PBMC_11_QC, PDAC_PBMC_12_QC, PDAC_PBMC_13_QC, PDAC_PBMC_14_QC, PDAC_PBMC_15_QC, PDAC_PBMC_16_QC,
   Healthy_PBMC_1_QC, Healthy_PBMC_2_QC, Healthy_PBMC_3_QC, Healthy_PBMC_4_QC)
gc()  # 强制垃圾回收


table(PBMC.merged$orig.ident)

# 使用更少的内存密集型设置
invisible(gc())

PBMC.merged <- NormalizeData(PBMC.merged,verbose = FALSE)
invisible(gc())

PBMC.merged <- FindVariableFeatures(PBMC.merged,verbose = FALSE)
invisible(gc())

PBMC.merged <- ScaleData(PBMC.merged,verbose = FALSE)
invisible(gc())
# 在运行PCA时使用更少的特征和细胞（如果需要）

#Sample.merged <- RunPCA(Sample.merged)
PBMC.merged <- RunPCA(PBMC.merged, 
                      npcs = 30,  # 只计算前30个主成分
                      verbose = FALSE)
invisible(gc())

PBMC.Integrated1 <- IntegrateLayers(object = PBMC.merged,
                                   method = HarmonyIntegration,
                                   orig.reduction = 'pca',
                                   new.reduction = "harmony",
                                   verbose = FALSE)

# 删除merged对象以释放内存

rm(PBMC.merged)
gc()


PBMC.Integrated1 <- FindNeighbors(PBMC.Integrated1,
                                 reduction = 'harmony',
                                 dims = 1:30,
                                 verbose = FALSE)

PBMC.Integrated1 <- FindClusters(PBMC.Integrated1,
                                resolution = 2,
                                cluster.name = "harmonyCluster",
                                verbose = FALSE)


PBMC.Integrated1 <- RunUMAP(PBMC.Integrated1,
                           reduction = "harmony",
                           dims = 1:30,
                           reduction.name = "harmonyUMAP",
                           verbose = FALSE)


# 保存harmony integration的UMAP图
umap_plot <- DimPlot(PBMC.Integrated1,
                     reduction = "harmonyUMAP",
                     group.by = "orig.ident",
                     label = FALSE,
                     repel = TRUE) +
  ggtitle("PBMC Samples - Harmony Integration") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# 保存UMAP图
umap_file_name <- "~/Desktop/基因癌症预测/UMAP图/PBMC_harmony_integration_umap.pdf"
ggsave(
  filename = umap_file_name,
  plot = umap_plot,
  width = 12,
  height = 8,
  device = "pdf",
  bg = "white"
)



prep_for_azimuth <- function(obj, filename = "~/Desktop/Tissue_For_Azimuth.rds") {
  library(Seurat)
  library(SeuratObject)
  
  # 检查 RNA assay
  if (!"RNA" %in% Assays(obj)) {
    stop("❌ 没有 RNA assay，Azimuth 无法识别。")
  }
  DefaultAssay(obj) <- "RNA"
  
  # 如果是 Assay5 (Seurat v5)
  if (inherits(obj[["RNA"]], "Assay5")) {
    message("🔄 检测到 Seurat v5 Assay5，正在转换为 v4 兼容 Assay...")
    
    # 合并 layers
    assay_v5 <- SeuratObject::JoinLayers(obj[["RNA"]])
    
    # 提取 counts / data
    counts_mat <- GetAssayData(assay_v5, layer = "counts")
    data_mat   <- GetAssayData(assay_v5, layer = "data")
    
    # 创建 v4 Assay（只能传 counts）
    assay_v4 <- CreateAssayObject(counts = counts_mat)
    
    # 如果 data 存在，直接赋值到 slot
    if (!is.null(data_mat) && nrow(data_mat) > 0) {
      assay_v4@data <- data_mat
    }
    
    # 替换
    obj[["RNA"]] <- assay_v4
  } else {
    message("✅ RNA assay 已经是 v4 Assay，无需转换")
  }
  
  # 瘦身
  message("📦 正在用 DietSeurat 瘦身对象 ...")
  obj <- DietSeurat(
    obj,
    assays = "RNA",
    dimreducs = "pca"
  )
  
  # 保存
  saveRDS(obj, filename)
  full_path <- normalizePath(filename)
  size_mb <- file.info(full_path)$size / (1024^2)
  message("💾 已保存到: ", full_path)
  message("📏 文件大小: ", round(size_mb, 2), " MB")
  
  invisible(obj)
}

PBMC.Integrated1 <- prep_for_azimuth(
  PBMC.Integrated1,
  filename = "~/Desktop/PBMC_For_Azimuth.rds"
)

#完成Azimuth后
predictions1 <- read.delim('~/Desktop/基因癌症预测/high_quality_PBMC.tsv', row.names = 1)
head(predictions1)
colnames(predictions1)

# 添加细胞注释到Seurat对象
PBMC.Integrated1 <- AddMetaData(
  object = PBMC.Integrated1,
  metadata = predictions1
)

# 检查添加注释后的情况
cat("总细胞数:", ncol(PBMC.Integrated1), "\n")
cat("有注释的细胞数:", sum(!is.na(PBMC.Integrated1$predicted.celltype.l1)), "\n")
cat("NA细胞数:", sum(is.na(PBMC.Integrated1$predicted.celltype.l1)), "\n")

# 关键步骤：只保留有高质量注释的细胞
cells_with_annotation <- colnames(PBMC.Integrated1)[!is.na(PBMC.Integrated1$predicted.celltype.l1)]
cat("筛选后保留的细胞数:", length(cells_with_annotation), "\n")

# 筛选Seurat对象，只保留有注释的细胞
PBMC.Integrated1_filtered <- subset(PBMC.Integrated1, cells = cells_with_annotation)

cat("筛选后的细胞数:", ncol(PBMC.Integrated1_filtered), "\n")

# 验证是否还有NA
cat("筛选后NA细胞数:", sum(is.na(PBMC.Integrated1_filtered$predicted.celltype.l1)), "\n")

# 更新对象名称
PBMC.Integrated1 <- PBMC.Integrated1_filtered
rm(PBMC.Integrated1_filtered)  # 释放内存
gc()

head(PBMC.Integrated1@meta.data)
colnames(PBMC.Integrated1@meta.data)

# 首先检查是否有harmonyUMAP降维结果
if (!"harmonyUMAP" %in% names(PBMC.Integrated1@reductions)) {
  cat("harmonyUMAP降维结果不存在，检查可用降维结果...\n")
  print(names(PBMC.Integrated1@reductions))
  
  # 如果没有harmonyUMAP，检查是否有其他UMAP
  umap_reductions <- grep("umap|UMAP", names(PBMC.Integrated1@reductions), ignore.case = TRUE, value = TRUE)
  
  if (length(umap_reductions) > 0) {
    cat("使用现有的UMAP降维结果:", umap_reductions[1], "\n")
    reduction_to_use <- umap_reductions[1]
  } else {
    cat("没有UMAP降维结果，需要创建...\n")
    
    # 检查是否有harmony降维结果
    if ("harmony" %in% names(PBMC.Integrated1@reductions)) {
      cat("基于harmony降维创建UMAP...\n")
      PBMC.Integrated1 <- RunUMAP(PBMC.Integrated1,
                                  reduction = "harmony",
                                  dims = 1:30,
                                  reduction.name = "harmonyUMAP",
                                  verbose = FALSE)
      reduction_to_use <- "harmonyUMAP"
    } else if ("pca" %in% names(PBMC.Integrated1@reductions)) {
      cat("基于PCA创建UMAP...\n")
      PBMC.Integrated1 <- RunUMAP(PBMC.Integrated1,
                                  reduction = "pca",
                                  dims = 1:30,
                                  reduction.name = "harmonyUMAP",
                                  verbose = FALSE)
      reduction_to_use <- "harmonyUMAP"
    } else {
      stop("没有找到可用的降维结果来创建UMAP")
    }
  }
} else {
  cat("harmonyUMAP降维结果存在，直接使用...\n")
  reduction_to_use <- "harmonyUMAP"
}

# 检查细胞类型分布
table(PBMC.Integrated1$predicted.celltype.l1)

#画图 - 现在不会有NA细胞了
umap_celltype_legend <- DimPlot(PBMC.Integrated1,
                                reduction = reduction_to_use,
                                group.by = "predicted.celltype.l1",
                                label = FALSE,  # 不显示标签
                                pt.size = 0.1) +
  ggtitle("PBMC Samples - Cell Type l1 (High Quality)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

umap_legend_file <- "~/Desktop/基因癌症预测/UMAP图/PBMC_azimuth_celltype_l1_filtered.pdf"
ggsave(
  filename = umap_legend_file,
  plot = umap_celltype_legend,
  width = 12,
  height = 8,
  device = "pdf",
  bg = "white"
)


# 设置细胞类型变量
cell_type <- "HSPC"  # 改为你需要的细胞类型

# DE分析
PBMC.Integrated1$celltype.condition <- paste(PBMC.Integrated1$predicted.celltype.l2, PBMC.Integrated1$condition, sep = "_")
table(PBMC.Integrated1$celltype.condition)
Idents(PBMC.Integrated1) <- "celltype.condition"

# 使用变量自动生成ident名称
ident_1 <- paste0(cell_type, "_Cancer")
ident_2 <- paste0(cell_type, "_Healthy")

mono.de <- FindMarkers(PBMC.Integrated1, ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE)
head(mono.de, n = 10)

library(dplyr)

# 数据处理
x_lim_min <- -10       # X轴最小值
x_lim_max <- 10       # X轴最大值
y_lim_max <- 20
df <- mono.de %>%
  mutate(gene = rownames(mono.de)) %>%
  mutate(log2fc = avg_log2FC,
         pval = p_val_adj,
         diffexpressed = case_when(
           log2fc > 0.6 & pval < 0.05 ~ "Upregulated",
           log2fc < -0.6 & pval < 0.05 ~ "Downregulated",
           TRUE ~ "Not significant"
         ),
         delabel = ifelse(diffexpressed == "Not significant", NA, gene))

# 选取上下调里最显著的前10个基因
df_sig <- df %>%
  filter(diffexpressed != "Not significant", pval < 0.05, abs(log2fc) > 0.6) %>%
  group_by(diffexpressed) %>%
  arrange(pval, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

# 使用变量生成标题
plot_title <- paste0("Differential Gene Expression in ", tools::toTitleCase(cell_type), " Cells: Cancer vs Healthy")

# 火山图
myvolcanoplot <- ggplot(
  data = df, 
  aes(x = log2fc, y = -log10(pval), col = diffexpressed)
) +
  geom_point(size = 2.5, alpha = 0.8, shape = 16) +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.8) +
  
  scale_color_manual(
    values = c("Downregulated" = "#00AFBB", 
               "Not significant" = "grey70",
               "Upregulated" = "#bb0c00"),
    labels = c("Downregulated", "Not significant", "Upregulated"),
    guide = guide_legend(override.aes = list(shape = 16, size = 4))
  ) +
  coord_cartesian(ylim = c(0, y_lim_max), xlim = c(x_lim_min, x_lim_max)) +
  labs(
    color = '',
    x = expression("log"[2]*"FC"),
    y = expression("-log"[10]*"P"),
    title = plot_title  # 使用变量标题
  ) +
  scale_x_continuous(breaks = seq(x_lim_min, x_lim_max, by = ifelse((x_lim_max - x_lim_min) > 20, 5, 2))) + 
  
  geom_text_repel(
    data = df_sig,
    aes(label = delabel),
    size = 6,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.5,
    segment.alpha = 0.8,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0.2,
    force = 10,
    direction = "both",
    nudge_x = ifelse(df_sig$log2fc > 0, 1, -1),
    nudge_y = 1
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black")
  )

print(myvolcanoplot)

# 显著基因
significant_genes <- df %>%
  filter(diffexpressed != "Not significant", pval < 0.05, abs(log2fc) > 0.6) %>%
  pull(gene)

#gene_file_name1 <- paste0("~/Desktop/基因癌症预测/火山图/PBMC_new/l2/", cell_type, "_sig_genes_comma_separated.txt")
#writeLines(paste(significant_genes, collapse = ", "), gene_file_name1)

#up
up_genes <- df %>%
  filter(diffexpressed != "Not significant", pval < 0.05, log2fc > 0.6) %>%
  pull(gene)

#gene_file_name2 <- paste0("~/Desktop/基因癌症预测/火山图/PBMC_new/l2/", cell_type, "_up_genes_comma_separated.txt")
#writeLines(paste(up_genes, collapse = ", "), gene_file_name2)

#down
down_genes <- df %>%
  filter(diffexpressed != "Not significant", pval < 0.05, log2fc < -0.6) %>%
  pull(gene)

#gene_file_name3 <- paste0("~/Desktop/基因癌症预测/火山图/PBMC_new/l2/", cell_type, "_down_genes_comma_separated.txt")
#writeLines(paste(down_genes, collapse = ", "), gene_file_name3)

# 打印所有基因名称
all_genes <- df$gene

# 使用变量生成文件名
#gene_file_name <- paste0("~/Desktop/基因癌症预测/火山图/PBMC_new/l2/", cell_type, "_all_genes_comma_separated.txt")
#writeLines(paste(all_genes, collapse = ", "), gene_file_name)

# 使用变量生成图片文件名
plot_file_name <- paste0("~/Desktop/基因癌症预测/火山图/PBMC_new/l2/", cell_type, "_cancer_vs_healthy_volcano_bold.pdf")

ggsave(
  filename = plot_file_name,
  plot = myvolcanoplot,
  width = 10,
  height = 8,
  device = "pdf",
  bg = "white"
)
