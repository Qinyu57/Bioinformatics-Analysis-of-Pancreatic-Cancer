# 安装monocle3 - 在终端安装完hdf5和cmake后运行此代码

# 步骤1：设置HDF5路径（Apple Silicon Mac用这个）
Sys.setenv(HDF5_DIR = "/opt/homebrew/opt/hdf5")
# 注意：如果是Intel Mac，改用这个路径：
# Sys.setenv(HDF5_DIR = "/usr/local/opt/hdf5")

# 步骤2：确保有devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# 步骤3：安装BPCells（monocle3的重要依赖）
cat("正在安装BPCells...\n")
devtools::install_github('bnprks/BPCells', upgrade = "never")

# 步骤4：安装monocle3
cat("\n正在安装monocle3...\n")
devtools::install_github('cole-trapnell-lab/monocle3', upgrade = "never")

# 步骤5：验证安装
cat("\n==================\n")
cat("验证安装结果：\n")
if (requireNamespace("monocle3", quietly = TRUE)) {
  library(monocle3)
  cat("✓ monocle3安装成功！\n")
  cat(paste("版本:", packageVersion("monocle3"), "\n"))
} else {
  cat("✗ 安装失败，请查看上面的错误信息\n")
}
cat("==================\n")

# ============================================================
# 轨迹分析 - 带细胞亚类标注的癌症vs健康对比版
# ============================================================

# 加载必需库
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(patchwork)

# ------------------------------------------------------------
# 参数设置区（可修改）
# ------------------------------------------------------------

# 目标细胞类型
target_celltype <- "other"

# 质量筛选阈值
confidence_threshold <- 0.7      
mapping_score_threshold <- 0.5   

# 轨迹分析参数
n_pcs <- 30
n_neighbors <- 15
min_dist <- 0.3

# 输出路径
output_dir <- paste0("~/Desktop/基因癌症预测/monocle3/", 
                     gsub(" ", "_", target_celltype), "_comparison_with_subtypes/")

# ------------------------------------------------------------
# 检查输入对象
# ------------------------------------------------------------

if (!exists("PBMC.Integrated1")) {
  stop("❌ 未找到 PBMC.Integrated1 对象，请先运行主分析脚本")
}

cat("✓ 检测到 PBMC.Integrated1 对象\n")
cat("总细胞数:", ncol(PBMC.Integrated1), "\n")
cat("可用细胞类型:\n")
print(table(PBMC.Integrated1$predicted.celltype.l1))
cat("\n")

# ------------------------------------------------------------
# 1. 筛选目标细胞
# ------------------------------------------------------------

cat("=== 筛选", target_celltype, "细胞 ===\n")

cells_of_interest <- colnames(PBMC.Integrated1)[
  PBMC.Integrated1$predicted.celltype.l1 == target_celltype
]
cat(target_celltype, "总数:", length(cells_of_interest), "\n")

has_confidence <- "predicted.celltype.l1.score" %in% colnames(PBMC.Integrated1@meta.data)
has_mapping <- "mapping.score" %in% colnames(PBMC.Integrated1@meta.data)

if (has_confidence && has_mapping) {
  cells_to_use <- colnames(PBMC.Integrated1)[
    PBMC.Integrated1$predicted.celltype.l1 == target_celltype &
      PBMC.Integrated1$predicted.celltype.l1.score > confidence_threshold &
      PBMC.Integrated1$mapping.score > mapping_score_threshold
  ]
  cat("筛选条件:\n")
  cat("  置信度 >", confidence_threshold, "\n")
  cat("  Mapping score >", mapping_score_threshold, "\n")
} else if (has_confidence) {
  cells_to_use <- colnames(PBMC.Integrated1)[
    PBMC.Integrated1$predicted.celltype.l1 == target_celltype &
      PBMC.Integrated1$predicted.celltype.l1.score > confidence_threshold
  ]
  cat("⚠️  仅使用置信度筛选\n")
} else {
  cells_to_use <- cells_of_interest
  cat("⚠️  无质量筛选\n")
}

cat("高质量细胞数:", length(cells_to_use), "\n")

if (length(cells_to_use) < 50) {
  stop("❌ 高质量细胞数量不足(<50)，请降低阈值或更换细胞类型")
}

seurat_subset <- subset(PBMC.Integrated1, cells = cells_to_use)
cat("\n条件分布:\n")
print(table(seurat_subset$condition))

# 检查并显示亚类信息
if ("predicted.celltype.l2" %in% colnames(seurat_subset@meta.data)) {
  cat("\n亚类分布:\n")
  print(table(seurat_subset$predicted.celltype.l2))
  has_subtype <- TRUE
} else {
  cat("\n⚠️  未找到 predicted.celltype.l2，将仅使用主类型\n")
  has_subtype <- FALSE
}
cat("\n")

# ------------------------------------------------------------
# 2. 转换为 CDS 对象
# ------------------------------------------------------------

cat("=== 创建 Monocle3 CDS 对象 ===\n")

DefaultAssay(seurat_subset) <- "RNA"

expression_matrix <- GetAssayData(seurat_subset, layer = "counts")
cell_metadata <- seurat_subset@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cat("✓ CDS 创建成功\n")
cat("细胞数:", ncol(cds), "| 基因数:", nrow(cds), "\n\n")

rm(expression_matrix, cell_metadata, gene_metadata)
gc()

# ------------------------------------------------------------
# 3. 预处理和轨迹学习
# ------------------------------------------------------------

cat("=== 预处理和降维 ===\n")

cat("1/4 预处理...\n")
cds <- preprocess_cds(cds, num_dim = n_pcs, method = "PCA")

cat("2/4 批次校正...\n")
cds <- align_cds(cds, alignment_group = "orig.ident")

cat("3/4 UMAP降维...\n")
cds <- reduce_dimension(cds,
                        reduction_method = "UMAP",
                        preprocess_method = "Aligned",
                        umap.n_neighbors = n_neighbors,
                        umap.min_dist = min_dist)

cat("4/4 聚类...\n")
cds <- cluster_cells(cds, resolution = 1e-3)

cat("✓ 预处理完成\n\n")

cat("=== 学习细胞轨迹 ===\n")
cds <- learn_graph(cds, use_partition = TRUE)
cat("✓ 轨迹学习完成\n\n")

# ------------------------------------------------------------
# 4. 癌症vs健康轨迹对比分析（带亚类）
# ------------------------------------------------------------

cat("=== 癌症vs健康轨迹对比分析（带亚类标注）===\n")

# 检查细胞数量
cancer_cells <- colnames(cds)[colData(cds)$condition == "Cancer"]
healthy_cells <- colnames(cds)[colData(cds)$condition == "Healthy"]

cat("Cancer 细胞数:", length(cancer_cells), "\n")
cat("Healthy 细胞数:", length(healthy_cells), "\n")

if (length(cancer_cells) < 30 | length(healthy_cells) < 30) {
  stop("❌ 至少一组细胞数不足30个，无法进行可靠分析")
}

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 分别处理癌症和健康组的函数
process_subset <- function(cds_subset, group_name) {
  cat("处理", group_name, "组...\n")
  
  cds_subset <- preprocess_cds(cds_subset, num_dim = n_pcs)
  cds_subset <- align_cds(cds_subset, alignment_group = "orig.ident")
  cds_subset <- reduce_dimension(cds_subset,
                                 umap.n_neighbors = n_neighbors,
                                 umap.min_dist = min_dist)
  cds_subset <- cluster_cells(cds_subset, resolution = 1e-3)
  cds_subset <- learn_graph(cds_subset, use_partition = TRUE)
  
  return(cds_subset)
}

# 处理癌症组
cds_cancer <- cds[, cancer_cells]
cds_cancer <- process_subset(cds_cancer, "Cancer")

# 处理健康组  
cds_healthy <- cds[, healthy_cells]
cds_healthy <- process_subset(cds_healthy, "Healthy")

# ------------------------------------------------------------
# 生成对比图（带亚类标注）
# ------------------------------------------------------------

# 如果有亚类信息，创建亚类对比图
if (has_subtype) {
  cat("\n=== 生成带亚类标注的轨迹图 ===\n")
  
  # 图1: 按亚类着色的并排对比
  p_cancer_subtype <- plot_cells(cds_cancer,
                                 color_cells_by = "predicted.celltype.l2",
                                 label_cell_groups = TRUE,
                                 label_leaves = TRUE,
                                 label_branch_points = TRUE,
                                 label_roots = TRUE,
                                 graph_label_size = 3,
                                 group_label_size = 3.5,
                                 cell_size = 1.2,
                                 alpha = 0.8) +
    ggtitle(paste(target_celltype, "Subtypes - Cancer")) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, color = "#E64B35"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
  
  p_healthy_subtype <- plot_cells(cds_healthy,
                                  color_cells_by = "predicted.celltype.l2",
                                  label_cell_groups = TRUE,
                                  label_leaves = TRUE,
                                  label_branch_points = TRUE,
                                  label_roots = TRUE,
                                  graph_label_size = 3,
                                  group_label_size = 3.5,
                                  cell_size = 1.2,
                                  alpha = 0.8) +
    ggtitle(paste(target_celltype, "Subtypes - Healthy")) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, color = "#4DBBD5"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
  
  # 并排亚类对比图
  p_subtype_side_by_side <- p_cancer_subtype | p_healthy_subtype +
    plot_annotation(
      title = paste(target_celltype, "Trajectory with Subtypes: Cancer vs Healthy"),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
    )
  
  ggsave(paste0(output_dir, "01_trajectory_subtypes_side_by_side.pdf"),
         p_subtype_side_by_side, width = 18, height = 8)
  cat("✓ 亚类并排对比图已保存\n")
  
  # 图2: 重叠对比图（按亚类着色）
  cat("生成亚类重叠对比图...\n")
  
  cds_combined_comparison <- cds
  colData(cds_combined_comparison)$comparison_group <- ifelse(
    colnames(cds_combined_comparison) %in% cancer_cells, "Cancer", "Healthy"
  )
  
  p_overlay_subtype <- plot_cells(cds_combined_comparison,
                                  color_cells_by = "predicted.celltype.l2",
                                  label_cell_groups = TRUE,
                                  label_leaves = TRUE,
                                  label_branch_points = TRUE,
                                  label_roots = TRUE,
                                  graph_label_size = 3.5,
                                  group_label_size = 3.5,
                                  cell_size = 1,
                                  cell_stroke = 0.3,
                                  alpha = 0.7) +
    ggtitle(paste(target_celltype, "Trajectory Overlay with Subtypes")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
  
  ggsave(paste0(output_dir, "02_trajectory_subtypes_overlay.pdf"),
         p_overlay_subtype, width = 14, height = 10)
  cat("✓ 亚类重叠对比图已保存\n")
  
  # 图3: 分面对比图（按亚类着色，按condition分面）
  cat("生成亚类分面对比图...\n")
  
  p_facet_subtype <- plot_cells(cds_combined_comparison,
                                color_cells_by = "predicted.celltype.l2",
                                label_cell_groups = TRUE,
                                label_leaves = FALSE,
                                label_branch_points = FALSE,
                                label_roots = FALSE,
                                graph_label_size = 0,
                                group_label_size = 3,
                                cell_size = 0.8,
                                alpha = 0.7) +
    facet_wrap(~comparison_group, ncol = 2) +
    ggtitle(paste(target_celltype, "Trajectory Facet with Subtypes")) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", size = 12),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
  
  ggsave(paste0(output_dir, "03_trajectory_subtypes_facet.pdf"),
         p_facet_subtype, width = 16, height = 7)
  cat("✓ 亚类分面对比图已保存\n")
  
  # 图4: 双重分面图（按亚类分面）
  cat("生成亚类双重分面图...\n")
  
  # 只保留最常见的几个亚类以保持清晰
  subtype_counts <- table(colData(cds_combined_comparison)$predicted.celltype.l2)
  top_subtypes <- names(sort(subtype_counts, decreasing = TRUE)[1:min(6, length(subtype_counts))])
  
  cds_top_subtypes <- cds_combined_comparison[, 
                                              colData(cds_combined_comparison)$predicted.celltype.l2 %in% top_subtypes]
  
  p_double_facet <- plot_cells(cds_top_subtypes,
                               color_cells_by = "comparison_group",
                               label_cell_groups = FALSE,
                               label_leaves = FALSE,
                               label_branch_points = FALSE,
                               label_roots = FALSE,
                               cell_size = 0.8,
                               alpha = 0.7) +
    facet_wrap(~predicted.celltype.l2, ncol = 3) +
    scale_color_manual(values = c("Cancer" = "#E64B35", "Healthy" = "#4DBBD5"),
                       name = "Condition") +
    ggtitle(paste(target_celltype, "Trajectory by Subtype (Top 6)")) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", size = 10),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
  
  ggsave(paste0(output_dir, "04_trajectory_by_subtype_facet.pdf"),
         p_double_facet, width = 15, height = 10)
  cat("✓ 亚类双重分面图已保存\n")
}

# 原始对比图（按cluster着色）
cat("\n=== 生成原始聚类对比图 ===\n")

p_cancer <- plot_cells(cds_cancer,
                       color_cells_by = "cluster",
                       label_cell_groups = TRUE,
                       label_leaves = TRUE,
                       label_branch_points = TRUE,
                       label_roots = TRUE,
                       graph_label_size = 3,
                       group_label_size = 4,
                       cell_size = 1.2,
                       alpha = 0.8) +
  ggtitle(paste(target_celltype, "- Cancer (Clusters)")) +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, color = "#E64B35"),
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.title = element_blank()
  )

p_healthy <- plot_cells(cds_healthy,
                        color_cells_by = "cluster",
                        label_cell_groups = TRUE,
                        label_leaves = TRUE,
                        label_branch_points = TRUE,
                        label_roots = TRUE,
                        graph_label_size = 3,
                        group_label_size = 4,
                        cell_size = 1.2,
                        alpha = 0.8) +
  ggtitle(paste(target_celltype, "- Healthy (Clusters)")) +
  scale_color_brewer(palette = "Set2", name = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, color = "#4DBBD5"),
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.title = element_blank()
  )

p_side_by_side <- p_cancer | p_healthy +
  plot_annotation(
    title = paste(target_celltype, "Trajectory Clusters: Cancer vs Healthy"),
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
  )

ggsave(paste0(output_dir, "05_trajectory_clusters_side_by_side.pdf"),
       p_side_by_side, width = 16, height = 7)
cat("✓ 聚类并排对比图已保存\n")
