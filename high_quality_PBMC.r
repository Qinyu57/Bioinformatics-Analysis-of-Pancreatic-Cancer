# 设置文件路径
input_file <- "~/Desktop/基因癌症预测/azimuth_PBMC.tsv"
output_dir <- "~/Desktop/基因癌症预测"
output_file <- file.path(output_dir, "high_quality_PBMC.tsv")

# 读取数据
df <- read.delim(input_file, sep = '\t')

# 定义筛选阈值
score_threshold <- 0.7
mapping_threshold <- 0.5

# 筛选高质量细胞 - 三个分类级别的置信指数都在0.7以上
high_quality_cells <- df[
  df$predicted.celltype.l1.score > score_threshold & 
    df$predicted.celltype.l2.score > score_threshold & 
    df$predicted.celltype.l3.score > score_threshold &
    df$mapping.score > mapping_threshold,
]

# 保存结果
write.table(high_quality_cells, output_file, sep = '\t', row.names = FALSE)

# 打印详细信息
cat("筛选完成！\n")
cat("原始细胞数量:", nrow(df), "\n")
cat("高质量细胞数量:", nrow(high_quality_cells), "\n")
cat("筛选比例:", round(nrow(high_quality_cells)/nrow(df)*100, 2), "%\n")
cat("结果已保存到:", output_file, "\n")