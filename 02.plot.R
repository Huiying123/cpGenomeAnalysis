#!/usr/bin/env Rscript

# =============================================================================
#         R SCRIPT FOR PLOTTING SLIDING WINDOW PI (WITH ANNOTATIONS)
# =============================================================================

# --- 1. 加载库并设置参数 ---
# ... (与之前脚本相同) ...
library(ggplot2)
library(scales)

# --- 用户配置区 ---
input_file <- "sliding_window_pi2.tsv"
# 新增：注释文件名
annotation_file <- "annotation.txt"
output_file <- "Sliding_Window_Pi_Plot_Annotated.pdf"
# ... (其他配置与之前相同) ...


# --- 2. 读取和处理数据 ---
# ... (读取滑动窗口数据，计算平均值，与之前脚本相同) ...
pi_data <- read.table(input_file, header = TRUE, sep = "\t", comment.char = "")
average_pi <- mean(pi_data$Pi, na.rm = TRUE)

# 新增：读取注释文件
anno_data <- read.table(annotation_file, header = TRUE, sep = "\t")
# 计算每个注释区域的中点，用于放置文本标签
anno_data$Midpoint <- (anno_data$Start_Pos + anno_data$End_Pos) / 2


# --- 3. 使用 ggplot2 绘制图表 (核心修改部分) ---
cat("--> Generating annotated plot...\n")

pi_plot <- ggplot(pi_data, aes(x = Window_Midpoint, y = Pi)) +
  
  # --- 新增：绘制基因/区域的色块 ---
  # 使用 geom_rect 绘制矩形来标记区域位置
  geom_rect(
    data = anno_data,
    aes(xmin = Start_Pos, xmax = End_Pos, ymin = -Inf, ymax = Inf, fill = Type),
    inherit.aes = FALSE,
    alpha = 0.2  # 设置透明度
  ) +
  
  # 自定义不同类型区域的填充颜色
  scale_fill_manual(values = c("Gene" = "lightgrey", "IGS" = "palegreen", "IR" = "lightblue", "Region" = "khaki")) +
  
  # --- 原有的滑动窗口曲线 ---
  geom_line(color = "dodgerblue", linewidth = 0.8) +
  geom_hline(yintercept = average_pi, linetype = "dashed", color = "red", linewidth = 0.6) +
  
  # --- 新增：添加基因/区域的文本标签 ---
  # 使用 geom_text 在每个区域下方添加名字
  geom_text(
    data = anno_data,
    aes(x = Midpoint, y = Label_Position_Y, label = Feature_Name),
    inherit.aes = FALSE,
    size = 3, # 调整字体大小
    angle = 90, # 将文字旋转90度，垂直显示
    hjust = 0   # 调整文字对齐方式
  ) +
  
  # --- 其他图表元素 (与之前脚本相同) ---
  labs(
    title = "Sliding Window Analysis of Nucleotide Diversity (π)",
    x = "Position along Genome (bp)",
    y = expression(paste("Nucleotide Diversity (", pi, ")")),
    fill = "Feature Type" # 为图例命名
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(min(anno_data$Label_Position_Y) * 1.1, NA)) +
  scale_x_continuous(labels = comma) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom" # 将图例放在底部
  )

# --- 4. 保存图表 ---
ggsave(output_file, plot = pi_plot, width = 15, height = 7, device = cairo_pdf)

cat("--> Annotated plot successfully saved to:", output_file, "\n")