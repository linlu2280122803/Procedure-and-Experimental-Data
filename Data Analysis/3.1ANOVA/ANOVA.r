# ==================== 三因素重复测量ANOVA完整分析脚本（修正版） ====================
# 实验设计：2 (警报类型：earcon vs auditory_icon) × 2 (噪音：低 vs 高) × 2 (任务：视觉 vs 听觉)
# 其中警报类型为被试间因素，噪音和任务为被试内因素。
# 数据文件：glmNDRT分析用表.csv
# =======================================================================

# 1. 加载必要的包 ----------------------------------------------------------
# 如果尚未安装，请先运行 install.packages(c("tidyverse", "afex", "emmeans"))
library(tidyverse)   # 数据操作和可视化
library(afex)        # 简化重复测量ANOVA
library(emmeans)     # 边际均值和事后比较
library(ggplot2)     # 绘图（tidyverse已包含）

# 2. 读取数据 --------------------------------------------------------------
# 请将文件路径修改为实际路径，或者将文件放在R工作目录下
df <- read.csv("glmNDRT分析用表.csv", stringsAsFactors = FALSE, fileEncoding = "GBK")

# 查看原始数据前几行
head(df)

# 3. 创建唯一被试ID --------------------------------------------------------
df$subject <- paste(df$警报类型, df$序号, sep = "_")

# 4. 转换为长格式并提取条件因子 ---------------------------------------------
df_long <- df %>%
  pivot_longer(
    cols = starts_with("任务绩效"),
    names_to = "condition",
    values_to = "performance"
  ) %>%
  mutate(
    noise = case_when(
      str_detect(condition, "低噪音") ~ "low",
      str_detect(condition, "高噪音") ~ "high"
    ),
    task = case_when(
      str_detect(condition, "视觉") ~ "visual",
      str_detect(condition, "听觉") ~ "auditory"
    ),
    alert = factor(警报类型, levels = c(1, 2), labels = c("earcon", "auditory_icon")),
    noise = factor(noise, levels = c("low", "high")),
    task  = factor(task,  levels = c("visual", "auditory"))
  ) %>%
  select(subject, alert, noise, task, performance) %>%
  filter(!is.na(performance))

# 检查转换后的数据
head(df_long)
summary(df_long)

# ==================== 5. 三因素重复测量ANOVA ==============================
model <- aov_ez(
  id       = "subject",
  dv       = "performance",
  data     = df_long,
  between  = "alert",
  within   = c("noise", "task")
)

# 打印ANOVA结果
print(model)

# ==================== 6. 提取边际均值（描述性统计表） =======================
# 获取所有8个条件的估计边际均值
emm_all <- emmeans(model, ~ alert * noise * task)

# 转换为数据框并重命名列
emm_desc <- as.data.frame(emm_all) %>%
  rename(
    Mean = emmean,
    SE = SE,
    df = df,
    CI_Lower = lower.CL,
    CI_Upper = upper.CL
  ) %>%
  arrange(alert, noise, task)

# 导出描述性表格（不含伪SD）
write.csv(emm_desc, "Estimated_Marginal_Means.csv", row.names = FALSE)

# 可选：在控制台打印
print(emm_desc)

# ==================== 7. 事后比较 ========================================

### 7.1 主效应的成对比较（Bonferroni校正，共3对）---------------------------
# 分别获取 task, noise, alert 的主效应边际均值
emm_task <- emmeans(model, ~ task)
emm_noise <- emmeans(model, ~ noise)
emm_alert <- emmeans(model, ~ alert)

# 成对比较（每个因素只有一对），并应用 Bonferroni 校正
pairs_task <- pairs(emm_task, adjust = "bonferroni")
pairs_noise <- pairs(emm_noise, adjust = "bonferroni")
pairs_alert <- pairs(emm_alert, adjust = "bonferroni")

# 查看结果
print("=== Task 主效应比较（视觉 vs 听觉）===")
print(pairs_task)
print("=== Noise 主效应比较（低 vs 高）===")
print(pairs_noise)
print("=== Alert 主效应比较（earcon vs auditory_icon）===")
print(pairs_alert)

# 合并三个主效应的比较结果为一个表格
main_effects <- rbind(
  as.data.frame(pairs_task),
  as.data.frame(pairs_noise),
  as.data.frame(pairs_alert)
)
write.csv(main_effects, "main_effects_pairwise_bonferroni.csv", row.names = FALSE)

### 7.2 特定对比（听觉图标+视觉 vs 耳标+听觉）在不同噪声水平下 --------------
# 注意：这个对比不是标准的主效应或简单效应，而是自定义对比
# 定义两个对比的系数（基于 emm_all 的顺序：earcon low visual, auditory_icon low visual, earcon high visual, auditory_icon high visual, earcon low auditory, auditory_icon low auditory, earcon high auditory, auditory_icon high auditory）
# 低噪声下：auditory_icon low visual (索引2) 减去 earcon low auditory (索引5)
# 高噪声下：auditory_icon high visual (索引4) 减去 earcon high auditory (索引7)
custom_contrasts <- list(
  "LowNoise: AudIcon_Visual minus Earcon_Auditory" = c(0, 1, 0, 0, -1, 0, 0, 0),
  "HighNoise: AudIcon_Visual minus Earcon_Auditory" = c(0, 0, 0, 1, 0, 0, -1, 0)
)

# 计算对比（仅这两个对比，使用Bonferroni校正，因为共2对）
custom_comp <- contrast(emm_all, custom_contrasts, adjust = "bonferroni")
custom_summary <- summary(custom_comp, infer = c(TRUE, TRUE))
print("=== 特定对比结果（听觉图标+视觉 vs 耳标+听觉）===")
print(custom_summary)

# 导出特定对比结果
write.csv(as.data.frame(custom_summary), "custom_contrasts_bonferroni.csv", row.names = FALSE)

# ==================== 8. 可视化 ===========================================
p <- ggplot(df_long, aes(x = noise, y = performance, color = alert, group = alert)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.8,
               position = position_dodge(0.1)) +
  facet_wrap(~ task, labeller = labeller(task = c(visual = "视觉任务", auditory = "听觉任务"))) +
  labs(x = "噪音水平", y = "任务绩效", color = "警报类型") +
  scale_color_manual(values = c("earcon" = "#7BB6E3",
                                "auditory_icon" = "#8EB052"),
                     labels = c("earcon" = "耳标", "auditory_icon" = "听觉图标")) +
  theme_minimal() +
  theme(legend.position = "top",
        strip.text = element_text(size = 12))

print(p)

# 保存图形（可选）
# ggsave("ANOVA_plot.pdf", plot = p, width = 8, height = 6)
# ggsave("ANOVA_plot.png", plot = p, width = 8, height = 6, dpi = 300)

# ==================== 9. 残差诊断 =========================================
residuals <- resid(model$lm)
qqnorm(residuals)
qqline(residuals, col = "red")

fitted <- fitted(model$lm)
plot(fitted, residuals)
abline(h = 0, lty = 2)

# ==================== 10. 结束 ===========================================
cat("\n分析完成！所有结果已导出。\n")