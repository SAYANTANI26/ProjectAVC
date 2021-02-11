rm(list = ls())

# Scripted by @ArpanKumarBasak
# email: arpankbasak@gmail.com

require(tidyverse)
require(cowplot)

# Script for analysis of COVID data and widlife

# Make theme
theme_set <- ggplot2::theme_bw() +
  ggplot2::theme(text = element_text(size = 16, hjust = 0.5, vjust = 0.5),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        legend.key=element_blank(),
        legend.text.align=0,
        legend.text=element_text(size = 10, hjust = 0.5, vjust = 0.5),
        legend.position="top",
        strip.text=element_text(face="bold", size = 16),
        axis.line=element_line(),
        axis.line.x=element_line(size = 1),
        axis.line.y=element_line(size = 1),
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        strip.background=element_rect(fill="transparent", colour=NA),
        strip.placement="outside")

# Read data
# df <- read.table("/klaster/work/abasak/project_SB/data/dataset.txt", 
#                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# # Table curation done after gathering the data
# df$species <- str_replace(df$species, pattern = "\312", replacement = "")
# df <- df %>% select(-10,-11)
# colnames(df)[6] <- "reported_time"
# df$species <- str_replace(df$species, "bird$", "birds")
# write.table(df, "./data/dataset_carpented.txt", sep = "\t", quote = FALSE)

df <- read.table("/klaster/work/abasak/project_SB/data/dataset_carpented.txt", 
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
filter(year %in% c(2019,2020), month != 12)

# df$year <- str_replace(as.character(df$year), "202$|220$", "2020")
# df$Gmina <- str_replace(as.character(df$Gmina), "202$|220$", "2020")
# df$Gmina <- str_replace(as.character(df$Gmina), "\312", "")
# df$species <- str_replace(as.character(df$species), "wild boar", "")

df <- df %>% mutate(year = as.numeric(year), 
  species = as.factor(species),
  month = as.numeric(month), 
  day = as.numeric(day),
  ID = paste(ID, Gmina, day, month, year, sep = "_"),
  group = as.factor(paste(Gmina, year, month, sep = "_")),
  lockdown = month %in% c(3, 4, 5)
)

# Factor levels and colours
time_stamps <- c("0000_0400", "0400_0800", "0800_1200", "1200_1600", "1600_2000", "2000_0000")
col.pal <- c("darkblue", "skyblue", "lightyellow", "lightcoral", "lavender", "darkviolet")
levs <- sapply(unique(df$Gmina), function(x) paste(x, time_stamps, sep = "_"))

df_wide <- df %>% 
mutate(reported_time = as.numeric(str_replace(as.character(reported_time), ":", ""))) %>%
filter(!is.na(reported_time)) %>%
mutate(time_range = factor(cut(reported_time, breaks = c(0, 400, 800, 1200, 1600, 2000, 2400),
  include.lowest = TRUE, labels = time_stamps), levels = time_stamps)) %>%
group_by(Gmina, year, month, time_range, species) %>% 
summarise(vals = n()) %>%
spread(key = species, value = "vals", fill = 0, convert = FALSE) %>%
mutate(group = as.factor(paste(Gmina, year, month, time_range, sep = "_")),
  Gmina_time = factor(paste(Gmina, time_range, sep = "_"), levels = levs)) %>%
data.frame(., stringsAsFactors = FALSE)

# Data summary
df_wide %>%
group_by(Gmina, year, month) %>%
summarise(n = n()) %>%
data.frame(.) %>%
write.table(., "./output/data_temporal_summary.txt", 
  sep = "\t", 
quote = FALSE, 
  row.names = FALSE)

trend_df <- df_wide %>% group_by(Gmina, year, month, time_range) %>% 
gather(key = "species", value = "vals", convert = FALSE, -Gmina, -year, -month, -group, -Gmina_time, -time_range) %>%
summarise(tot_hwc = sum(vals)) %>%
spread(key = time_range, value = "tot_hwc", fill = 0, convert = FALSE) %>%
mutate(group = as.factor(paste(Gmina, year, month, sep = "_")),
  lockdown = as.factor(c(month %in% c(3, 4, 5) & year == 2020))) %>%
data.frame(., stringsAsFactors = FALSE)

mat_time <- apply((trend_df %>% select(-Gmina, -year, -month, -group, -lockdown)), 2, function(x) x/sum(x))
row.names(mat_time) <- as.character(trend_df$group)
row.names(trend_df) <- as.character(trend_df$group)

# d <- 1 - cor(t(mat_time))
d <- vegan::vegdist(mat_time, method = "jaccard")
hc <- hclust(as.dist(d), "ward.D2")
data_sorted_time <- row.names(mat_time)[hc$order]
mds_obj <- cmdscale(as.dist(d), eig = TRUE, k = 3)
variance <- round(100*mds_obj$eig/sum(mds_obj$eig), 2) 

d <- 1 - cor(mat_time)
hc <- hclust(as.dist(d), "ward.D2")
time_sorted <- row.names(t(mat_time))[hc$order]

# Time span
df_trendline <- df_wide %>% group_by(Gmina, year, month, time_range) %>% 
gather(key = "species", value = "vals", convert = FALSE, -Gmina, -year, -month, -group, -Gmina_time, -time_range) %>%
summarise(tot_hwc = sum(vals)) %>%
data.frame(.) %>%
separate(time_range, into = c("t1", "t2")) %>%
mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = as.factor(c(month %in% c(3, 4, 5) & year == 2020)),
        t1 = as.numeric(t1)
  ) %>% 
data.frame(.)

# Effect of lockdown on the trend
fit <- aov(glm(tot_hwc ~ 0 + Gmina:lockdown, data = df_trendline))
mod <- broom::tidy(TukeyHSD(fit)) %>%
mutate(FDR = p.adjust(adj.p.value, "fdr")) %>%
mutate(significance = FDR < 0.05) %>%
arrange(desc(significance)) %>%
data.frame(.) 

write.table(mod, "./statistics/model_comparison_trendline.txt", quote = FALSE, sep = "\t")

df_trendline %>%
ggplot(aes(x = t1, y = tot_hwc)) +
geom_point(aes(colour = lockdown, size = lockdown, alpha = lockdown), position = "jitter") +
geom_smooth(method = "loess", se = FALSE, aes(colour = lockdown, linetype = lockdown)) +
facet_grid(Gmina ~., scale = "fixed", space = "free") +
scale_x_continuous(breaks = c(0, 400, 800, 1200, 1600, 2000, 2400)) +
scale_colour_manual(values = c(`FALSE` = "lightgrey", `TRUE` = "black"), guide = FALSE) +
scale_size_manual(values = c(`FALSE` = 0.3, `TRUE` = 1), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.6), guide = FALSE) +
scale_linetype_manual(values = c(`FALSE` = "dashed", `TRUE` = "solid"), guide = FALSE) +
labs(x = "Time (hrs)", y = "Total AVC") +
theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14)) +
  ggsave("/klaster/work/abasak/project_SB/figures/trendline_temporal.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

# Plot NMDS -- temporal trend
(plot_mds12_time <- cbind.data.frame(trend_df[,c(1,2,3)], mat_time, 
                 mds1 = mds_obj$points[,1], 
                 mds2 = mds_obj$points[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = mds1, y = mds2, group = lockdown)) +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(colour = lockdown, alpha = lockdown, fill = Gmina, shape = Gmina), size = 3) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  # scale_colour_manual(values = col.idx) +
  scale_fill_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_colour_manual(values = c(`FALSE` = "lightgrey", `TRUE` = "black"), guide = FALSE) +
  scale_shape_manual(values = c(21,22,23)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  labs(x = paste0("MDS1 - ", variance[1], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS12_analysis_temporal_trend.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

(plot_mds23_time <- cbind.data.frame(trend_df[,c(1,2,3)], mat_time, 
                 mds3 = mds_obj$points[,3], 
                 mds2 = mds_obj$points[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = mds3, y = mds2, group = lockdown)) +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(colour = lockdown, alpha = lockdown, fill = Gmina, shape = Gmina), size = 3) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  # scale_colour_manual(values = col.idx) +
  scale_fill_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_colour_manual(values = c(`FALSE` = "lightgrey", `TRUE` = "black"), guide = FALSE) +
  scale_shape_manual(values = c(21,22,23)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  labs(x = paste0("MDS3 - ", variance[3], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS23_analysis_temporal_trend.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

cbind.data.frame(trend_df[,c(1,2,3)], mat_time, 
                 mds1 = mds_obj$points[,1], 
                 mds2 = mds_obj$points[,2]) %>%
mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>%
write.table(., file  = "./statistics/mds_trendine_dataframe.txt", sep = "\t", quote = FALSE)

# Temporal multivariate analysis
set.seed(1)
cca_obj <- vegan::capscale(formula = mat_time ~ lockdown:Gmina + Condition(month + year), 
  data = trend_df, perm = 1000, dist = "jaccard", sqrt.dist = TRUE)
anova_obj <- vegan::anova.cca(cca_obj)

# Fetch stats
variance <- round(100*cca_obj$CCA$eig/sum(cca_obj$CCA$eig), 2)
pval <- anova_obj$`Pr(>F)`[1]
con_var <- round(100*cca_obj$CCA$tot.chi/cca_obj$tot.chi, 2)
chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                       row.names = c("total", "constrianed", "unconstrained"))

ti <- paste0("(Time_range_HWC) ~ Lockdown:Region + Condition(month + year), Variance explained = ", con_var, " % ; p.value = ", round(pval, 3), "*")

# CCA analysis
(cca_plot_trend <- cbind.data.frame(trend_df[,c(1,2,3)], mat_time, 
                 cca1 = cca_obj$CCA$wa[,1], 
                 cca2 = cca_obj$CCA$wa[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
         lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = cca1, y = cca2)) +
  ggtitle(ti) +
  # facet_wrap(year ~ season, scale = "free", switch = "x") +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(colour = lockdown, fill = Gmina, shape = Gmina, alpha = lockdown, size = lockdown)) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black")) +
  scale_shape_manual(values = c(21,22,23)) +
  scale_fill_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 3), guide = FALSE) +
  labs(x = paste0("CCA1 - ", variance[1], " %"), 
       y = paste0("CCA2 - ", variance[2], " %"), 
       colour = "",
       shape = "") +
  theme_set + 
  theme(title = element_text(size = 4),
        legend.text = element_text(size = 4),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/cca_analysis_time_trend.png", 
         dpi = 600, units = "in", 
         width = 5, height = 5, limitsize = FALSE, bg = "transparent",
         device = "png")

# Save stats
write.table(variable, file = "/klaster/work/abasak/project_SB/statistics/temporal_trendline_cca_summary.txt",
           quote = F, sep = "\t")

# Spatio temporal analysis based on time range
mat <- apply((df_wide %>% select(-Gmina, -year, -month, -group, -Gmina_time, -time_range)), 2, function(x) x/sum(x))
row.names(mat) <- as.character(df_wide$group)
row.names(df_wide) <- as.character(df_wide$group)

# d <- 1 - cor(t(mat))
d <- vegan::vegdist(mat, method = "jaccard")
hc <- hclust(as.dist(d), "ward.D2")
data_sorted <- row.names(mat)[hc$order]
mds_obj <- cmdscale(as.dist(d), eig = TRUE, k = 3)
variance <- round(100*mds_obj$eig/sum(mds_obj$eig), 2) 

d <- 1 - cor(mat)
hc <- hclust(as.dist(d), "ward.D2")
animals_sorted <- row.names(t(mat))[hc$order]

# Plot NMDS -- lockdown
(plot_mds12 <- cbind.data.frame(df_wide[,c(1,2,3,4)], mat, 
                 mds1 = mds_obj$points[,1], 
                 mds2 = mds_obj$points[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = mds1, y = mds2, group = time_range)) +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(colour = lockdown, alpha = lockdown, fill = time_range, shape = Gmina), size = 3) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  # scale_colour_manual(values = col.idx) +
  scale_fill_manual(values = col.pal) +
  scale_colour_manual(values = c(`FALSE` = "lightgrey", `TRUE` = "black"), guide = FALSE) +
  scale_shape_manual(values = c(21,22,23)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  labs(x = paste0("MDS1 - ", variance[1], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS12_analysis_temporal.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

(plot_mds23 <- cbind.data.frame(df_wide[,c(1,2,3,4)], mat, 
                 mds3 = mds_obj$points[,3], 
                 mds2 = mds_obj$points[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = mds3, y = mds2, group = lockdown)) +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(colour = lockdown, alpha = lockdown, fill = time_range, shape = Gmina), size = 3) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  # scale_colour_manual(values = col.idx) +
  scale_fill_manual(values = col.pal) +
  scale_colour_manual(values = c(`FALSE` = "lightgrey", `TRUE` = "black"), guide = FALSE) +
  scale_shape_manual(values = c(21,22,23)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  labs(x = paste0("MDS3 - ", variance[3], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS23_analysis_temporal.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

# Plot NMDS -- location
(plot_mds12_l <- cbind.data.frame(df_wide[,c(1,2,3,4)], mat, 
                 mds1 = mds_obj$points[,1], 
                 mds2 = mds_obj$points[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = mds1, y = mds2, group = Gmina)) +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(alpha = lockdown, size = lockdown, colour = Gmina, shape = time_range)) +
  stat_ellipse(aes(group = Gmina, colour = Gmina), 
               lty = "dashed", size = 1, 
               type = "norm") +
  scale_colour_manual(values = c("tomato", "olivedrab", "slateblue")) +
  # scale_fill_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_shape_manual(values = c(0,1,2,5,6,19)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 3), guide = FALSE) +
  labs(x = paste0("MDS1 - ", variance[1], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS12_loc_analysis_temporal.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

(plot_mds23_l <- cbind.data.frame(df_wide[,c(1,2,3,4)], mat, 
                 mds3 = mds_obj$points[,3], 
                 mds2 = mds_obj$points[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = mds3, y = mds2, group = Gmina)) +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(alpha = lockdown, size = lockdown, colour = Gmina, shape = time_range)) +
  stat_ellipse(aes(group = Gmina, colour = Gmina), 
               lty = "dashed", size = 1, 
               type = "norm") +
  scale_colour_manual(values = c("tomato", "olivedrab", "slateblue")) +
  # scale_fill_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_shape_manual(values = c(0,1,2,5,6,19)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 3), guide = FALSE) +
  labs(x = paste0("MDS3 - ", variance[3], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS23_loc_analysis_temporal.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")


# Write the MDS output
cbind.data.frame(df_wide[,c(1,2,3,4)], mat, 
                 mds1 = mds_obj$points[,1], 
                 mds2 = mds_obj$points[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% write.table(., "./statistics/mds_analysis_temporal.txt", sep = "\t", quote = FALSE)

# Cannonical analysis
# cca_obj <- vegan::capscale(formula = mat ~ month + Condition(Gmina + year), 
#   data = df_wide)
# aov_cca <- vegan::anova.cca(cca_obj)

# cca_obj <- vegan::capscale(formula = mat ~ year + Condition(Gmina + month), 
#   data = df_wide)
# aov_cca <- vegan::anova.cca(cca_obj)
set.seed(1)
df_wide$lockdown <- as.factor(c(df_wide$month %in% c(3, 4, 5) & df_wide$year == 2020))
cca_obj <- vegan::capscale(formula = mat ~ time_range:Gmina:lockdown + Condition(month + year), 
  data = df_wide, perm = 1000, dist = "jaccard", sqrt.dist = TRUE)
anova_obj <- vegan::anova.cca(cca_obj)

# d <- as.dist(1-cor(t(mat)))
# maov <- vegan::adonis(formula =  d ~ time_range + Gmina + month + year, data = df_wide,
#   perm = 1000)

# Fetch stats
variance <- round(100*cca_obj$CCA$eig/sum(cca_obj$CCA$eig), 2)
pval <- anova_obj$`Pr(>F)`[1]
con_var <- round(100*cca_obj$CCA$tot.chi/cca_obj$tot.chi, 2)
chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                       row.names = c("total", "constrianed", "unconstrained"))

ti <- paste0("(AVC) ~ Time:Gmina:Lockdown + Condition(month + year), Variance explained = ", con_var, " % ; p.value = ", round(pval, 3), "*")

# Pairwise PERMANOVA for the factors
df_wide$group <- as.factor(paste(df_wide$time_range, df_wide$Gmina, df_wide$lockdown, sep = "_"))
pwm_mod <- RVAideMemoire::pairwise.perm.manova(vegan::vegdist(mat, "jaccard"), df_wide$group, 
  test = "Wilks", 
  nperm = 1000, 
  p.method = "fdr")

# Obtain the pairwise statistics for the comparison
(pw_stats <- as.data.frame(pwm_mod$p.value) %>%
add_column(A = row.names(.), .before = 1) %>%
gather(key = "B", value = "FDR", convert = FALSE, -A) %>%
na.omit(.) %>%
separate(A, into = c("A_time_1", "A_time_2", "A_location", "A_lockdown"), sep = "_", convert = FALSE) %>%
separate(B, into = c("B_time_1", "B_time_2", "B_location", "B_lockdown"), sep = "_", convert = FALSE) %>%
mutate(A_time = paste(A_time_1, A_time_2, sep = "_"),
  B_time = paste(B_time_1, B_time_2, sep = "_"),
  significance = FDR < 0.05) %>%
arrange(desc(significance)) %>%
select(-A_time_1, -A_time_2, -B_time_1, -B_time_2)) %>%
write.table(., file = "./statistics/pw_perm_temporal.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot airwise permanova
pw_stats %>%
mutate(A_location = as.factor(A_location),
A_lockdown = as.factor(A_lockdown),
B_location = as.factor(B_location),
B_lockdown = as.factor(B_lockdown), 
A_time = as.factor(A_time),
B_time = as.factor(B_time), 
significance = as.factor(significance),
FDR = -log10(FDR)) %>%
ggplot(aes(x = A_time, y = B_time)) +
geom_raster(aes(fill = FDR)) +
geom_tile(aes(colour = significance), 
  fill = "transparent", 
  size = 0.9,
  width = 0.95,
  height = 0.95) +
facet_grid(A_lockdown + A_location ~ B_lockdown + B_location, scale = "free", switch = "both", space = "free") +
scale_fill_gradient(low = "white", 
                      high = "darkred", 
                      na.value = "darkgrey") +
scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "darkgrey"), guide = FALSE) +
theme_set +
theme(panel.spacing = unit(0.1, "lines"),
      axis.text.x = element_text(angle = 90, size = 8, hjust = 0.8, vjust = 0.5),
      axis.text.y = element_text(size = 8), 
      strip.text.y = element_text(size = 8), 
      strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(hjust = 0.5, vjust = 0.8, size = 5)
      ) +
labs(x = "", fill = "", y = "") +
ggsave("/klaster/work/abasak/project_SB/figures/pw_pnova_stats_temporal.png",
       dpi = 600, units = "in",
       width = 8, height = 9, limitsize = FALSE, 
       bg = "transparent",
       device = "png")


# CCA analysis
(cca_plot <- cbind.data.frame(df_wide[,c(1,2,3,4)], mat, 
                 cca1 = cca_obj$CCA$wa[,1], 
                 cca2 = cca_obj$CCA$wa[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = cca1, y = cca2)) +
  ggtitle(ti) +
  # facet_wrap(year ~ season, scale = "free", switch = "x") +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(colour = lockdown, shape = Gmina, 
    fill = time_range,
    alpha = lockdown, 
    size = lockdown)) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"), guide = FALSE) +
  scale_shape_manual(values = c(21,22,23)) +
  scale_fill_manual(values = col.pal) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 3), guide = FALSE) +
  labs(x = paste0("CCA1 - ", variance[1], " %"), 
       y = paste0("CCA2 - ", variance[2], " %"), 
       colour = "",
       shape = "") +
  theme_set + 
  theme(title = element_text(size = 4),
        legend.text = element_text(size = 4),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/cca_analysis_time.png", 
         dpi = 600, units = "in", 
         width = 5, height = 5, limitsize = FALSE, bg = "transparent",
         device = "png")

# CCA highlighting the location density
(cca_plot_l <- cbind.data.frame(df_wide[,c(1,2,3,4)], mat, 
                 cca1 = cca_obj$CCA$wa[,1], 
                 cca2 = cca_obj$CCA$wa[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = cca1, y = cca2)) +
  ggtitle(ti) +
  # facet_wrap(year ~ season, scale = "free", switch = "x") +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(shape = time_range, 
    alpha = lockdown,
    size = lockdown, colour = Gmina)) +
  stat_ellipse(aes(group = Gmina, colour = Gmina), 
               lty = "dashed", size = 1, 
               type = "norm") +
  scale_colour_manual(values = c("tomato", "olivedrab", "slateblue")) +
  # scale_fill_manual(values = c("tomato", "olivedrab", "slateblue")) +
  # scale_shape_manual(values = c(1,19)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 3), guide = FALSE) +
  labs(x = paste0("CCA1 - ", variance[1], " %"), 
       y = paste0("CCA2 - ", variance[2], " %"), 
       colour = "",
       shape = "") +
  theme_set + 
  theme(title = element_text(size = 4),
        legend.text = element_text(size = 4),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/cca_analysis_loc_time.png", 
         dpi = 600, units = "in", 
         width = 5, height = 5, limitsize = FALSE, bg = "transparent",
         device = "png")

# Relative abundance of the animals observed
mat_df <- cbind.data.frame(df_wide[,c(1,2,3,4)], mat)

df_mat_long <- mat_df %>% 
gather(key = "species", value = "count", 
  convert = FALSE, -Gmina, -year, -month, -time_range) %>%
mutate(group = as.factor(paste(Gmina, year, month, time_range, sep = "_")),
  lockdown = as.factor(c(month %in% c(3, 4, 5) & year == 2020))) %>%
data.frame(., stringsAsFactors = FALSE)

animals <- colnames(mat)
res <- lapply(animals, function(x){

    # x = "badger"
    temp <- df_mat_long %>% filter(species == x) %>%
    mutate(year = as.factor(year), month = as.factor(month))

    fit <- aov(lm(count ~ 0 + lockdown:Gmina, data = temp))
    mod <- broom::tidy(TukeyHSD(fit)) %>% 
    mutate(animal = x) %>%
    data.frame(., stringsAsFactors = FALSE)

  })

# Satistics data carpent
stat_df <- do.call(rbind.data.frame, res) %>%
mutate(FDR = p.adjust(adj.p.value, "fdr"))%>%
mutate(significance = FDR <= 0.05) %>%
separate(contrast, into = c("A", "B"), convert = FALSE, sep = "-") %>%
separate(A, into = c("A_lockdown", "A_region"), convert = FALSE, sep = ":") %>%
separate(B, into = c("B_lockdown", "B_region"), convert = FALSE, sep = ":") %>%
arrange(desc(significance))

# stat_df %>% filter(A_lockdown == TRUE)

# Plot of relative abundance
(hmap1a <- df_mat_long %>% 
  mutate(species = factor(species, levels = animals_sorted, labels = str_replace_all(animals_sorted, "\\.", " ")),
        id = factor(group, levels = data_sorted),
        Gmina = as.factor(Gmina),
        month = as.factor(month),
        year = as.factor(year),
        lockdown = as.factor(ifelse(month %in% c(3, 4, 5) & year == 2020, "Locked month", ""))) %>%
  group_by(species, Gmina, lockdown, time_range) %>%
  summarise(count = mean(count)) %>%
  mutate(logra = ifelse(!is.finite(log10(count)), NA, log10(count))) %>% 
  ggplot(aes(x = time_range, y = species)) +
  geom_raster(alpha = 1, aes(fill = logra), show.legend = T) +
  facet_grid(.~ Gmina + lockdown, space = "fixed", scale = "free") +
  scale_fill_gradient(low = "black", high = "yellow", na.value = "black",
                      breaks = c(-2, -1.5, -1, -0.5, 0), 
                      labels = c(paste(c(0.01, 0.1, 1, 10, 100), sep = ""))
                      ) +
  theme_set +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 0.8, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        # strip.text.y = element_text(angle = 180), 
        strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.8)
        ) +
  labs(x = "", fill = "% meanHWC", y = "")) +
  ggsave("/klaster/work/abasak/project_SB/figures/abundance_HWC_temporal.png",
         dpi = 600, units = "in",
         width = 7, height = 10, limitsize = FALSE, 
         device = "png")


# Plot the summary statistics
(hmap1b <- stat_df %>% 
  filter(A_lockdown == TRUE, B_region == "KRAKOW") %>%
  mutate(animal = factor(animal, levels = animals_sorted, labels = str_replace_all(animals_sorted, "\\.", " ")), 
    # comparison = as.factor(comparison), 
    significance = as.factor(significance)
    # B_time = factor(B_time, levels = time_stamps),
    # A_time = factor(A_time, levels = time_stamps)) %>%
    # term = as.factor(term)
    ) %>% 
  na.omit(.) %>%  
  ggplot(aes(x = A_region, y = animal)) +
  geom_raster(alpha = 1, aes(fill = -estimate), show.legend = T) +
  geom_tile(alpha = 1, aes(colour = significance), fill = NA, 
    width = 0.9, height = 0.9, size = 1) +
  facet_grid(.~ B_lockdown, space = "free", scale = "free") +
  scale_fill_gradient2(low = "darkblue", 
                      high = "darkred", mid = "white",
                      na.value = "darkgrey"
                      # breaks = c(-4, -2, 0, 2, 4), 
                      # labels = c(paste(c(-4, -2, 0, 2, 4), sep = ""))
                      ) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA), guide = FALSE) +
  theme_set +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 0.8, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        # strip.text.y = element_text(angle = 180), 
        strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.8, size = 5)
        ) +
  labs(x = "", fill = "Mean estimate", y = "")) +
  ggsave("/klaster/work/abasak/project_SB/figures/stats_HWC_gmina_temporal.png",
         dpi = 600, units = "in",
         width = 10, height = 9, limitsize = FALSE, 
         device = "png")

# Make boxplots for which are consistently signinficnant on absolute counts
sig_animals <- na.omit(unique(stat_df$animal[stat_df$significance == TRUE])) %>% as.character()
parallel::mclapply(sig_animals, function(x){

  # Plots for bxp
  # x = "wild.boar"
  temp <- df_wide
  colnames(temp)[colnames(temp) == x] <- "id"
  
  temp %>% select(Gmina, year, month, time_range, id) %>%
  separate(time_range, into = c("t1", "t2"), sep = "_") %>%
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020),
        id = 100*id/sum(id), 
        t1 = as.numeric(t1)
  ) %>%
  ggplot(aes(x = t1, y = id)) +
  geom_hline(yintercept = mean(100*temp$id/sum(temp$id)), colour = "darkred", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = c(600,1800), colour = "lightgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(stat = "identity", 
              position = "jitter", 
              alpha = 0.6, 
              aes(shape = year, alpha = lockdown, size = lockdown)) +
  # geom_boxplot(outlier.alpha = 0, aes(colour = time_range), 
  #              size = 2, 
  #              fill = NA, 
  #              alpha = 0.6, 
  #              position = position_dodge2()) +
  geom_smooth(method = "loess", se = FALSE, aes(colour = lockdown, linetype = lockdown)) +
  facet_grid(Gmina ~., scale = "fixed", space = "free") +
  scale_shape_manual(values = c(1,19)) +
  scale_x_continuous(breaks = c(0, 400, 800, 1200, 1600, 2000, 2400)) +
  scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"), guide = FALSE) +
  scale_size_manual(values = c(`FALSE` = 0.3, `TRUE` = 1), guide = FALSE) +
  scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.6), guide = FALSE) +
  scale_linetype_manual(values = c(`FALSE` = "dashed", `TRUE` = "solid"), guide = FALSE) +
  labs(x = "Time (hrs)", y = paste0(x, " AVC")) +
  theme_set +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 12), 
        # strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.5, size = 5)) +
  ggsave(paste0("/klaster/work/abasak/project_SB/figures/boxplots_temporal/trendline_", x, ".png"),
         dpi = 600, units = "in",
         width = 4, height = 4.5, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

}, mc.cores = 2)

# ======================

# Filter the main dataset by city == Krakow
df_krk <- read.table("/klaster/work/abasak/project_SB/data/dataset_carpented.txt", 
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
filter(year %in% c(2019,2020), month != 12, Gmina == "KRAKOW")

df_krk_w <- df_krk %>% 
mutate(reported_time = as.numeric(str_replace(as.character(reported_time), ":", ""))) %>%
filter(!is.na(reported_time)) %>%
mutate(time_range = factor(cut(reported_time, breaks = c(0, 400, 800, 1200, 1600, 2000, 2400),
  include.lowest = TRUE, labels = time_stamps), levels = time_stamps)) %>%
group_by(Gmina, year, month, day, time_range, species) %>% 
summarise(vals = n()) %>%
spread(key = species, value = "vals", fill = 0, convert = FALSE) %>%
mutate(group = as.factor(paste(Gmina, year, month, day, time_range, sep = "_"))) %>%
data.frame(., stringsAsFactors = FALSE)

# Compute GLM for the filtered data
mat <- apply((df_krk_w %>% select(-Gmina, -year, -month, -time_range, -group, -day)), 2, function(x) x/sum(x))
row.names(mat) <- df_krk_w$group
row.names(df_krk_w) <- df_krk_w$group

# Compute statistics G-linear model
mat_df <- cbind.data.frame(df_krk_w[,c(1,2,3,4,5,23)], mat)

df_mat_long <- mat_df %>% 
gather(key = "species", value = "count", 
  convert = FALSE, -Gmina, -year, -month, -time_range, -day, -group) %>%
# mutate(group = as.factor(paste(Gmina, year, month, sep = "_"))) %>%
data.frame(., stringsAsFactors = FALSE)

animals <- colnames(mat)
res <- parallel::mclapply(animals, function(x){

    # x = "badger"
    temp <- df_mat_long %>% 
    filter(species == x) %>%
    mutate(year = as.factor(year), month = as.factor(month), 
      time_range = as.factor(time_range),
      count = as.numeric(count)) %>%
    mutate(group = as.factor(paste(year, month, sep = "_")))

    fit <- aov(lm(count ~ 0 + year:month:time_range, data = temp))
    mod <- broom::tidy(TukeyHSD(fit)) %>% 
    mutate(animal = x) %>%
    data.frame(., stringsAsFactors = FALSE)

}, mc.cores = 8)

stat_df <- do.call(rbind.data.frame, res) %>%
# mutate(FDR = p.adjust(adj.p.value, "fdr"))%>%
separate(contrast, into = c("A", "B"), sep = "-", remove = FALSE) %>%
separate(A, into = c("A_year", "A_month", "A_time_range"), sep = ":") %>%
separate(B, into = c("B_year", "B_month", "B_time_range"), sep = ":") %>%
filter(A_year == "2020", A_month == B_month, A_time_range == B_time_range) %>%
mutate(significance = adj.p.value <= 0.05) %>%
arrange(desc(significance))

# Make individual plot of within Krakow AVC

animal_mat <- stat_df %>% select(contrast, estimate, animal) %>%
spread(key = animal, value = estimate, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = F)

# Compute correlation between the animal variations within comparison
d <- 1 - cor(animal_mat[,-1])
hc <- hclust(as.dist(d), "ward.D2")
animal_sorted <- colnames(animal_mat[,-1])[hc$order]

# Plot on the basis of mean estimate
(hmap_2a <- stat_df %>% 
  # filter(A_lockdown == TRUE, B_region == "KRAKOW") %>%
  mutate(animal = factor(animal, levels = animal_sorted, 
    labels = str_replace_all(animal_sorted, "\\.", " ")), 
    # comparison = as.factor(comparison), 
    significance = as.factor(significance)
    # B_time = factor(B_time, levels = time_stamps),
    # A_time = factor(A_time, levels = time_stamps)) %>%
    # term = as.factor(term)
    ) %>% 
  na.omit(.) %>%  
  ggplot(aes(y = A_time_range, x = animal)) +
  geom_raster(alpha = 1, aes(fill = estimate), show.legend = T) +
  geom_tile(alpha = 1, aes(colour = significance), fill = NA, 
    width = 0.9, height = 0.9, size = 1) +
  facet_grid(A_month ~., switch = "y", space = "free", 
    scale = "free") +
  scale_fill_gradient2(low = "darkblue", 
                      high = "darkred", mid = "white",
                      na.value = "darkgrey"
                      # breaks = c(-4, -2, 0, 2, 4), 
                      # labels = c(paste(c(-4, -2, 0, 2, 4), sep = ""))
                      ) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA), guide = FALSE) +
  theme_set +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 0.8, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        # strip.text.y = element_text(angle = 180), 
        strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.8, size = 5)
        ) +
  labs(x = "", fill = "Mean difference", y = "")) +
  ggsave("/klaster/work/abasak/project_SB/figures/stats_krakow_avc.png",
         dpi = 600, units = "in",
         width = 10, height = 9, limitsize = FALSE, 
         device = "png")


# Compute total AVC
df_tot <- df_mat_long %>%
group_by(year, month, day, time_range) %>%
summarise(tot_avc = sum(count)) %>%
ungroup(.) %>%
mutate(year = as.factor(year),
  month = as.factor(month),
  time_range = as.factor(time_range)) %>%
data.frame(.)

# Make a matrix for integrative calculation
mat_tot <- df_tot %>%
spread(key = time_range, value = tot_avc, fill = 0, convert = FALSE) %>%
data.frame(.)

row.names(mat_tot) <- paste(mat_tot$year, mat_tot$month, mat_tot$day, sep = "_")
mat_tot <- mat_tot[,-c(1, 2, 3)]

xdis <- vegan::vegdist(as.matrix(mat_tot), method = "jaccard")


fit_tot <- aov(glm(tot_avc ~ 0 + year:month:time_range, data = df_tot))
mod_tot <- broom::tidy(TukeyHSD(fit_tot)) %>% 
separate(contrast, into = c("A", "B"), sep = "-", remove = FALSE) %>%
separate(A, into = c("A_year", "A_month", "A_time_range"), sep = ":") %>%
separate(B, into = c("B_year", "B_month", "B_time_range"), sep = ":") %>%
filter(A_year == "2020", A_month == B_month, A_time_range == B_time_range) %>%
mutate(significance = adj.p.value <= 0.05) %>%
arrange(desc(significance)) %>%
data.frame(.)

(hmap_2b <- mod_tot %>% 
  # filter(A_lockdown == TRUE, B_region == "KRAKOW") %>%
  mutate(
    # animal = factor(animal, levels = animal_sorted, 
    # labels = str_replace_all(animal_sorted, "\\.", " ")), 
    # # comparison = as.factor(comparison), 
    significance = as.factor(significance)
    # B_time = factor(B_time, levels = time_stamps),
    # A_time = factor(A_time, levels = time_stamps)) %>%
    # term = as.factor(term)
    ) %>% 
  na.omit(.) %>%  
  ggplot(aes(y = A_time_range, x = "")) +
  geom_raster(alpha = 1, aes(fill = estimate), show.legend = T) +
  geom_tile(alpha = 1, aes(colour = significance), fill = NA, 
    width = 0.9, height = 0.9, size = 1) +
  facet_grid(A_month ~., switch = "y", space = "free", 
    scale = "free") +
  scale_fill_gradient2(low = "darkblue", 
                      high = "darkred", mid = "white",
                      na.value = "darkgrey"
                      # breaks = c(-4, -2, 0, 2, 4), 
                      # labels = c(paste(c(-4, -2, 0, 2, 4), sep = ""))
                      ) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA), guide = FALSE) +
  theme_set +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 0.8, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        # strip.text.y = element_text(angle = 180), 
        strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.8, size = 5)
        ) +
  labs(x = "", fill = "Mean difference", y = "")) +
  ggsave("/klaster/work/abasak/project_SB/figures/stats_krakow_tot_avc.png",
         dpi = 600, units = "in",
         width = 5, height = 9, limitsize = FALSE, 
         device = "png")


# Analysis of the traffic intensity
df_traffic <- read.table("/klaster/work/abasak/project_SB/data/dataset_traffic.txt", 
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
separate(DATA, into = c("day", "month", "year"), sep = "/", convert = FALSE) %>%
filter(month %in% as.character(c("01","02","03","04","05","06"))) %>%
mutate(HOUR = as.numeric(str_replace_all(as.character(HOUR), ":", ""))/100) %>%
mutate(time_range = factor(cut(HOUR, breaks = c(0, 400, 800, 1200, 1600, 2000, 2400),
  include.lowest = TRUE, labels = time_stamps), levels = time_stamps),
month = str_replace_all(month, "0", "")) %>%
mutate(year = str_replace_all(year, "20", "2020")) %>%
mutate(year = str_replace_all(year, "19", "2019")) %>%
data.frame(.)

# Carpent data and stratify at time level
df_traffic_l <- df_traffic %>% 
select(-HOUR) %>%
gather(key = "street", value = "intensity", -year, -month, -day, -time_range) %>%
group_by(year, month, day, time_range) %>%
summarise(intensity = (sum(intensity))) %>%
mutate(group  = paste(year, month, time_range, sep = ":")) %>%
data.frame(.)

mat_traffic <- df_traffic_l %>%
select(-group) %>%
spread(key = time_range, value = intensity, fill = 0, convert = FALSE) %>%
mutate(day = str_replace_all(day, "^0", "")) %>%
add_column(group = paste(.$year, .$month, .$day, sep = "_"), .before = 1) %>%
data.frame(.)

row.names(mat_traffic) <- mat_traffic$group
mat_traffic <- mat_traffic[,-c(1:4)]

idx <- which(row.names(mat_traffic) %in% row.names(mat_tot))

y_dis <- vegan::vegdist(as.matrix(mat_traffic[idx,]), method = "jaccard")

set.seed(1)
mantel_obj <- vegan::mantel(xdis, y_dis, method = "pearson", permutations = 1000)

fit <- aov(glm(log10(intensity) ~ 0 + group, data = df_traffic_l))
mod <- broom::tidy(TukeyHSD(fit)) %>% 
separate(contrast, into = c("A", "B"), sep = "-", remove = FALSE) %>%
separate(A, into = c("A_year", "A_month", "A_time"), sep = ":") %>%
separate(B, into = c("B_year", "B_month", "B_time"), sep = ":") %>%
filter(B_year == 2019, A_time == B_time, B_month == A_month) %>%
data.frame(.)

idx <- match(stat_df$contrast, mod$contrast)
int_df <- stat_df %>%
mutate(traffic_estimates = mod$estimate[idx],
  traffic_pvals = mod$adj.p.value[idx],
  lockdown = B_month %in% c("3", "4", "5")) %>%
data.frame(.)

idx <- match(mod_tot$contrast, mod$contrast)
int_df_tot <- mod_tot %>%
mutate(traffic_estimates = mod$estimate[idx],
  traffic_pvals = mod$adj.p.value[idx],
  lockdown = B_month %in% c("3", "4", "5")) %>%
data.frame(.)

# Heatmap for traffic intensity
(hmap_2c <- mod[idx,] %>% 
  # filter(A_lockdown == TRUE, B_region == "KRAKOW") %>%
  mutate(
    # animal = factor(animal, levels = animal_sorted, 
    # labels = str_replace_all(animal_sorted, "\\.", " ")), 
    # # comparison = as.factor(comparison), 
    significance = as.factor(adj.p.value <= 0.05)
    # B_time = factor(B_time, levels = time_stamps),
    # A_time = factor(A_time, levels = time_stamps)) %>%
    # term = as.factor(term)
    ) %>% 
  na.omit(.) %>%  
  ggplot(aes(y = A_time, x = "")) +
  geom_raster(alpha = 1, aes(fill = estimate), show.legend = T) +
  geom_tile(alpha = 1, aes(colour = significance), fill = NA, 
    width = 0.9, height = 0.9, size = 1) +
  facet_grid(A_month ~., switch = "y", space = "free", 
    scale = "free") +
  scale_fill_gradient2(low = "darkblue", 
                      high = "darkred", mid = "white",
                      na.value = "darkgrey"
                      # breaks = c(-4, -2, 0, 2, 4), 
                      # labels = c(paste(c(-4, -2, 0, 2, 4), sep = ""))
                      ) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA), guide = FALSE) +
  theme_set +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 0.8, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        # strip.text.y = element_text(angle = 180), 
        strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.8, size = 5)
        ) +
  labs(x = "", fill = "Mean difference", y = "")) +
  ggsave("/klaster/work/abasak/project_SB/figures/stats_krakow_traffic.png",
         dpi = 600, units = "in",
         width = 5, height = 9, limitsize = FALSE, 
         device = "png")

mark_x <- quantile(abs(int_df_tot$estimate[!is.na(int_df_tot$estimate)]), 0.1)
mark_y <- quantile(abs(int_df_tot$traffic_estimates[!is.na(int_df_tot$traffic_estimates)]), 0.01)

#
idx <- which(int_df_tot$lockdown == TRUE)
int_cor <- cor.test(int_df_tot$traffic_estimates[idx], int_df_tot$estimate[idx], method = "spearman")


# Correlation between the mean difference AVC 
int_df_tot %>%
na.omit(.) %>%
mutate(mark = abs(estimate) > mark_x) %>%
ggplot(aes(x = estimate, y = traffic_estimates)) +
geom_vline(xintercept = c(-mark_x, mark_x), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_hline(yintercept = c(-mark_y, mark_y), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_point(aes(fill = A_time_range, colour = lockdown, size = mark, alpha = lockdown), shape = 23) +
# facet_grid(.~ A_month, 
#   space = "free", 
#   switch = "x", 
#   scale = "free") +
scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"), guide = FALSE) +
scale_fill_manual(values = col.pal, guide = FALSE) +
scale_size_manual(values = c(`FALSE` = 0.5, `TRUE` = 2), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.8), guide = FALSE) +
labs(x = "MeanDifference AVC vs2019", y = "MeanDifference Traffic vs2019") +
  theme_set +
  theme(axis.text.x = element_text(size = 14, angle = 90), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 12), 
        # strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.5, size = 5)) +
  ggsave(paste0("/klaster/work/abasak/project_SB/figures/integrated_traffic_avc_estimates.png"),
         dpi = 600, units = "in",
         width = 5, height = 5.5, limitsize = FALSE, 
         bg = "transparent",
         device = "png")





# Mantel correlation between AVC in krakow and traffic intensity

# Matrix for distance calculation

# Concatenate the heatmaps
hmap2 <- cowplot::plot_grid(
  hmap_2a + theme_void() + 
  theme(
    strip.text = element_text(size = 0), 
    legend.position = "none", panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_text(size = 12)), NULL,
  hmap_2b + theme_void() + 
  theme(strip.text = element_text(size = 0),
    legend.position = "none", panel.spacing = unit(0.1, "lines")), NULL,
  hmap_2c + theme_void() + 
  theme(strip.text = element_text(size = 0),
    legend.position = "none", panel.spacing = unit(0.1, "lines")),
  nrow = 1, 
  ncol = 5, 
  axis = "tblr", 
  rel_widths = c(0.25, 0.0025, 0.01, 0.0025, 0.01), 
  rel_height = c(1, 1, 1, 1, 1)
)

ggsave(hmap2, 
  filename = "/klaster/work/abasak/project_SB/figures/hmap_summary_krakow.png",
         dpi = 600, units = "in",
         width = 5.2, height = 5, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

# Multi-variate analysis
mds_x <- cmdscale(xdis, k = 2, eig = TRUE)
mds_y <- cmdscale(y_dis, k = 2, eig = TRUE)

var_traffic <- round(100*(mds_y$eig/sum(mds_y$eig)), 2)
var_avc <- round(100*(mds_x$eig/sum(mds_x$eig)), 2)
idx <- match(row.names(mds_x$points), row.names(mds_x$points))

# MDS-MDS plot
mds_obj <- cbind.data.frame(avc_comp1 = mds_x$points[,1], 
  avc_comp2 = mds_x$points[,2], 
  traffic_comp1 = mds_y$points[idx,1],
  traffic_comp2 = mds_y$points[idx,2]) %>%
add_column(id = row.names(.), .before = 1) %>%
separate(id, into = c("year", "month", "day"), sep = "_", convert = FALSE) %>%
mutate(year = as.factor(year),
  month = as.factor(month),
  day = as.factor(day),
  lockdown = as.factor(month %in% c("3", "4", "5") & year == "2020")) %>%
data.frame(.)

(mds_c1 <- mds_obj %>%
ggplot(aes(x = avc_comp1, y = traffic_comp1)) +
geom_vline(xintercept = c(0), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_hline(yintercept = c(0), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_point(aes(fill = year, 
  colour = lockdown, 
  size = lockdown, 
  alpha = lockdown), shape = 22) +
stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"), guide = FALSE) +
scale_fill_manual(values = c(`2019` = "darkgrey", `2020` = "black"), guide = FALSE) +
scale_size_manual(values = c(`FALSE` = 0.5, `TRUE` = 2), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.8), guide = FALSE) +
labs(x = paste0("AVC component 1: ", var_avc[1], "%"), 
  y = paste0("Traffic component 1: ", var_traffic[1], "%")) +
  theme_set +
  theme(axis.text.x = element_text(size = 14, angle = 90), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 12), 
        # strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.5, size = 5))) +
  ggsave(paste0("/klaster/work/abasak/project_SB/figures/integrated_traffic_avc_MDS1.png"),
         dpi = 600, units = "in",
         width = 5, height = 5.5, limitsize = FALSE, 
         bg = "transparent",
         device = "png")


(mds_c2 <- mds_obj %>%
ggplot(aes(y = avc_comp2, x = traffic_comp2)) +
geom_vline(xintercept = c(0), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_hline(yintercept = c(0), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_point(aes(fill = year, 
  colour = lockdown, 
  size = lockdown, 
  alpha = lockdown), shape = 22) +
stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"), guide = FALSE) +
scale_fill_manual(values = c(`2019` = "darkgrey", `2020` = "black"), guide = FALSE) +
scale_size_manual(values = c(`FALSE` = 0.5, `TRUE` = 2), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.8), guide = FALSE) +
labs(y = paste0("AVC component 2: ", var_avc[2], "%"), 
  x = paste0("Traffic component 2: ", var_traffic[2], "%")) +
  theme_set +
  theme(axis.text.x = element_text(size = 14, angle = 90), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 12), 
        # strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.5, size = 5))) +
  ggsave(paste0("/klaster/work/abasak/project_SB/figures/integrated_traffic_avc_MDS2.png"),
         dpi = 600, units = "in",
         width = 5, height = 5.5, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

(mds_avc <- mds_obj %>%
ggplot(aes(y = avc_comp2, x = avc_comp1)) +
geom_vline(xintercept = c(0), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_hline(yintercept = c(0), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_point(aes(fill = year, 
  colour = lockdown, 
  size = lockdown, 
  alpha = lockdown), shape = 22) +
stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"), guide = FALSE) +
scale_fill_manual(values = c(`2019` = "darkgrey", `2020` = "black"), guide = FALSE) +
scale_size_manual(values = c(`FALSE` = 0.5, `TRUE` = 2), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.8), guide = FALSE) +
labs(x = paste0("AVC component 1: ", var_avc[1], "%"), 
  y = paste0("AVC component 2: ", var_avc[2], "%")) +
  theme_set +
  theme(axis.text.x = element_text(size = 14, angle = 90), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 12), 
        # strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.5, size = 5))) +
  ggsave(paste0("/klaster/work/abasak/project_SB/figures/kakow_avc_MDS.png"),
         dpi = 600, units = "in",
         width = 5, height = 5.5, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

(mds_traffic <- mds_obj %>%
ggplot(aes(x = traffic_comp2, y = traffic_comp1)) +
geom_vline(xintercept = c(0), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_hline(yintercept = c(0), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_point(aes(fill = year, 
  colour = lockdown, 
  size = lockdown, 
  alpha = lockdown), shape = 22) +
stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"), guide = FALSE) +
scale_fill_manual(values = c(`2019` = "darkgrey", `2020` = "black"), guide = FALSE) +
scale_size_manual(values = c(`FALSE` = 0.5, `TRUE` = 2), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.8), guide = FALSE) +
labs(y = paste0("Traffic component 1: ", var_traffic[1], "%"), 
  x = paste0("Traffic component 2: ", var_traffic[2], "%")) +
  theme_set +
  theme(axis.text.x = element_text(size = 14, angle = 90), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 12), 
        # strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.5, size = 5))) +
  ggsave(paste0("/klaster/work/abasak/project_SB/figures/kakow_traffic_MDS.png"),
         dpi = 600, units = "in",
         width = 5, height = 5.5, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

mds_composite_int <- cowplot::plot_grid(
  mds_traffic,
  mds_c1,
  mds_c2,
  mds_avc,
  nrow = 2, ncol = 2, "tblr",
  rel_widths = c(1, 1), 
  rel_height = c(1, 1),
  scale = c(1, 1),
  greedy = FALSE
)

ggsave(mds_composite_int, 
  filename = "/klaster/work/abasak/project_SB/figures/mds_summary_krakow.png",
         dpi = 600, units = "in",
         width = 8, height = 8, limitsize = FALSE, 
         bg = "transparent",
         device = "png")



# ============
# TEmporal trend
mds_composite_time <- cowplot::plot_grid(
  plot_mds12_time,
  # plot_mds23_time,
  cca_plot_trend + ggtitle(""),
  nrow = 1, ncol = 2, "tblr",
  rel_widths = c(1, 1), 
  rel_height = c(1, 1),
  scale = c(1, 1),
  greedy = FALSE
)

ggsave(mds_composite_time, 
  filename = "/klaster/work/abasak/project_SB/figures/mds_summary_temporal_trend.png",
         dpi = 600, units = "in",
         width = 8, height = 4, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

# Spatio tempora for individual AVC
mds_composite_l <- cowplot::plot_grid(
  plot_mds12_l,
  # plot_mds23_l,
  cca_plot_l + ggtitle(""),
  nrow = 1, ncol = 2, "tblr",
  rel_widths = c(1, 1), 
  rel_height = c(1, 1),
  scale = c(1, 1),
  greedy = FALSE
)

ggsave(mds_composite_l, 
  filename = "/klaster/work/abasak/project_SB/figures/mds_loc_summary_temporal.png",
         dpi = 600, units = "in",
         width = 8, height = 4, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

mds_composite <- cowplot::plot_grid(
  plot_mds12,
  # plot_mds23,
  cca_plot + ggtitle(""),
  nrow = 1, ncol = 2, "tblr",
  rel_widths = c(1, 1), 
  rel_height = c(1, 1),
  scale = c(1, 1),
  greedy = FALSE
)

ggsave(mds_composite, 
  filename = "/klaster/work/abasak/project_SB/figures/mds_summary_temporal.png",
         dpi = 600, units = "in",
         width = 8, height = 4, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

# Concatenate the heatmaps
hmap1 <- cowplot::plot_grid(
  hmap1a + theme_void() + 
  theme(
    strip.text = element_text(size = 0), 
    legend.position = "none", panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_text(size = 12)), NULL,
  hmap1b + theme_void() + 
  theme(strip.text = element_text(size = 0),
    legend.position = "none", panel.spacing = unit(0.1, "lines")),
  nrow = 1, 
  ncol = 3, 
  axis = "tblr", 
  rel_widths = c(1, 0.025, 0.075), 
  rel_height = c(1, 1, 1)
)

ggsave(hmap1, 
  filename = "/klaster/work/abasak/project_SB/figures/hmap_summary_temporal.png",
         dpi = 600, units = "in",
         width = 10, height = 4, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

# Save tables
write.table(stat_df, file = "/klaster/work/abasak/project_SB/statistics/temporal_thsd_summary.txt", 
           quote = F, sep = "\t")
write.table(variable, file = "/klaster/work/abasak/project_SB/statistics/temporal_cca_summary.txt",
           quote = F, sep = "\t")


# END OF SCRIPT
sessionInfo()