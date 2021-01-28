rm(list = ls())
## Add GG repel
# Scripted by @ArpanKumarBasak
# email: arpankbasak@gmail.com

require(tidyverse)
require(cowplot)
require(ggrepel)

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
filter(year %in% c(2019, 2020), month != 12, Gmina == "KRAKOW")

# Read mobility data
df_mob <- read.csv2("/klaster/work/abasak/project_SB/data/raw/apple_mobility_data.csv", 
                 header = TRUE, sep = ",", stringsAsFactors = FALSE) %>% 
filter(sub.region %in% c("Lesser Poland Province")) %>%
select(-geo_type, -region, -alternative_name, -country) %>%
gather(key = "key", value = "vals", convert = FALSE, -transportation_type, -sub.region) %>%
separate(key, into = c("year", "month", "day"), convert = FALSE, sep = "\\.", remove = FALSE) %>%
mutate(year = as.numeric(str_replace_all(year, "X", "")),
  month = as.numeric(month),
  day = as.numeric(day), 
  vals = as.numeric(vals)
  ) %>%
filter(year %in% c(2019, 2020), month %in% c(1, 2, 3, 4, 5, 6))

# Segregate and carpent the month dataset for integration
month_df <- df_mob %>%
na.omit(.) %>%
group_by(month, transportation_type) %>%
summarise(mobility = mean(vals)) %>%
spread(key = transportation_type, value = mobility, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = FALSE)

# Data carpenting
df <- df %>% mutate(year = as.numeric(year), 
  species = as.factor(species),
  month = as.numeric(month), 
  day = as.numeric(day),
  ID = paste(ID, Gmina, day, month, year, sep = "_"),
  group = as.factor(paste(Gmina, year, month, sep = "_")),
  lockdown = month %in% c(3, 4, 5)
)

# Merge temperature data here
df_wide <- df %>% group_by(Gmina, year, month, day, species) %>% 
summarise(vals = n()) %>%
spread(key = species, value = "vals", fill = 0, convert = FALSE) %>%
mutate(group = as.factor(paste(Gmina, year, month, day, sep = "_"))) %>%
data.frame(., stringsAsFactors = FALSE)

mat <- apply((df_wide %>% select(-Gmina, -year, -month, -day, -group)), 2, function(x) x/sum(x))
row.names(mat) <- df_wide$group
row.names(df_wide) <- df_wide$group

# Compute statistics G-linear model
mat_df <- cbind.data.frame(df_wide[,c(1,2,3,4,22)], mat)

df_mat_long <- mat_df %>% 
gather(key = "species", value = "count", 
  convert = FALSE, -Gmina, -year, -month, -day, -group) %>%
# mutate(group = as.factor(paste(Gmina, year, month, sep = "_"))) %>%
data.frame(., stringsAsFactors = FALSE)

animals <- colnames(mat)
res <- lapply(animals, function(x){

    # x = "badger"
    temp <- df_mat_long %>% 
    filter(species == x) %>%
    mutate(year = as.factor(year), month = as.factor(month), 
      # lockdown = factor(ifelse(month %in% c(5,3,4), "Not lockdown", "lockdown"), 
        # levels = c("Not lockdown", "lockdown"),
      count = as.numeric(count)) %>%
    mutate(group = as.factor(paste(year, month, sep = "_")))

    fit <- aov(lm(count ~ 0 + year:month, data = temp))
    mod <- broom::tidy(TukeyHSD(fit)) %>% 
    mutate(animal = x) %>%
    data.frame(., stringsAsFactors = FALSE)

})

stat_df <- do.call(rbind.data.frame, res) %>%
# mutate(FDR = p.adjust(adj.p.value, "fdr"))%>%
separate(contrast, into = c("A", "B"), sep = "-", remove = FALSE) %>%
separate(A, into = c("A_year", "A_month"), sep = ":") %>%
separate(B, into = c("B_year", "B_month"), sep = ":") %>%
filter(A_year == "2020", A_month == B_month) %>%
mutate(significance = adj.p.value <= 0.05) %>%
arrange(desc(significance))

animal_mat <- stat_df %>% select(contrast, estimate, animal) %>%
spread(key = animal, value = estimate, fill = 0, convert = FALSE) %>%
data.frame(., stringsAsFactors = F)

# Compute correlation between the animal variations within comparison
d <- 1 - cor(animal_mat[,-1])
hc <- hclust(as.dist(d), "ward.D2")
animal_sorted <- colnames(animal_mat[,-1])[hc$order]

idx <- match(stat_df$A_month, month_df$month)
stat_df <- stat_df %>% mutate(driving = month_df$driving[idx], 
  walking = month_df$walking[idx])


# heatmap for year month effect within Krakow
(hmap1a <- stat_df %>% 
  mutate(
    animal = factor(animal, levels = animal_sorted, labels = str_replace_all(animal_sorted, "\\.", " ")), 
    significance = as.factor(significance), 
    A_month = as.factor(A_month)) %>%
  ggplot(aes(x = A_month, y = animal)) +
  geom_raster(alpha = 1, aes(fill = estimate), show.legend = T) +
  geom_tile(alpha = 1, aes(colour = significance), fill = NA, 
    width = 0.9, height = 0.9, size = 1) +
  # facet_grid(.~ B_month, space = "free", scale = "free") +
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
  ggsave("/klaster/work/abasak/project_SB/figures/stats_krakow_lockdown.png",
         dpi = 600, units = "in",
         width = 4, height = 8, limitsize = FALSE, 
         device = "png")

fdr_line <- min(-log10(stat_df$adj.p.value[stat_df$adj.p.value < 0.05]))

# Comutation vs animal effect size --add GGrepel for animals
(com_sig <- stat_df %>% 
  gather(key = "commutation", value = "RI", driving, walking) %>%
  mutate(animal = factor(animal, levels = animal_sorted, labels = str_replace_all(animal_sorted, "\\.", " ")), 
    significance = as.factor(abs(estimate) > 0.01), 
    lockdown = as.factor(A_month %in% c(3, 4, 5)),
    lockd = as.factor(B_month %in% c(3, 4, 5)),
    trend = sign(estimate)) %>%
  mutate(animal_text = ifelse(significance == "TRUE", 
    as.character(animal), 
    "")) %>%
  ggplot(aes(y = log10(RI), x = abs(estimate), label = animal_text)) +
  geom_vline(xintercept = c(0, 0.0001, 0.001, 0.01), lty = "solid", colour = "darkgrey", lwd = 1) +
  geom_point(aes(colour = lockdown, size = significance, alpha = significance, fill = commutation), shape = 21) +
  geom_text_repel(colour = "black", size = 2) +
  facet_grid(lockd ~., space = "free_x", scale = "free_x", switch = "y") +
  # scale_shape_manual(values = c(`walking` = 21, `driving` = 22), guide = FALSE) +
  scale_size_manual(values = c(`TRUE` = 1, `FALSE` = .5), guide = FALSE) +
  scale_alpha_manual(values = c(`TRUE` = 0.6, `FALSE` = 0.3), guide = FALSE) +
  scale_fill_manual(values = c(`walking` = "darkmagenta", `driving` = "darkgreen"), guide = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "darkgrey"), guide = FALSE) +
  theme_set +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 0.8, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        # strip.text.y = element_text(angle = 180), 
        strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.8, size = 5)
        ) +
  labs(x = "", y = "", colour = "", size = "", shape = "", alpha = "")) +
  ggsave("/klaster/work/abasak/project_SB/figures/stats_commutation.png",
         dpi = 600, units = "in",
         width = 4, height = 5, limitsize = FALSE, 
         device = "png")

# Comutation for animals for the month
(strip_traffic <- month_df %>% 
  # filter(month != 4) %>%
  gather(key = "key", value = "value", convert = FALSE, -month) %>%
  mutate(key = as.factor(key),
    month = as.factor(month),
    value = log10(value)) %>%
  ggplot(aes(x = month, y = key)) +
  # geom_raster(alpha = 1, aes(fill = value), show.legend = T) +
  geom_tile(alpha = 1, colour = "black", aes(fill = value)) +
  scale_fill_gradient(low = "white", 
                      high = "darkgrey", 
                      breaks = c(0, 1, 1.5, 2, 2.5), 
                      labels = c(paste(c(0, 1, 1.5, 2, 2.5), sep = ""))
                      ) +
  theme_set +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 0.8, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        # strip.text.y = element_text(angle = 180), 
        strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.8, size = 5)
        ) +
  labs(x = "", fill = "Total traffic intensity", y = "")) +
  ggsave("/klaster/work/abasak/project_SB/figures/strip_krakow_lockdown.png",
         dpi = 600, units = "in",
         width = 4, height = 2, limitsize = FALSE, 
         device = "png")


# Multi variate
d <- 1-cor(t(mat), method = "spearman")
mds_obj <- cmdscale(d, k = 3, eig = TRUE)
variance <- round(100*(mds_obj$eig/sum(mds_obj$eig)), 2)


# Plot MDS
(plot_mds12 <- cbind.data.frame(df_wide[,c(1,2,3)], mat, 
                 mds1 = mds_obj$points[,1], 
                 mds2 = mds_obj$points[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  ggplot(aes(x = mds1, y = mds2, group = lockdown)) +
  geom_hline(yintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_vline(xintercept = 0.0, colour = "darkgrey", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_point(aes(
    fill = year, 
    alpha = lockdown, 
    colour = lockdown, 
    shape = year), size = 3) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  # scale_colour_manual(values = col.idx) +
  scale_fill_manual(values = c(`2019` = "darkblue",`2020` = "indianred")) +
  scale_colour_manual(values = c(`FALSE` = "lightgrey", `TRUE` = "black"), guide = FALSE) +
  scale_shape_manual(values = c(22,23)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  labs(x = paste0("MDS1 - ", variance[1], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS12_analysis_krakow.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")


  (plot_mds23 <- cbind.data.frame(df_wide[,c(1,2,3)], mat, 
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
  geom_point(aes(
    fill = year, 
    alpha = lockdown, 
    colour = lockdown, 
    shape = year), size = 3) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  # scale_colour_manual(values = col.idx) +
  scale_fill_manual(values = c(`2019` = "darkblue",`2020` = "indianred")) +
  scale_colour_manual(values = c(`FALSE` = "lightgrey", `TRUE` = "black"), guide = FALSE) +
  scale_shape_manual(values = c(22,23)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  labs(x = paste0("MDS3 - ", variance[3], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS23_analysis_krakow.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

set.seed(1)
df_wide$lockdown <- as.factor(df_wide$month %in% c(3, 4, 5) & df_wide$year == 2020)

cca_obj <- vegan::capscale(formula = mat ~ year + month + Condition(day), 
  data = df_wide, dist = "jaccard", perm = 1000, sqrt.dist = TRUE)
anova_obj <- vegan::anova.cca(cca_obj)

# Fetch stats
variance <- round(100*cca_obj$CCA$eig/sum(cca_obj$CCA$eig), 2)
pval <- anova_obj$`Pr(>F)`[1]
con_var <- round(100*cca_obj$CCA$tot.chi/cca_obj$tot.chi, 2)
chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                       row.names = c("total", "constrianed", "unconstrained"))

ti <- paste0("(HWC) ~ Year:Month + Condition(days), Variance explained = ", con_var, " % ; p.value = ", round(pval, 3), "*")

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
  geom_point(aes(
    fill = year, 
    alpha = lockdown, 
    colour = lockdown, 
    shape = year), size = 3) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  # scale_colour_manual(values = col.idx) +
  scale_fill_manual(values = c(`2019` = "darkblue",`2020` = "indianred")) +
  scale_colour_manual(values = c(`FALSE` = "lightgrey", `TRUE` = "black"), guide = FALSE) +
  scale_shape_manual(values = c(22,23)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  labs(x = paste0("CCA1 - ", variance[1], " %"), 
       y = paste0("CCA2 - ", variance[2], " %"), 
       colour = "",
       shape = "") +
  theme_set + 
  theme(title = element_text(size = 4),
        legend.text = element_text(size = 4),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/cca_analysis_mobility.png", 
         dpi = 600, units = "in", 
         width = 5, height = 5, limitsize = FALSE, bg = "transparent",
         device = "png")


# Concatenate the heatmaps
# hmap1 <- cowplot::plot_grid(
#   hmap1a + theme_void() + theme(legend.position = "none", panel.spacing = unit(0.1, "lines"),
#     axis.text.y = element_text(size = 12),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank()), NULL,
#   strip_traffic + theme_void() + theme(legend.position = "none", panel.spacing = unit(0.1, "lines"),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank()),
#   nrow = 3, 
#   ncol = 1, 
#   axis = "tblr", 
#   rel_height = c(1, 0.025, 0.075), 
#   rel_width = c(1, 1, 1)
# )

hmap1a + theme_void() + 
theme(legend.position = "none", panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()) +
ggsave( filename = "/klaster/work/abasak/project_SB/figures/hmap_mobility1.png",
         dpi = 600, units = "in",
         width = 3, height = 6, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

strip_traffic + theme_void() + 
theme(legend.position = "none", panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()) +
ggsave( filename = "/klaster/work/abasak/project_SB/figures/hmap_strip.png",
         dpi = 600, units = "in",
         width = 2.5, height = 1, limitsize = FALSE, 
         bg = "transparent",
         device = "png")


# Save tables
# write.table(stat_df, file = "/klaster/work/abasak/project_SB/statistics/thsd_summary.txt", 
#            quote = F, sep = "\t")
# write.table(variable, file = "/klaster/work/abasak/project_SB/statistics/cca_summary.txt",
#            quote = F, sep = "\t")


# END OF SCRIPT
sessionInfo()