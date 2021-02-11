rm(list = ls())
setwd("/klaster/work/abasak/project_SB/")
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
# df$species <- str_replace(df$species, pattern = "\302\240", replacement = "")
# # df$Gmina <- str_replace(df$Gmina, pattern = "\303\263", replacement = "o")
# # df$Gmina <- str_replace(df$Gmina, pattern = "\305\201", replacement = "o")
# # df$species <- str_replace(df$species, pattern = "\\<U+00A0>", replacement = "")

# # df <- df %>% select(-10,-11)
# colnames(df)[6] <- "reported_time"
# df$species <- str_replace(df$species, "bird$", "birds")
# write.table(df, "./data/dataset_carpented.txt", sep = "\t", quote = FALSE)

df <- read.table("/klaster/work/abasak/project_SB/data/dataset_carpented.txt", 
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
filter(year %in% c(2019,2020), month != 12)

# Read temperature data
# df_temp <- read.table("/klaster/work/abasak/project_SB/data/dataset_temperature.txt", 
#                  header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
# gather(key = "year", convert = FALSE, value = "temperature", -month, -day) %>%
# mutate(year = as.numeric(str_replace_all(year, "X", ""))) %>%
# mutate(group = paste(day, month, year, sep = "_")) %>% 
# na.omit(.) %>%
# # group_by(group) %>%
# # summarise(mtemp = mean(temperature)) %>%
# data.frame(., stringsAsFactors = FALSE)

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
# mutate(temperature = df_temp$temperature[match(.$temp_id, df_temp$group)])

# Merge temperature data here

df_wide <- df %>% group_by(Gmina, year, month, species) %>% 
summarise(vals = n()) %>%
spread(key = species, value = "vals", fill = 0, convert = FALSE) %>%
mutate(group = as.factor(paste(Gmina, year, month, sep = "_"))) %>%
data.frame(., stringsAsFactors = FALSE)

# Data summary
df_wide %>%
group_by(Gmina, year) %>%
summarise(n = n()) %>%
data.frame(.) %>%
write.table(., "./output/data_spatial_summary.txt", 
  sep = "\t", 
quote = FALSE, 
  row.names = FALSE)

mat <- apply((df_wide %>% select(-Gmina, -year, -month, -group)), 2, function(x) x/sum(x))
row.names(mat) <- paste(df_wide[,1], df_wide[,2], df_wide[,3], sep = "_")
row.names(df_wide) <- paste(df_wide[,1], df_wide[,2], df_wide[,3], sep = "_")

# d <- 1 - cor(t(mat))
d <- vegan::vegdist(mat, "jaccard")
hc <- hclust(as.dist(d), "ward.D2")
data_sorted <- row.names(mat)[hc$order]
mds_obj <- cmdscale(as.dist(d), eig = TRUE, k = 3)
variance <- round(100*mds_obj$eig/sum(mds_obj$eig), 2) 

d <- 1 - cor(mat)
hc <- hclust(as.dist(d), "ward.D2")
animals_sorted <- row.names(t(mat))[hc$order]

# Plot NMDS -- lockdown
(plot_mds12 <- cbind.data.frame(df_wide[,c(1,2,3)], mat, 
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
  geom_point(aes(colour = lockdown, alpha = lockdown, fill = year, shape = Gmina), size = 3) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  # scale_colour_manual(values = col.idx) +
  scale_fill_manual(values = c("darkblue","indianred")) +
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
  ggsave("/klaster/work/abasak/project_SB/figures/MDS12_analysis.png", 
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
  geom_point(aes(colour = lockdown, alpha = lockdown, fill = year, shape = Gmina), size = 3) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  # scale_colour_manual(values = col.idx) +
  scale_fill_manual(values = c("darkblue", "indianred")) +
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
  ggsave("/klaster/work/abasak/project_SB/figures/MDS23_analysis.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

# Plot NMDS -- location
(plot_mds12_l <- cbind.data.frame(df_wide[,c(1,2,3)], mat, 
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
  geom_point(aes(alpha = lockdown, size = lockdown, colour = Gmina, shape = year)) +
  stat_ellipse(aes(group = Gmina, colour = Gmina), 
               lty = "dashed", size = 1, 
               type = "norm") +
  scale_colour_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_fill_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_shape_manual(values = c(1,19)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 3), guide = FALSE) +
  labs(x = paste0("MDS1 - ", variance[1], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS12_loc_analysis.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

(plot_mds23_l <- cbind.data.frame(df_wide[,c(1,2,3)], mat, 
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
  geom_point(aes(alpha = lockdown, size = lockdown, colour = Gmina, shape = year)) +
  stat_ellipse(aes(group = Gmina, colour = Gmina), 
               lty = "dashed", size = 1, 
               type = "norm") +
  scale_colour_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_fill_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_shape_manual(values = c(1,19)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 3), guide = FALSE) +
  labs(x = paste0("MDS3 - ", variance[3], " %"), 
       y = paste0("MDS2 - ", variance[2], " %"), colour = "", shape = "", fill = "") +
  theme_set +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 3),
        axis.text = element_text(size = 12),
        axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 14))) +
  ggsave("/klaster/work/abasak/project_SB/figures/MDS23_loc_analysis.png", 
         dpi = 600, units = "in", 
         width = 4, height = 4, limitsize = FALSE, bg = "transparent",
         device = "png")

# Write out
cbind.data.frame(df_wide[,c(1,2,3)], mat, 
                 mds1 = mds_obj$points[,1], 
                 mds2 = mds_obj$points[,2]) %>% 
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020)
  ) %>% 
  write.table(., file = "./statistics/mds_analysis_spatio.txt", quote = FALSE, sep = "\t")
# Cannonical analysis -- Region and lockdown

# cca_obj <- vegan::capscale(formula = mat ~ year + Condition(Gmina + month), 
#   data = df_wide)
# aov_cca <- vegan::anova.cca(cca_obj)
set.seed(1)
cca_obj <- vegan::capscale(formula = mat ~ Gmina + Condition(month + year), 
  data = df_wide, dist = "jaccard", perm = 1000, sqrt.dist = TRUE)
anova_obj <- vegan::anova.cca(cca_obj)

# Fetch stats
variance <- round(100*cca_obj$CCA$eig/sum(cca_obj$CCA$eig), 2)
pval <- anova_obj$`Pr(>F)`[1]
con_var <- round(100*cca_obj$CCA$tot.chi/cca_obj$tot.chi, 2)
chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                       row.names = c("total", "constrianed", "unconstrained"))

ti <- paste0("(AVC) ~ Region + Condition(month + year), Variance explained = ", con_var, " % ; p.value = ", round(pval, 3), "*")


# CCA analysis
(cca_plot <- cbind.data.frame(df_wide[,c(1,2,3)], mat, 
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
  geom_point(aes(colour = Gmina, shape = year, alpha = lockdown, size = lockdown)) +
  stat_ellipse(aes(group = Gmina, colour = Gmina), 
               lty = "dashed", size = 1, 
               type = "norm") +
  scale_colour_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_shape_manual(values = c(1,19)) +
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
  ggsave("/klaster/work/abasak/project_SB/figures/cca_analysis.png", 
         dpi = 600, units = "in", 
         width = 5, height = 5, limitsize = FALSE, bg = "transparent",
         device = "png")

# Effect of lockdown
df_wide$lockdown <- factor(ifelse(df_wide$month %in% c(3, 4, 5) & df_wide$year == 2020, "lockdown", ""), 
  levels = c("", "lockdown"))
cca_obj <- vegan::capscale(formula = mat ~ lockdown:Gmina + Condition(Gmina + year + month), 
  data = df_wide, dist = "jaccard", sqrt.dist = TRUE, perm = 1000)
anova_obj <- vegan::anova.cca(cca_obj)

variance <- round(100*cca_obj$CCA$eig/sum(cca_obj$CCA$eig), 2)
pval <- anova_obj$`Pr(>F)`[1]
con_var <- round(100*cca_obj$CCA$tot.chi/cca_obj$tot.chi, 2)
chis <- c(cca_obj$tot.chi, cca_obj$CCA$tot.chi, cca_obj$CA$tot.chi)
variable <- data.frame(inertia = chis, proportion = chis/chis[1], 
                       row.names = c("total", "constrianed", "unconstrained"))

ti <- paste0("(AVC) ~ lockdown + Condition(month + year), Variance explained = ", con_var, " % ; p.value = ", round(pval, 3), "*")

# Pairwise PERMANOVA for the factors
df_wide$group <- as.factor(paste(df_wide$Gmina, df_wide$lockdown, sep = "_"))
pwm_mod <- RVAideMemoire::pairwise.perm.manova(vegan::vegdist(mat, "jaccard"), df_wide$group, 
  test = "Wilks", 
  nperm = 1000, 
  p.method = "fdr")

# Obtain the pairwise statistics for the comparison
pw_stats <- as.data.frame(pwm_mod$p.value) %>%
add_column(A = row.names(.), .before = 1) %>%
gather(key = "B", value = "FDR", convert = FALSE, -A) %>%
na.omit(.) %>%
separate(A, into = c("A_location", "A_lockdown"), sep = "_", convert = FALSE) %>%
separate(B, into = c("B_location", "B_lockdown"), sep = "_", convert = FALSE) %>%
mutate(significance = FDR < 0.05) %>%
arrange(desc(significance))

write.table(pw_stats, file = "./statistics/pw_perm_spatio.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot the significance as heatmap
pw_stats %>%
mutate(A_location = as.factor(A_location),
A_lockdown = as.factor(A_lockdown),
B_location = as.factor(B_location),
B_lockdown = as.factor(B_lockdown), 
significance = as.factor(significance),
FDR = -log10(FDR)) %>%
ggplot(aes(x = A_location, y = B_location)) +
geom_raster(aes(fill = FDR)) +
geom_tile(aes(colour = significance), 
  fill = "transparent", 
  size = 0.9,
  width = 0.95,
  height = 0.95) +
facet_grid(A_lockdown ~ B_lockdown, scale = "free", switch = "both", space = "free") +
scale_fill_gradient(low = "white", 
                      high = "darkred", 
                      na.value = "darkgrey") +
scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "darkgrey"), guide = FALSE) +
theme_set +
theme(panel.spacing = unit(0.1, "lines"),
      axis.text.x = element_text(angle = 90, size = 8, hjust = 0.8, vjust = 0.5),
      axis.text.y = element_text(size = 8), 
      # strip.text.y = element_text(angle = 180), 
      strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(hjust = 0.5, vjust = 0.8, size = 5)
      ) +
labs(x = "", fill = "", y = "") +
ggsave("/klaster/work/abasak/project_SB/figures/pw_pnova_stats_spatio.png",
       dpi = 600, units = "in",
       width = 4, height = 5, limitsize = FALSE, 
       bg = "transparent",
       device = "png")

# CCA analysis
(cca_plot_2 <- cbind.data.frame(df_wide[,c(1,2,3)], mat, 
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
  	alpha = lockdown, size = lockdown, fill = Gmina)) +
  stat_ellipse(aes(group = lockdown, colour = lockdown), 
               lty = "dashed", size = 1, 
               type = "norm") +
  scale_fill_manual(values = c("tomato", "olivedrab", "slateblue")) +
  scale_shape_manual(values = c(21,22,23)) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "lightgrey", `TRUE` = "black"), guide = FALSE) +
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
  ggsave("/klaster/work/abasak/project_SB/figures/cca_analysis_lockdown.png", 
         dpi = 600, units = "in", 
         width = 5, height = 5, limitsize = FALSE, bg = "transparent",
         device = "png")

# Relative abundance of the animals observed
mat_df <- cbind.data.frame(df_wide[,c(1,2,3,25)], mat)

df_mat_long <- mat_df %>% 
gather(key = "species", value = "count", 
  convert = FALSE, -Gmina, -year, -month, -group) %>%
# mutate(group = as.factor(paste(Gmina, year, month, sep = "_"))) %>%
data.frame(., stringsAsFactors = FALSE)

animals <- colnames(mat)
res <- lapply(animals, function(x){

    # x = "badger"
    temp <- df_mat_long %>% filter(species == x) %>%
    mutate(year = as.factor(year), month = as.factor(month), 
      group = factor(ifelse(month %in% c(3, 4, 5) & year == 2020, "lockdown", ""), 
  levels = c("", "lockdown")))

    fit <- aov(lm(count ~ 0 + group + Gmina, data = temp))
    mod <- broom::tidy(TukeyHSD(fit)) %>% 
    mutate(animal = x) %>%
    data.frame(., stringsAsFactors = FALSE)

})

stat_df <- do.call(rbind.data.frame, res) %>%
mutate(FDR = p.adjust(adj.p.value, "fdr"))%>%
mutate(significance = FDR <= 0.05) %>%
arrange(desc(significance))


# Plot of relative abundance
(hmap1a <- df_mat_long %>% 
  mutate(species = factor(species, levels = animals_sorted, labels = str_replace_all(animals_sorted, "\\.", " ")),
        id = factor(group, levels = data_sorted),
        Gmina = as.factor(Gmina),
        month = as.factor(month),
        year = as.factor(year),
        lockdown = as.factor(ifelse(month %in% c(3, 4, 5) & year == 2020, "Locked month", "")),
        logra = ifelse(!is.finite(log10(count)), NA, log10(count))) %>% 
  ggplot(aes(x = month, y = species)) +
  geom_raster(alpha = 1, aes(fill = logra), show.legend = T) +
  facet_grid(.~ Gmina + lockdown + year, space = "fixed", scale = "free") +
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
  labs(x = "", fill = "% AVC", y = "")) +
  ggsave("/klaster/work/abasak/project_SB/figures/abundance_HWC.png",
         dpi = 600, units = "in",
         width = 12, height = 8, limitsize = FALSE, 
         device = "png")


# Plot the summary statistics
(hmap1b <- stat_df %>% 
  filter(term == "Gmina") %>%
  mutate(animal = factor(animal, levels = animals_sorted, labels = str_replace_all(animals_sorted, "\\.", " ")), 
    contrast = as.factor(contrast), 
    significance = as.factor(significance)) %>%
    # term = as.factor(term)) %>% 
  ggplot(aes(x = contrast, y = animal)) +
  geom_raster(alpha = 1, aes(fill = -estimate), show.legend = T) +
  geom_tile(alpha = 1, aes(colour = significance), fill = NA, 
    width = 0.9, height = 0.9, size = 1) +
  # facet_grid(.~ , space = "free", scale = "free") +
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
  ggsave("/klaster/work/abasak/project_SB/figures/stats_HWC_gmina.png",
         dpi = 600, units = "in",
         width = 4, height = 9, limitsize = FALSE, 
         device = "png")


(hmap1c <- stat_df %>% 
  filter(term != "Gmina") %>%
  # separate(comparison, into = c("A", "B"), sep = "-") %>%
  # separate(A, into = c("A_year", "A_month"), sep = ":") %>%
  # separate(B, into = c("B_year", "B_month"), sep = ":") %>%
  # filter(B_year == 2020, B_month %in% c(3, 4, 5)) %>%
  mutate(animal = factor(animal, levels = animals_sorted, labels = str_replace_all(animals_sorted, "\\.", " ")), 
    # A_year = as.factor(A_year), 
    # A_month = as.factor(A_month),
    significance = as.factor(significance)) %>%
    # term = as.factor(term)) %>% 
  na.omit(.) %>%
  ggplot(aes(x = contrast, y = animal)) +
  geom_raster(alpha = 1, aes(fill = estimate), show.legend = T) +
  geom_tile(alpha = 1, aes(colour = significance), fill = NA, 
    width = 0.9, height = 0.9, size = 1) +
  # facet_grid(.~ A_year + B_month, space = "free", scale = "free") +
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
  ggsave("/klaster/work/abasak/project_SB/figures/stats_HWC_lockdown.png",
         dpi = 600, units = "in",
         width = 3, height = 9, limitsize = FALSE, 
         device = "png")

# Make boxplots for which are consistently signinficnant on absolute counts
sig_animals <- na.omit(unique(stat_df$animal[stat_df$significance == TRUE])) %>% as.character()
parallel::mclapply(sig_animals, function(x){

  # Plots for bxp
  # x = "wild.boar"
  temp <- df_wide
  colnames(temp)[colnames(temp) == x] <- "id"
  
  temp %>% select(Gmina, year, month, id) %>%
  mutate(year = as.factor(year),
         month = as.factor(month), 
         Gmina = as.factor(Gmina),
        lockdown = c(month %in% c(3, 4, 5) & year == 2020),
        id = 100*id/sum(id)
  ) %>%
  ggplot(aes(x = Gmina, y = id)) +
  geom_hline(yintercept = mean(100*temp$id/sum(temp$id)), colour = "darkred", 
             alpha = 0.5, 
             lty = "solid", size = 0.8) +
  geom_boxplot(outlier.alpha = 0, aes(colour = Gmina), 
               size = 2, 
               fill = NA, 
               alpha = 0.6, 
               position = position_dodge2()) +
  geom_point(stat = "identity", 
              position = "jitter", 
              alpha = 0.6, 
              aes(shape = year, alpha = lockdown, size = lockdown), guide = FALSE) +
  scale_colour_manual(values = c("tomato", "olivedrab", "slateblue"), guide = FALSE) +
  scale_shape_manual(values = c(1,19), guide = FALSE) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.8), guide = FALSE) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 3), guide = FALSE) +
  labs(x = "", y = "%", colour = "", shape = "") +
  theme_set +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16), 
        axis.title.y = element_text(size = 12), 
        # strip.text.x = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(hjust = 0.5, vjust = 0.5, size = 5)) +
  ggsave(paste0("/klaster/work/abasak/project_SB/figures/boxplots/boxplots_", x, ".png"),
         dpi = 600, units = "in",
         width = 4, height = 4.5, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

}, mc.cores = 2)

# Concatenate the multidimensional analysis
mds_composite <- cowplot::plot_grid(
  plot_mds12,
  plot_mds23,
  nrow = 1, ncol = 2, "tblr",
  rel_widths = c(1, 1), 
  rel_height = c(1, 1),
  scale = c(1, 1),
  greedy = FALSE
)

ggsave(mds_composite, 
  filename = "/klaster/work/abasak/project_SB/figures/mds_summary.png",
         dpi = 600, units = "in",
         width = 8, height = 4, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

mds_composite_l <- cowplot::plot_grid(
  plot_mds12_l,
  plot_mds23_l,
  cca_plot + ggtitle(""),
  cca_plot_2 + ggtitle(""),
  nrow = 1, ncol = 4, "tblr",
  rel_widths = c(1, 1, 1, 1), 
  rel_height = c(1, 1, 1, 1),
  scale = c(1, 1, 1, 1),
  greedy = FALSE
)

ggsave(mds_composite_l, 
  filename = "/klaster/work/abasak/project_SB/figures/mds_loc_summary.png",
         dpi = 600, units = "in",
         width = 16, height = 4, limitsize = FALSE, 
         bg = "transparent",
         device = "png")


# PQ figre composite
mds_composite_fig1a <- cowplot::plot_grid(
  # plot_mds12,
  plot_mds12_l,
  cca_plot + ggtitle(""),
  # cca_plot_2 + ggtitle(""),
  nrow = 1, ncol = 2, "tblr",
  rel_widths = c(1, 1), 
  rel_height = c(1, 1),
  scale = c(1, 1),
  greedy = FALSE
)

mds_composite_fig1b <- cowplot::plot_grid(
  plot_mds12,
  # plot_mds12_l,
  # cca_plot + ggtitle(""),
  cca_plot_2 + ggtitle(""),
  nrow = 1, ncol = 2, "tblr",
  rel_widths = c(1, 1), 
  rel_height = c(1, 1),
  scale = c(1, 1),
  greedy = FALSE
)

ggsave(mds_composite_fig1a, 
  filename = "/klaster/work/abasak/project_SB/figures/fig1a_mds_summary.png",
         dpi = 600, units = "in",
         width = 8, height = 4, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

ggsave(mds_composite_fig1b, 
  filename = "/klaster/work/abasak/project_SB/figures/fig1b_mds_summary.png",
         dpi = 600, units = "in",
         width = 8, height = 4, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

# Concatenate the heatmaps
hmap1 <- cowplot::plot_grid(
  hmap1a + theme_void() + theme(legend.position = "none", panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 0)), NULL,
  hmap1b + theme_void() + theme(legend.position = "none", panel.spacing = unit(0.1, "lines"),
    strip.text = element_text(size = 0)), NULL,
  hmap1c + theme_void() + theme(legend.position = "none", panel.spacing = unit(0.1, "lines"),
    strip.text = element_text(size = 0)), 
  nrow = 1, 
  ncol = 5, 
  axis = "tblr", 
  rel_widths = c(1, 0.025, 0.075, 0.035, 0.035), 
  rel_height = c(1, 1, 1, 1, 1)
)

ggsave(hmap1, 
  filename = "/klaster/work/abasak/project_SB/figures/hmap_summary.png",
         dpi = 600, units = "in",
         width = 12, height = 4, limitsize = FALSE, 
         bg = "transparent",
         device = "png")

# Save tables
write.table(stat_df, file = "/klaster/work/abasak/project_SB/statistics/thsd_summary.txt", 
           quote = F, sep = "\t")
write.table(variable, file = "/klaster/work/abasak/project_SB/statistics/cca_summary.txt",
           quote = F, sep = "\t")


# END OF SCRIPT
sessionInfo()