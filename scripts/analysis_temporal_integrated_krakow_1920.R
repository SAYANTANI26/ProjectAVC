rm(list = ls())

# Scripted by @ArpanKumarBasak
# email: arpankbasak@gmail.com; akumarbasak@mpipz.mpg.de; arpan.kumar.basak@doctoral.uj.edu.pl

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


# Factor levels and colours
time_stamps <- c("0000_0400", "0400_0800", "0800_1200", "1200_1600", "1600_2000", "2000_0000")
col.pal <- c("darkblue", "skyblue", "lightyellow", "lightcoral", "lavender", "darkviolet")
# levs <- sapply(unique(df$Gmina), function(x) paste(x, time_stamps, sep = "_"))


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

    set.seed(1)
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

# set.seed(1)
# mantel_obj <- vegan::mantel(xdis, y_dis, method = "pearson", permutations = 1000)

fit <- aov(glm(intensity/1e6 ~ 0 + group, data = df_traffic_l))
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
int_cor <- broom::tidy(cor.test(int_df_tot$traffic_estimates, int_df_tot$estimate, method = "spearman", alternative = "greater")) %>% mutate(test = "overall")
int_cor_idx <- broom::tidy(cor.test(int_df_tot$traffic_estimates[idx], int_df_tot$estimate[idx], method = "spearman", alternative = "less")) %>% mutate(test = "lockdown")

rbind.data.frame(int_cor, int_cor_idx) %>% write.table("./statistics/integrated_correlation.txt", 
quote = FALSE, 
row.names = FALSE, 
sep = "\t")

# Correlation between the mean difference AVC 
int_df_tot %>%
na.omit(.) %>%
mutate(mark = abs(estimate) > mark_x) %>%
ggplot(aes(y = estimate, x = traffic_estimates)) +
geom_vline(xintercept = c(-mark_y, mark_y), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_hline(yintercept = c(-mark_x, mark_x), colour = "darkgrey", lty = "solid", lwd = 0.8) +
geom_abline(intercept = -1, colour = "darkred", alpha = 0.6, lty = "dashed", lwd = 0.8) +
geom_point(aes(fill = A_time_range, colour = lockdown, size = mark, alpha = lockdown), shape = 23) +
# facet_grid(.~ A_month, 
#   space = "free", 
#   switch = "x", 
#   scale = "free") +
scale_colour_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"), guide = FALSE) +
scale_fill_manual(values = col.pal, guide = FALSE) +
scale_size_manual(values = c(`FALSE` = 0.5, `TRUE` = 2), guide = FALSE) +
scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.8), guide = FALSE) +
labs(y = "MeanDifference AVC vs2019", x = "MeanDifference Traffic vs2019") +
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
# Save tables
# write.table(stat_df, file = "/klaster/work/abasak/project_SB/statistics/temporal_thsd_summary.txt", 
#            quote = F, sep = "\t")
# write.table(variable, file = "/klaster/work/abasak/project_SB/statistics/temporal_cca_summary.txt",
#            quote = F, sep = "\t")


# END OF SCRIPT
sessionInfo()