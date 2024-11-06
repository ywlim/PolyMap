# PolyMap analysis of COVID antibody library and antigen library interaction mapping

Author: Yoong Wearn Lim and Cindy Yeh

Date: 11/05/2024

## Goal

To analyze sequencing data from PolyMap experiments to map antibody library and antigen library interaction.

Antibody libraries used:

* Antibody library from COVID convalescent donors from 2020 (antibody_lib_2020_donors_rep1 and antibody_lib_2020_donors_rep2)
* Antibody library from COVID convalescent and vaccinated donors from 2023 (antibody_lib_2023_donors_rep1 and antibody_lib_2023_donors_rep2)

Antigen library used: 

* COVID spike antigen library consisting of 18 unique sequences representing phylogenetically distinct clades (and non-target control antigens)

## Analysis

### Generating antigen barcode map

The first part of the analysis is performed in Python.

The antigen, cell barcode, and UMI for each sequencing read can be retrieved using the script `Barcodes/dropseq_custombc_cellbc_UMI.py`. This script searches for known barcodes linked to antigens used in the Polymap experiment. Antigen/barcode file is `ab_ag_barcodes/spike_2024_bcs.txt`. It allows for up to two mismatches in the barcode and requires SciPy to run.

```
# If SciPy is not installed
python3 -m pip install scipy

# To generate the antigen-barcode-UMI file
cd Barcodes
# python3 dropseq_custombc_cellbc_UMI.py <R1.fastq with antigen barcode> <R2.fastq with UMI/cell barcode> <antigen/barcode file> <output file>
python3 dropseq_custombc_cellbc_UMI.py R1_antigen.fastq R2_umi_cell_bc.fastq spike_2024_bcs.txt sample_bc_map.txt

```

### Generating antibody barcode map

The antibody CDR3H, cell barcode, and UMI for each sequencing read can be retrieved using the script `denovo_cdr3_cellbc_UMI.py`. This script extracts the CDR3H amino acid sequence and filters out reads with a high frequency of low-quality nucleotides. CDR3H sequences can be generated as previously described ([Adler 2017]([https://www.tandfonline.com/doi/full/10.1080/19420862.2017.1371386#d1e1138)).

```
# To generate the antibody-barcode-UMI file
cd Barcodes
# python3 denovo_cdr3_cellbc_UMI.py <R1.fastq with UMI/cell barcode> <CDR3 file> <output file>
python3 denovo_cdr3_cellbc_UMI.py R1_umi_cell_bc.fastq sample_CDR3s.txt sample_bc_map.txt 
```


### Antigen


All following analyses are performed in R.


```r
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(cowplot)
library(kneedle)

filter = dplyr::filter
select = dplyr::select
summarise = dplyr::summarise
arrange = dplyr::arrange
rename = dplyr::rename

custom_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"))

theme_set(custom_theme)
```



Here we process the antigen bc_map files to generate a data frame of cell barcode and antigen identity. We only consider reads with unique UMI, and we only consider cell barcodes that pass the knee test.


```r
info <- read_delim("info.txt", delim = "\t")
info <- info %>% 
  mutate(sample = as.character(sample))

info_ant <- info %>% filter(type == "antigen")

# pdf(file = "antigen_knee_plots.pdf", width = 5, height = 4)

df_antigen <- lapply(info_ant$sample, function(samplex) {
  mini_ant <- read_delim(paste0("Barcodes/", samplex, "_bc_map.txt"), delim = "\t") %>% 
    transmute(antigen = ag_variant, 
              cell_barcode,
              antigen_UMI = UMI) %>% 
    mutate(cell_barcode_umi = paste0(cell_barcode, "_", antigen_UMI)) %>% 
    filter(!duplicated(cell_barcode_umi)) 
  
  # arrange cell bc in descending order
  # calculate cumulative fraction read 
  bc_ant <- mini_ant %>%
    group_by(cell_barcode) %>%
    tally() %>% 
    mutate(fraction = n / sum(n)) %>%
    arrange(desc(fraction))
  bc_ant$cumsum <- cumsum(bc_ant$fraction)
  bc_ant$cell_id <- 1:nrow(bc_ant)
  
  # find knee
  knee_ant <- kneedle(bc_ant$cell_id, y = bc_ant$cumsum)
  
  # round down the knee
  knee_ant <- floor(knee_ant[1])

  print(ggplot(bc_ant, aes(x = cell_id, y = cumsum)) +
    geom_point() +
    geom_vline(xintercept = knee_ant) +
    annotate(geom = "text", label = paste0("knee=", knee_ant), 
             x = -Inf, y = Inf, hjust = -0.5, vjust = 2) +
   labs(title = paste0("Antigen (", samplex, ")"), x = "Cell", y = "Cumulative fraction of reads"))

  # remove bad cells as determined by knee test
  good_bc_ant <- bc_ant$cell_barcode[bc_ant$cell_id <= knee_ant]

  mini_ant <- mini_ant %>% 
    filter(cell_barcode %in% good_bc_ant) %>% 
    select(antigen, cell_barcode)
  
  # get antigen ID
  mini_ant <- mini_ant %>% 
    group_by(cell_barcode, antigen) %>% 
    tally() 
  colnames(mini_ant)[3] <- "read"
 
  tot <- mini_ant %>% 
    ungroup() %>% 
    group_by(cell_barcode) %>% 
    summarise(total = sum(read))

  mini_ant <- left_join(mini_ant, tot)
  
  mini_ant <- mini_ant %>% 
    mutate(percent = read / total * 100)
  mini_ant$sample <- samplex
  
  return(mini_ant)
})
# dev.off()

df_antigen <- do.call(rbind, df_antigen)

saveRDS(df_antigen, file = "df_antigen.RDS")
```

#### Top antigen

Here we determine the antigen identity for each cell barcode, requiring each cell barcode to be associated with ≥10 antigen reads and ≥90% single antigen reads to filter out droplets containing multiple cells. 

```r
antigen_name <- read_delim("antigen_names.txt")

antigens <- c("CoV2_WT",
"CoV2_20E_EU1",
"CoV2_20H_betaV2",
"CoV2_20I_alphaV1",
"CoV2_20J_gammaV3",
"CoV2_21A_delta",
"CoV2_21B_kappa",
"CoV2_21C_epsilon",
"CoV2_21D_eta",
"CoV2_21F_iota",
"CoV2_21G_lambda",
"CoV2_21H_mu",
"CoV2_21K_omicron",
"CoV2_21M_omicron",
"CoV2_22A_omicron",
"CoV2_22C_omicron",
"CoV2_22D_BA.2.75",
"CoV2_23B_XBB.1.16",
"CoV1_WT",
"PD1",
"CTLA4",
"BFP")

df_antigen_top <- df_antigen %>% 
  filter(antigen %in% antigens) %>% 
  group_by(sample, cell_barcode) %>% 
  top_n(n = 1, wt = percent) 

df_antigen_top <- left_join(df_antigen_top, antigen_name)

# Filter for antigens with >= 10 reads and >= 90% of total antigen reads per cell
df_antigen_top2 <- df_antigen_top %>% 
  filter(total >= 10, percent >= 90)
```

### Antibody

Here we process the antibody bc_map files to generate a data frame of cell barcode and antibody identity. We only consider reads with unique UMI. We don't run the knee test on antibody because we expect some cell barcodes to not have lots of antibody reads (negative target antigens).

```r
info_ab <- info %>% filter(type == "antibody")

df_ab <- lapply(info_ab$sample, function(samplex) {
  mini_ab <- read_delim(paste0("Barcodes/", samplex, "_bc_map.txt"), delim = "\t") 

  mini_ab <- mini_ab %>% 
    transmute(antibody = cdr3, 
              cell_barcode,
              ab_UMI = UMI) %>% 
    mutate(cell_barcode_umi = paste0(cell_barcode, "_", ab_UMI)) %>% 
    filter(!duplicated(cell_barcode_umi)) %>% 
    select(antibody, cell_barcode)
  
  mini_ab <- mini_ab %>% 
    filter(!grepl("\\*", antibody),
           !grepl("!", antibody))
  
  # get ab ID
  mini_ab <- mini_ab %>% 
    group_by(cell_barcode, antibody) %>% 
    tally() 
  colnames(mini_ab)[3] <- "read"
  tot <- mini_ab %>% 
    ungroup() %>% 
    group_by(cell_barcode) %>% 
    summarise(total = sum(read))

  mini_ab <- left_join(mini_ab, tot)
  mini_ab <- mini_ab %>% 
    mutate(percent = read / total * 100)
  mini_ab$sample <- samplex
  
  return(mini_ab)
})
df_ab <- do.call(rbind, df_ab)
```

### Merging antibody and antigen

Merging the antibody and antigen data frame by cell barcode.


```r
info_ab <- info_ab %>% select(sample, name)
info_ant <- info_ant %>% select(sample, name)

df_antigen_top3 <- df_antigen_top2 %>% 
  ungroup() %>% 
  rename(read_antigen = read,
         total_read_antigen = total,
         percent_antigen = percent) %>% 
  left_join(info_ant) %>% 
  select(-sample) 

df_ab <- df_ab %>% 
  filter(sample %in% info_ab$sample) %>% 
  ungroup() %>% 
  rename(read_antibody = read,
         total_read_antibody = total,
         percent_antibody = percent) %>% 
  left_join(info_ab) %>% 
  select(-sample) 

df <- inner_join(df_ab, df_antigen_top3, by = c("name", "cell_barcode"))

write.table(df, file = "count_antibody_antigen.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

#### Number of  cells

Good cell = cell with antigen ID and antibody reads.

```r
temp <- df %>% 
  select(name, cell_barcode, antigen2) %>% 
  unique()
n <- as.data.frame(table(temp$name, temp$antigen2))
colnames(n) <- c("name", "antigen", "n")

tot_cell <- n %>% group_by(name) %>% summarise(sum = sum(n))
tot_cell$name <- factor(tot_cell$name, levels = unique(info$name))

n$antigen <- factor(n$antigen, levels = antigen_name$antigen2)

pdf(file = "num_good_cells.pdf", width = 11, height = 6)
ggplot(n, aes(x = antigen, y = n, fill = antigen, label = n)) +
  geom_col() +
  geom_text() +
  geom_text(data = tot_cell, aes(x = Inf, y = Inf, label = paste0("n=", sum)), inherit.aes = FALSE, vjust = 2, hjust = 2) +
  facet_wrap(~ name, ncol= 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(fill = "none") +
  labs(title = "Number of cells", 
       subtitle = "Filtered for cells >= 10 antigen reads, top antigen >= 90% of reads")
dev.off()
```

Starting from this point, we will focus on antigens with >= 10 cells in both replicates.

```r
n_filtered <- n %>% 
  filter(n >= 10)

# which antigen passes the filter for all 4 samples?
n_num_exp <- data.frame(table(n_filtered$antigen))
good_antigens <- n_num_exp$Var1[n_num_exp$Freq == 4]

antigen_name_good <- antigen_name %>% 
  filter(antigen2 %in% good_antigens)

df <- df %>% 
  filter(antigen2 %in% antigen_name_good$antigen2)
```

### Normalized across antigens

For each antibody, we compute the total antibody reads across all cells. We then calculate the percentage of antibody reads associated with each antigen. We divide the percent antibody reads for each antigen by the number of cells expressing that antigen, resulting in normalized percent reads reflecting each antibody’s relative binding across different antigen variants. 

```r
# num of ab reads for a given antigen
ab_read <- df %>% 
  ungroup() %>% 
  group_by(name, antibody, antigen2) %>% 
  summarise(read_ab = sum(read_antibody)) %>% 
  rename(antigen = antigen2)

# total # reads for each ab
ab_read_total <- ab_read %>% group_by(name, antibody) %>% 
  summarise(total = sum(read_ab))

# for each ab, what is the relative binding to each antigen
ab_read <- left_join(ab_read, ab_read_total) 
ab_read <- ab_read %>% mutate(percent = read_ab / total * 100)

# normalize by number of cells per antigen
ab_read_all <- ab_read %>% ungroup() %>% 
  left_join(n)
ab_read_all <- ab_read_all %>%  mutate(percent_scaled = percent / n)
saveRDS(ab_read_all, file = "ab_read_all.RDS")
```

### Replicate analysis

Here we want to identify antibodies with sufficient read counts to provide reliable and reproducible binding profiles across experimental replicates. Steps:

* To establish a minimum read count threshold, iteratively vary the minimum total antibody read count threshold from 0 to 5000
* For each threshold value, filter the antibodies to include only those with total reads greater than or equal to the threshold, then calculate the Pearson correlation coefficient of the normalized percent reads for all antibody-antigen pairs between two experimental replicates
* Identify the smallest threshold value at which the Pearson correlation coefficient reached at least 0.8
* Exclude antibodies with total reads below this threshold 


```r
# join by replicates to create min_ab_total column
ab_read_rep1 <- ab_read_all %>%
  filter(grepl("rep1", name)) %>%
  mutate(name2 = gsub("_rep1", "", name)) %>%
  select(name2, antibody, antigen, total, percent_scaled) %>%
  rename(name = name2, R1_ab_total = total, percent_scaled_R1 = percent_scaled)

# divide by total num
ab_read_rep1_sum <- ab_read_rep1 %>%
  group_by(name) %>%
  summarise(sum_total = sum(R1_ab_total))

ab_read_rep1 <- ab_read_rep1 %>%
  left_join(ab_read_rep1_sum) %>%
  mutate(R1_ab_frac = R1_ab_total / sum_total) %>%
  select(-sum_total)

ab_read_rep2 <- ab_read_all %>%
  filter(grepl("rep2", name)) %>%
  mutate(name2 = gsub("_rep2", "", name)) %>%
  select(name2, antibody, antigen, total, percent_scaled) %>%
  rename(name = name2, R2_ab_total = total, percent_scaled_R2 = percent_scaled)

# divide by total num
ab_read_rep2_sum <- ab_read_rep2 %>%
  group_by(name) %>%
  summarise(sum_total = sum(R2_ab_total))
ab_read_rep2 <- ab_read_rep2 %>%
  left_join(ab_read_rep2_sum) %>%
  mutate(R2_ab_frac = R2_ab_total / sum_total) %>%
  select(-sum_total)

ab_read_all_merge <- ab_read_rep1 %>%
  inner_join(ab_read_rep2) %>%
  mutate(min_ab_total = pmin(R1_ab_total, R2_ab_total),
         min_ab_frac = pmin(R1_ab_frac, R2_ab_frac))
```

#### Correlation

```r
thresholds = seq(0, 5000, by = 1) 

r2_df <- lapply(thresholds, function(thresholdx) {
  # filter by min_ab_total, which is the smaller of R1_ab_total and R2_ab_total
  temp <- ab_read_all_merge %>%
    filter(min_ab_total >= thresholdx)
  
  corr <- temp %>%
    group_by(name) %>%
    summarise(correlation = cor(percent_scaled_R1, percent_scaled_R2, use = "pairwise.complete.obs"),
              unique_abs = length(unique(antibody)))
  
  corr$threshold <- thresholdx 

  return(corr)
})
r2_df <- do.call(rbind, r2_df)
```

```r
r2_df2 <- r2_df %>%
  filter(correlation >= 0.8) %>%
  group_by(name) %>%
  slice_min(correlation, with_ties = FALSE)

# color by unique antibodies
pdf(file = "correlation_by_read_threshold.pdf", width = 8, height = 3)
ggplot(r2_df, aes(x = threshold, y = correlation, color = unique_abs)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red",
                       limits = c(0, 200),
                       oob = scales::squish,
                       breaks = c(0, 50, 100, 150, 200),
                       labels = c("0", "50","100", "150",">=200")) +
  geom_vline(data = r2_df2, aes(xintercept = threshold)) + 
  geom_text(data = r2_df2, y = -Inf, aes(x = plot_nums[10], label = paste0("# unique abs @ Rsq=0.8: ", unique_abs)), color = "black", size = 4, hjust = 1, vjust = -1) +
  facet_wrap(~ name, scales = "free", nrow = 1) + 
  coord_cartesian(xlim = c(0, plot_nums[10])) + 
  labs(title = "Correlation by minimum read threshold") +
  ylab("Correlation coefficient (R-squared)")
dev.off()
```

#### Good antibodies

Filter for good antibodies that pass the total antibody read threshold determined above.

```r
threshold1 <- r2_df2 %>% 
  filter(name == "antibody_lib_2020_donors") %>% 
  pull(threshold)

df1 <- ab_read_all_merge %>% 
  filter(name == "antibody_lib_2020_donors",
         min_ab_total >= threshold1)
  
threshold2 <- r2_df2 %>% 
  filter(name == "antibody_lib_2023_donors") %>% 
  pull(threshold)

df2 <- ab_read_all_merge %>% 
  filter(name == "antibody_lib_2023_donors",
         min_ab_total >= threshold2)

saveRDS(df1, file = "df_antibody_lib_2020_donors.RDS")
saveRDS(df2, file = "df_antibody_lib_2023_donors.RDS")
```

#### Scatterplots

```r
library(viridis)

correlation_final <- function(dfx, x_pos, y_pos) {
  rsq <- round(cor(dfx$percent_scaled_R1, dfx$percent_scaled_R2), 3)
  model <- lm(dfx$percent_scaled_R2 ~ dfx$percent_scaled_R1) 
  slope = round(coef(model)[2], 2)
  intercept = round(coef(model)[1], 2)

  ggplot(dfx, aes(x = percent_scaled_R1, y = percent_scaled_R2)) + 
    geom_pointdensity() +
    scale_color_viridis() +
    xlab("Polymap score (replicate 1)") + 
    ylab("Polymap score (replicate 2)") +
    geom_abline(slope = coef(model)[2], intercept = coef(model)[1], color = "grey38") + 
    annotate(geom = "text", label = paste("Pearson\ncorrelation=\n", rsq), 
             x = x_pos, y = y_pos, color = "black", hjust = 0.5, vjust = 1) +
    labs(color = "Point\ndensity")
}

pdf(file = "replicate_correlation_antibody_lib_2020_donors.pdf", width = 4, height = 3)
correlation_final(df1, 0.2, 0.5)
dev.off()

pdf(file = "replicate_correlation_antibody_lib_2023_donors.pdf", width = 4, height = 3)
correlation_final(df2, 0.1, 0.6)
dev.off()
```

### Heatmap (hierarchical clustering)

Heatmaps to visualize antibody-antigen polymap scores.

```r
antigens <- antigen_name_good$antigen2

hierarchical_clustering <- function(dfx) {
  mini <- dfx
  mini_wide <- mini %>% 
    filter(!antigen %in% c("BFP", "CTLA-4", "PD-1")) %>% 
    select(antibody, antigen, percent_scaled_R1) %>% 
    spread(antigen, percent_scaled_R1)
  mini_wide <- column_to_rownames(mini_wide, var = "antibody") %>% 
    as.matrix()
  mini_wide[is.na(mini_wide)] <- 0
  mini_wide <- mini_wide[, antigens[!antigens %in% c("BFP", "CTLA-4", "PD-1")]]
  
  # Compute distance matrices
  row_dist <- dist(mini_wide)

  # Perform hierarchical clustering
  row_hclust <- hclust(row_dist)
  
  # Reorder the matrix based on clustering results
  mini_ordered <- mini_wide[row_hclust$order, ]
  
  mini$antibody <- factor(mini$antibody, levels = rev(rownames(mini_ordered)))
  
  mini$antigen <- factor(mini$antigen, levels = antigens)
  
  return(mini)
}
```

```r
theme_heatmap <- theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_blank())
        
make_heatmap <- function(dfx, fill_upper_limit_R1, fill_upper_limit_R2, y_axis_text_size) {
  
  if (is.na(fill_upper_limit_R1)) {
    fill_upper_limit_R1 <- max(dfx$percent_scaled_R1)
  }
  if (is.na(fill_upper_limit_R2)) {
    fill_upper_limit_R2 <- max(dfx$percent_scaled_R2)
  }
  
  p1 <- ggplot(dfx, aes(x = antigen, y = antibody, fill = percent_scaled_R1)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue",
                        limits = c(0, fill_upper_limit_R1),
                        oob = scales::squish) +
    theme_heatmap +
    theme(axis.text.y = element_text(size = ifelse(is.null(y_axis_text_size), element_text()$size, y_axis_text_size))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(y = "Antibody", x = "Antigen", title = "Replicate 1", fill = "Normalized\n% reads")
  
  p2 <- ggplot(dfx, aes(x = antigen, y = antibody, fill = percent_scaled_R2)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue",
                        limits = c(0, fill_upper_limit_R2),
                        oob = scales::squish) +
    theme_heatmap +
    theme(axis.text.y = element_text(size = ifelse(is.null(y_axis_text_size), element_text()$size, y_axis_text_size))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(y = "Antibody", x = "Antigen", title = "Replicate 2", fill = "Normalized\n% reads")
  p <- plot_grid(p1, p2, ncol = 2,  align = "h", axis = "bt", rel_widths = c(1, 1))
  return(p)
}
```

```r
# all antibodies
df1_ordered <- hierarchical_clustering(df1)
df2_ordered <- hierarchical_clustering(df2)

# top 40 antibodies
# no need to make one for the 2023 library because there are only/exactly 40 antibodies in total
temp <- df1 %>% 
  select(antibody, R1_ab_total) %>% 
  unique() %>% 
  arrange(desc(R1_ab_total)) %>% 
  head(40)
df1_top <- df1 %>% filter(antibody %in% temp$antibody)

df1_top_ordered <- hierarchical_clustering(df1_top)
```

```r
pdf(file = "heatmap_2020_lib_all.pdf", width = 12, height = 11)
make_heatmap(df1_ordered, NA, NA, 6)
dev.off()

pdf(file = "heatmap_2023_lib_all.pdf", width = 12, height = 6.5)
make_heatmap(df2_ordered, 0.3, 0.39, 9)
dev.off()

pdf(file = "heatmap_2020_lib_top40.pdf", width = 12, height = 6.5)
make_heatmap(df1_top_ordered, 0.45, 0.4, 9)
dev.off()
```


