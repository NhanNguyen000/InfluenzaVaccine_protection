rm(list = ls())

library(tidyverse)
library(venn)

get.DAs <- function(resLimma) {
  # Aim: extract the DAPs/DAMs with p.value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_DA <- resLimma$p.value %>% as.data.frame %>% 
    select(matches("log2")) %>% filter(. <0.05)
  
  return(res_DA)
}

get.tstat <- function(resLimma) {
  # Aim: extract the t-statistic value from linnear comparision comparison
  # following the get.limmaRes result
  
  res_tstat <- resLimma$t %>% as.data.frame %>% 
    select(matches("log2"))
  
  return(res_tstat)
}

get.plotDat_clusterRow <- function(plotDat, colName, valColumn) {
  # Aim: make the heatmap data with dendrogram clustering order, so can have clustered heatmap with ggplot
  # this function comvert the plot data (long format) to wide format 
  # with column name from the given "colName" column, value from given "valName" colum, and rowname from fixed "valName" column
  
  # outcome: plotDat_order - plot data long format with order in valName, so it will show clustering with ggplot heatmap
  
  plotDat_wide <- plotDat %>% select(c("valName", colName, valColumn)) %>% 
    pivot_wider(names_from = colName, values_from = valColumn) %>% 
    column_to_rownames("valName")
  
  plotDat_dendrogram <- as.dendrogram(hclust(d = dist(plotDat_wide)))
  plotDat_denOrder <- order.dendrogram(plotDat_dendrogram)
  
  plotDat_order <- plotDat %>%
    mutate(valName = factor(valName, levels = valName[plotDat_denOrder], ordered = TRUE))
  
  return(plotDat_order)
}

# load data --------------------------------------------------------
load("resPro_abT1_abFC.RData")

# significant proteins / metabolites ---------------------------------------------------------
## DAs data ------------------
DAs <- resPro_abT1_abFC %>% 
  lapply(function(x) x %>% lapply(function(y) y %>% lapply(function(z) get.DAs(z))))

# DAs for abT1
DAs_abT1 <- DAs %>%
  lapply(function(x) x %>% 
           lapply(function(y) y$abT1 %>% as.data.frame %>% rownames_to_column("Assay")) %>% 
           imap(~mutate(.x, strain = .y)) %>% purrr::reduce(full_join)) %>%
  imap(~mutate(.x, season = .y)) %>% purrr::reduce(full_join)

# DAs for abFC
DAs_abFC <- DAs %>%
  lapply(function(x) x %>% 
           lapply(function(y) y$abFC %>% as.data.frame %>% rownames_to_column("Assay")) %>% 
           imap(~mutate(.x, strain = .y)) %>% purrr::reduce(full_join)) %>%
  imap(~mutate(.x, season = .y)) %>% purrr::reduce(full_join)

## tstat data ------------------
tstat_dat <- resPro_abT1_abFC %>%
  lapply(function(x) x %>% lapply(function(y) y %>% lapply(function(z) get.tstat(z))))

# tstat for abT1
tstat_abT1 <- tstat_dat %>%
  lapply(function(x) x %>% 
           lapply(function(y) y$abT1 %>% as.data.frame %>% rownames_to_column("Assay")) %>% 
           imap(~mutate(.x, strain = .y)) %>% purrr::reduce(full_join)) %>%
  imap(~mutate(.x, season = .y)) %>% purrr::reduce(full_join)

# tstat for abFC
tstat_abFC <- tstat_dat %>%
  lapply(function(x) x %>% 
           lapply(function(y) y$abFC %>% as.data.frame %>% rownames_to_column("Assay")) %>% 
           imap(~mutate(.x, strain = .y)) %>% purrr::reduce(full_join)) %>%
  imap(~mutate(.x, season = .y)) %>% purrr::reduce(full_join)

## merge DAs and tstat data together --------
abFC <- DAs_abFC %>% rename("pvalue" = "abFC_log2") %>%
  full_join(tstat_abFC %>% rename("tstat" = "abFC_log2"))

abT1 <- DAs_abT1 %>% rename("pvalue" = "T1_log2") %>%
  full_join(tstat_abT1 %>% rename("tstat" = "T1_log2"))

dat_total <- abFC %>% mutate(type = "abFC") %>%
  full_join(abT1 %>% mutate(type = "abT1")) %>% 
  mutate(group = paste0(type, "_", season, "_", strain)) %>%
  mutate(group = gsub("iMED_|ZirFlu_", "", group))

# check certain proteins ---------------------------
load("selected_DAPs.RData")
selected_DAs
selected_DAs <- c("CD83", "CXCL8")
tstat_abT1_DAPs <- tstat_abT1 %>% filter(Assay %in% selected_DAs)
tstat_abFC_DAPs <- tstat_abFC %>% filter(Assay %in% selected_DAs)

tstat_abFC_DAPs_v2 <- tstat_abFC_DAPs %>% 
  mutate(abFC = paste0(season, "_", strain)) %>% 
  rename("value" = "abFC_log2") %>% select(-strain, -season) %>% 
  relocate(value, .after = abFC)

circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(10, 10))

col_mat <- ifelse(tstat_abFC_DAPs_v2$value>0, "red", "blue")
chordDiagram(
  x = tstat_abFC_DAPs_v2,
  col = col_mat,
  transparency = 0.5,
  directional = 1,
  annotationTrack = "grid",
  link.largest.ontop = TRUE
)

circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 7, 
      labels = sector.index, 
      facing = "reverse.clockwise", 
      cex = 0.8
    )
  }
)

# both abFC and abT1 -----------
selected_DAs <- c("CD83", "CXCL8", 
                  #"TNFRSF11A", "LSP1", 
                   # "TNFSF13C", "VEGFA", "PGF", "KRT19", "MLN", "SIGLEC1", "CCL7",
                   # "OSCAR", "SERPINB8", "KLRB1", "CD4", "TNFSF10",
                   # "IL6", "CD160"
                  )
selected_DAs <- c("CD83", "CXCL8", "TNFRSF11A", "LSP1")
dat_temp <- dat_total %>% filter(Assay %in% selected_DAs)

dat_temp_v2 <- dat_temp %>% 
  select(Assay, group, tstat) 

col_df <- dat_temp %>% 
  select(Assay, group, tstat, pvalue) %>%
  mutate(link_color = ifelse(is.na(pvalue), "gray",
                             ifelse(tstat > 0, "red", "blue"))) %>%
  select(-tstat, -pvalue)

circos.clear()
circos.par(start.degree = 0, gap.degree = 2, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar=c(0,0,0,0))

chordDiagram(
  x = dat_temp_v2,
  col = col_df$link_color,
  transparency = 0.5,
  directional = 1,
  annotationTrack = "grid",
  link.largest.ontop = TRUE
)

circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 5, 
      labels = sector.index, 
      facing = "reverse.clockwise", 
      cex = 0.6
    )
  }
)

circos.clear()
circos.par(start.degree = 180, gap.degree = 2, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar=c(0,0,0,0))
dat_temp_v3 <- dat_temp %>% 
  select(Assay, group, tstat) 

sector <- sort(unique(c(dat_temp_v3$Assay, dat_temp_v3$group)))
df.groups <- structure(c(rep("abFC", 14), rep("abT1", 14), rep("protein", 4)), names = sector)
chordDiagram(
  x = dat_temp_v3,
  order = sector,
  col = col_df$link_color,
  transparency = 0.5,
  directional = 1,
  group = df.groups,
  annotationTrack = "grid",
  link.largest.ontop = TRUE
)


circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3, 
      labels = gsub("abFC_|abT1_", "", sector.index), 
      facing = "reverse.clockwise", 
      cex = 0.5
    )
  }
)
for (i in unique(df.groups)) {
  highlight.sector(grep(i, get.all.sector.index(), value = T), track.index = 1, 
                   #col = Type_Cols[match(i, labels(Type_Cols))], 
                   text = i, 
                   text.vjust = 2, padding = c(-.2, 0, -.5, 0))
}

# check this test code later -------------------


mycolor2 <- viridis(2, alpha = 1, begin = 0, end = 1, option = "D")
mycolor2 <- mycolor2[sample(1:2)]

chordDiagram(
  x = tstat_abFC_DAPs_v2, 
 # grid.col = mycolor2,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)


circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
      minor.ticks = 1, 
      major.tick.percentage = 0.5,
      labels.niceFacing = FALSE)
  }
)


# Load the circlize library
library(circlize)

# test code --------
data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/13_AdjacencyDirectedWeighted.csv", header=TRUE)

# short names
colnames(data) <- c("Africa", "East Asia", "Europe", "Latin Ame.",   "North Ame.",   "Oceania", "South Asia", "South East Asia", "Soviet Union", "West.Asia")
rownames(data) <- colnames(data)

# I need a long format
data_long <- data %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname)
# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor <- viridis(10, alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- mycolor[sample(1:10)]

chordDiagram(
  x = data_long, 
  grid.col = mycolor,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
      minor.ticks = 1, 
      major.tick.percentage = 0.5,
      labels.niceFacing = FALSE)
  }
)
