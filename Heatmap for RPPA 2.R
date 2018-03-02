# Creates Large Heatmap Demonstrating the logFC of proteins over time for each Cell Line and Treatment

library(limma)
library(lumi)
library(gplots)
library(readxl)
library(tidyverse)
library(pvclust)
library(lmdme)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(plyr)
library(ggrepel)
source("FunctionsFromBrian.R")

setwd("~/Desktop/RPPA PLS-DR")

# crossref <- read.csv("paths.csv")

# read in RPPA data saved as csv from excel
# CHANGE THIS TO RUN LOCALLY
# allData <- read.csv("~/BTSync/MASS_SPEC/CrazyRPPA/RPPA_EGFR_IGF1R.csv")

allData <- read.csv("dasat_bms_rppa.csv", stringsAsFactors = FALSE)

# grab only "targets" data. Limma uses the language of "targets" and "expression matrix". "targets"
# refers to annotation data that describes the experimental design. The "expression matrix" refers
# to the expression data for which the linear model will be fit. 

targets <- allData[,c(1:5)]
head(allData[,c(6:151)])


# limma relies heavily on factors for experimental design, so levels must be explictly set
library(reshape2)
# grabbing expression matrix data here
rppa.raw <- as.matrix(allData[,c(7:151)])
# limma expects an expression matrix to have genes as rows and samples as columns, so we take the transverse
rppa.trnsv <- t(rppa.raw)
class(rppa.trnsv) <- "numeric"

# this will induce NAs ("ND" is not able to be cast as a numeric)

# "ND"  means unsure of value 
# "0" means the analyte is truly not there0

# limma expects log transformed data
rppa.BC <- log2(rppa.trnsv)

# turn matrix into object that limma needs
row.names(rppa.BC)[145] <- paste(row.names(rppa.BC)[150], "ALT", " ")


# Rather than trying to run all analysis at once, we are running two independent paired analyses
# One to extract the consistent response to stimulation
# One to extract the consistent response to treatment with either OSI or BMS

# Make a targets data frame with data for the treatment effects
# treat_targets <- targets %>%
#   # filter() assumes "AND" behavior when you use ","
#   # Any data we grab must satisfy ALL filters specified in the command
#   filter(Stimulus %in% c("IGF"), # grab stimulated data only 
#          `Inhibitor - EGFR` == "0", # grab data that wasn't treated with EGFR inhibitor
#          `Inhibitor - IGF1R` %in% c("0", "BMS", "OSI")) %>% # grab IGF1R inhibitor data
#   # As mentioned earlier, limma needs factors to behave properly, doubly so for paired analysis
#   mutate(Treatment = factor(`Inhibitor - IGF1R`, levels = c("0", "BMS", "OSI")),
#          `Cell Line` = factor(`Cell Line`, levels = unique(`Cell Line`))
#   )
# 
# treat_targets <- targets  %>% filter(Time.Point %in% c("15 min"))

# treat_targets <- targets %>% mutate('Cell Line' = factor('Cell.Line', levels = unique('Cell.Line')), Treatment = factor('Treatment', levels = c("BMS5", "BMS7", "BMS5+7", "BMS7+Dasat", "Control", "Dasat")))
treat_rppa <- rppa.BC
targets$Target <- make.names(with(targets, paste(Cell.Line, gsub(" ", "", Treatment, fixed = TRUE), gsub(" ", "", Time.Point, fixed = TRUE), sep = ".")))
# targets$Target <- with(targets, paste(Cell.Line, Treatment, Time.Point, sep = "."))

lev <- levels(as.factor(targets$Target))
design <- model.matrix(~ targets$Target)
colnames(design) <- lev
fit <-lmFit(treat_rppa, design)

### END LOADING DESIGN and TARGETS

### BEGIN REORG LARGE HEAT MAP
celllines <- levels(as.factor(make.names(targets$Cell.Line)))
treatments <- make.names(levels(as.factor(gsub(" ", "", targets[which(targets$Treatment != "Control"),]$Treatment, fixed = TRUE))))
timepoints <- levels(as.factor(gsub(" ", "", targets[which(targets$Time.Point != "15 min"),]$Time.Point, fixed = TRUE)))

# storage.df <- data.frame(T1HR = numeric(), T24HR = numeric(), T3HR = numeric(), T8HR = numeric(), AveExpr = numeric(), F = numeric(), P.Value = numeric(), adj.P.Val = numeric(), Treatment = character(), Cell.Line = character(), stringsAsFactors = FALSE)
storage.df <- data.frame()
newcontrast <- vector()
for(trt in treatments)
{
  for(cls in celllines)
  {
    newcontrast <- vector()
    for(tp in timepoints)
    {
      nexttime <- paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
      newcontrast <- c(newcontrast, nexttime)
    }
    cont.dif <- makeContrasts(contrasts = newcontrast, levels = design)
    fit2 <- contrasts.fit(fit, cont.dif)
    fit2 <- eBayes(fit2)
    stats <- topTableF(fit2, number = Inf)
    colnames(stats) <- c("1HR", "24HR", "3HR", "8HR", "AveExpr", "F", "P.Value", "adj.P.Val")
    stats$Treatment <- rep(trt, nrow(stats))
    stats$Cell.Line <- rep(cls, nrow(stats))
    storage.df <- merge(storage.df, stats %>% rownames_to_column("Protein"), all = TRUE)
  }
}


allStatsData <- storage.df
allStatsData <- allStatsData %>% melt(id = c("Protein", "AveExpr", "F", "P.Value", "adj.P.Val", "Treatment", "Cell.Line"), variable.name = "Time.Point", value.name = "logFC")
allStatsData$Time.Point<- factor(allStatsData$Time.Point, levels = c("1HR", "3HR", "8HR", "24HR"))
# allStatsData$Treatment <- factor(allStatsData$Treatment, levels = c("BMS5996226", "BMS754807", "Dasatinib", "BMS754807.BMS599626", "BMS754807.Dasatinib"))

# 
# prot.toplot <-allStatsData %>% 
#   filter(adj.P.Val <= 0.05)
# prots <- unique(prot.toplot$Protein)

plot.made <- allStatsData %>%
  filter(adj.P.Val <= 0.05, Cell.Line %in% c("Cal27", "FaDu", "OSC.19", "SCC25", "SCC9")) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(Cell.Line, reorder(Protein, -logFC))) +
  geom_tile(aes(fill = logFC), color = "white") + 
  scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red3", mid = "white") + 
  theme_bw() + 
  scale_y_discrete(name = "") +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  facet_grid(Treatment ~ Time.Point, scales = "free", space = "free") + 
  xlab("Time (hrs)") +
  geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.01), vjust = 0.75, hjust = 0.5)

# Saves Image
ggsave(plot = plot.made, filename = "RPPA2 Heatmap Reorg All.png", width = 11, height = 27, units = "in", dpi = 300, device = "png")

### END LARGE REORG HEAT MAP


### BEGIN LARGE HEAT MAP
celllines <- levels(as.factor(make.names(targets$Cell.Line)))
treatments <- make.names(levels(as.factor(gsub(" ", "", targets[which(targets$Treatment != "Control"),]$Treatment, fixed = TRUE))))
timepoints <- levels(as.factor(gsub(" ", "", targets[which(targets$Time.Point != "15 min"),]$Time.Point, fixed = TRUE)))

# storage.df <- data.frame(T1HR = numeric(), T24HR = numeric(), T3HR = numeric(), T8HR = numeric(), AveExpr = numeric(), F = numeric(), P.Value = numeric(), adj.P.Val = numeric(), Treatment = character(), Cell.Line = character(), stringsAsFactors = FALSE)
storage.df <- data.frame()
newcontrast <- vector()
for(trt in treatments)
{
  for(cls in celllines)
  {
    newcontrast <- vector()
    for(tp in timepoints)
    {
      nexttime <- paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
      newcontrast <- c(newcontrast, nexttime)
    }
    cont.dif <- makeContrasts(contrasts = newcontrast, levels = design)
    fit2 <- contrasts.fit(fit, cont.dif)
    fit2 <- eBayes(fit2)
    stats <- topTableF(fit2, number = Inf)
    colnames(stats) <- c("1HR", "24HR", "3HR", "8HR", "AveExpr", "F", "P.Value", "adj.P.Val")
    stats$Treatment <- rep(trt, nrow(stats))
    stats$Cell.Line <- rep(cls, nrow(stats))
    storage.df <- merge(storage.df, stats %>% rownames_to_column("Protein"), all = TRUE)
  }
}


allStatsData <- storage.df
allStatsData <- allStatsData %>% melt(id = c("Protein", "AveExpr", "F", "P.Value", "adj.P.Val", "Treatment", "Cell.Line"), variable.name = "Time.Point", value.name = "logFC")
allStatsData$Time.Point<- factor(allStatsData$Time.Point, levels = c("1HR", "3HR", "8HR", "24HR"))
# allStatsData$Treatment <- factor(allStatsData$Treatment, levels = c("BMS5996226", "BMS754807", "Dasatinib", "BMS754807.BMS599626", "BMS754807.Dasatinib"))

# 
# prot.toplot <-allStatsData %>% 
#   filter(adj.P.Val <= 0.05)
# prots <- unique(prot.toplot$Protein)

plot.made <- allStatsData %>%
  filter(adj.P.Val <= 0.05, Cell.Line %in% c("Cal27", "FaDu", "OSC.19", "SCC25", "SCC9")) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(Time.Point, fct_rev(Protein))) +
  geom_tile(aes(fill = logFC), color = "white") + 
  scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red3", mid = "white") + 
  theme_bw() + 
  scale_y_discrete(name = "") +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  facet_grid(Cell.Line ~ Treatment, scales = "free", space = "free", switch = "y") + 
  xlab("Time (hrs)") +
  geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.01), vjust = 0.75, hjust = 0.5)

# Saves Image
ggsave(plot = plot.made, filename = "RPPA2 Heatmap All.png", width = 11, height = 27, units = "in", dpi = 300, device = "png")

### END LARGE HEAT MAP

### BEGIN Cluster HEAT MAP
celllines <- levels(as.factor(make.names(targets$Cell.Line)))
treatments <- make.names(levels(as.factor(gsub(" ", "", targets[which(targets$Treatment != "Control"),]$Treatment, fixed = TRUE))))
timepoints <- levels(as.factor(gsub(" ", "", targets[which(targets$Time.Point != "15 min"),]$Time.Point, fixed = TRUE)))

# storage.df <- data.frame(T1HR = numeric(), T24HR = numeric(), T3HR = numeric(), T8HR = numeric(), AveExpr = numeric(), F = numeric(), P.Value = numeric(), adj.P.Val = numeric(), Treatment = character(), Cell.Line = character(), stringsAsFactors = FALSE)
storage.df <- data.frame()
newcontrast <- vector()
for(trt in treatments)
{
  for(cls in celllines)
  {
    newcontrast <- vector()
    for(tp in timepoints)
    {
      nexttime <- paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
      newcontrast <- c(newcontrast, nexttime)
    }
    cont.dif <- makeContrasts(contrasts = newcontrast, levels = design)
    fit2 <- contrasts.fit(fit, cont.dif)
    fit2 <- eBayes(fit2)
    stats <- topTableF(fit2, number = Inf)
    colnames(stats) <- c("1HR", "24HR", "3HR", "8HR", "AveExpr", "F", "P.Value", "adj.P.Val")
    stats$Treatment <- rep(trt, nrow(stats))
    stats$Cell.Line <- rep(cls, nrow(stats))
    storage.df <- merge(storage.df, stats %>% rownames_to_column("Protein"), all = TRUE)
  }
}


allStatsData <- storage.df

distfunc <- function(x) daisy(x,metric="gower")
df_expression <- na.omit(allStatsData)
dat <- df_expression[,2:5]  # numerical columns
row.order <- hclust(distfunc(dat))$order # clustering
col.order <- hclust(dist(t(dat)))$order
dat2 <- cbind(df_expression[,1], dat, df_expression[,6:11], row.order)
dat_new <- dat2[row.order, col.order]
colnames(dat2)[1] <- "Protein"




allStatsData <- dat2 %>% melt(id = c("Protein", "AveExpr", "F", "P.Value", "adj.P.Val", "Treatment", "Cell.Line", "row.order"), variable.name = "Time.Point", value.name = "logFC")
allStatsData$Time.Point<- factor(allStatsData$Time.Point, levels = c("1HR", "3HR", "8HR", "24HR"))
# allStatsData$Treatment <- factor(allStatsData$Treatment, levels = c("BMS5996226", "BMS754807", "Dasatinib", "BMS754807.BMS599626", "BMS754807.Dasatinib"))

# 
# prot.toplot <-allStatsData %>% 
#   filter(adj.P.Val <= 0.05)
# prots <- unique(prot.toplot$Protein)

plot.made <- allStatsData %>%
  filter(adj.P.Val <= 0.05, Cell.Line %in% c("Cal27", "FaDu", "OSC.19", "SCC25", "SCC9"), Treatment %in% c("BMS754807", "Dasatinib", "BMS754807.Dasatinib")) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(Time.Point, reorder(Protein, -logFC))) +
  geom_tile(aes(fill = logFC), color = "white") + 
  scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red3", mid = "white") + 
  theme_bw() + 
  scale_y_discrete(name = "") +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  facet_grid(Cell.Line ~ Treatment, scales = "free", space = "free", switch = "y") + 
  xlab("Time (hrs)") +
  geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.01, Cell.Line %in% c("Cal27", "FaDu", "OSC.19", "SCC25", "SCC9"), Treatment %in% c("BMS754807", "Dasatinib", "BMS754807.Dasatinib")), vjust = 0.75, hjust = 0.5)

# Saves Image
ggsave(plot = plot.made, filename = "RPPA2 Heatmap Cluster.png", width = 11, height = 27, units = "in", dpi = 300, device = "png")

### END LARGE Cluster HEAT MAP

# plot.made <- allStatsData %>%
#   filter(Protein %in% prots, Cell.Line %in% c("Cal27", "FaDu", "OSC.19", "SCC25", "SCC9")) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
#   # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
#   ggplot(aes(Time.Point, fct_rev(Protein))) +
#   geom_tile(aes(fill = logFC), color = "white") + 
#   scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red", mid = "white") + 
#   theme_bw() + 
#   scale_y_discrete(name = "") +
#   scale_x_discrete(position = "bottom") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
#   facet_grid(Cell.Line ~ Treatment, scales = "free", space = "free", switch = "y") + 
#   xlab("Time (hrs)") +
#   geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.05), vjust = 0.75, hjust = 0.5)
# 
# # Saves Image
# ggsave(plot = plot.made, filename = "RPPA2 Heatmap All Prots.png", width = 11, limitsize = FALSE, height = 50, units = "in", dpi = 300, device = "png")

### Volcano Plot
volc.trt <- "BMS754807.Dasatinib"
volc.tp <- "8HR"
volc.cls <- c("Cal27", "FaDu", "OSC.19", "SCC25", "SCC9")
times <- c("1HR", "3HR", "8HR", "24HR")
for(volc.trt in treatments)
{
for(volc.tp in times)
{
plot <- allStatsData %>% 
  filter(Cell.Line %in% volc.cls,Treatment %in% volc.trt, Time.Point %in% volc.tp) %>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val))) + geom_point(aes(group = Protein, color = Cell.Line)) + theme_bw() + xlim(-8, 8) + 
  geom_hline(linetype = 2, yintercept = -log10(0.05), alpha = 0.5) + 
  geom_vline(linetype = 1, xintercept = -2, alpha = 0.5, color = "red") + 
  geom_vline(linetype = 1, xintercept = 2, alpha = 0.5, color = "red") + 
  annotate(geom="text", label="FDR 5% Threshold", alpha = 0.8, x=5, y=-log(0.05), vjust=-1) + 
  geom_text_repel(aes(label = Protein), color = "black", size = 3.5, segment.color = "grey", data = allStatsData %>% filter(Cell.Line %in% volc.cls,Treatment %in% volc.trt, Time.Point %in% volc.tp, adj.P.Val <= 0.05, logFC > 2 | logFC < -2)) + 
  ggtitle(paste("Treatment:", volc.trt, " @", volc.tp, sep = "")) + 
  ylim(0,12) +
  theme(legend.position = "bottom")
  
  # geom_text(aes(x = logFC, y = -log10(adj.P.Val), label = Protein), position=position_jitter(width=1,height=1), data = allStatsData %>% filter(Cell.Line == "Cal27", Treatment == "BMS754807.Dasatinib", adj.P.Val <= 0.001))
ggsave(plot = plot, path = "VolcanoPlots", filename = paste("Volcano",volc.trt,volc.tp,"png", sep = "."), width = 11, height = 8.5, units = "in", dpi = 300, device = "png")
}
}

### End Volcano Plot

plot.made2 <- allStatsData %>%
  filter(adj.P.Val <= 0.05, Cell.Line %in% c("Cal27", "FaDu", "OSC.19", "SCC25", "SCC9")) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(x = Time.Point, y = logFC, group = 1)) +
  geom_line(stat = "identity") + 
  theme_bw() + 
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  facet_grid(Cell.Line + Protein ~ Treatment, scales = "free", space = "free") + 
  xlab("Time (hrs)")

# Saves Image
ggsave(plot = plot.made2, filename = "RPPA2 Line All.png", width = 11, height = 49, units = "in", dpi = 300, device = "png")


# Synergism Scripts
# Formula: (Combo1hr - Combo15min) - 0.5 (DrugA1hr - DrugA15 min) - 0.5 (DrugB1hr - DrugB15 min) - 0.5 (Control1hr - Control15 min)
# nexttime <- paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
#
#



celllines <- levels(as.factor(make.names(targets$Cell.Line)))
treatments <- make.names(levels(as.factor(gsub(" ", "", targets[which(targets$Treatment != "Control"),]$Treatment, fixed = TRUE))))
timepoints <- levels(as.factor(gsub(" ", "", targets[which(targets$Time.Point != "15 min"),]$Time.Point, fixed = TRUE)))

newcontrast.bms5bms7 <- vector()
newcontrast.bms7dasat <- vector()
storage.bms5bms7 <- data.frame()
storage.bms7dasat <- data.frame()

for(cls in celllines)
{
  newcontrast.bms5bms7 <- vector()
  newcontrast.bms7dasat <- vector()
  for(tp in timepoints)
  {
    nexttime.bms5bms7 <- paste(cls, "bmsbms", tp, " = (", cls, ".BMS754807.BMS599626.", tp, " - ",cls,".BMS754807.BMS599626.15min",") - 0.5*(",cls,".Control.",tp," - ",cls,".Control.15min) - 0.5*(",cls,".BMS599626.",tp," - ",cls,".BMS599626.15min) - 0.5*(",cls,".BMS754807.",tp," - ",cls,".BMS754807.15min)", sep = "")
    nexttime.bms7dasat <- paste(cls, "bmsdasat", tp, " = (", cls, ".BMS754807.Dasatinib.", tp, " - ",cls,".BMS754807.Dasatinib.15min",") - 0.5*(",cls,".Control.",tp," - ",cls,".Control.15min) - 0.5*(",cls,".Dasatinib.",tp," - ",cls,".Dasatinib.15min) - 0.5*(",cls,".BMS754807.",tp," - ",cls,".BMS754807.15min)", sep = "")
    newcontrast.bms5bms7 <- c(newcontrast.bms5bms7, nexttime.bms5bms7)
    newcontrast.bms7dasat <- c(newcontrast.bms7dasat, nexttime.bms7dasat)
  }
  cont.dif.bms5bms7 <- makeContrasts(contrasts = newcontrast.bms5bms7, levels = design)
  fit2.bms5bms7 <- contrasts.fit(fit, cont.dif.bms5bms7)
  fit2.bms5bms7 <- eBayes(fit2.bms5bms7)
  stats <- topTableF(fit2.bms5bms7, number = Inf)
  colnames(stats) <- c("1HR", "24HR", "3HR", "8HR", "AveExpr", "F", "P.Value", "adj.P.Val")
  stats$Treatment <- rep("BMS754807.BMS599626", nrow(stats))
  stats$Cell.Line <- rep(cls, nrow(stats))
  storage.bms5bms7 <- merge(storage.bms5bms7, stats %>% rownames_to_column("Protein"), all = TRUE)
    
    
  cont.dif.bms7dasat <- makeContrasts(contrasts = newcontrast.bms7dasat, levels = design)
  fit2.bms7dasat <- contrasts.fit(fit, cont.dif.bms7dasat)
  fit2.bms7dasat <- eBayes(fit2.bms7dasat)
  stats <- topTableF(fit2.bms7dasat, number = Inf)
  colnames(stats) <- c("1HR", "24HR", "3HR", "8HR", "AveExpr", "F", "P.Value", "adj.P.Val")
  stats$Treatment <- rep("BMS754807.Dasatinib", nrow(stats))
  stats$Cell.Line <- rep(cls, nrow(stats))
  storage.bms7dasat <- merge(storage.bms7dasat, stats %>% rownames_to_column("Protein"), all = TRUE)
}
allStatsData <- merge(storage.bms5bms7, storage.bms7dasat, all = TRUE)
allStatsData <- allStatsData %>% melt(id = c("Protein", "AveExpr", "F", "P.Value", "adj.P.Val", "Treatment", "Cell.Line"), variable.name = "Time.Point", value.name = "logFC")
allStatsData$Time.Point<- factor(allStatsData$Time.Point, levels = c("1HR", "3HR", "8HR", "24HR"))

plot.made <- allStatsData %>%
  filter(adj.P.Val <= 0.05, Cell.Line %in% c("Cal27", "FaDu", "OSC.19", "SCC25", "SCC9")) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(Cell.Line, reorder(Protein, logFC))) +
  geom_tile(aes(fill = logFC), color = "white") + 
  scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red", mid = "white") + 
  theme_bw() + 
  scale_y_discrete(name = "") +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  facet_grid(Treatment ~ Time.Point, scales = "free", space = "free", switch = "y") + 
  xlab("Time (hrs)") +
  geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.01), vjust = 0.75, hjust = 0.5) +
  ggtitle("Synergism: Combo - 0.5 Drug A - 0.5 Drug B - 0.5 Control")

# Saves Image
ggsave(plot = plot.made, filename = "RPPA2 Heatmap All Synergy2.png", width = 7, height = 14, units = "in", dpi = 300, device = "png")
