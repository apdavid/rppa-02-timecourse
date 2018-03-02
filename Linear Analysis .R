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
library(magicfor)   
source("FunctionsFromBrian.R")
# crossref <- read.csv("paths.csv")

# read in RPPA data saved as csv from excel
# CHANGE THIS TO RUN LOCALLY
# allData <- read.csv("~/BTSync/MASS_SPEC/CrazyRPPA/RPPA_EGFR_IGF1R.csv")

allData <- read.csv("dasat_bms_rppa.csv", stringsAsFactors = FALSE)

# grab only "targets" data. Limma uses the language of "targets" and "expression matrix". "targets"
# refers to annotation data that describes the experimental design. The "expression matrix" refers
# to the expression data for which the linear model will be fit. 

targets <- allData[,c(1:6)]
head(allData[,c(7:151)])


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
cont.dif <- makeContrasts(SCC9.BMS5.Diff.1hr = (SCC9.BMS599626.1hr - SCC9.BMS599626.15min) - (SCC9.Control.1hr - SCC9.Control.15min), 
                          SCC9.BMS5.Diff.3hr = (SCC9.BMS599626.3hr - SCC9.BMS599626.15min) - (SCC9.Control.3hr - SCC9.Control.15min), 
                          SCC9.BMS5.Diff.8hr = (SCC9.BMS599626.8hr - SCC9.BMS599626.15min) - (SCC9.Control.8hr - SCC9.Control.15min), 
                          SCC9.BMS5.Diff.24hr = (SCC9.BMS599626.24hr - SCC9.BMS599626.15min) - (SCC9.Control.24hr - SCC9.Control.15min), 
                          SCC9.BMS7.Diff.1hr = (SCC9.BMS754807.1hr - SCC9.BMS754807.15min) - (SCC9.Control.1hr - SCC9.Control.15min), 
                          SCC9.BMS7.Diff.3hr = (SCC9.BMS754807.3hr - SCC9.BMS754807.15min) - (SCC9.Control.3hr - SCC9.Control.15min), 
                          SCC9.BMS7.Diff.8hr = (SCC9.BMS754807.8hr - SCC9.BMS754807.15min) - (SCC9.Control.8hr - SCC9.Control.15min), 
                          SCC9.BMS7.Diff.24hr = (SCC9.BMS754807.24hr - SCC9.BMS754807.15min) - (SCC9.Control.24hr - SCC9.Control.15min), 
                          SCC9.BMS57.Diff.1hr = (SCC9.BMS754807.BMS599626.1hr - SCC9.BMS754807.BMS599626.15min) - (SCC9.Control.1hr - SCC9.Control.15min), 
                          SCC9.BMS57.Diff.3hr = (SCC9.BMS754807.BMS599626.3hr - SCC9.BMS754807.BMS599626.15min) - (SCC9.Control.3hr - SCC9.Control.15min), 
                          SCC9.BMS57.Diff.8hr = (SCC9.BMS754807.BMS599626.8hr - SCC9.BMS754807.BMS599626.15min) - (SCC9.Control.8hr - SCC9.Control.15min), 
                          SCC9.BMS57.Diff.24hr = (SCC9.BMS754807.BMS599626.24hr - SCC9.BMS754807.BMS599626.15min) - (SCC9.Control.24hr - SCC9.Control.15min), 
                          levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
write.csv(topTableF(fit2, adjust="BH", number = 100, p.value = 0.05), file = "SCC9.BMS599626.Diff.csv")


cont.dif <- makeContrasts(Cal27.BMS5.Diff.1hr = (Cal27.BMS599626.1hr - Cal27.BMS599626.15min) - (Cal27.Control.1hr - Cal27.Control.15min), 
                          Cal27.BMS5.Diff.3hr = (Cal27.BMS599626.3hr - Cal27.BMS599626.15min) - (Cal27.Control.3hr - Cal27.Control.15min), 
                          Cal27.BMS5.Diff.8hr = (Cal27.BMS599626.8hr - Cal27.BMS599626.15min) - (Cal27.Control.8hr - Cal27.Control.15min), 
                          Cal27.BMS5.Diff.24hr = (Cal27.BMS599626.24hr - Cal27.BMS599626.15min) - (Cal27.Control.24hr - Cal27.Control.15min), 
                          Cal27.BMS7.Diff.1hr = (Cal27.BMS754807.1hr - Cal27.BMS754807.15min) - (Cal27.Control.1hr - Cal27.Control.15min), 
                          Cal27.BMS7.Diff.3hr = (Cal27.BMS754807.3hr - Cal27.BMS754807.15min) - (Cal27.Control.3hr - Cal27.Control.15min), 
                          Cal27.BMS7.Diff.8hr = (Cal27.BMS754807.8hr - Cal27.BMS754807.15min) - (Cal27.Control.8hr - Cal27.Control.15min), 
                          Cal27.BMS7.Diff.24hr = (Cal27.BMS754807.24hr - Cal27.BMS754807.15min) - (Cal27.Control.24hr - Cal27.Control.15min), 
                          Cal27.BMS57.Diff.1hr = (Cal27.BMS754807.BMS599626.1hr - Cal27.BMS754807.BMS599626.15min) - (Cal27.Control.1hr - Cal27.Control.15min), 
                          Cal27.BMS57.Diff.3hr = (Cal27.BMS754807.BMS599626.3hr - Cal27.BMS754807.BMS599626.15min) - (Cal27.Control.3hr - Cal27.Control.15min), 
                          Cal27.BMS57.Diff.8hr = (Cal27.BMS754807.BMS599626.8hr - Cal27.BMS754807.BMS599626.15min) - (Cal27.Control.8hr - Cal27.Control.15min), 
                          Cal27.BMS57.Diff.24hr = (Cal27.BMS754807.BMS599626.24hr - Cal27.BMS754807.BMS599626.15min) - (Cal27.Control.24hr - Cal27.Control.15min), 
                          levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
write.csv(topTableF(fit2, adjust="BH", number = 100, p.value = 0.05), file = "Cal27.BMS599626.Diff.csv")


cont.dif <- makeContrasts(SCC9.BMS57.Diff.1hr = (SCC9.BMS754807.BMS599626.1hr - SCC9.BMS754807.BMS599626.15min) - (SCC9.Control.1hr - SCC9.Control.15min), 
                          SCC9.BMS57.Diff.3hr = (SCC9.BMS754807.BMS599626.3hr - SCC9.BMS754807.BMS599626.15min) - (SCC9.Control.3hr - SCC9.Control.15min), 
                          SCC9.BMS57.Diff.8hr = (SCC9.BMS754807.BMS599626.8hr - SCC9.BMS754807.BMS599626.15min) - (SCC9.Control.8hr - SCC9.Control.15min), 
                          SCC9.BMS57.Diff.24hr = (SCC9.BMS754807.BMS599626.24hr - SCC9.BMS754807.BMS599626.15min) - (SCC9.Control.24hr - SCC9.Control.15min),
                          Cal27.BMS57.Diff.1hr = (Cal27.BMS754807.BMS599626.1hr - Cal27.BMS754807.BMS599626.15min) - (Cal27.Control.1hr - Cal27.Control.15min), 
                          Cal27.BMS57.Diff.3hr = (Cal27.BMS754807.BMS599626.3hr - Cal27.BMS754807.BMS599626.15min) - (Cal27.Control.3hr - Cal27.Control.15min), 
                          Cal27.BMS57.Diff.8hr = (Cal27.BMS754807.BMS599626.8hr - Cal27.BMS754807.BMS599626.15min) - (Cal27.Control.8hr - Cal27.Control.15min), 
                          Cal27.BMS57.Diff.24hr = (Cal27.BMS754807.BMS599626.24hr - Cal27.BMS754807.BMS599626.15min) - (Cal27.Control.24hr - Cal27.Control.15min), 
                          levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
write.csv(topTableF(fit2, adjust="BH", number = 100, p.value = 0.05), file = "Cal27SCC9.BMS599626.Diff.csv")

celllines <- levels(as.factor(make.names(targets$Cell.Line)))
treatments <- make.names(levels(as.factor(gsub(" ", "", targets[which(targets$Treatment != "Control"),]$Treatment, fixed = TRUE))))
timepoints <- levels(as.factor(gsub(" ", "", targets[which(targets$Time.Point != "15 min"),]$Time.Point, fixed = TRUE)))

newcontrast <- vector()
heatmaparray <- data.frame()
for(cls in celllines)
{
  for(trt in treatments)
  {
    newcontrast <- vector()
    for(tp in timepoints)
    {
      # print(paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
      nexttime <- paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
      newcontrast <- c(newcontrast, nexttime)
    }
    cont.dif <- makeContrasts(contrasts = newcontrast, levels = design)
    fit2 <- contrasts.fit(fit, cont.dif)
    fit2 <- eBayes(fit2)
    write.csv(topTableF(fit2, adjust="BH", number = 10), file = paste(cls, trt, ".csv"))
    
  }
}

cls <- "Cal27"
trt <- "BMS599626"
newcontrast <- vector()
for(tp in timepoints)
{
  # print(paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
  nexttime <- paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
  newcontrast <- c(newcontrast, nexttime)
}
cont.dif <- makeContrasts(contrasts = newcontrast, levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
stats <- topTableF(fit2, adjust="BH", number = Inf)
colnames(stats) <- c("1HR", "24HR", "3HR", "8HR", "AveExpr", "F", "P.Value", "adj.P.Val")
## Create Heatmap DF

allStatsData <- bind_rows(stats %>% rownames_to_column("Protein")) %>% melt(id = c("Protein", "AveExpr", "F", "P.Value", "adj.P.Val"), variable.name = "Time.Point", value.name = "logFC")
allStatsData$Time.Point<- factor(allStatsData$Time.Point, levels = c("1HR", "3HR", "8HR", "24HR"))
allStatsData %>%
  filter(adj.P.Val <=0.05) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(Time.Point, fct_rev(Protein))) +
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red", mid = "white") + 
  theme_bw() + 
  scale_y_discrete(name = "") +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.01), vjust = 0.75, hjust = 0.5) + ggtitle(paste("Cell Line: ", cls, " \nTreatment: ", trt))

allStatsData %>%
  filter(adj.P.Val <=0.05) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(x = Time.Point, y = logFC, group = 1)) +
  geom_line(stat = "identity") + 
  theme_bw() + 
  xlab("Time Point") +
  ylab(expression("log"[2] * "(Fold Change)")) +
  facet_wrap( ~ Protein, ncol = 1) + 
  ggtitle(paste("Cell Line: ", cls, "\nTreatment: ", trt))

### COMBINED PLOTS


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

plot.made <- allStatsData %>%
  filter(adj.P.Val <= 0.05, Cell.Line %in% c("Cal27", "FaDu", "OSC.19", "SCC25", "SCC9")) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(Time.Point, fct_rev(Protein))) +
  geom_tile(aes(fill = logFC), color = "gray20") + 
  scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red", mid = "white") + 
  theme_bw() + 
  scale_y_discrete(name = "") +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  facet_grid(Cell.Line ~ Treatment, scales = "free", space = "free") + 
  xlab("Time (hrs)") +
  geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.01), vjust = 0.75, hjust = 0.5)
ggsave(plot = plot.made, filename = "All.png", width = 10, height = 24, dpi = 300, device = "png")


cls <- "Cal27"
trt <- "BMS599626"
newcontrast <- vector()
for(tp in timepoints)
{
  # print(paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
  nexttime <- paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
  newcontrast <- c(newcontrast, nexttime)
}
cont.dif <- makeContrasts(contrasts = newcontrast, levels = design)
fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
stats <- topTableF(fit2, adjust="BH", number = Inf)
colnames(stats) <- c("1HR", "24HR", "3HR", "8HR", "AveExpr", "F", "P.Value", "adj.P.Val")
## Create Heatmap DF

allStatsData <- bind_rows(stats %>% rownames_to_column("Protein")) %>% melt(id = c("Protein", "AveExpr", "F", "P.Value", "adj.P.Val"), variable.name = "Time.Point", value.name = "logFC")
allStatsData$Time.Point<- factor(allStatsData$Time.Point, levels = c("1HR", "3HR", "8HR", "24HR"))
allStatsData %>%
  filter(adj.P.Val <=0.05) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(Time.Point, fct_rev(Protein))) +
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red", mid = "white") + 
  theme_bw() + 
  scale_y_discrete(name = "") +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.01), vjust = 0.75, hjust = 0.5) + ggtitle(paste("Cell Line: ", cls, " \nTreatment: ", trt))

allStatsData %>%
  filter(adj.P.Val <=0.05) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(x = Time.Point, y = logFC, group = 1)) +
  geom_line(stat = "identity") + 
  theme_bw() + 
  xlab("Time Point") +
  ylab(expression("log"[2] * "(Fold Change)")) +
  facet_wrap( ~ Protein, ncol = 1) + 
  ggtitle(paste("Cell Line: ", cls, "\nTreatment: ", trt))





# ggsave saves the last plot made, but you can pass it saved plots (i.e. p <- ggplot(...))
# ggsave(p, fileName, etc...)

with(allStatsData %>% filter(adj.P.Val <=0.05), table(Protein))

# ggsave("2017-Dec-01-OSI-Or-BMS-Vs-IGF-Stim-Paired-By-Cell-Line-Heat-Map APD 5-PCT New.pdf",
#        width = 8.5,
#        height = 11)
ggsave("Version 3 bottom.png",
       width = 5.5,
       height = 11, dpi = 300, device = "png")


# paste(cls,trt,tp," = (",cls,".",trt,".",tp," - ",cls,".",trt,".15min",") - (",cls,".Control.",tp," - ",cls,".Control.15min)", sep = "")
# Which genes respond differently over time in the BMS5 relative to control
# Make a targets data frame with data for the stimulation effect, much as the same as above
# stim_targets <- targets %>%
#   filter(`Cell Line` != "SCC25GR1",
#          Stimulus %in% c("0", "IGF"),
#          `Inhibitor - EGFR` == "0",
#          `Inhibitor - IGF1R` == "0") %>%
#   mutate(Stimulus = factor(Stimulus, levels = c("0", "IGF")),
#          `Cell Line` = factor(`Cell Line`, levels = unique(`Cell Line`))
#   )

# subset appropriate experimental data for each analysis, which returns the index that satisfy the constraints

# 
# stim_rppa <- rppa.BC[,which((targets$`Cell Line` != "SCC25GR1" &
#                                targets$Stimulus %in% c("0", "IGF") &
#                                targets$`Inhibitor - EGFR` == "0" &
#                                targets$`Inhibitor - IGF1R` == "0") == T)]

# Paired analysis for limma (see limma user manual section 9.4 for details)
# stim_design <- model.matrix(~stim_targets$`Cell Line` + stim_targets$Stimulus)
# stim_fit <- lmFit(stim_rppa, stim_design)
# stim_fit2 <- eBayes(stim_fit)
# stim_stats <- topTable(stim_fit2, coef = "stim_targets$StimulusIGF", number = Inf)


# WHERE PLS MATRICES COME IN
# treat_design <- model.matrix(~treat_targets$Cell.Line + treat_targets$Treatment + treat_targets$Time.Point)
# treat_fit <- lmFit(treat_rppa, treat_design)
# treat_fit2 <- eBayes(treat_fit)

# grab BMS results, then OSI stats. Usually we do this with a loop for contrasts, but that
# is overkill for just two contrasts
# BMS_stats <- topTable(treat_fit2, coef = "treat_targets$TreatmentBMS", number = Inf)
# OSI_stats <- topTable(treat_fit2, coef = "treat_targets$TreatmentOSI", number = Inf)

BMS59_stats <- topTable(treat_fit2, coef = "treat_targets$TreatmentBMS599626", number = Inf)
BMS75_stats <- topTable(treat_fit2, coef = "treat_targets$TreatmentBMS754807", number = Inf)
BMS59BMS75_stats <- topTable(treat_fit2, coef = "treat_targets$TreatmentBMS754807 + BMS599626", number = Inf)
BMS75Dasat_stats <- topTable(treat_fit2, coef = "treat_targets$TreatmentBMS754807 + Dasatinib", number = Inf)
Dasat_stats <- topTable(treat_fit2, coef = "treat_targets$TreatmentDasatinib", number = Inf)

allStatsData <- bind_rows(BMS75_stats %>% rownames_to_column("Protein"),
                          BMS59BMS75_stats %>% rownames_to_column("Protein"),
                          BMS75Dasat_stats %>% rownames_to_column("Protein"),
                          Dasat_stats %>% rownames_to_column("Protein"),
                          .id = "Treatment") %>%
  left_join(data.frame("Treatment" = c("1", "2", "3", "4"),
                       "Condition" = c("BMS75",
                                       "BMS59 + BMS75",
                                       "BMS75 + Dasat",
                                       "Dasat")))


# combine the _stats results together for plotting
# allStatsData <- bind_rows(stim_stats %>% rownames_to_column("Protein"),
#                           BMS_stats %>% rownames_to_column("Protein"),
#                           OSI_stats %>% rownames_to_column("Protein"),
#                           .id = "Treatment") %>%
#   left_join(data.frame("Treatment" = c("1", "2", "3"),
#                        "Condition" = c("StimVsUnstim",
#                                        "StimBMS_Vs_Stim",
#                                        "StimOSI_Vs_Stim")))
# 
# allStatsData <- bind_rows(stim_stats %>% rownames_to_column("Protein"),
#                           BMS_stats %>% rownames_to_column("Protein"),
#                           OSI_stats %>% rownames_to_column("Protein"),
#                           .id = "Treatment") %>%
#   left_join(data.frame("Treatment" = c("1", "2", "3"),
#                        "Condition" = c("Uninhibited",
#                                        "BMS",
#                                        "OSI")))
# allStatsData$Condition <- relevel(allStatsData$Condition, "Uninhibited")
# # filter to significant genes, then pass into ggplot for heatmap
# allStatsData %>% filter(adj.P.Val <=0.05) %>% write_csv("FDR 5 Percent.csv")
# # allStatsData <- allStatsData %>% filter(adj.P.Val <=0.05)
# allStatsData <- merge(allStatsData, crossref, by = "Protein", all.x = TRUE)


allStatsData %>%
  filter(adj.P.Val <=0.05) %>% 
  ggplot(aes(Time.Point, logFC)) + geom_point()

allStatsData %>%
  filter(adj.P.Val <=0.05) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(Condition, fct_rev(Protein))) +
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red", mid = "white") + 
  theme_bw() + 
  scale_y_discrete(name = "") +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.01), vjust = 0.75, hjust = 0.5)# ggsave saves the last plot made, but you can pass it saved plots (i.e. p <- ggplot(...))
# ggsave(p, fileName, etc...)

# ggsave("2017-Dec-01-OSI-Or-BMS-Vs-IGF-Stim-Paired-By-Cell-Line-Heat-Map APD 5-PCT New.pdf",
#        width = 8.5,
#        height = 11)
ggsave("Version 3 bottom.png",
       width = 5.5,
       height = 11, dpi = 300, device = "png")




### 10

# ggsave("~/BTSync/MASS_SPEC/CrazyRPPA/2017-Dec-01-OSI-Or-BMS-Vs-IGF-Stim-Paired-By-Cell-Line-Heat-Map APD.pdf",
#        width = 8.5,
#        height = 11)

# save stats to csv file, openable by excel and R
# allStatsData %>% write_csv("~/BTSync/MASS_SPEC/CrazyRPPA/2017-Dec-01-OSI-Or-BMS-Vs-IGF-Stim-Paired-By-Cell-Line-Results-Table.csv")

allStatsData %>% write_csv("2017-Dec-01-OSI-Or-BMS-Vs-IGF-Stim-Paired-By-Cell-Line-Results-Table.csv")

allStatsData %>% filter(adj.P.Val <=0.05) %>% write_csv("FDR 5 Percent.csv")

source("FunctionsFromBrian.R")
myClustResult <- pvclust(allStatsData, method.hclust="average", method.dist="correlation", parallel=TRUE)

plot(myClustResult) # dendogram with p values
#add rectangles around groups highly supported by the data
pvrect(myClustResult, alpha=.1)