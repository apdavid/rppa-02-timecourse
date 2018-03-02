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
rppa.melt <- melt(allData, id.vars = c("GMU.ID", "Cell.Line", "Date", "Time.Point", "Treatment", "Y"))
head(rppa.melt)
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
treat_targets <- targets

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
treat_rppa <- rppa.BC
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

#### LMDME -- 2nd attempt
library(lmdme)
design <- treat_targets
design$Time.Point <- as.factor(design$Time.Point)
design$Time.Point = factor(design$Time.Point,levels(design$Time.Point)[c(2,1,4,5,3)])


design$Treatment <- as.factor(design$Treatment)
design$Cell.Line <- as.factor(design$Cell.Line)

fit <- lmdme(model = ~ Time.Point * Treatment, data = treat_rppa, design = design, Bayes = TRUE, verbose = TRUE)
id.fit <- F.p.values(fit, term = "Time.Point:Treatment") < 0.05
pls.hits <- rownames(treat_rppa[id.fit,])
sum(id.fit) # 26

decomposition(fit, decomposition = "plsr", term = "Time.Point:Treatment", scale = "row", type = "coefficient", subset = id.fit)
# biplot(fit, term = "Time.Point:Treatment", xlabs = ".", expand = 0.9, subset = id.fit, scale = "row", which = "x")
# plsr.coefficients <- components(fit)
# 
# fit@p.values %>% unlist() %>% hist()
# fit@F.p.values %>% unlist() %>% hist()
################ SAVES THE P-VALUES
# test <- as.data.frame(paste(treat_rppa, F.p.values(fit, term = "Time.Point:Treatment")))
# test <- cbind(rownames(treat_rppa), F.p.values(fit, term = "Time.Point:Treatment"))[id.fit,]
# write.csv(test, file = "PLS Hits P-values.csv")


# >>Time Course Plots
# rppa.melt$Time.Point <- as.factor(rppa.melt$Time.Point)
# rppa.melt$Time.Point = factor(rppa.melt$Time.Point,levels(rppa.melt$Time.Point)[c(2,1,4,5,3)])
rppa.melt$Treatment <- as.factor(rppa.melt$Treatment)
rppa.melt$Treatment = factor(rppa.melt$Treatment,levels(rppa.melt$Treatment)[c(5,6,1,2,3,4)])
rppa.melt$Time.Point <- revalue(rppa.melt$Time.Point, c("15 min" = "0.25", "1 hr" = "1", "3 hr" = "3", "8 hr" = "8", "24 hr" = "24"))
rppa.melt$Time.Point <- as.numeric(rppa.melt$Time.Point)

## Stratify by Protein
protein <- "S6RibosomalProtein.S235_236"
for(protein in pls.hits) {
plot <- rppa.melt %>%
  filter(variable == protein) %>%
  ggplot(aes(x = Time.Point, y = log2(as.numeric(value)))) + geom_point() + facet_grid(Cell.Line ~ Treatment) + stat_summary(fun.y = mean, color = "blue", size = 1, geom = "line") + ggtitle(paste("Protein", protein, sep = ": ")) + theme_bw() + ylab("log2(Intensity)") + xlab("Time (hrs)")

  png(filename = paste(protein, "png", sep = "."), width = 11, height = 8.5, units = "in", res = 300)
  print(plot)
  dev.off()
}

## Stratify by Treatment


condition <- "Control"
for(condition in levels(rppa.melt$Treatment)) {
  plot <- rppa.melt %>%
    filter(Treatment == condition, variable %in% pls.hits) %>%
    ggplot(aes(x = Time.Point, y = log2(as.numeric(value)))) + geom_point() + facet_grid(variable ~ Cell.Line) + stat_summary(fun.y = mean, color = "blue", size = 1, geom = "line") + ggtitle(paste("Treatment", condition, sep = ": ")) + theme_bw() + ylab("log2(Intensity)") + xlab("Time (hrs)")
  
  png(filename = paste(condition, "png", sep = "."), width = 10, height = 30, units = "in", res = 300)
  print(plot)
  dev.off()
}

for(condition in levels(rppa.melt$Treatment)) {
  print(condition)
}

id2.fit <- F.p.values(fit, term = "Time.Point:Treatment") < 0.1
pls2.hits <- rownames(treat_rppa[id2.fit,])
condition <- "Control"
for(condition in levels(rppa.melt$Treatment)) {
plot <- rppa.melt %>%
  filter(Treatment == condition, variable %in% pls2.hits) %>%
  ggplot(aes(x = Time.Point, y = log2(as.numeric(value)))) + geom_point() + facet_grid(~ Cell.Line) + stat_summary(fun.y = mean, color = "blue", size = 1, geom = "line") + ggtitle(paste("Treatment", condition, sep = ": ")) + theme_bw() + ylab("log2(Intensity)") + xlab("Time (hrs)")
png(filename = paste("All Proteins", condition, "png", sep = "."), width = 11, height = 8.5, units = "in", res = 300)
print(plot)
dev.off()
}

# all proteins and treatments , FDR < 0.1
id2.fit <- F.p.values(fit, term = "Time.Point:Treatment") < 0.1
pls2.hits <- rownames(treat_rppa[id2.fit,])
plot <- rppa.melt %>%
  filter(variable %in% pls2.hits, Time.Point != 24) %>%
  ggplot(aes(x = Time.Point, y = log2(as.numeric(value)))) + geom_point() + facet_grid(Treatment ~ Cell.Line) + stat_summary(fun.y = mean, color = "blue", size = 1, geom = "line") + theme_bw() + ylab("log2(Intensity)") + xlab("Time (hrs)")
png(filename = paste("All Proteins.NarrowTime", "AllTreatments", "png", sep = "."), width = 11, height = 17, units = "in", res = 300)
print(plot)
dev.off()

png(filename = paste(condition, "png", sep = "."), width = 10, height = 30, units = "in", res = 300)
print(plot)
dev.off()

# LMDME -- First attempt
library(lmdme)
design <- targets
design$Time.Point <- as.factor(design$Time.Point)
design$Time.Point = factor(design$Time.Point,levels(design$Time.Point)[c(2,1,4,5,3)])
design$Treatment <- as.factor(design$Treatment)
design$Cell.Line <- as.factor(design$Cell.Line)

fit <- lmdme(model = ~ Cell.Line + Time.Point + Treatment + Time.Point:Treatment, data = treat_rppa, design = design)
decomposition(fit, decomposition="plsr", type="coefficient", term="Time.Point:Treatment", subset=id, scale="row")
id<-F.p.values(fit, term="Time.Point:Treatment")<0.05
rownames(treat_rppa[id,])
write.csv(rownames(treat_rppa[id,]), file = "file.csv")

decomposition(fit, decomposition="plsr", type="coefficient", term="Time.Point:Treatment", subset=id, scale="row")
biplot(fit, which="loadings", xlabs="o", ylabs=colnames(coefficients(fit, term="Time.Point:Treatment")), var.axes=TRUE)

decomposition(fit, decomposition="pca", type="coefficient", term="Time.Point:Treatment", subset=id, scale="row")

png("myplot.png", width=12, height=8, units="in", res=300)
loadingplot(fit, term.x = "Time.Point", term.y = "Treatment")
dev.off()
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