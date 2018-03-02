library(limma)
library(lumi)
library(tidyverse)
library(gplots)
library(reshape2)
library(readxl)
source("FunctionsFromBrian.R")


# # loads workbook into R as an "workbook"
# wb <- readxl::read_excel("/Users/Cl5pn/Box Sync/Jameson_Lab/RPPA/R analysis/DG Array1 Summary File-FINAL 070914.xlsx",sheet=2)
# 
# # converts excel file into "data.frame"
# allData <- wb
# rm(wb)

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
# limma relies heavily on factors for experimental design, so levels must be explictly set
targets$`Time.Point` <- gsub("\ ", "", targets$`Time.Point`) # spaces can be troublesome
targets$`Time.Point`<- factor(targets$`Time.Point`,
                              levels = c("15min", "1hr", "3hr", "8hr", "24hr"))
# if you must change the order, do as seen here
targets$Treatment <- gsub("\ \\+\ ", "_", targets$Treatment) # special chars are not allowed
targets$Treatment <- factor(targets$Treatment,
                            levels = c("Control", "BMS754807", "BMS599626",  "Dasatinib", "BMS754807_BMS599626", "BMS754807_Dasatinib"))

targets$`Cell.Line` <- gsub("\\-", "", targets$`Cell.Line`)


rppa.raw <- as.matrix(allData[,c(6:151)])
# limma expects an expression matrix to have genes as rows and samples as columns, so we take the transverse
rppa.trnsv <- t(rppa.raw[,-1])
class(rppa.trnsv) <- "numeric"

# this will induce NAs ("ND" is not able to be cast as a numeric)

# "ND"  means unsure of value 
# "0" means the analyte is truly not there

# limma expects log transformed data
rppa.BC <- log2(rppa.trnsv)

# turn matrix into object that limma needs
lumi.Q <- new("ExpressionSet", exprs = rppa.BC)
design <- (targets %>%
  transmute(design = paste(`Cell.Line`, `Time.Point`, Treatment, sep = "_")))[[1]] %>%
  as.factor()
  
design <- model.matrix(~0+design)
fit <- lmFit(lumi.Q, design)

# fill in the blank
cont.matrix <- makeOneWayAnovaContrastMatrix(design)
contMatrixRowIDs <- colsplit(row.names(cont.matrix), "\\_", c("Cell.Line", "Time.Point", "Treatment"))
contMatrixRowIDs <- contMatrixRowIDs %>% 
  mutate(`Cell.Line` = gsub("design", "", `Cell.Line`))

# Cell.Lines need to be the same
# Time.Points need to be the same
# one treatment must be control

validContrastIndex <- vector(length = ncol(cont.matrix))
humanReadableColNames <- vector(length = ncol(cont.matrix))
for ( i in 1:ncol(cont.matrix)){
  contrastColumn <- cont.matrix[,i]
  firstCellLine <- contMatrixRowIDs$`Cell.Line`[cont.matrix[,i] == 1]
  secondCellLine <- contMatrixRowIDs$`Cell.Line`[cont.matrix[,i] == -1]
  firstTimePoint <- contMatrixRowIDs$`Time.Point`[cont.matrix[,i] == 1]
  secondTimePoint <- contMatrixRowIDs$`Time.Point`[cont.matrix[,i] == -1]
  firstTreatment <- contMatrixRowIDs$Treatment[cont.matrix[,i] == 1]
  secondTreatment <- contMatrixRowIDs$Treatment[cont.matrix[,i] == -1]
  
  if ((any(c(firstTreatment, secondTreatment) == "Control") == T &
      all(c(firstTreatment, secondTreatment) == "Control") == F) &
      firstCellLine == secondCellLine  &
      firstTimePoint == secondTimePoint){
      validContrastIndex[i] <- T
      humanReadableColNames[i] <- paste(firstCellLine, firstTimePoint, firstTreatment, secondTreatment, sep = ".")
  }else{validContrastIndex[i] <- F}
}

cont.matrix <- cont.matrix[,validContrastIndex]
colnames(cont.matrix) <- humanReadableColNames[validContrastIndex]
for (cell in unique(contMatrixRowIDs$`Cell.Line`)){
  for(time in unique(contMatrixRowIDs$`Time.Point`)){
    for(combo in unique(contMatrixRowIDs$Treatment[grep("\\_", contMatrixRowIDs$Treatment)])){
      nextcontvector <- rep(0, nrow(cont.matrix))
      nextcontvector[contMatrixRowIDs$`Cell.Line`==cell &
                       contMatrixRowIDs$`Time.Point`==time &
                       contMatrixRowIDs$Treatment== combo] <- 1
      nextcontvector[contMatrixRowIDs$`Cell.Line`==cell &
                       contMatrixRowIDs$`Time.Point`==time &
                       contMatrixRowIDs$Treatment %in% c("BMS599626", "Dasatinib")] <- -0.5
    cont.matrix <- cbind(cont.matrix, nextcontvector)
    colnames(cont.matrix)[ncol(cont.matrix)] <- paste(cell, time, combo, "synergy", sep=".")
      }
  }
}

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allStatsData <- data.frame()
for ( i in 1:length(colnames(fit2))){
  top <- topTable(fit2, coef = i, number = Inf)
  top$Source <- rep(colnames(fit2)[i], nrow(top))
  top$Protein <- row.names(top)
  row.names(top) <- NULL
  allStatsData <- rbind(allStatsData, top)
}

allStatsDataTargets <- colsplit(allStatsData$Source, "\\.", c("Cell.Line", "Time.Point", "Treatment 1", "Treatment 2"))

allStatsData$logFC[allStatsDataTargets$`Treatment 1` == "Control"] <- allStatsData$logFC[allStatsDataTargets$`Treatment 1` == "Control"]*(-1) 

allStatsDataTargets$`Treatment 1`[allStatsDataTargets$`Treatment 1` == "Control"] <- allStatsDataTargets$`Treatment 2`[allStatsDataTargets$`Treatment 1` == "Control"] 
allStatsDataTargets$`Treatment 2`[allStatsDataTargets$`Treatment 1` == allStatsDataTargets$`Treatment 2`] <- "Control" 

allStatsData <- bind_cols(allStatsData, allStatsDataTargets) %>%
  select(-Source)

save(allStatsData, file="RPPA_contrasts.Rdata")

allStatsData$`Treatment 1`[allStatsData$`Treatment 2` == "synergy"] <- paste0(allStatsData$`Treatment 1`[allStatsData$`Treatment 2` == "synergy"], "SYN")
allStatsData <- allStatsData %>%
  mutate(`Time.Point` = factor(`Time.Point`, levels = c("15min", "1hr", "3hr", "8hr", "24hr")),
         `Cell.Line` = factor(`Cell.Line`, levels = c("Cal27", "OSC19", "SCC25", "FaDu", "SCC9")),
         `Treatment 1` = factor(`Treatment 1`, levels = c("BMS754807", "BMS599626", "Dasatinib", "BMS754807_BMS599626", "BMS754807_Dasatinib", "BMS754807_BMS599626SYN", "BMS754807_DasatinibSYN")))

###if we want to plot only select epitopes
#FilterList <- unique(allStatsData$Protein)[grep("Jak2 Y1007|Src Y527|Src Family Y416|RON Y1353|Pyk2 Y402|PRAS40 T246|Paxillin Y118|FAK Y576_577|ErbB3???HER3 Y1289|AcetylCoA Carboxylase S79", unique(allStatsData$Protein))]

ggplot(allStatsData %>% 
         filter(adj.P.Val < 0.05,
                `Treatment 1` %in% c("BMS754807_DasatinibSYN")),
                #Protein %in% FilterList),
                #`Time.Point` == "24hr"),
                #`Cell.Line` == "FaDu"),
       aes(`Cell.Line`, Protein, fill = logFC)) +
  geom_tile() +
  facet_grid(`Treatment 1` ~ `Time.Point`) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5.5)
        ) + 
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0) 


####diff heat map instead of geom_tile setup####


DataToBeHeatMapped<- allStatsData %>%
  select(logFC, Protein:`Treatment 2`) %>%
  spread(Protein,
         logFC)

FDRMatrix <- allStatsData %>%
  select(adj.P.Val, Protein:`Treatment 2`) %>%
  spread(Protein,
         adj.P.Val)
# apply same filter and select call on matrix 

FDRMatrix <- FDRMatrix %>%
  as.matrix()

FDRCutoff <- 0.05
FDRMatrix[FDRMatrix <= FDRCutoff] <- "*"
FDRMatrix[FDRMatrix > FDRCutoff] <- ""

#DataToBeHeatMapped <- DataToBeHeatMapped[,c(12,13,11)]

source("FunctionsFromBrian.R")

heatmap.2(as.matrix(DataToBeHeatMapped %>% 
                      #filter() %>% # put in filter criteria here to grab specific comparisons of interest
                      select(-(`Cell.Line`:`Treatment 2`))) %>%
            t(), # must be a numeric matrix of values only, row names should be proteins, column names should be samples
          Rowv = T, # cluster rows (proteins)?
          Colv = F, # cluster columns (samples)?
          col = bluered(50),
          #labRow = row.names(allStatsData), # labels for rows (proteins)? By default, set to row.names(x)
          #labCol = colnames(allStatsData), #c("Cal27", "FaDu", "OSC19", "SCC25", "scc9"),#c("599 v ctrl", "754 v ctrl", "combo"),
          # labels for columns (samples)? By default, set to colnames(x)
          key=T, # plot key?
          srtCol = 25, cexCol = 1.25,
          #main= "Dasatinib Combo 1hr",
          #density.info = "none", # controls if key has histogram overlay
          scale="none", # "none" is used for fold changes, "row" is used for expression values as it'll scale the colors by Z score for that protein.
          symbreaks=T, # ensures 0 fold change or median expression value is always white
          trace="none", # I have no idea why this option exists, keep it on "none"
          dendrogram = "row", # draw dendrogram for "row", "col" or "both" (rows and columns),
          #cellnote = FDRMatrix,
          margins=c(4.5,8), lhei=c(1.5,8.5), lwid=c(2,5), keysize=1, key.par = list(cex=0.5),
          cexRow = 0.5)
          #notecol = "black")

DataToBeHeatMapped %>% 
  select(`Cell.Line`:`Treatment 2`) %>%
  ggplot()

#ggplot(allStatsData, aes(`Cell.Line`, Protein, fill = -log10(adj.P.Val))) +
# geom_tile() +
# facet_grid(`Time.Point` ~ `Treatment 1`) + 
# theme_bw() + 
# theme(axis.text.x = element_text(angle = 90, hjust = 1),
#      axis.text.y = element_blank()) +
#  scale_fill_gradient2(low = "blue", high = "red", midpoint = -log10(0.10))

