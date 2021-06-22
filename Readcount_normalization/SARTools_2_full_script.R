#SARTools DESeq2 pipeline
#source("http://bioconductor.org/biocLite.R")
#BiocInstaller::biocLite(c("DESeq2", "edgeR", "genefilter"))
#install.packages("devtools")

### added:
#install.packages("ghit")
#devtools::install_github("r-lib/remotes")
#remotes::install_github("PF2-pasteur-fr/SARTools")                 # working


library(devtools)
library(DESeq2)
library(SARTools)


################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### April 20th, 2015
### designed to be executed with SARTools 1.1.0
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

loadTargetFile <- function(targetFile, varInt, condRef, batch){
  target <- read.table(targetFile, header=TRUE, sep="\t")
  if (!I(varInt %in% names(target))) stop(paste("The factor of interest", varInt, "is not in the target file"))
  if (!is.null(batch) && !I(batch %in% names(target))) stop(paste("The batch effect", batch, "is not in the target file")) 
  target[,varInt] <- as.factor(target[,varInt])
  if (!I(condRef %in% as.character(target[,varInt]))) stop(paste("The reference level", condRef, "is not a level of the factor of interest"))
  target[,varInt] <- relevel(target[,varInt],ref=condRef)
  target <- target[order(target[,varInt]),]
  rownames(target) <- as.character(target[,1])
  # check if varInt contains replicates
  if (min(table(target[,varInt]))<1) stop(paste("The factor of interest", varInt, "has a level without replicates"))
  # warning message if batch is numeric
  if (!is.null(batch) && is.numeric(target[,batch])) warning(paste("The", batch, "variable is numeric. Use factor() or rename the levels with letters to convert it into a factor"))
  cat("Target file:\n")
  print(target)
  return(target)
}

workDir <- ""                                                     # working directory for the R session
setwd(workDir)

projectName <- "Readcount normalization for all datasets"         # name of the project
author <- "Teresa Coimbra"                                        # author of the statistical analysis/report


targetFile <- "target_file.txt"                                   # path to the design/target file

rawDir <- ""                                                      # path to the directory containing raw counts files
featuresToRemove <- c("alignment_not_unique",                     # names of the features to be removed
                      "ambiguous", "no_feature",                  # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")

varInt <- "detail"                                                # factor of interest
condRef <-  "SRR5CD2"                                             # reference biological condition
batch <- NULL                                                     # blocking factor: NULL (default) or "batch" for example

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.01                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors


#colors <- c("cornflowerblue","green1", "palegreen3", "forestgreen", "brown3", "coral1", "red", "blue","#00B8E5", "blueviolet", "lightcoral","red3", "palegreen3", "forestgreen", "orangered", "orangered4")              # vector of colors of each biological condition on the plots
colors <- c("#F8766D", "#E68613", "#CD9600", "#ABA300", "#7CAE00", "#0CB702", "#00BE67", "#00C19A", "#00BFC4", "#00B8E7", "#00A9FF", "#8494FF", "#C77CFF","#AE87FF", "#ED68ED", "#FF61CC", "#FF68A1")              # vector of colors of each biological condition on the plots


################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file

target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)
# 
# # description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)
# 
# # analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                          locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                          cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
# 
# # PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)
# 
# # summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)
# 
# # save image of the R session
save.image(file=paste0(projectName, ".RData"))
# 
# # generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)

################################################################################
###                             post processing                              ###
################################################################################



