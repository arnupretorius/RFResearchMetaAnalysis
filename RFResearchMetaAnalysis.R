##################################################################
# A Meta-Analysis of Research in Random Forests for Classification
##################################################################

# Check for missing packages and install if missing
list.of.packages <- c("MASS","dplyr","latex2exp", "mlbench", "ggplot2", "caret", "doSNOW", "lattice",
                      "obliqueRF", "stargazer", "rotationForest", "randomForest",
                      "scmamp", "surv2sampleComp", "ElemStatLearn", "hmeasure")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load required packages
load <- lapply(list.of.packages, require, character.only = TRUE)

# required packages from bioconductor for scmamp package
source("https://bioconductor.org/biocLite.R")
biocLite("graph")
n
biocLite("Rgraphviz")
n

# ---------------------------------------------------------------------------
############################################################################
# Performance estimation method used in the papers considered in 
# the meta-analysis.
############################################################################
# ---------------------------------------------------------------------------
# load data
data <- read.csv("RFComparisonsData.csv")
evalMeth <- data %>% select(paper_title, evaluation)
evalMeth <- unique(evalMeth)
evalMeth <- as.data.frame(evalMeth %>% count(evaluation))
evalMeth <- evalMeth[order(evalMeth$n, decreasing = TRUE),]

# mark which estimation methods are not "reliable"
notIndex <- c(4, 9, 22, 27)
grp <- rep("Reliable", 27)
grp[notIndex] <- "Not reliable"
evalMeth$grp <- grp

# plot estimation methods used
ggplot(evalMeth, aes(x=evaluation, y=n, fill=grp)) + geom_bar(stat="identity") +
    scale_x_discrete(limits=evalMeth$evaluation) +
    scale_fill_manual(name="" ,values=c("darkgreen", "skyblue")) +
    theme_bw() + ylab("#Papers") + xlab("Estimation method") +
    theme(legend.position = c(0.9,0.7), axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
# ---------------------------------------------------------------------------
#########################################################################
# Reported error rates for Breiman's Forest-RI on the top ten 
# most popular data sets used in papers.
#########################################################################
# ---------------------------------------------------------------------------
# load data
data <- read.csv("RFComparisonsData.csv")

# compute top used data sets across papers
ds <- data %>% select(paper_title, dataset)
allDS <- unique(ds$dataset)
allDSCount <- rep(0, length(allDS))
dsSplit <- split(ds, ds$paper_title)
for(i in 1:length(dsSplit)){
    dsPerPaper <-  unique(dsSplit[[i]][,2])
    for(j in 1:length(dsPerPaper)){
        index <- which(allDS == dsPerPaper[j])
        allDSCount[index] <- allDSCount[index] + 1
    }
}
dsFrame <- data.frame(dataset=allDS, freqUsed=allDSCount)
dsFrame <- dsFrame[order(dsFrame$freqUsed, decreasing = TRUE),]

# get dataset characteristics
dataChar <- data %>% select(dataset, dataset_size, num_inputs, classes) %>% distinct(.keep_all=TRUE)
dataChar <- merge(dsFrame, dataChar)
dataChar <- dataChar[order(dataChar$freqUsed, decreasing = TRUE),]

# plot variatability of RF on above data sets
rfDataSets <- factor(unique(dataChar$dataset)[1:10])
keepIndex <- NULL
count <- 1
for(i in 1:nrow(data)){
    if(data[i,]$method == "rf" && data[i,]$dataset %in% rfDataSets){
        keepIndex[count] <- i
        count <- count + 1
    }
}
rfCompData <- data[keepIndex,] %>% select(paper_title, dataset, method, error)
rfCompData <- rfCompData[order(rfCompData$dataset),]
# remove other lymphoma
lympRemove <- NULL
count <- 1
for(i in 1:nrow(rfCompData)){
    if(rfCompData[i,]$dataset == "lymphoma" && rfCompData[i,]$error > 3){
        lympRemove[count] <- i
        count <- count + 1
    }
}
rfCompData <- rfCompData[-lympRemove,]
rfCompData <- rfCompData[-101, ] # remove gross outlier in: On extreme pruning (possibly reported error instead of acc)
# plot errors
ggplot(rfCompData, aes(y=error, x=dataset)) + geom_boxplot(fill="skyblue", outlier.colour = "red")+
    theme_bw() + ylab("Reported error rates") + xlab("Benchmark data set")+
    scale_x_discrete(limits=unique(dataChar$dataset)[1:10])+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

# outlier breast cancer: On extreme pruning (Possibly not the same breast cancer dataset)
# outlier glass: Tripoli et al. paper
# outlier sonar: On extreme pruning (possibly reported error instead of acc)
# ---------------------------------------------------------------------------
############################################################
# Omnibus p-val for Forest-RI over different papers
############################################################
# ---------------------------------------------------------------------------
# test between random forests from different papers
lop <- arrange(data, paper_title) %>% select(paper_title, dataset, method, error)
lop <- filter(lop, lop$method == "rf")
lop$acc <- round(100-lop$error,3)

# top 10 used data sets for Forest-RI
topDataset <- factor(dsFrame[1:15,1])
papers <- factor(unique(lop$paper_title))

# compute all combinations
allcomb <- combn(1:15, 10, simplify = FALSE)
candidateList <- list()
canL <- NULL
# find combinations of 15 choose 10 such that largest number of rfs can be compared
for (k in 1:length(allcomb)) {
    datasetTop10 <- topDataset[allcomb[[k]]]
    lop10 <- filter(lop, lop$dataset %in% datasetTop10)
    splitlop <- split(lop10, factor(lop10$dataset))
    candidatePapers <- NULL
    check <- 0
    count <- 1
    for(i in 1:length(papers)){
        for(j in 1:length(splitlop)){
            papersPresent <- factor(splitlop[[j]][,1])
            if(!(papers[i] %in% papersPresent)){
                check <- 1
            }
        }
        if(check == 0){
            candidatePapers[count] <- as.character(papers[i])
            count <- count + 1
        } else {
            check <- 0
        }
    }
    candidateList[[k]] <- candidatePapers
    canL[[k]] <- length(candidatePapers)
}

# find combination sharing the maximum number of datasets
maxIndex <- which(canL == max(canL))[1]
# list of papers
candidateList[[maxIndex]]

# create compare matrix
dsets <- factor(topDataset[allcomb[[maxIndex]]])
filterIndex <- sapply(lop$dataset, function(x){ifelse(x %in% dsets, TRUE, FALSE)})
lop <- lop[filterIndex,]
lopSplit <- split(lop, factor(lop$dataset))
compareMat <- matrix(0, nrow=length(dsets), ncol=length(papers))
rownames(compareMat) <- names(lopSplit)
colnames(compareMat) <- papers
for(i in 1:length(lopSplit)){
    for(j in 1:nrow(lopSplit[[i]])){
        compareMat[i, which(papers == as.character(lopSplit[[i]]$paper_title[j]))] <- lopSplit[[i]]$acc[j]
    }
}
# prune to include papers containing all top ten data sets
keepIndex <- apply(compareMat, 2, function(x){
    ifelse(length(which(x == 0)) == 0, TRUE, FALSE)
})
rCompareMat <- compareMat[,keepIndex]

# compute omnibus tests (#algorithm < 5)
imanDavenportTest(rCompareMat)
friedmanAlignedRanksTest(rCompareMat)
quadeTest(rCompareMat)
# ---------------------------------------------------------------------------
########################################################################
# Methods used to compare different algorithms over multiple 
# data sets in the papers considered for the meta-analysis.
########################################################################
# ---------------------------------------------------------------------------
# plot evaluation method used
evalsData <- data %>% select(paper_title, comparison) %>% distinct(.keep_all=TRUE)
lims <- unique(evalsData$comparison)[c(1,2,4,3,5,6,7,8)]
ggplot(evalsData, aes(x=comparison)) + geom_bar(fill="darkgreen") +
    xlab("Comparison method") + ylab("#Papers")+
    theme_bw()+
    scale_x_discrete(limits=lims)+
    scale_y_continuous(breaks = seq(0, 20, by = 2))+
    theme(axis.text.x = element_text(angle = 10, vjust = 1, hjust = 1))
# ---------------------------------------------------------------------------
########################################################################
# Methods used to compare different algorithms over multiple 
# data sets in the papers considered for the meta-analysis.
########################################################################
# ---------------------------------------------------------------------------
# Redo analyses using omnibus and post-hoc tests
dataSplit <- split(data, factor(data$paper_title))
pvals <- NULL
checkPhTest <- NULL
phTests <- list()

# for each paper build a compare matrix and compute the Ivan-Davenport test p-value
for(k in 1:length(dataSplit)){
    lop <- arrange(dataSplit[[k]], dataset) %>% select(dataset, method, error)
    lop <- as.data.frame(summarise(group_by(lop, dataset, method), mean(error)))
    lop$acc <- round(100-lop$`mean(error)`,3)
    
    # create compare matrix
    lopSplit <- split(lop, factor(lop$dataset))
    dsets <- unique(lop$dataset)
    methods <- unique(lop$method)
    compareMat <- matrix(0, nrow=length(dsets), ncol=length(methods))
    rownames(compareMat) <- dsets
    colnames(compareMat) <- methods
    for(i in 1:length(lopSplit)){
        for(j in 1:nrow(lopSplit[[i]])){
            compareMat[i, which(methods == as.character(lopSplit[[i]]$method[j]))] <- lopSplit[[i]]$acc[j]
        }
    }
    pvals[k] <- round(imanDavenportTest(compareMat)[[3]], 3)
    if(!is.na(pvals[k]) && pvals[k] < 0.05 && "rf" %in% colnames(compareMat)){
        phTests[[k]] <- postHocTest(compareMat, test = "friedman",
                                    control = "rf", correct = "finner", alpha = 0.05)
        checkPhTest[k] <- ifelse(length(which(phTests[[k]]$corrected.pval < 0.05)) == 0, TRUE, FALSE)
    }
}

# paper where no significant result was found
omnibusFailed <- which(pvals > 0.05)
PHvsRFFailed <- which(checkPhTest)
nonSigResult <- c(omnibusFailed, PHvsRFFailed)

# plot pvals
pvals[PHvsRFFailed] <- 0.05
grp <- rep("Omnibus: Iman-Davenport", length(pvals))
grp[PHvsRFFailed] <- "Post-hoc: Finner (Forest-RI as control)"
grp <- factor(grp)
pvalData <- data.frame(pv = pvals[-29], Test=grp[-29])
ggplot(pvalData, aes(x=1:nrow(pvalData), y=pv, fill=Test)) + geom_bar(stat="identity") +
    theme_bw() + xlab("'Papers' (no names given)") +
    ylab("p-value") + geom_hline(yintercept = 0.05, col="red", linetype="dashed") +
    scale_fill_manual(values=c("darkgreen", "skyblue"))+
    scale_x_continuous(breaks = seq(from=1,to=34, by=1))+
    annotate("text", x=0.5, y=0.09, label = "alpha==0.05 ",parse = TRUE)+
    theme(legend.position=c(0.2,0.8))
# ---------------------------------------------------------------------------