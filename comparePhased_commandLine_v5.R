

# this version is for phased vs unphased; for looking at variable rec or BGS, use different version
# Instead of plotting the RAE (relative abs error) for each POD,
# I'm going to try outputting the MRAE (median of the errors, summarizing all PODs) for each* parameter, so that I can look at effects on individual parameters
# So I want to output the effect on each parameter, also outputting the name of the parameter.
# Although a slight pain in the ass, I can type out and save some commands for looking all declines across all models, or migration starts, or whatever




args = commandArgs(trailingOnly=TRUE)
print(args)
print("\nWARNING: If comparing sim data to PODs simuated with a different (even if similar) command, first compare headers from each file\n")
### params
paramColumns = args[1]
threads = as.numeric(args[2])
simDataPath = args[3]
obsDataPath = args[4]
outputName = args[5]
statIndicesWeDontWant = c(13:18, 25:30, 47:48, 57:58, 59:88) # excluding ZnS because leads to missing data which is a pain to deal with 	





### read in obs
myObs = read.table(obsDataPath, header = T)

### read in simulation data
data = read.table(simDataPath, header = T)

### filter some stats that we don't want	
firstStatIndex = match( "s_average_segs_1", names(data) )
filteredData = data[, -(statIndicesWeDontWant + (firstStatIndex-1))]
filteredData = na.omit(filteredData)
firstObsStatIndex = match( "s_average_segs_1", names(myObs) ) # this should allow the comparison of data without -Ne to pods simulated with -Ne as the last param
filteredObs = myObs[, -(statIndicesWeDontWant + (firstObsStatIndex-1))]
filteredObs = na.omit(filteredObs)
numPhasedFilteredObs = dim(filteredObs)[1]

### separate priors and sum stats
priors = filteredData[,c(1:firstStatIndex-1)]
stats = filteredData[,c(firstStatIndex:dim(filteredData)[2])]
params = filteredObs[,c(1:firstObsStatIndex-1)]
obsStats = filteredObs[,c(firstObsStatIndex:dim(filteredObs)[2])]

### get param columns, if not specified on command line
if (paramColumns == "all")
{
	newParams = c()
	library(matrixStats)
	for (p in 1:dim(priors)[2])
	{
		if (sd(priors[,p]) != 0)
		{
			newParams = c(newParams,p)
		}
	}
	paramColumns = newParams
} else {
     paramColumns = as.numeric(strsplit(paramColumns, split=",")[[1]])
}
print(paramColumns)



### random forests
library(abcrf)
modindex = rep("theMod", dim(stats)[1]) # placeholder only
myList = list(factor(modindex), data.frame(priors), data.frame(stats))
names(myList) = c("modindex", "priors", "stats")
sumsta <- myList$stats[modindex == "theMod",]
obsList = obsStats
predictionErrorsP = c()
MRAEs = c()
print("random forests, phased data")
for (p in paramColumns)
{
	print(p)
	r <- myList$priors[, p]
	df <- data.frame(r, sumsta)
	model.rf.r <- regAbcrf(r~., df, ntree=100, ncores = threads, paral=T) # can take a VERY long time, days if default settings
	meds1 = predict(model.rf.r, obsList, df, ncores = threads, paral = T)$med # also takes a while (doesn't use "ntrees")
	RAE = abs((params[,p] - meds1) / params[,p])
	predictionErrorsP = cbind( predictionErrorsP, c(names(myObs)[p], RAE) )
	MRAEs = rbind(MRAEs, c("phased", names(myObs)[p], median(RAE)))
}



### repeated the above procedure with un-phased data
phasedStats = c(25:42, 53:56) 
removeStats = c(phasedStats, statIndicesWeDontWant)
filteredData = data[, -(removeStats + (firstStatIndex-1))]
filteredData = na.omit(filteredData)
filteredObs = myObs[, -(removeStats + (firstObsStatIndex-1))]
filteredObs = na.omit(filteredObs)
numUnhasedFilteredObs = dim(filteredObs)[1]

if (numUnhasedFilteredObs == numPhasedFilteredObs)
{
### separate priors and sum stats
priors = filteredData[,c(1:firstStatIndex-1)]
stats = filteredData[,c(firstStatIndex:dim(filteredData)[2])]
params = filteredObs[,c(1:firstObsStatIndex-1)]
obsStats = filteredObs[,c(firstObsStatIndex:dim(filteredObs)[2])]

### random forests
myList = list(factor(modindex), data.frame(priors), data.frame(stats))
names(myList) = c("modindex", "priors", "stats")
sumsta <- myList$stats[modindex == "theMod",]
obsList = obsStats
predictionErrorsU = c()
print("random forets, unphased data")
for (p in paramColumns)
{
	print(p)
	r <- myList$priors[, p]
	df <- data.frame(r, sumsta)
	model.rf.r <- regAbcrf(r~., df, ntree=100, ncores = threads, paral=T) # can take a VERY long time, days if default settings
	meds2 = predict(model.rf.r, obsList, df, ncores = threads, paral = T)$med # also takes a while (doesn't use "ntrees")
	RAE = abs((params[,p] - meds2) / params[,p])
        predictionErrorsU = cbind( predictionErrorsU, c(names(myObs)[p], RAE) )
	MRAEs = rbind(MRAEs, c("unphased", names(myObs)[p], median(RAE)))
}


### output
outputFile = paste(outputName, "_mrad", ".txt", sep = "")
write.table(MRAEs, outputFile, sep = "\t", row.names = F, col.names = F, quote = F)
outputFile = paste(outputName, "_phasedError", ".txt", sep = "")
write.table(predictionErrorsP, outputFile, sep = "\t", row.names = F, col.names = F, quote = F)
outputFile = paste(outputName, "_unphasedError", ".txt", sep = "")
write.table(predictionErrorsU, outputFile, sep = "\t", row.names = F, col.names = F, quote = F)



} else {
print("phased-stat PODs retained post filtering contain NAs")
}




