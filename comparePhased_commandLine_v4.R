


# this version will not be comparing phased and unphased stats at all, but rather comparing two treatments with jsut phased (or just unphased, if you change it) stats
# so it takes an extra set of sims


args = commandArgs(trailingOnly=TRUE)
print(args)
print("\nWARNING: If comparing sim data to PODs simuated with a different (even if similar) command, first compare headers from each file\n")
### params
paramColumns = args[1]
threads = as.numeric(args[2])
simDataPath = args[3]
neutralDataPath = args[4]
obsDataPath = args[5]
outputName = args[6]




#statIndicesWeDontWant = c(13:18, 25:30, 47:48, 57:58, 59:88)
statIndicesWeDontWant = c(13:18, 25:42, 47:48, 53:58, 59:88) # phasedStats 25:42, 53:56





### read in obs
myObs = read.table(obsDataPath, header = T)

### read in simulation data
data = read.table(simDataPath, header = T)
dataN = read.table(neutralDataPath, header = T)
### filter some stats that we don't want	
firstStatIndex = match( "s_average_segs_1", names(data) )
filteredData = data[, -(statIndicesWeDontWant + (firstStatIndex-1))]
filteredDataN = dataN[, -(statIndicesWeDontWant + (firstStatIndex-1))]
filteredData = na.omit(filteredData)
filteredDataN = na.omit(filteredDataN)
firstObsStatIndex = match( "s_average_segs_1", names(myObs) ) # this should allow the comparison of data without -Ne to pods simulated with -Ne as the last param
filteredObs = myObs[, -(statIndicesWeDontWant + (firstObsStatIndex-1))]
filteredObs = na.omit(filteredObs)
numPhasedFilteredObs = dim(filteredObs)[1]



### reduce both datasets to equal numbers of rows / simulations
if ( dim(filteredData)[1] > dim(filteredDataN)[1] )
{
	filteredData = filteredData[1:dim(filteredDataN)[1], ] 	
}
if ( dim(filteredData)[1] < dim(filteredDataN)[1] )
{
	filteredDataN = filteredDataN[1:dim(filteredData)[1], ]
}



### separate priors and sum stats
priors = filteredData[,c(1:firstStatIndex-1)]
priorsN = filteredDataN[,c(1:firstStatIndex-1)]
stats = filteredData[,c(firstStatIndex:dim(filteredData)[2])]
statsN = filteredDataN[,c(firstStatIndex:dim(filteredDataN)[2])]
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
modindexN = rep("theMod", dim(statsN)[1]) # placeholder only
myList = list(factor(modindex), data.frame(priors), data.frame(stats))
myListN = list(factor(modindexN), data.frame(priorsN), data.frame(statsN))
names(myList) = c("modindex", "priors", "stats")
names(myListN) = c("modindex", "priors", "stats")
sumsta <- myList$stats[modindex == "theMod",]
sumstaN <- myListN$stats[modindex == "theMod",]
obsList = obsStats
MRAEs = c()
predictionErrors1 = c()
predictionErrors2 = c()
print("random forests")
for (p in paramColumns)
{
	print(p)
	r <- myList$priors[, p]
	df <- data.frame(r, sumsta)
	model.rf.r <- regAbcrf(r~., df, ntree=100, ncores = threads, paral=T) # can take a VERY long time, days if default settings
	meds1 = predict(model.rf.r, obsList, df, ncores = threads, paral = T)$med # also takes a while (doesn't use "ntrees")
	RAE = abs((params[,p] - meds1) / params[,p]) # relative absolute error
	predictionErrors1 = cbind(predictionErrors1, c(names(myObs)[p], RAE))
	MRAEs = rbind(MRAEs, c("control", names(myObs)[p], median(RAE)))

	# again with neutral sims
	r <- myListN$priors[, p]
	df <- data.frame(r, sumstaN)
	model.rf.r <- regAbcrf(r~., df, ntree=100, ncores = threads, paral=T) # can take a VERY long time, days if default settings
	meds1 = predict(model.rf.r, obsList, df, ncores = threads, paral = T)$med # also takes a while (doesn't use "ntrees")
	RAE = abs((params[,p] - meds1) / params[,p]) # relative absolute error
	predictionErrors2 = cbind(predictionErrors2, c(names(myObs)[p], RAE))
	MRAEs = rbind(MRAEs, c("ignore", names(myObs)[p], median(RAE)))
}
#phasedRAE = rowSums(predictionErrors1) 
#phasedRAE_N = rowSums(predictionErrors2)









######### mean +/- SD plot
# outputFile = paste(outputName, "_meanSD", ".jpeg", sep = "")
# jpeg(outputFile)
# u_mean=mean(phasedRAE); p_mean=mean(phasedRAE_N); u_sd=sd(phasedRAE); p_sd=sd(phasedRAE_N)
# plot(0, ylim = c(0, max(c(u_sd+u_mean, p_sd+p_mean))*1.25 ), xlim = c(0.45,.55), xaxt='n', xlab = "", ylab = "prediction error", col = "white")
# segments(0.5-0.05, u_mean, 0.5-0.005, u_mean, lwd = 3, col = "red") # mean lines
# segments(0.5+0.005, p_mean, 0.5+0.05, p_mean, lwd = 3, col = "blue")
# segments(0.5-0.05, u_mean+u_sd, 0.5-0.005, u_mean+u_sd, lwd = 1, col = "red") # upper sd lines
# segments(0.5+0.005, p_mean+p_sd, 0.5+0.05, p_mean+p_sd, lwd = 1, col = "blue")
# segments(0.5-0.05, u_mean-u_sd, 0.5-0.005, u_mean-u_sd, lwd = 1, col = "red") # lower sd lines
# segments(0.5+0.005, p_mean-p_sd, 0.5+0.05, p_mean-p_sd, lwd = 1, col = "blue")
# segments(0.4725, u_mean-u_sd, 0.4725, u_mean+u_sd, lwd = 1, lty = "dashed", col = "red") # vertical dashed lines
# segments(0.5275, p_mean-p_sd, 0.5275, p_mean+p_sd, lwd = 1, lty = "dashed", col = "blue")
# dev.off()



# ######### box and whiskers
# outputFile = paste(outputName, "_boxNwhis", ".jpeg", sep = "")
# jpeg(outputFile)
# boxplot( log(df$predictionError) ~ df$phase )
# dev.off()



### output
outputFile = paste(outputName, "_mrad", ".txt", sep = "")
write.table(MRAEs, outputFile, sep = "\t", row.names = F, col.names = F, quote = F)
outputFile = paste(outputName, "_control", ".txt", sep = "")
write.table(predictionErrors1, outputFile, sep = "\t", row.names = F, col.names = F, quote = F)
outputFile = paste(outputName, "_ignore", ".txt", sep = "")
write.table(predictionErrors2, outputFile, sep = "\t", row.names = F, col.names = F, quote = F)

	
	




