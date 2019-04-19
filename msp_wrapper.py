


# NOTE(s):
#          1. events must be specified in increasing order, beginning with 0, going further back in time
#          2. currently only works with fixed theta, although it wouldn't be too hard to add in a theta-prior
#          3. requires nsites to be specified, even if rho = 0
#          4. assuming locus length is the same for each population for the same locus (in locfile)
#          5. currently assuming multiple populations

import sys, numpy, os, math
mspmsVersion = "~/.local/bin/mspms"
command = [mspmsVersion] + sys.argv[1:]
priors = {"-U" : 2,
          "-R" : 2} # flags for prior ranges, and the number of parameters they take: currently uniform and log uniform



### function for random numbers
def uniform(low, high):
    r = numpy.random.uniform(0, 1, 1)
    return float((high - low)*r + low)


def loguniform(low, high):
    r = numpy.random.uniform(0, 1, 1)
    ta = math.log(low)
    tb = math.log(high)
    tr = ((tb - ta)*r + ta)
    return math.exp(tr)



### get output info
oInd = command.index("-output")
outInfo = command[oInd + 1].split("/")
outName = outInfo[-1] # the end of the path
outPath = "/".join(outInfo[:-1]) # everything but the end of the path
if outPath != "":
    outPath += "/" # if there is a path or directory, add a slash to the end; otherwise, don't do a slash or it will try to write to root
command = command[:oInd] + command[oInd+2:]



### read in loc file
locInfo = []
locPath = None
currentLocusNum = 0
if "-locfile" in command:
    lInd = command.index("-locfile")
    locPath = command[lInd + 1]
    command = command[:lInd] + command[lInd+2:]
    with open(locPath, "r") as infile:
        head = infile.readline()
        for line in infile:
            newline = line.strip().split()
            locusNum = int(newline[0].split("f")[1])
            if locusNum != currentLocusNum:
                locInfo.append( [map(int,newline[1:3])] )
                currentLocusNum += 1
            else:
                locInfo[locusNum-1].append( map(int,newline[1:3]) )


### organize prior information
priorRanges = []
priorCounter = 0
remainingPriors = True
while remainingPriors == True:
    remainingPriors = False
    for field in range(len(command)):
        if command[field] in priors:
            infoCounter = priors[command[field]]
            priorRanges.append( [command[field]] + map(float,command[field+1:field+infoCounter+1]) )
            command[field] = "XPRIORX" # placeholder for priors
            command = command[0:field+1] + command[field+infoCounter+1:]
            remainingPriors = True
            break



### grab a few more parameters from the command line
numSamples = command[1]
numSims = int(command[2])
command[2] = len(locInfo) # replacing the number of "full simulations" with the number of loci to simulate
tind = command.index("-t")
theta = float(command[tind+1])
rind = command.index("-r")
rho = float(command[rind+1])
nsites = float(command[rind+2])
iind = command.index("-I")
numPops = int(command[iind+1])
priorHeader = []
for p in range(len(priorRanges)):
    priorHeader.append("p"+str(p+1))
priorHeader = "\t".join(priorHeader)
priorHeader += "\n"
msABCcommand = "msABC NUMSEQS 1 -I NUMPOPS --frag-begin --finp LOCFILE --frag-end --obs POLDATA > " + outPath + "temp_stats_" + outName
msABCcommand = msABCcommand.replace("NUMSEQS", numSamples)
dashIinfo = [str(numPops)]
for pop in range(numPops):
    dashIinfo.append( str( command[iind+2+pop] ))
msABCcommand = msABCcommand.replace("NUMPOPS", " ".join(dashIinfo))
msABCcommand = msABCcommand.replace("LOCFILE", locPath)
msABCcommand = msABCcommand.replace("POLDATA", outPath + "temp_polData_" + outName)









### do simulations
firstSim = True
for sim in range(numSims):

    # go through priors and draw random parameters 
    params = []
    priorCounter = 0 
    priorDrawsForOutput = []
    for field in range(len(command)):
        if command[field] == "XPRIORX":
            prior = priorRanges[priorCounter]
            if prior[0] == "-U": # uniform
                param = uniform(prior[1], prior[2])
                params.append([field, param])
            elif prior[0] == "-R": # log uniform
                param = loguniform(prior[1], prior[2])
                params.append([field, param])
            priorDrawsForOutput.append(param)
            priorCounter += 1

    # simulate data
    indivCommand = list(command)
    for p in range(len(params)):
        indivCommand[params[p][0]] = params[p][1]
    indivCommand = " ".join(map(str,indivCommand)) + " | tail -n +3 > " + outPath + "temp_polData_" + outName # the tail -n +3 is shaves off the mspms command 
    os.system( indivCommand )

    # calculate sum stats on the set of simulated loci
    os.system( str(msABCcommand) )

    # organize output
    if firstSim == True: # if first sim, output header
        with open(outPath + "temp_priors_" + outName, "w") as priorFile:
            priorFile.write(priorHeader)
        os.system("head -3 " + outPath + "temp_stats_" + outName + " | tail -1 > " + outPath + "temp_statHeader_" + outName )
        os.system("paste " + outPath + "temp_priors_" + outName + " " + outPath + "temp_statHeader_" + outName + " >" + outPath + "finalOut_" + outName)
        firstSim = False
    with open(outPath + "temp_priors_" + outName, "w") as summaryFile:
        summaryFile.write("\t".join(map(str,priorDrawsForOutput)))
    os.system("tail -1 " + outPath + "temp_stats_" + outName + " | tail -1 > " + outPath + "temp_statsTail_" + outName )
    os.system("paste " + outPath + "temp_priors_" + outName + " " + outPath + "temp_statsTail_" + outName + " >>" + outPath + "finalOut_" + outName)









# rm temp files
os.system("rm " + outPath + "temp_polData_" + outName)
os.system("rm " + outPath + "temp_priors_" + outName)
os.system("rm " + outPath + "temp_statHeader_" + outName)
os.system("rm " + outPath + "temp_statsTail_" + outName)
os.system("rm " + outPath + "temp_stats_" + outName)







