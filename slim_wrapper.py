


# Usage: python slim_wrapper.py out/path/out_name num_sims path/to/locfile slim arguments (arguments meaning... we'll come back to this)

# NOTE(s):
#          1. all priors must be integers or else will be rounded, except for migration rates
#          2. the string "migRate" must be part of migration rate parameters, or else will be rounded to integer (e.g. 0)
#          3. assumes a single locus length is used 
#          4. currently assuming 2 populations

import sys, numpy, os, math
slimVersion = "/projects/chriscs/Software/SLiM/slim"
out_name = sys.argv[1] # important for writing temp files that don't get overwritten by other by other jobs
num_sims = int(sys.argv[2])
loc_path = sys.argv[3]
# skip the next argv, which = "slim"  
command = [slimVersion] + sys.argv[5:]
priors = {"-U" : 2,
          "-R" : 2} # flags for prior ranges, and the number of parameters they take: currently uniform and log uniform
out_path = "/".join(out_name.split("/")[0:-1]) 
if out_path != []:
    out_path += "/"
out_name = out_name.split("/")[-1]




### function for random numbers
def uniform(low, high):
    r = numpy.random.uniform(0, 1, 1)
    return float( ((high - low)*r) + low ) # floating because it comes out as an "array" at first


def loguniform(low, high):
    r = numpy.random.uniform(0, 1, 1)
    ta = math.log(low)
    tb = math.log(high)
    tr = ((tb - ta)*r + ta)
    return math.exp(tr) 





### organize prior information
priorRanges = []
priorCounter = 0
remainingPriors = True
while remainingPriors == True:
    remainingPriors = False
    for field in range(len(command)):
        if "=" in command[field]:
            if command[field].split("=")[1] in priors:
                infoCounter = priors[command[field].split("=")[1]]
                priorRanges.append( [command[field].split("=")[1]] + map(float,command[field+1:field+infoCounter+1]) )
                command[field] = command[field].split("=")[0] +"=XPRIORX" # placeholder for priors
                command = command[0:field+1] + command[field+infoCounter+1:]
                remainingPriors = True
                break
priorHeader = []
for p in range(len(priorRanges)):
    priorHeader.append("p"+str(p+1))
priorHeader = "\t".join(priorHeader)
priorHeader += "\n"




### grab number of loci, locus length, and num samples from loc file
with open(loc_path) as infile:
    head = infile.readline()
    locusCounter = 0
    for line in infile:
        locusCounter += 1
num_samples = line.strip().split()[1]
num_loci = locusCounter / 2
locus_length = line.strip().split()[2]
firstHalf, confFile = command[:-1], command[-1]
command = firstHalf + ["-d numSamples="+str(num_samples), "-d numLoci="+str(num_loci), "-d locLength="+str(locus_length)] + [confFile]





### msABC command
msABCcommand = "msABC NUMSEQS 1 -I NUMPOPS --frag-begin --finp LOCFILE --frag-end --obs " 
msABCcommand = msABCcommand.replace("NUMSEQS", str(int(num_samples)*2)) 
iind = msABCcommand.index("-I")
dashIinfo = [str(2)]
for pop in range(2):   
    dashIinfo.append( str( num_samples ))
msABCcommand = msABCcommand.replace("NUMPOPS", " ".join(dashIinfo))
msABCcommand = msABCcommand.replace("LOCFILE", loc_path) 



### do simulations
firstSim = True
for sim in range(num_sims):

    # go through priors and draw random parameters 
    params = []
    priorCounter = 0 
    priorDrawsForOutput = []
    for field in range(len(command)):
        if "=" in command[field]:
            parameterPrior = command[field].split("=")
            if parameterPrior[1] == "XPRIORX":
                prior = priorRanges[priorCounter]
                if prior[0] == "-U": # uniform
                    param = uniform(prior[1], prior[2])
                elif prior[0] == "-R": # log uniform
                    param = loguniform(prior[1], prior[2])
                if "migRate" not in parameterPrior[0]:
                    param = int(round(param))
                params.append([field, param])
                priorDrawsForOutput.append(param)
                priorCounter += 1

    # simulate data
    indivCommand = list(command)
    for p in range(len(params)):
        indivCommand[params[p][0]] = indivCommand[params[p][0]].split("=")[0] + "=" + str(params[p][1])
    indivCommand = " ".join(map(str,indivCommand)) + " > " + out_path + "temp_polData_" + out_name
    print indivCommand
    os.system( indivCommand )

    # convert output to ms format
    convertCommand = "python /projects/chriscs/Scripts/slim2ms.py " + out_path + "temp_polData_" + out_name + " > " + out_path + "temp_polData_" + out_name + ".ms"
    os.system(convertCommand)

    # calculate sum stats on the set of simulated loci
    os.system( str(msABCcommand) + out_path + "temp_polData_" + out_name + ".ms > " + out_path + "temp_stats_" + out_name  )

    # organize output
    if firstSim == True: # if first sim, output header
        with open(out_path + "temp_priors_" + out_name, "w") as priorFile:
            priorFile.write(priorHeader)
        os.system("head -3 " + out_path + "temp_stats_" + out_name + " | tail -1 > " + out_path + "temp_statHeader_" + out_name )
        os.system("paste " + out_path + "temp_priors_" + out_name + " " + out_path + "temp_statHeader_" + out_name + " >" + out_path + "finalOut_" + out_name)
        firstSim = False
    with open(out_path + "temp_priors_" + out_name, "w") as summaryFile:
        summaryFile.write("\t".join(map(str,priorDrawsForOutput)))
    os.system("tail -1 " + out_path + "temp_stats_" + out_name + " | tail -1 > " + out_path + "temp_statsTail_" + out_name )
    os.system("paste " + out_path + "temp_priors_" + out_name + " " + out_path + "temp_statsTail_" + out_name + " >>" + out_path + "finalOut_" + out_name)





### rm temp files
os.system("rm " + out_path + "temp_polData_" + out_name)
os.system("rm " + out_path + "temp_polData_" + out_name + ".ms")
os.system("rm " + out_path + "temp_priors_" + out_name)
os.system("rm " + out_path + "temp_statHeader_" + out_name)
os.system("rm " + out_path + "temp_statsTail_" + out_name)
os.system("rm " + out_path + "temp_stats_" + out_name)







