

# convert slim output to ms format

import sys

firstLocus = True
polys = False
seqs = []
covered = []
genomicPositions = {}
end = 0
locusNum = 1
locusDict = {}
locusDict[1] = []
with open(sys.argv[1]) as infile:
    for line in infile:
        if line.strip()[:25] == "initializeGenomicElement(":
            start = int(line.strip().split()[1].split(",")[0])
            if firstLocus == True:
                first = int(start) # the genomic start position of the first locus
                firstLocus = False
                for i in range(start): # fill in blanks in "covered"
                    covered.append(0)
            else:
                if (start-end) > 1:
                    for i in range(end+1, start): # fill in blanks in "covered"
                        covered.append(0)
                    locusNum += 1 
                    locusDict[locusNum] = []
            end = int(line.strip().split()[2].split(")")[0]) # genomic end position of the final locus  
            for i in range(start, end+1):
                covered.append(locusNum)
            genomicPositions[locusNum] = [start, end] # start and end genomic positions for this locus
        if polys == False:
            if line.strip().split() != []:
                if line.strip().split()[0] == "positions:": 
                    positions = line.strip().split()[1:]
                    polys = True
        else:
            seqs.append( line.strip() )
print covered
 


# length of the total simulated dataset. To loop through all sites: for i in range(first, end+1): 
length = int(end) - int(first) + 1

# add each snp to the corresponding locus
for snp in range(len(positions)):
    bp = int(round( ( float(positions[snp]) * end )))
    locusNumber = int(covered[bp])
    locusDict[locusNumber].append([snp, bp])

# now, break up the polymorphism data by locus
for locus in locusDict:
    print # blank
    print "//" + "\t" + str(locus)
    snps = locusDict[locus]
    print "segsites:", len(snps)
    if snps != []:
        outline = ["positions:"]
        for snp in locusDict[locus]: # need to convert each genomic position to locus position
            locusLength = genomicPositions[locus][1]+1 - genomicPositions[locus][0] 
            locusBp = snp[1] - genomicPositions[locus][0]
            locusPosition = float(locusBp) / float(locusLength)
            outline.append(locusPosition)
        print " ".join(map(str,outline))
    for seq in seqs:
        outline = ""
        for snp in locusDict[locus]:
            outline += seq[snp[0]]
        print outline


