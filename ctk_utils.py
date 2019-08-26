import numpy as np
from pathlib import Path

def createCountDictFromBedgraph(bedgraphFileName):

    #count_bedgraph = np.loadtxt(fname=bedgraphFileName, dtype="str")
    #count_bedgraph = open(bedgraphFileName)

    with open(bedgraphFileName,'r') as count_bedgraph:
        count_dict = {}
        for line in count_bedgraph:
            if (line != "\n" and not line.startswith("track")):
                count_data = line.split("\t")
                start = int(count_data[1])
                stop = int(count_data[2])

                for i in range(start, stop):
                    count_dict[i] = int(count_data[3].split("\n")[0])

    return count_dict

def convert2rpmBedgraph(bedgraphFileName,totalread):
    # we adjust rpm_count into scale (10^5). This will give us a more lower number on count but preserve the sense of rpm count
    p = Path(bedgraphFileName)
    inputFolder = Path(p.parent)
    fileName = p.stem
    suffix = p.suffix

    saveFileName = fileName + "_rpm" + suffix
    bedgraphRPM = inputFolder.joinpath(saveFileName)
    bedgraphRPMFile = open(bedgraphRPM, 'w')

    rpm_factor = totalread/1000000

    with open(bedgraphFileName,'r') as count_bedgraph:
        for line in count_bedgraph:
            if (line != "\n" and not line.startswith("track")):
                count_data = line.split("\t")
                chr = count_data[0]
                start = int(count_data[1])
                stop = int(count_data[2])
                count = int(count_data[3].split("\n")[0])
                rpm_count = round((count/rpm_factor)/100000,2) # adjust rpm count to 10^5 scale
                bedgraphline = chr + "\t" + str(start) + "\t" + str(stop) + "\t" + str(rpm_count) + "\n"
                bedgraphRPMFile.write(bedgraphline)
            else:
                bedgraphRPMFile.write(line)

    bedgraphRPMFile.close()

