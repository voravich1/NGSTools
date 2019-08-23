class utils:

    def createCountDictFromBedgraph(self,bedgraphFileName):

        count_bedgraph = np.loadtxt(fname=bedgraphFileName, dtype="str")

        count_dict = {}
        for count_data in count_bedgraph:

            if (count_data[0] != "" and count_data[0] != "track"):
                start = int(count_data[1])
                stop = int(count_data[2])

                for i in range(start, stop + 1):
                    count_dict[i] = count_data[3]

        return count_dict