#Author: Jose Lugo-Martinez
#File: runHypergraphletsRegressionSuite.py
#Date: October 23, 2017
#Advisor: Prof. Predrag Radivojac

#Sample call: python runHypergraphletsRegressionSuite.py  

import os, sys, glob, copy, math, random

def runHypergraphletsRegressionSuite():
    dataPATH = './data/'
    outputPATH = './results/'
    testCaseName = 'example'
    node2nodeID = {'R':'0', 'A':'1', 'B':'2', 'C':'3'}
    MAX_NODE_LABELS = 'NNNNN' 
    MAX_HYPEREDGE_LABELS = 'EEEEEEEEEEEEE' 
    HYPERGRAPHLETS_TYPES = 472
    SAMPLES = 50

    FLAG_rijk = '\t0\t1\t2\t3\t4'
    FLAG_ijk = '\t1\t2\t3\t4'
    SUBGRAPH_TEST = False
    
    infilename = 'hypergraphlets_description.txt'
    try:
        #Open input file 
        infile = open(infilename, "r")
    except IOError as e:
        print ("<<ERROR>> Unable to open the file", infilename, "\nThis program will be quiting now.", e)
        sys.exit()

    hypergraph = {}
    rijkHypergraphlets = []
    ijkHypergraphlets = []
    rijk_ijkHypergraphlets = []
    #Iterate over file
    for line in infile:
        #Eliminate any unnecessary white spaces
        line = line.strip()
        hgType = int(line.split(':')[0])
        hyperedges_desc = line.split(':')[1]
        hyperedgesDesc = hyperedges_desc.split('|')
        totalHyperedges = 0
        temp = []
        if hgType == 1:
            totalNodes = 2
        elif hgType > 1 and  hgType < 11:
            totalNodes = 3
        else:
            totalNodes = 4
        for currHyperedge in hyperedgesDesc:
            hyperedge = ''
            for currNode in currHyperedge.split(','):
                nodeID = node2nodeID[currNode]
                hyperedge += '\t' + nodeID
            temp.append(hyperedge)
            totalHyperedges += 1
            if currHyperedge == 'R,A,B,C' or currHyperedge == 'A,B,C':
                totalNodes = 5
                if currHyperedge == 'R,A,B,C':
                    temp.append(FLAG_rijk)
                    totalHyperedges += 1
                    rijkHypergraphlets.append(hgType)
                else:
                    temp.append(FLAG_ijk)
                    totalHyperedges += 1
                    ijkHypergraphlets.append(hgType)
        if (hgType in rijkHypergraphlets) and (hgType in ijkHypergraphlets):
            rijk_ijkHypergraphlets.append(hgType)
        hypergraph[hgType] = copy.copy(temp)
        #Output test case for current hypergraphlet type
        #Output hypergraph file
        outfilenameHG = dataPATH + testCaseName + str(hgType) + '.hypergraph'
        outfileHG = open(outfilenameHG, "w")
        hyperedges = hypergraph[hgType]
        random.shuffle(hyperedges)
        for i in range(0, len(hyperedges)):
            temp = hyperedges[i].strip().split('\t')
            random.shuffle(temp)
            permutedHyperedge = ''
            for nodeID in temp:
                permutedHyperedge += '\t' + nodeID
            outlineHG = str(i) + permutedHyperedge + '\n'
            outfileHG.writelines(outlineHG)
        outfileHG.close()

        #Output node labels file
        outfilenameNL = dataPATH + testCaseName + str(hgType) + '.nlabels'
        outfileNL = open(outfilenameNL, "w")            
        outlineNL = MAX_NODE_LABELS[0:totalNodes] + '\n'
        outfileNL.writelines(outlineNL)
        outfileNL.close()

        #Output hyperedge labels file
        outfilenameEL = dataPATH + testCaseName + str(hgType) + '.elabels'
        outfileEL = open(outfilenameEL, "w")
        outlineEL = MAX_HYPEREDGE_LABELS[0:totalHyperedges] + '\n'
        outfileEL.writelines(outlineEL)
        outfileEL.close()
    #Close input file
    infile.close()

    print (len(rijkHypergraphlets), len(ijkHypergraphlets), len(rijk_ijkHypergraphlets))
    
    casesFailed = 0
    notFound = 0
    for hgType in range(1, HYPERGRAPHLETS_TYPES):
        failed = 0
        for sampleRun in range(0, SAMPLES):
            #Output hypergraph file
            outfilenameHG = dataPATH + testCaseName + str(hgType) + '.hypergraph'
            outfileHG = open(outfilenameHG, "w")
            hyperedges = hypergraph[hgType]
            random.shuffle(hyperedges)
            for i in range(0, len(hyperedges)):
                temp = hyperedges[i].strip().split('\t')
                random.shuffle(temp)
                permutedHyperedge = ''
                for nodeID in temp:
                    permutedHyperedge += '\t' + nodeID
                outlineHG = str(i) + permutedHyperedge + '\n'
                outfileHG.writelines(outlineHG)
            outfileHG.close()

            #Run hyperkernel code
            command = './run_hyperkernel -p examples.pos -n examples.neg -g ' + dataPATH + testCaseName + str(hgType) + ' -l ' + dataPATH + testCaseName + str(hgType) +  ' -e ' + dataPATH + testCaseName + str(hgType) + ' -t 2 -z 0 -s ' + outputPATH + testCaseName + str(hgType) + '_shgk.svml '
            os.system(command)

            #Open SVML results file
            resultFilename = outputPATH + testCaseName + str(hgType) + '_shgk.svml'
            infileSVML = open(resultFilename, "r")
            temp = infileSVML.readline().strip()
            pairs = temp.split(' ')
            typeFound = False
            for i in range(1, (len(pairs) - 1)):
                EXPECTED_COUNT = 1.0
                countType = int(pairs[i].split(':')[0])
                count = float(pairs[i].split(':')[1])
#                print ("\t", sampleRun, hgType, countType, count)
                if countType == hgType:
                    typeFound = True
                    if count != EXPECTED_COUNT and not SUBGRAPH_TEST:
                        failed += 1
                        print ('Type ', hgType, ' ... MISMATCH', pairs)
                    break
            if (typeFound == False):
                print ('Type ', hgType, ' ... NOT FOUND', countType, count)
                failed += 1
                notFound += 1
                print (pairs)
                return
            infileSVML.close()

        if failed > 0:
            print ('Type ', hgType, ' ... FAILED on', (float(failed)/float(SAMPLES) * 100.0), '% of inputs.')
            casesFailed += 1
        else:
            print ('Type ', hgType, ' ... PASSED on', SAMPLES, 'different inputs.')

    print ("FAILED test cases ... ", casesFailed)
    print ("Hypergraphlets NOT FOUND ... ", notFound)

    return

if __name__ == '__main__':
    runHypergraphletsRegressionSuite() #Alternate call, if you are using IDLE.
##    runHypergraphletsRegressionSuite(sys.argv)
