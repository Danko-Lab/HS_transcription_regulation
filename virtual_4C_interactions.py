# Virtual 4C plot - show interactions with a given locus

import os, sys
import timing

import numpy as np
import getopt
from collections import defaultdict

# Initialize globals
TSSDict = {}
interactionsDict = defaultdict(dict)
cumulativeDict = defaultdict(dict)
binRange = {}
windowSize = 100000
binSize = 5000

# Is the read being examined within the window surrounding the TSS 
def InTSS(tssRowNumber, chromo, position, strand):
    global TSSDict
    global binSize

    # If the end of the TSS file has been reached, return False
    #print 'InTSS:', chromo, position, strand
    if tssRowNumber not in TSSDict: return False

    tss = TSSDict[tssRowNumber]
    tssChromo = tss[0]
    tssPosition = tss[1]
    tssStrand = tss[5]

    #print 'InTSS:', tssChromo, chromo, tssPosition, position, tssStrand, strand

    # If chromosomes or strands are different, return False
    if tssChromo != chromo: return False
    #if tssStrand != strand: return False

    # If read not within the window surrounding the TSS return False
    if position < tssPosition - (binSize/2) or position > tssPosition + (binSize/2): return False

    #print 'InTSS: return True'

    # Otherwise, the read is within the window surrounding the TSS, so return True
    return True

# Return the bin number for the read
def GetBin(tssRowNumber, chromo, position, strand):
    global TSSDict
    global windowSize
    global binSize
    global binRange

    binNum = -1

    # If the end of the TSS file has been reached, return an invalid bin
    if tssRowNumber not in TSSDict: return -1

    tss = TSSDict[tssRowNumber]
    tssChromo = tss[0]
    tssPosition = tss[1]
    tssStrand = tss[5]

    # If chromosomes are different, return an invalid bin
    if tssChromo != chromo: return -1

    # Figure out the bin number
    distance = position - tssPosition
    #print 'Distance: position - tssPosition:', position, ' - ', tssPosition, ' = ', distance
    for binTemp in range(windowSize / binSize):
        if distance in binRange[binTemp]: 
            binNum = binTemp
            break

    return binNum

# Return the cumulative bin number for the read
def GetCumulativeBin(tssRowNumber, chromo, position, strand):
    global TSSDict
    global windowSize
    global binSize

    binNum = -1

    # If the end of the TSS file has been reached, return an invalid bin
    if tssRowNumber not in TSSDict: return -1

    tss = TSSDict[tssRowNumber]
    tssChromo = tss[0]
    tssPosition = tss[1]
    tssStrand = tss[5]

    # If chromosomes are different, return an invalid bin
    if tssChromo != chromo: return -1

    # Figure out the bin number
    binNum = position - (position % binSize)

    return binNum

def main(argv):
    global TSSDict
    global interactionsDict
    global cumulativeDict
    global windowSize
    global binSize
    global binRange

    inputTSSFile = 'HSP70Aa_gene_promoter.bed'
    inputPairedContactsFile = 'aligned_InSitu_HS_1K/merged_nodups.txt'
    outputFile = 'aligned_InSitu_HS_1K/4C_like_output_InSitu_HS_HSP70Aa_test'
    qualityThreshold = 0
    testMode = False
    delimiterForPairedContacts = ' '

    try:
        opts, args = getopt.getopt(argv,"hts:r:w:b:q:o:g:",["test", "ofile="])
    except getopt.GetoptError:
        print '4C_like_plot_cumulative_for_juicer.py -s <input TSS file> -r <input paired contacts file> -w <window size> -b <bin size> -q <quality threshold> -o <outputfile> -t'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '4C_like_plot_cumulative_for_juicer.py -s <input TSS file> -r <input paired contacts file> -w <window size> -b <bin size> -q <quality threshold> -o <outputfile> -t'
            sys.exit()
        elif opt in ("-s"):
            inputTSSFile = arg
        elif opt in ("-r"):
            inputPairedContactsFile = arg
        elif opt in ("-w"):
            windowSize = int(arg)
        elif opt in ("-b"):
            binSize = int(arg)
        elif opt in ("-q"):
            qualityThreshold = int(arg)
        elif opt in ("-o", "--ofile"):
            outputFile = arg
        elif opt in ("-t", "--test"):
            testMode = True
    
    print '\nInput TSS file is:', inputTSSFile
    print 'Input paired contacts file is:', inputPairedContactsFile
    print 'Output file is:', outputFile
    print 'Window size is:', str(windowSize)
    print 'Bin size is:', str(binSize)
    print 'Quality threshold is:', str(qualityThreshold)
    print 'Test mode is:', str(testMode)
    print 'Opts:', opts, '\n'
    
    # Set testing flags
    testingNumberOfRows = 10000000

    # Initialize
    output_folder  = ''
    outputFile     = output_folder + outputFile

    # Read TSS file
    #chr2L	87381	87382	CG11450-RB	1052.35161509999	+
    rowNumber = -1
    with open(inputTSSFile, 'r') as inf:
        for line in inf:
            line = line.strip('\n').strip('\r').split('\t')

            rowNumber += 1
            if testMode:
                if rowNumber >= 2: break
        
            chromo = line[0]
            if chromo[:3] == 'chr': chromo = chromo[3:]
            start = int(line[1])
            end = int(line[2])
            name = line[3]
            score = line[4]
            strand = line[5]

            dataList = [chromo, start, end, name, score, strand]

            # Store data
            TSSDict[rowNumber] = dataList
            print 'TSSDict[', rowNumber, ']:', TSSDict[rowNumber]

            # Set up interaction counts
            for binNum in range(windowSize / binSize): 
                interactionsDict[rowNumber][binNum] = 0
                #print 'interactionsDict[', rowNumber, '][', binNum, ']:', interactionsDict[rowNumber][binNum]

            # Set up cumulative counts
            if chromo not in cumulativeDict: cumulativeDict[chromo] = defaultdict(dict)

    # Set up bins
    binStartList = []
    for binNum in range(windowSize / binSize): 
        binStart = (0 - (windowSize / 2)) + (binNum * binSize)
        binEnd = (0 - (windowSize / 2)) + (binNum * binSize) + binSize
        binRange[binNum] = range(binStart, binEnd)
        binStartList.append(binStart)
        #if testMode and binNum == 0: print 'binNum:', binNum, ':', binRange[binNum], '\n'

    # Process paired contacts file
    #0 chr2L 11 0 0 chr2L 48 0 0 
    print 'Processing paired contacts file:', inputPairedContactsFile
    rowNumber = 0
    tssRowNumber = 0
    with open(inputPairedContactsFile, 'r') as inr:
        for line in inr:

            rowNumber += 1
            if testMode:
                if rowNumber%10000000 == 0: print 'Row number:', str(rowNumber)
                if rowNumber > testingNumberOfRows: break
        
            line = line.strip('\n').strip('\r').split(delimiterForPairedContacts)
            if len(line) < 4: continue

            #print 'line:', line
            strand1 = int(line[0])
            chromo1 = line[1]
            position1 = int(line[2])
            strand2 = int(line[4])
            chromo2 = line[5]
            position2 = int(line[6])
            quality1 = int(line[8])
            quality2 = int(line[11])

            if strand1 == '0': strand1 = '+'
            else: strand1 = '-'
            if strand2 == '0': strand2 = '+'
            else: strand2 = '-'

            if quality1 >= qualityThreshold and quality2 >= qualityThreshold:
                # Loop thru start sites and count paired contacts
                for tssRowNumber in TSSDict:
                    if InTSS(tssRowNumber, chromo1, position1, strand1):
                        binNum = GetBin(tssRowNumber, chromo2, position2, strand2)
                        cumulativeBinNum = GetCumulativeBin(tssRowNumber, chromo2, position2, strand2)
                        #print 'rowNumber:', str(tssRowNumber), 'binNumber:', str(binNum)
                        if binNum in interactionsDict[tssRowNumber]: interactionsDict[tssRowNumber][binNum] += 1
                        if cumulativeBinNum > -1: 
                            #print 'chromo2:', chromo, 'cumulativeBinNum:', str(cumulativeBinNum)
                            if cumulativeBinNum in cumulativeDict[chromo2]: cumulativeDict[chromo2][cumulativeBinNum] += 1
                            else: cumulativeDict[chromo2][cumulativeBinNum] = 1
                    #print 'chromo2:', chromo2, ' position2:', position2
                    if InTSS(tssRowNumber, chromo2, position2, strand2):
                        binNum = GetBin(tssRowNumber, chromo1, position1, strand1)
                        cumulativeBinNum = GetCumulativeBin(tssRowNumber, chromo1, position1, strand1)
                        #print 'rowNumber:', str(tssRowNumber), 'binNumber:', str(binNum)
                        if binNum in interactionsDict[tssRowNumber]: interactionsDict[tssRowNumber][binNum] += 1
                        if cumulativeBinNum > -1: 
                            #print 'chromo1:', chromo, 'cumulativeBinNum:', str(cumulativeBinNum)
                            if cumulativeBinNum in cumulativeDict[chromo1]: cumulativeDict[chromo1][cumulativeBinNum] += 1
                            else: cumulativeDict[chromo1][cumulativeBinNum] = 1

    # Write data to output file
    outputFileCumulative = outputFile + '.bedgraph'
    print 'Writing to cumulative output file:', outputFileCumulative
    with open(outputFileCumulative, 'w') as outf:
        header = 'track type=bedGraph name="BedGraph Format" description="BedGraph format" priority=20'
        #outf.write(header + '\n')
        for chromo in sorted(cumulativeDict):
            for binNum, data in sorted(cumulativeDict[chromo].items()):
                outf.write(chromo + '\t' + str(binNum) + '\t' + str(binNum + (binSize - 1)) + '\t' + str(data) + '\n')
                #print 'chromo:', chromo, 'binNum:', str(binNum), ', count:', str(data)

    # Plot results
    print 'Writing to individual output files'
    for rowNumber in sorted(interactionsDict):
        totalCount = 0
        binList = []
        tss = TSSDict[rowNumber]
        tssChromo = tss[0]
        tssPosition = tss[1]
        tssName = tss[3]
        tssStrand = tss[5]
        outputFile_temp = outputFile + '_' + tssChromo + '_' + str(tssPosition) + '_' + str(tssStrand) + '_' + tssName + '.bedgraph'

        with open(outputFile_temp, 'w') as outf:
            for binNumber, data in sorted(interactionsDict[rowNumber].items()):
                binList.append(data)
                totalCount += data
                binStart = (tssPosition - (windowSize / 2)) + (binNumber * binSize)
                if binStart > 0: outf.write(tssChromo + '\t' + str(binStart) + '\t' + str(binStart + (binSize - 1)) + '\t' + str(data) + '\n')
                if testMode and data > 0: print 'rowNumber:', str(rowNumber), 'binNumber:', str(binNumber), ', count:', str(data)

    # Done
    print ("\nDone!\nTotal count:", totalCount)
    
if __name__ == "__main__":
    main(sys.argv[1:])

