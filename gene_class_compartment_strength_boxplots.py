# Boxplot for gene classes

import numpy as np
import scipy.stats
import getopt
import random
import math

from collections import defaultdict
from sklearn.utils import resample

import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt

import os, sys
import timing

# Set up globals
useLogs = False

def plotData(data_to_plot, colorList, xAxisLabels, yAxis, title, outputFile):
    # Create a figure instance
    fig = plt.figure(1, figsize=(9, 6))

    # Create an axes instance
    ax = fig.add_subplot(111)

    # Create the boxplot
    bp = ax.boxplot(data_to_plot)

    ## add patch_artist=True option to ax.boxplot() 
    ## to get fill color
    bp = ax.boxplot(data_to_plot, patch_artist=True)

    ## change outline color, fill color and linewidth of the boxes
    boxNumber = 0
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = colorList[boxNumber] )
        boxNumber += 1
        #box.set( facecolor = '#1b9e77' )

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)

    ## Custom x-axis labels
    ax.set_xticklabels(xAxisLabels, fontsize=14)

    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    fig.suptitle(title, fontsize=14)
    #plt.xlabel('xlabel', fontsize=14)
    plt.ylabel(yAxis, fontsize=14)
    plt.ylim([-0.15,0.15])
    #plt.ylim([-0.04,0.04])
    plt.tick_params(axis='y', which='major', labelsize=14)

    # draw horizontal line at zero
    plt.plot([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], lw=1, linestyle='--', color='0.25')

    # Save the figure
    #fig.savefig(outputFile, bbox_inches='tight', format='pdf')
    fig.savefig(outputFile+'.pdf', bbox_inches='tight', format='pdf')
    fig.savefig(outputFile+'.png', bbox_inches='tight', format='png')
    plt.show()

def main(argv):

    inputFile1 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_HSF-dependent.txt'
    inputFile2 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_upreg.txt'
    inputFile3 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_downreg.txt'
    inputFile4 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_unreg.txt'
    inputFile5 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_HSF-dependent.txt'
    inputFile6 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_upreg.txt'
    inputFile7 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_downreg.txt'
    inputFile8 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_unreg.txt'
    inputFile9 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_downreg.txt'
    inputFile10 = '../juicer/aligned_S2_inSitu_HS/S2_inSitu_HS_KRnorm_chr3R_50K_eigen_unreg.txt'
    xAxis = ',HSF-dependent,,HSF-indpendent,,Down Regulated,,Unregulated,,Silent'
    yAxis = 'Compartment Strength'
    title = 'Comparison of Heat Shock Genes in S2 Cells'
    outputFile = 'S2_gene_class_boxplot'
    testMode = False

    try:
        opts, args = getopt.getopt(argv,"htx:y:l:a:b:c:d:e:f:g:i:j:k:o:",["test"])
    except getopt.GetoptError:
        print 'gene_class_boxplots_eight_boxes.py -a <inputfile 1> -b <inputfile 2> -c <inputfile 3> -d <inputfile 4> -e <inputfile 5> -f <inputfile 6> -g <inputfile 7> -i <inputfile 8> -j <inputfile 9> -k <inputfile 10> -o <outputfile> -l <title> -x <x-axis label> -y <y-axis label> -t'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'gene_class_boxplots_eight_boxes.py -a <inputfile 1> -b <inputfile 2> -c <inputfile 3> -d <inputfile 4> -e <inputfile 5> -f <inputfile 6> -g <inputfile 7> -i <inputfile 8> -j <inputfile 9> -k <inputfile 10> -o <outputfile> -l <title> -x <x-axis label> -y <y-axis label> -t'
            sys.exit()
        elif opt in ("-a"):
            inputFile1 = arg
        elif opt in ("-b"):
            inputFile2 = arg
        elif opt in ("-c"):
            inputFile3 = arg
        elif opt in ("-d"):
            inputFile4 = arg
        elif opt in ("-e"):
            inputFile5 = arg
        elif opt in ("-f"):
            inputFile6 = arg
        elif opt in ("-g"):
            inputFile7 = arg
        elif opt in ("-i"):
            inputFile8 = arg
        elif opt in ("-j"):
            inputFile9 = arg
        elif opt in ("-k"):
            inputFile10 = arg
        elif opt in ("-o"):
            outputFile = arg
        elif opt in ("-l"):
            title = arg
        elif opt in ("-x"):
            xAxis = arg
        elif opt in ("-y"):
            yAxis = arg
        elif opt in ("-t", "--test"):
            testMode = True

    # Set testing flags
    testingNumberOfRows = 10
    loadTestRows = 1000000
    rowNumber = 0

    # Initialize
    if useLogs: sampleSize = 10000
    else: sampleSize = 100000000
    work_folder    = ''
    output_folder  = ''
    if useLogs: outputFile = outputFile + '_log'
    outputFile = output_folder + outputFile

    print 'Input file a is:', inputFile1
    print 'Input file b is:', inputFile2
    print 'Input file c is:', inputFile3
    print 'Input file d is:', inputFile4
    print 'Input file e is:', inputFile5
    print 'Input file f is:', inputFile6
    print 'Input file g is:', inputFile7
    print 'Input file i is:', inputFile8
    print 'Input file j is:', inputFile9
    print 'Input file k is:', inputFile10
    print 'Output graph is:', outputFile
    print 'Test mode is:', str(testMode)
    print ''
    
    # Read input files into lists
    list1 = []
    rowNumber = 0
    with open(inputFile1) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list1.append(val)

    print 'First data set loaded'
    print 'Len of list1:', len(list1)

    list2 = []
    rowNumber = 0
    with open(inputFile2) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list2.append(val)

    print 'Second data set loaded'
    print 'Len of list2:', len(list2)

    list3 = []
    rowNumber = 0
    with open(inputFile3) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list3.append(val)

    print 'Third data set loaded'
    print 'Len of list3:', len(list3)

    list4 = []
    rowNumber = 0
    with open(inputFile4) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list4.append(val)

    print 'Fourth data set loaded'
    print 'Len of list4:', len(list4)

    list5 = []
    rowNumber = 0
    with open(inputFile5) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list5.append(val)

    print 'Fifth data set loaded'
    print 'Len of list5:', len(list5)

    list6 = []
    rowNumber = 0
    with open(inputFile6) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list6.append(val)

    print 'Sixth data set loaded'
    print 'Len of list6:', len(list6)

    list7 = []
    rowNumber = 0
    with open(inputFile7) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list7.append(val)

    print 'Seventh data set loaded'
    print 'Len of list7:', len(list7)

    list8 = []
    rowNumber = 0
    with open(inputFile8) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list8.append(val)

    print 'Eighth data set loaded'
    print 'Len of list8:', len(list8)

    list9 = []
    rowNumber = 0
    with open(inputFile9) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list9.append(val)

    print 'Nineth data set loaded'
    print 'Len of list9:', len(list9)

    list10 = []
    rowNumber = 0
    with open(inputFile10) as f:
        for line in f:
            rowNumber += 1
            if testMode and rowNumber > loadTestRows: break
            line = line.strip('\n').strip('\r')
            if line == 'NaN' or line == 'nan': line = '0.0'
            val = float(line)
            list10.append(val)

    print 'Tenth data set loaded'
    print 'Len of list10:', len(list10)

    ## Create test data
    #np.random.seed(10)
    #collectn_1 = np.random.normal(100, 10, 200)
    #collectn_2 = np.random.normal(80, 30, 200)
    #collectn_3 = np.random.normal(90, 20, 200)
    #collectn_4 = np.random.normal(70, 25, 200)

    ## combine these different collections into a list    
    data_to_plot = [list1, list2, list3, list4, list5, list6, list7, list8, list9, list10]
    colorList = ['darkred', 'darkred', 'red', 'red', 'blue', 'blue', 'darkorange', 'darkorange', 'darkgrey', 'darkgrey']

    # Plot
    xAxisLabels = xAxis.split(',')
    plotData(data_to_plot, colorList, xAxisLabels, yAxis, title, outputFile)

    print ("\nDone!")
    
if __name__ == "__main__":
    main(sys.argv[1:])


