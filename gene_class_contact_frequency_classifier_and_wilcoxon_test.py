# Process table of contact frequency data for multiple gene classes

import os, sys
import timing

import numpy as np
import math
import getopt
from collections import defaultdict

import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
from sklearn.utils.fixes import signature

# Load various scikit classifiers
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import auc

# Load pandas
import pandas as pd

# Load stats from scipy
from scipy import stats

# Enable R script calling
import rpy2.robjects as robjects

# Set random seed
np.random.seed(0)

# Initialize globals
interactionsDict = defaultdict(dict)

geneClassList = ['HSF-dep', 'HSF-indep', 'Down reg', 'Unreg']

columnList = ['CF to nearest transcribed', 'CF to nearest non-transcribed', 'CF to before bin transcribed', 'CF to before bin non-transcribed', 'CF to after bin transcribed', 'CF to after bin non-transcribed']

# Process data file
def ProcessInteractionsFile(inputDataFile, featureNames, targetNames, testMode):
    global interactionsDict
    global geneClassList
    global columnList

    rowNumber = -1
    nRows = 0
    with open(inputDataFile, 'r') as inf:
        for line in inf:
            line = line.strip('\n').strip('\r').split('\t')

            rowNumber += 1
        
            if testMode: 
                print 'line:', line
                if rowNumber >= 10: break

            # Set up non-counted columns
            geneName = line[0]
            interactionsDict[geneName] = defaultdict(dict)

            interactionsDict[geneName]["geneName"] = line[0]
            interactionsDict[geneName]["HSF-dep"] = line[1]
            interactionsDict[geneName]["HSF-indep"] = line[2]
            interactionsDict[geneName]["Down reg"] = line[3]
            interactionsDict[geneName]["Unreg"] = line[4]
            interactionsDict[geneName]["chromo"] = line[5]
            interactionsDict[geneName]["start"] = int(line[6])
            interactionsDict[geneName]["end"] = int(line[7])
            interactionsDict[geneName]["score"] = float(line[8])
            interactionsDict[geneName]["strand"] = line[9]
            interactionsDict[geneName]["geneClass"] = line[10]
            interactionsDict[geneName]["distance transcribed"] = int(line[11])
            interactionsDict[geneName]["distance non-transcribed"] = int(line[12])
            interactionsDict[geneName]["CF to nearest transcribed"] = int(line[13])
            interactionsDict[geneName]["CF to nearest non-transcribed"] = int(line[14])
            interactionsDict[geneName]["CF to before bin transcribed"] = int(line[15])
            interactionsDict[geneName]["CF to before bin non-transcribed"] = int(line[16])
            interactionsDict[geneName]["CF to after bin transcribed"] = int(line[17])
            interactionsDict[geneName]["CF to after bin non-transcribed"] = int(line[18])
            interactionsDict[geneName]["siteWindowTranscribed"]["CF to all transcribed"] = line[19]
            interactionsDict[geneName]["siteWindowNonTranscribed"]["CF to all non-transcribed"] = line[20]
            interactionsDict[geneName]["siteBinTranscribedFoldchange"] = float(line[21])
            interactionsDict[geneName]["siteBinNonTranscribedFoldchange"] = float(line[22])
            interactionsDict[geneName]["siteWindowTranscribed"]["foldchange"] = line[23]
            interactionsDict[geneName]["siteWindowNonTranscribed"]["foldchange"] = line[24]
            interactionsDict[geneName]["siteWindowTranscribed"]["distance"] = line[25]
            interactionsDict[geneName]["siteWindowNonTranscribed"]["distance"] = line[26]
            interactionsDict[geneName]["polII"] = int(line[27])

            # Custom features
            interactionsDict[geneName]["random"] = 0.0
            interactionsDict[geneName]["transcribedDosage"] = interactionsDict[geneName]["CF to nearest transcribed"] * interactionsDict[geneName]["siteBinTranscribedFoldchange"]
            interactionsDict[geneName]["modifiedTranscribedDosageForClosest"] = 0.0
            interactionsDict[geneName]["modifiedTranscribedDosageForAllPeaks"] = 0.0

            if interactionsDict[geneName]["distance transcribed"] <= interactionsDict[geneName]["distance non-transcribed"]:
                interactionsDict[geneName]["CF to nearest"] = interactionsDict[geneName]["CF to nearest transcribed"]
                interactionsDict[geneName]["siteBinFoldchange"] = interactionsDict[geneName]["siteBinTranscribedFoldchange"]
            else:
                interactionsDict[geneName]["CF to nearest"] = interactionsDict[geneName]["CF to nearest non-transcribed"]
                interactionsDict[geneName]["siteBinFoldchange"] = interactionsDict[geneName]["siteBinNonTranscribedFoldchange"]

            for tarNam in targetNames:
                if interactionsDict[geneName][tarNam] == '1': nRows += 1

    print 'Number of genes for', inputDataFile, 'is:', str(rowNumber + 1)

    nCols = len(featureNames)
    interactionsData = np.ndarray(shape=(nRows, nCols), dtype=float)
    target = np.ndarray(shape=(nRows), dtype=int)
    geneNamesArray = np.ndarray(shape=(nRows), dtype='|S35')

    # Load features and target values
    rowNumber = -1
    for geneName in interactionsDict:
        for itemNumber in range(len(targetNames)):
            tarNam = targetNames[itemNumber]
            #if interactionsDict[geneName][tarNam] == '1':
            if interactionsDict[geneName]['geneClass'] == tarNam:
                rowNumber += 1
                target[rowNumber] = itemNumber
                geneNamesArray[rowNumber] = geneName

                for featNumber in range(len(featureNames)):
                    featNam = featureNames[featNumber]
                    interactionsData[rowNumber, featNumber] = interactionsDict[geneName][featNam]

    print 'Number of target rows for', inputDataFile, 'is:', str(rowNumber + 1)

    return interactionsData, target, geneNamesArray

def plotData(data_to_plot, figureNum, colorList, alphaList, xAxisLabels, yAxis, title, outputFile, barCount):
    # Create a figure instance
    figWidth = barCount * 1.25
    fig = plt.figure(figureNum, figsize=(figWidth, 6))

    # Create an axes instance
    ax = fig.add_subplot(111)

    ## add patch_artist=True option to ax.boxplot() 
    ## to get fill color
    bp = ax.boxplot(data_to_plot, patch_artist=True)

    ## change outline color, fill color and linewidth of the boxes
    boxNumber = 0
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = colorList[boxNumber], alpha = alphaList[boxNumber] )
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
    #plt.ylim([-0.15,0.15])
    plt.tick_params(axis='y', which='major', labelsize=14)

    # draw horizontal line at zero
    xLineList = []
    yLineList = []
    for i in range(barCount + 2):
        xLineList.append(i)
        yLineList.append(0)
    plt.plot(xLineList, yLineList, lw=1, linestyle='--', color='0.25')
    #plt.plot([0, 1, 2, 3], [0, 0, 0, 0], lw=1, linestyle='--', color='0.25')

    # Save the figure
    #fig.savefig(outputFile, bbox_inches='tight', format='pdf')
    fig.savefig(outputFile+'.pdf', bbox_inches='tight', format='pdf')
    fig.savefig(outputFile+'.png', bbox_inches='tight', format='png')
    plt.show()

def main(argv):
    global interactionsDict
    global geneClassList
    global columnList

    programName = 'gene_class_contact_frequency_random_forest_wilcoxon'
    inputDataFile = 'K562_gene_class_contact_frequency_3K_w_PROseq_counts_and_rpkm_minus_HS_table.txt'
    outputFile = 'K562_gene_class_contact_frequency_3K_w_PROseq_counts_and_rpkm_minus_HS_wilcoxon_boxplot_1000'
    outputScoresFile = 'classifier_results_reversed_labels.txt'
    testMode = False
    delimiterForPairedContacts = ' '

    try:
        opts, args = getopt.getopt(argv,"hr:q:o:",["test", "ofile=", "data="])
    except getopt.GetoptError:
        print 'gene_class_contact_frequency_random_forest_wilcoxon.py --data <input data file>-o <outputfile prefix> --test'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'gene_class_contact_frequency_random_forest_wilcoxon.py --data <input data file> -o <outputfile prefix> --test'
            sys.exit()
        elif opt in ("--data"):
            inputDataFile = arg
        elif opt in ("-o", "--ofile"):
            outputFile = arg
        elif opt in ("-t", "--test"):
            testMode = True
    
    print ''
    print 'Input data file is:', inputDataFile
    print 'Output file prefix is:', outputFile
    print 'Test mode is:', str(testMode)
    print ''
    
    # Set testing flags
    testingNumberOfRows = 1000000

    # Initialize
    output_folder  = ''
    outputFile     = output_folder + outputFile

    plotTitle = 'PR AUC'
    barTitleList = ['control', 'mod all', 'mod close', 'cf/ps', 'cf/ps txs', 'cf/ps/dist txs']
    barColorList = ['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue']
    alphaList = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
    barCount = 6
    barLists = [[] for _ in range(barCount)]
    data_to_plot = []
    numLoops = 1000

    # Set up multiple trials
    trials = []
    trials = [['random'], ['modifiedTranscribedDosageForAllPeaks'], ['modifiedTranscribedDosageForClosest'], ['CF to nearest', 'siteBinFoldchange'], ['CF to nearest transcribed', 'CF to nearest non-transcribed', 'siteBinTranscribedFoldchange', 'siteBinNonTranscribedFoldchange'],['distance transcribed', 'distance non-transcribed', 'CF to nearest transcribed', 'CF to nearest non-transcribed', 'siteBinTranscribedFoldchange', 'siteBinNonTranscribedFoldchange']]

    for trialNum in range(len(trials)):
        featureNames = trials[trialNum]

        # Target names must be sorted alphbetically for pd.factorize to work correctly
        targetNames = np.sort(np.array(['HSF-dep', 'HSF-indep']))[::-1]

        # Read Data file
        interactionsData, target, geneNamesArray = ProcessInteractionsFile(inputDataFile, featureNames, targetNames, testMode)

        #print interactionsData
        print featureNames
        print targetNames
        print 'Target:', target[0:10]
        print 'geneNamesArray:', geneNamesArray[0:10]
    
        # Create a dataframe with the feature variables
        df = pd.DataFrame(interactionsData, columns=featureNames)

        # View the top rows
        #print df.head(n=10)

        # Add a new column with the target names, this is what we are going to try to predict
        df['gene_class'] = pd.Categorical.from_codes(target, targetNames)

        # Add gene names column
        geneNamesDf = pd.DataFrame(geneNamesArray)
        df = pd.concat([df, geneNamesDf], axis=1)

        # View the top 5 rows
        print df.head(n=3)

        # Create a new column that for each row, generates a random number between 0 and 1, and
        # if that value is less than or equal to .8, then sets the value of that cell as True
        # and false otherwise. This is a quick and dirty way of randomly assigning some rows to
        # be used as the training data and some as the test data.
        for loopNum in range(numLoops):
            df['is_train'] = np.random.uniform(0, 1, len(df)) <= .8

            # View the top 5 rows
            #print df.head(n=5)

            # Write taining flags and call R script to produce values for a, b, and z
            if featureNames[0] == 'modifiedTranscribedDosageForClosest' or featureNames[0] == 'modifiedTranscribedDosageForAllPeaks':
                df['is_train'].to_csv(r'is_train.txt', header=False, index=None, sep='\t', mode='w')
                r=robjects.r
                r.source("get_abz.R")
                with open('tmp.abz.out', 'r') as infabz:
                    for values in infabz:
                        values = values.strip('\n').strip('\r').split('\t')
                        if values[0] == 'b': continue
                        b = float(values[0])
                        a = float(values[1])
                        z = float(values[2])
                        print 'a:', str(a), ' b:', str(b), ' z:', str(z)

                r.source("get_dose.R")
                doseDict = {}
                with open('tmp.dose.out', 'r') as infdose:
                    for values in infdose:
                        values = values.strip('\n').strip('\r').split('\t')
                        if values[0] == 'geneID': continue
                        doseDict[values[0]] = values[1]
                        #print 'gene:', values[0], ' dose:', values[1]

                indexNum = -1
                for geneName in geneNamesArray:
                    indexNum += 1
                    # Scale dosage using sigmoid curve
                    txcf = interactionsDict[geneName]["CF to nearest transcribed"]
                    txps = interactionsDict[geneName]["siteBinTranscribedFoldchange"]
                    ntxcf = interactionsDict[geneName]["CF to nearest non-transcribed"]
                    ntxps = interactionsDict[geneName]["siteBinNonTranscribedFoldchange"]

                    txAdjusted = (1 / (1 + math.exp(-a * (txcf - b)))) * txps
                    ntxAdjusted = z * (1 / (1 + math.exp(-a * (ntxcf - b)))) * ntxps

                    #if indexNum < 5: print 'index:', str(indexNum), ' geneName:', geneName, ' txAdjusted:', str(txAdjusted + ntxAdjusted)
                    if featureNames[0] == 'modifiedTranscribedDosageForClosest':
                        df.loc[indexNum, 'modifiedTranscribedDosageForClosest'] = txAdjusted + ntxAdjusted
                    if featureNames[0] == 'modifiedTranscribedDosageForAllPeaks':
                        df.loc[indexNum, 'modifiedTranscribedDosageForAllPeaks'] = doseDict[geneName]

            # View the top 5 rows
            #print df.head(n=5)

            # Create two new dataframes, one with the training rows, one with the test rows
            train, test = df[df['is_train']==True], df[df['is_train']==False]

            # Show the number of observations for the test and training dataframes
            #print('Number of observations in the training data:', len(train))
            #print('Number of observations in the test data:',len(test))

            # Create a list of the feature column's names
            nCols = len(featureNames)
            features = df.columns[:nCols]

            # View features
            #print 'Features:', features

            # train['gene_class'] contains the actual gene class names. Before we can use it,
            # we need to convert each gene class name into a digit. So, in this case there
            # are two gene class's, which have been coded as 0 or 1.
            #print train['gene_class'].head()
            y = pd.factorize(train['gene_class'], sort=True)[0]

            # View target
            #print train['gene_class'][0:30]
            #print 'Target:', y[0:30]

            # Create a random forest Classifier. By convention, clf means 'Classifier' 
            #clf = RandomForestClassifier(n_jobs=2, random_state=0, n_estimators=100)
            #clf = KNeighborsClassifier(2)
            #clf = SVC(kernel="linear", C=0.025)
            #clf = MLPClassifier(alpha=1)
            clf = AdaBoostClassifier(n_estimators=100)
            scores = cross_val_score(clf, interactionsData, target, cv=5)
            #print 'Scores:', scores

            # Train the Classifier to take the training features and learn how they relate
            # to the training y (the gene class)
            #print 'clf:', clf.fit(train[features], y)
            clf.fit(train[features], y)

            y_test = pd.factorize(test['gene_class'], sort=True)[0]
            accuracy = clf.score(test[features], y_test)
            #print "Accuracy: ", accuracy 

            #print y_test[0:20]
            #print clf.predict(test[features])[0:20]

            # Apply the Classifier we trained to the test data (which, remember, it has never seen before)
            #print 'targetNames:', targetNames
            #print test[features].head()
            if featureNames[0] == 'random':
                y_score = clf.predict(test[features])
                HSFdepCutoff = float(227)/1172
                #print 'HSFdepCutoff:',HSFdepCutoff
                for i in range(len(y_score)):
                    randNum = np.random.uniform()
                    if randNum <= HSFdepCutoff: y_score[i] = 1
                    else: y_score[i] = 0
            else:
                y_score = clf.predict(test[features])

            #print 'y_test: ', y_test[0:30]
            #print 'y_score:', y_score[0:30]

            # View the predicted probabilities of the first 10 observations
            y_probs = clf.predict_proba(test[features])
            #print 'Probs:', y_probs[0:10]

            # Create actual english names for the classes for each predicted gene class
            preds = targetNames[y_score]

            # View the PREDICTED species for the first five observations
            #print 'Predictions:', preds[0:10]

            # View the ACTUAL species for the first five observations
            #print 'Actual:', test['gene_class'][0:10]

            # Create confusion matrix
            #print 'Confusion matrix:'
            #print pd.crosstab(test['gene_class'], preds, rownames=['Actual Gene Class'], colnames=['Predicted Gene Class'])
    
            # View a list of the features and their importance scores
            #print 'Importance:', list(zip(train[features], clf.feature_importances_))

            precision, recall, thresholds = precision_recall_curve(y_test, y_score)
            average_precision = average_precision_score(y_test, y_score)

            #precision, recall, thresholds = precision_recall_curve(y_true[:, c], y_pred[:, c])
            pr_auc = auc(recall, precision)
            print 'pr_auc:', pr_auc
            #print 'average_precision:', average_precision

            barLists[trialNum].append(pr_auc)

        data_to_plot.append(barLists[trialNum])

    # Build dataframe from plot data add perform wilcoxon test
    auprcDF = pd.DataFrame(data_to_plot).T
    auprcDF.columns = barTitleList

    # Wilcoxon test on on first column compared to each of the others
    wilcoxonResults = []
    for col in range(1,len(barTitleList)):
        wilcoxonResults.append(stats.wilcoxon(auprcDF[barTitleList[0]], auprcDF[barTitleList[col]]))

    #result = stats.wilcoxon(auprcDF[barTitleList[0]], auprcDF[barTitleList[1]])
    print 'results:'
    print wilcoxonResults

    # Write columns to text file
    auprcDF.to_csv(outputFile+'.txt', header=True, index=None, sep='\t', mode='w')

    # Write Wilcoxon results
    with open(outputFile+'_stats.txt', 'w') as outStats:
        outStats.write('\t'.join(str(val) for val in wilcoxonResults) + '\n')

    # Boxplot
    yAxis = 'Area under PR curve'
    xAxisLabels = barTitleList
    plotData(data_to_plot, 0, barColorList, alphaList, xAxisLabels, yAxis, plotTitle, outputFile, barCount)

    # Done
    print ''
    print 'Done!'
    print ''

if __name__ == "__main__":
    main(sys.argv[1:])

