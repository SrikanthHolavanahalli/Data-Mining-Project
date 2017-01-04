import random
import os
import math
import operator
#import matplotlib.pyplot as plt
#from sklearn import metrics



muproEssentialData = {}
muproEssentialDataPositive = {}
muproEssentialDataNegative = {}



def readFromMupro(filename):


    chainId = ""
    linecount = 1
    totalItemCount = 0
    checkDupData ={}

    AnalysisProtMap = {}
    AnalysisProtChainMap = {}


    fST = open('seqPlusStrData.txt' ,'w')

    fpos = open('muproPos.txt' ,'w')
    fneg = open('muproNeg.txt' ,'w')



    sequence = ""

    with open(filename) as f:
        for line in f:

            if linecount % 10 == 1:
                mutData = line.split(" ")
                chainId =mutData[0][-2]
            if linecount % 10 == 3:
                sequence = line[:-1]

            if linecount % 10 == 9 :
                mutData = line.split(" ")
                protId = mutData[0]
                wildtype = mutData[1]
                position = int(mutData[2])
                mutation = mutData[3]
                SA = mutData[4]
                pH = mutData[5]
                temp = mutData[6]
                energy_change = mutData[7][:-1]
                sign = "+"
                if float(energy_change) <= 0:
                    sign = "-"


                reqVals = []
                uniqID = protId + chainId  + wildtype + str(position) + mutation + pH + temp  + sign

                #protID + chainID + position + wildType + mutation

                if uniqID not in checkDupData.keys():


                # Analysis code of no. of postive and negative examples per protien
                    if uniqID not in checkDupData.keys():
                        if protId not in AnalysisProtMap.keys():
                            if sign == "+":
                                AnalysisProtMap[protId] = [1,0]
                            else:
                                AnalysisProtMap[protId] = [0,1]
                            totalItemCount +=1
                        else:
                            if sign == "+":
                                AnalysisProtMap[protId][0] += 1
                            else:
                                AnalysisProtMap[protId][1] += 1
                            totalItemCount +=1



                # Analysis code of no. of postive and negative examples per protien chain
                    protPlusChainID = protId + chainId

                    if uniqID not in checkDupData.keys():
                        if protPlusChainID not in AnalysisProtChainMap.keys():
                            if sign == "+":
                                AnalysisProtChainMap[protPlusChainID] = [1,0]
                            else:
                                AnalysisProtChainMap[protPlusChainID] = [0,1]

                        else:
                            if sign == "+":
                                AnalysisProtChainMap[protPlusChainID][0] += 1
                            else:
                                AnalysisProtChainMap[protPlusChainID][1] += 1


                    window = ""
                    if position - 3 < 0:
                        w = sequence[:position + 4]
                        q = 0

                        while(q < 7- len(w)):
                            window += "-"
                            q += 1

                        window += w

                    elif position + 4 > len(sequence):
                        window = sequence[position-3:]

                        while(len(window) < 7):
                            window += "-"

                    else:
                        window = sequence[position-3:position+4]


                    data = sign + "\t" + protId + "\t" + chainId  + "\t" +  wildtype + "\t"
                    data += str(position) + "\t" +  mutation +  "\t" + pH + "\t"  + temp + "\t" + window +  "\n"


                    if sign == "+":
                        fpos.write(data)
                    else:
                        fneg.write(data)
                    fST.write(data)
                    checkDupData[uniqID] = 1
            else:
                pass
            linecount += 1



    analysis_file_per_protein_str = open('analysis_file_per_protein_str.txt', 'w')

    analysis_file_per_protein_Chain_str = open('analysis_file_per_protein__chain_str.txt', 'w')

    
    for k in AnalysisProtMap.keys():
        fileStr = ""
        fileStr = str(k) + " " + str(AnalysisProtMap[k][0]) + " " + str(AnalysisProtMap[k][1]) + '\n'
        analysis_file_per_protein_str.write(fileStr)

    
    for k in AnalysisProtChainMap.keys():
        fileStr = ""
        fileStr = str(k) +  " " + str(AnalysisProtChainMap[k][0])+  " " + str(AnalysisProtChainMap[k][1]) + '\n'
        analysis_file_per_protein_Chain_str.write(fileStr)


    analysis_file_per_protein_str.close()
    analysis_file_per_protein_Chain_str.close()
    fST.close()
    fpos.close()
    fneg.close()
    x = 1



def convertToSVMFormat(filename):
    #f = open(filename,'r')


    fw = open('myfile.txt','w')


    protID = ""
    protMap = {}

    protPlusChainID =""
    protPlusChainMap = {}

    lineCount = 1
    protIndex = 0
    protPlusChainIndex = 0
    with open(filename) as f:
        for line in f:
            aminoAcidMap = {}
            seqRow = line.split('\t')
            seq = seqRow[5]

            protID = seqRow[0]
            protPlusChainID = seqRow[0] + seqRow[1]

            if protPlusChainID not in protPlusChainMap.keys():
                protPlusChainIndex +=1
                protPlusChainMap[protPlusChainID] = protPlusChainIndex

            if protID not in protMap.keys():
                protIndex += 1
                protMap[protID] = protIndex
            svmStr = ""
            classVal = ""

            if lineCount < len(muproEssentialDataPositive):
                classVal = "1"
            else:
                classVal = "-1"

            lineCount += 1

            svmStr += classVal + " "

            for s in seq:
                if s not in aminoAcidMap.keys():
                    aminoAcidMap[s] = 1
                else:
                    aminoAcidMap[s] += 1

            #def_amino_acid_map ={}
            #for i in range(1,27):
            #    def_amino_acid_map[str(unichr(i + 64))] = i


            neighbor_list=[]

            for k in aminoAcidMap.keys():
                neighbor_list.append(k)

            neighbor_list.sort()

            for i in range(len(neighbor_list)):
                aminoResidue = ord(neighbor_list[i]) - 64
                count = aminoAcidMap[neighbor_list[i]]
                svmStr += str(aminoResidue) + ":" + str(count) + " "

            wild_type_feature_number = ord(seqRow[3]) + 26 - 64
            mutation_feature_number  = ord(seqRow[4]) + 52 - 64



            svmStr += str(wild_type_feature_number) + ":" + "-1" + " " + str(mutation_feature_number) + ":" +"1"
            #svmStr +=  "53:"+ pH+ " " "54:" + temp + " " + "55:" + SA           ########------bottleneck---------------#########
            svmStr += '\n'

            fw.write(svmStr)

    return protMap,protPlusChainMap



# following 4 variables used to create aux data
protID = ""
protMap = {}

protPlusChainID =""
protPlusChainMap = {}


def createAuxDataForStrData(filename, protMap, protPlusChainMap ):

    #Structure data aux files
    fpch = open('StrAuxPerProt.txt','w')
    fpp  = open('StrAuxPerChain.txt','w')


    with open(filename) as f:
        for line in f:
            seqRow = line.split('\t')
            protID = seqRow[0]
            protPlusChainID = seqRow[0] + seqRow[1]

            #aux data per protien
            fpch.write(str(protMap[protID]) + '\n')

            #aux data per chain
            fpp.write( str( protPlusChainMap[protPlusChainID]) + '\n')

    fpch.close()
    fpp.close()



cluster_Protein_map = {}
def createAuxDataForStrDataPerCluster(filename, protMap ,auxfile):


    rowLen = 0
    cluster_id = -1
    protID = ""
    with open(filename) as f:
        for line in f:
            rowLen = len(line)
            if rowLen < 13:
                cluster_id = line[line.index("r") + 2: -1]
            else:
                ind = line.index('>')
                protID = line[ind+1:ind+5]
                cluster_Protein_map[protID] = cluster_id


    fpch = open(auxfile,'w')

    filename = 's1496_pdb.features'

    with open(filename) as f:
        for line in f:
            seqRow = line.split('\t')
            protID = seqRow[0]

            #aux data per protien
            fpch.write(str(cluster_Protein_map[protID]) + '\n')

nFold = 10


strDataHashMap = {}




def readFromStr(filename):
    protIndex = 0
    protPlusChainIndex = 0

    f = open(filename, 'r')
    fseq = open('mupro_1496.txt','r')
    fsvm = open('strSVMfileIP.txt' , 'w')

    for data in f:
        rowVal = data.split('\t')

        seqDataRow = fseq.readline().split('\t')

        protID = rowVal[0]
        chainID = rowVal[1]
        position = rowVal[2]
        wildType = rowVal[3]
        mutation = rowVal[4]
        neighbors = rowVal[5]

        protPlusChainID= protID + chainID

        if protID not in protMap.keys():
            protMap[protID] = protIndex
            protIndex += 1

        if protPlusChainID not in protPlusChainMap.keys():
            protPlusChainMap[protPlusChainID] = protPlusChainIndex
            protPlusChainIndex += 1



        window = seqDataRow[-1][:-1]

        i = 0
        strSVM = ""

        sign = seqDataRow[0]

        if sign == "+":
            strSVM += "+1 "
        else:
            strSVM += "-1 "

        for w in window:
            if w != "-":
                if w != wildType:
                    strSVM += str(ord(w) - 64 + i*26) + ":1 "
                    i += 1
            else:
                i += 1

        if ord(mutation) > ord(wildType):
            strSVM += str(ord(wildType) - 64 + i*26) + ":-1 " + str(ord(mutation) - 64 + i*26) + ":1 "
        else:
            strSVM += str(ord(mutation) - 64 + i*26)+ ":1 " + str(ord(wildType) - 64 + i*26) + ":-1 "

        i += 1

        neighborsMap = {}
        neighborsList = []
        for n in neighbors:
            if n in neighborsMap.keys():
                neighborsMap[n] += 1
            else:
                neighborsMap[n] = 1
            neighborsList.append(n)

        neighborsList = list(set(neighborsList))
        neighborsList.sort()

        for k in range(len(neighborsList)):
            strSVM += str(ord(neighborsList[k]) - 64 + i*26) + ":" + str(neighborsMap[neighborsList[k]]) + " "

        label = -1
        if sign == "+":
            label = 1
        else:
            label = -1

        strSVM += "#" + str(label) + '\n'

        fsvm.write(strSVM)

    fsvm.close()
    x = 1


def nValid(mapFile,auxiliaryFile,dataFile):
    length=len(mapFile)
    sizeOfOneSplit=length/nFold
    left=length-(nFold*sizeOfOneSplit)
    left=nFold-left
    s=range(1,length+1)
    testList=[]
    trainingFile = open('trainData.txt','w')
    testFile = open('testData.txt','w')
    auxiliaryData = open(auxiliaryFile, 'r')
    featureData = open(dataFile, 'r')
    dataList=[]
    auxiliaryList=[]
    for data in featureData:
        dataList.append(data.strip())
    for data in auxiliaryData:
        auxiliaryList.append(data.strip())
    temp = open('tempfile.txt','w')
    temp.close()
    testFinalFile = open('testFinalData.txt','w')
    outputFile = open('finalOutputFile.txt','w')
    for i in range(0,nFold):
        trainingFile = open('trainData.txt','w')
        testFile = open('testData.txt','w')
        for j in range(0,sizeOfOneSplit):
            var=random.choice(list(s))
            testList.append(var)
            s.remove(var)
        if i == left:
            left+=1
            var=random.choice(list(s))
            testList.append(var)
            s.remove(var)
        for traverse in range(0,len(auxiliaryList)):
            number=int(auxiliaryList[traverse])
            if number in testList:
                testFile.write(str(dataList[traverse])+"\n")
                testFinalFile.write(str(dataList[traverse])+"\n")
            else:
                trainingFile.write(str(dataList[traverse])+"\n")
        trainingFile.close()
        testFile.close()
        testList=[]

        commandLearn="svm_learn.exe -z c -t 2 -g 0.1 -j 2 -c 5 trainData.txt model_file"
        os.system(commandLearn)


        commandClassify="svm_classify.exe testData.txt model_file output.txt >> tempfile.txt"

        os.system(commandClassify)
        outputTemp = open('output.txt', 'r')

        for out in outputTemp:
            outputFile.write(out)
        outputTemp.close()
    testFinalFile.close()
    outputFile.close()
    analysis()


def sigmod(value):
    return 1/(1+ math.exp(-value))


totalPrecision = 0
totalRecall = 0
totalF1Measure = 0
totalAccuracy = 0
totalCorelationCofficient = 0




def analysis():
    outputFinalFile = open('finalOutputFile.txt', 'r')
    testFinalFile = open('testFinalData.txt', 'r')
    truePositive=0
    falsePositive=0
    trueNegative=0
    falseNegative=0
    length=testFinalFile.read().split('\n')
    testFinalFile.close()
    testFinalFile = open('testFinalData.txt', 'r')


    rocList = []

    for data in range(0,len(length)-1):
        output=outputFinalFile.readline()
        classLabel=int(testFinalFile.readline().split('#')[1].strip())
        prediction=-1
        output=sigmod(float(output))
        if output >= 0.5:
            prediction=+1
        rocList.append([float(output),int(classLabel)])

        if classLabel == 1 and prediction == 1:
            truePositive+=1
        elif classLabel == 1 and prediction == -1:
            falseNegative+=1
        elif classLabel == -1 and prediction == -1:
            trueNegative+=1
        elif classLabel == -1 and prediction == 1:
            falsePositive+=1
    precision=float(truePositive)/(truePositive+falsePositive)*100
    recall=truePositive/float(truePositive+falseNegative)*100
    f1Measure=float(2*truePositive)/(2*truePositive+falsePositive+falseNegative)*100
    accuracy=float(truePositive+trueNegative)/(truePositive+trueNegative+falseNegative+falsePositive)*100

    corelationCofficient=((trueNegative*truePositive)-(falsePositive*falseNegative))/(math.sqrt((truePositive+falseNegative)*(truePositive+falsePositive)*(trueNegative+falseNegative)*(trueNegative+falsePositive)))

    tpr=recall
    fpr=falsePositive/float(trueNegative+falsePositive)*100
    print "\n\n"


    #print "Precision: "+str(precision)+"%"
    #print "Recall: "+str(recall)+"%"
    #print "f1-Measure: "+str(f1Measure)+"%"
    #print "Accuracy: "+str(accuracy)+"%"
    #print "Corelation Cofficient: "+str(corelationCofficient)

    global totalPrecision
    global totalRecall
    global totalF1Measure
    global totalAccuracy
    global totalCorelationCofficient

    totalPrecision += precision
    totalRecall += recall
    totalF1Measure += f1Measure
    totalAccuracy += accuracy
    totalCorelationCofficient += corelationCofficient

    x = 1
    x += 1
    #print "True Positive Rate: "+str(tpr)
    #print "False Positive Rate: "+str(fpr)
    #print "check"

    plotROC(rocList)



def plotROC(rocList):
    #list1 = sorted(rocList, key=operator.itemgetter(0), reverse=True)
    sortedList = sorted(rocList, key=operator.itemgetter(0))
    #######################################
    trueClass=[]
    logSigValue=[]
    for data in sortedList:
        logSigValue.append(data[0])
        trueClass.append(data[1])
    fpr, tpr, thresholds = metrics.roc_curve(trueClass, logSigValue)
    roc_auc= metrics.auc(fpr, tpr)
    print fpr
    print tpr
    print thresholds
    ###################################################
    """
    tpr=[]
    fpr=[]
    #for i in increment:
    for i in sortedList:
        truePositive=0
        falsePositive=0
        trueNegative=0
        falseNegative=0
        for data in sortedList:
            prediction = float(data[0])
            dataSign=1
            if prediction < i[0]:
                dataSign=-1
            classLabel=int(data[1])
            if classLabel == 1 and dataSign == 1:
                truePositive+=1
            elif classLabel == 1 and dataSign == -1:
                falseNegative+=1
            elif classLabel == -1 and dataSign == -1:
                trueNegative+=1
            elif classLabel == -1 and dataSign == 1:
                falsePositive+=1
        tpr.append(truePositive/float(truePositive+falseNegative))
        fpr.append(falsePositive/float(trueNegative+falsePositive))
    """
    plt.figure()
    plt.plot(fpr, tpr,'g', label='AUC = %0.2f'% roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic for Structure only per protein chain')
    plt.legend(loc="lower right")
    
    plt.savefig('ROC for Structure Only per Protein chain')
    plt.show()
    #print "check"



def mergePosNeg(posFile, negFile):

    f = open('mupro_1496.txt' , 'w')

    fPos = open(posFile,  'r')
    fNeg = open(negFile,  'r')

    for data in fPos:
        f.write(data)

    for data in fNeg:
        f.write(data)

    f.close()
    fPos.close()
    fNeg.close()

def main():

    #this call will give me 1496 unique values from mupro.data and divides them to postives , negatives energy changes hashmaps
    #readFromMupro('mupro.data')


    #divide data to pos and neg files
    readFromMupro('mupro.data')

    #merge mupro pos and neg file
    mergePosNeg('muproPos.txt' , 'muproNeg.txt' )

    #read str data and create feature file
    readFromStr('s1496_pdb.features')

    #per protien and per chain aux file created in same function
    createAuxDataForStrData('s1496_pdb.features', protMap, protPlusChainMap)

   #per protein

    CVIterations = 5

    #while(CVIterations > 0):
    #    nValid(protMap, 'StrAuxPerProt.txt',  'strPlusSeqSVMfileIP.txt' )
    #    CVIterations -= 1

    #per chain
    #nValid(protPlusChainMap, 'StrAuxPerChain.txt',  'strPlusSeqSVMfileIP.txt' )

    while(CVIterations > 0):
        nValid(protMap, 'StrAuxPerChain.txt',  'strPlusSeqSVMfileIP.txt' )
        CVIterations -= 1

    global totalPrecision
    global totalRecall
    global totalF1Measure
    global totalAccuracy
    global totalCorelationCofficient

    print "Precision: "+str(totalPrecision/5)+"%"
    print "Recall: "+str(totalRecall/5)+"%"
    print "f1-Measure: "+str(totalF1Measure/5)+"%"
    print "Accuracy: "+str(totalAccuracy/5)+"%"
    print "Corelation Cofficient: "+str(totalCorelationCofficient/5)


if __name__ == "__main__":
    main()
