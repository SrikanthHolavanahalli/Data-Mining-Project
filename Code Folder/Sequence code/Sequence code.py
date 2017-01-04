import random
#import svmlight
from asynchat import async_chat
import os
import math


pathForBaseFile='E:\Semester 2\Data Mining\Project\Stability\Stability\\mupro.data'
pathForResultFile="C:\Users\Vikrant Kaushal\PycharmProjects\DM_Project\\tempfile.txt"
allLengths=[]
proteinIdDict={}
proteinChainIdDict={}
nFold=10

def readFile(path):
    f = open(path, 'r')
    required=[3,9]
    count=0
    toSave=[]
    #result = open('myfile.txt','w')
    features182 = open('features182.txt','w')
    totalCOunt={}
    global proteinChainIdDict
    global proteinIdDict
    extra=0
    duplicateRemoval={}
    proteinChainID=1
    proteinID=1

    for data in f:



        duplicateRemovalList=[]
        if data in '\n':
            toSave=[]
            continue
        count+=1
        #row=data.split('\t').rstrip()
        row=data.rstrip()
        if count==1:
            pdbCode=data.rstrip()
            if not proteinChainIdDict.has_key(pdbCode):
                proteinChainIdDict.update({pdbCode:proteinChainID})
                proteinChainID+=1


        if count in required:
            toSave.append(row)


        if count is 9:


            words=data.split(" ")
            if not proteinIdDict.has_key(words[0]):
                proteinIdDict.update({words[0]:proteinID})
                proteinID+=1
            sign="+"
            if float(words[7].rstrip())<=0:
                sign="-"


            string=words[1]+words[2]+words[3]+words[5]+words[6]+sign
            key=pdbCode+string
            if duplicateRemoval.has_key(key):
                count=0
                continue
            else:
                duplicateRemoval.update({key:1})
            if float(words[7].rstrip())<=0:
                features182.write("-1")
                sign="-"
            else:
                features182.write("1")
            index=int(words[2])-1
            countAminoAcids={}
            multiplication=0

            for i in range(index-3,index+4):
                #print ord(toSave[0][i])-ord('A')+1
                #print toSave[0][i]
                if i == index or i < 0:
                    #multiplication+=1
                    continue
                if i > len(toSave[0])-1 :
                    break

                features182.write(" "+str(((ord(toSave[0][i]))-ord('A')+1)+(26*multiplication))+":"+str(1))
                multiplication+=1





            count=0


            wildType=ord(words[1])-ord('A')+1
            wildType=6*26+wildType
            mutation=ord(words[3])-ord('A')+1

            mutation=6*26+mutation

            if wildType < mutation:
                features182.write(" "+str(wildType)+":-1 "+str(mutation)+":1")
            else:
                features182.write(" "+str(mutation)+":1 "+str(wildType)+":-1")

            if float(words[7].rstrip())<=0:
                features182.write(" #-1")
            else:
                features182.write(" #1")
            features182.write("\n")



    #result.close()
    features182.close()

def createMap(path):
    f = open(path, 'r')
    required=[3,9]
    count=0
    global proteinChainIdDict
    global proteinIdDict
    extra=0
    toSave=[]
    proteinChainID=1
    proteinID=1
    for data in f:
        if data in '\n':
            toSave=[]
            continue
        count+=1
        row=data.rstrip()
        if count==1:
            pdbCode=data.rstrip()
            if not proteinChainIdDict.has_key(pdbCode):
                proteinChainIdDict.update({pdbCode:proteinChainID})
                proteinChainID+=1
        if count in required:
            toSave.append(row)
        if count is 9:
            words=data.split(" ")
            if not proteinIdDict.has_key(words[0]):
                proteinIdDict.update({words[0]:proteinID})
                proteinID+=1
            count=0



def auxiliaryFile(path):
    f = open(path, 'r')
    required=[3,9]
    count=0
    toSave=[]

    proteinIDFile = open('ProteinId.txt','w')
    proteinChainIDFile = open('proteinChainId.txt','w')

    global proteinChainIdDict
    global proteinIdDict
    extra=0
    duplicateRemoval={}

    for data in f:
        if data in '\n':
            toSave=[]
            continue
        count+=1
        #row=data.split('\t').rstrip()
        row=data.rstrip()
        if count==1:
            pdbCode=data.rstrip()


        if count in required:
            toSave.append(row)


        if count is 9:

            words=data.split(" ")
            sign="+"
            if float(words[7].rstrip())<=0:
                sign="-"
            string=words[1]+words[2]+words[3]+words[5]+words[6]+sign
            key=pdbCode+string
            if duplicateRemoval.has_key(key):
                count=0
                continue
            else:
                duplicateRemoval.update({key:1})
            index=int(words[2])-1

            proteinChainIDFile.write(str(proteinChainIdDict.get(pdbCode))+"\n")
            proteinIDFile.write(str(proteinIdDict.get(words[0]))+"\n")

            count=0


    proteinIDFile.close()
    proteinChainIDFile.close()


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

        #commandLearn="svm_learn.exe -z c -t 0 trainData.txt model_file"
        #for all parameters considered
        commandLearn="svm_learn.exe -z c -t 2 -g 0.1 -j 1 -c 5 trainData.txt model_file"
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
    rocList=[]
    testList=[]
    for data in range(0,len(length)-1):
        output=outputFinalFile.readline()
        classLabel=int(testFinalFile.readline().split('#')[1].strip())
        prediction=-1
        output=sigmod(float(output))
        if output >= 0.5:
            prediction=+1
        rocList.append([float(output),int(classLabel)])
        testList.append([int(prediction),int(classLabel)])


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
    predictionList=[]
    trueList=[]
    for data in testList:
        predictionList.append(data[0])
        trueList.append(data[1])
    print metrics.accuracy_score(trueList,predictionList)
    print "\n\n"
    print "Precision: "+str(precision)+"%"
    print "Recall: "+str(recall)+"%"
    print "f1-Measure: "+str(f1Measure)+"%"
    print "Accuracy: "+str(accuracy)+"%"
    print "Corelation Cofficient: "+str(corelationCofficient)
    #print "True Positive Rate: "+str(tpr)
    #print "False Positive Rate: "+str(fpr)

    plotROC(rocList)



    print "check"

import operator
import matplotlib.pyplot as plt
from sklearn import metrics



#increment=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]

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
    #print fpr
    #print tpr
    #print thresholds
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
    plt.plot(fpr, tpr,'b', label='AUC = %0.2f'% roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic for Sequence only and per Protein chain')
    plt.legend(loc="lower right")
    plt.show()

    print "check"








def readResult(path):
    f = open(path, 'r')

    totalAverage=0
    for data in f:
        if data[0:8]== "Accuracy":
            string=""
            for i in data[21:29]:
                if i == "%":
                    break;
                if i.isdigit() or i=='.':
                    string+=i

            totalAverage+=float(string)
    if method==1:
        string="proteinID"
    else:
        string="protein + chainID"
    nameofFile=string+"_result.txt"
    result = open(nameofFile,'w')
    result.write(str(totalAverage/nFold))
    result.close()

#createMap(pathForBaseFile)
readFile(pathForBaseFile)
auxiliaryFile(pathForBaseFile)

pathForAuxiliaryProteinIDFile="C:\Users\Vikrant Kaushal\PycharmProjects\DM_Project\\ProteinId.txt"
pathForAuxiliaryProteinChainIDFile="C:\Users\Vikrant Kaushal\PycharmProjects\DM_Project\\proteinChainId.txt"

pathForDataFile="C:\Users\Vikrant Kaushal\PycharmProjects\DM_Project\\features182.txt"
method=int(raw_input("1 for proteinID\n2 for proteinID+ chaninID"))
if method==1:
    nValid(proteinIdDict,pathForAuxiliaryProteinIDFile,pathForDataFile)
else:
    nValid(proteinChainIdDict,pathForAuxiliaryProteinChainIDFile,pathForDataFile)
readResult(pathForResultFile)
print "done"