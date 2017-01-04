import os
import random
import operator
import matplotlib.pyplot as plt
import math
from sklearn import metrics

precision = 0
accuracy = 0
f1Measure = 0
corelationCofficient = 0
tpr = 0
fpr = 0
recall = 0

length = 4
cluster_files = ["protein_40.data.clstr", "protein_50.data.clstr", "protein_70.data.clstr", "protein_90.data.clstr", "protein_chain_40.data.clstr", "protein_chain_50.data.clstr", "protein_chain_70.data.clstr", "protein_chain_90.data.clstr"]
auxilary = {}
mupro_unique = {}
nFold=10

def read_mupro():
	no_of_lines = 0
	mupro = open("mupro.data")
	for line in mupro:
		no_of_lines = no_of_lines+1
		if(no_of_lines%10==1):
			unique_key = line.strip()
		elif(no_of_lines%10==9):
			line = line.split()
			unique_key = unique_key+line[1]+line[2]+line[3]+line[5]+line[6]
			if(float(line[7])<=0):
				unique_key = unique_key+'-'
			else:
				unique_key = unique_key+'+'

			if(not mupro_unique.has_key(unique_key)):
				mupro_unique[unique_key] = 1

def sigmod(value):
	return 1/(1+ math.exp(-value))



def create_auxilary_file(file_name):
	auxilary_file = open(file_name.split('.')[0]+"_auxilary.data", "w")
	mupro_file = open("mupro.data")
	no_of_lines = 0
	for line in mupro_file:
		no_of_lines = no_of_lines+1
		if(no_of_lines%10==1):
			auxilary_key = (line.strip())[0:length]
			unique_key = line.strip()
		elif(no_of_lines%10==9):
			line = line.split()
			unique_key = unique_key+line[1]+line[2]+line[3]+line[5]+line[6]
			if(float(line[7])<=0):
				unique_key = unique_key+'-'
			else:
				unique_key = unique_key+'+'

			if(unique_key in mupro_unique):
				del mupro_unique[unique_key]
				auxilary_file.write(str(auxilary[auxilary_key])+'\n')
	auxilary_file.close()

def nValid(mapFile,auxiliaryFile,dataFile, cluster_type, penalty_ratio):
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
      commandLearn="./svm_learn -z c -t 2 -g 0.1 -j %s -c 5 trainData.txt model_file > /dev/null 2>&1"%(str(penalty_ratio))
      os.system(commandLearn)


      commandClassify="./svm_classify testData.txt model_file output.txt > tempfile.txt"

      os.system(commandClassify)
      outputTemp = open('output.txt', 'r')

      for out in outputTemp:
          outputFile.write(out)
      outputTemp.close()
  testFinalFile.close()
  outputFile.close()
  analysis(cluster_type, penalty_ratio)

def create_auxilary():
	global length
	read_mupro()
	file_no = 0
	for cluster_file_name in cluster_files:
		file_no = file_no+1

		if(file_no>=5):
			length = 5

		read_mupro()
		cluster = -1
		cluster_file = open(cluster_file_name)
		for line in cluster_file:
			if line[0:1] == '>':
				cluster = cluster+1
			else:
				unique_key = line.split('>')[1][0:length]
				auxilary[unique_key] = cluster
		create_auxilary_file(cluster_file_name)
		for i in range(0,5):
			nValid(auxilary, cluster_file_name.split('.')[0]+"_auxilary.data", "sequence_features.txt" ,"sequence_"+cluster_file_name, 1)
		display_analysis()
		for i in range(0,5):
			nValid(auxilary, cluster_file_name.split('.')[0]+"_auxilary.data", "structure_features.txt" ,"structure_"+cluster_file_name, 2)
		display_analysis()
		for i in range(0,5):
			nValid(auxilary, cluster_file_name.split('.')[0]+"_auxilary.data", "structure_plus_sequence.txt" ,"structure+sequence_"+cluster_file_name, 2)
		display_analysis()
		auxilary.clear()	
		cluster_file.close()

def display_analysis():
	global accuracy
	global f1Measure
	global recall
	global tpr
	global fpr
	global corelationCofficient
	global precision


	print "\n\n"
	print "Precision: "+str(float(precision)/5)+"%"
	print "Recall: "+str(float(recall)/5)+"%"
	print "f1-Measure: "+str(float(f1Measure)/5)+"%"
	print "Accuracy: "+str(float(accuracy)/5)+"%"
	print "Corelation Cofficient: "+str(float(corelationCofficient)/5)
	print "True Positive Rate: "+str(float(tpr)/5)
	print "False Positive Rate: "+str(float(fpr)/5)

	precision = 0
	accuracy = 0
	f1Measure = 0
	corelationCofficient = 0
	tpr = 0
	fpr = 0
	recall = 0



def analysis(cluster_type, penalty_ratio):
	global accuracy
	global f1Measure
	global recall
	global tpr
	global fpr
	global corelationCofficient
	global precision

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

	precision += float(truePositive)/(truePositive+falsePositive)*100
	recall += truePositive/float(truePositive+falseNegative)*100
	f1Measure += float(2*truePositive)/(2*truePositive+falsePositive+falseNegative)*100
	accuracy += float(truePositive+trueNegative)/(truePositive+trueNegative+falseNegative+falsePositive)*100
	corelationCofficient += ((trueNegative*truePositive)-(falsePositive*falseNegative))/(math.sqrt((truePositive+falseNegative)*(truePositive+falsePositive)*(trueNegative+falseNegative)*(trueNegative+falsePositive)))
	tpr += recall
	fpr += falsePositive/float(trueNegative+falsePositive)*100

	plotROC(rocList, cluster_type)

increment=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]

def find_name(cluster_type):
	cluster_type = cluster_type.split('.')[0]
	cluster_type = cluster_type.split('_')
	type_c = cluster_type[0]+" "
	protein = cluster_type[1]+" "
	if("chain" in cluster_type[2]):
		protein += "chain "
		similarity = "Similarity "+cluster_type[3]+"%"
	else:
		similarity = "cluster similarity "+cluster_type[2]+"%"
	print type_c+protein+similarity
	return type_c, protein, similarity


def plotROC(rocList, cluster_type):
    type_c, protein, similarity = find_name(cluster_type)		
    sortedList = sorted(rocList, key=operator.itemgetter(0))
    #######################################
    trueClass=[]
    logSigValue=[]
    for data in sortedList:
        logSigValue.append(data[0])
        trueClass.append(data[1])
    fpr, tpr, thresholds = metrics.roc_curve(trueClass, logSigValue)
    roc_auc= metrics.auc(fpr, tpr)
 
    plt.figure()
    plt.plot(fpr, tpr,'g', label='AUC = %0.2f'% roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC for '+type_c+protein+similarity)
    plt.legend(loc="lower right")
    plt.savefig(cluster_type+'.png')
    plt.close()

    print "check"



create_auxilary()
read_mupro()
