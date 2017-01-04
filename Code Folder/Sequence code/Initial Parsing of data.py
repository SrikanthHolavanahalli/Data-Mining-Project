from asynchat import async_chat

pathForBaseFile='E:\Semester 2\Data Mining\Project\Stability\Stability\mupro.data'
allLengths=[]



def readFile(path):
    f = open(path, 'r')
    required=[3,9]
    count=0
    toSave=[]
    result = open('myfile.txt','w')
    features140 = open('features140.txt','w')
    totalCOunt={}
    extra=0
    duplicateRemoval={}
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

        if count in required:
            toSave.append(row)


        if count is 9:
            checking=str(toSave[0])
            for i in range(0,len(checking)):
                if totalCOunt.has_key((checking[i])):
                        totalCOunt[checking[i]]=totalCOunt.get(checking[i])+1
                else:
                    totalCOunt.update({checking[i]:1})
            length=len(totalCOunt)
            if length is 21:
                extra+=1
                #print "catch"
            if length not in allLengths:
                allLengths.append(length)
            totalCOunt={}

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
            countAminoAcids={}
            for i in range(index-3,index+4):
                #print ord(toSave[0][i])-ord('A')+1
                #print toSave[0][i]
                if i is index:
                    continue
                if i >len(toSave[0])-1:
                    break
                if countAminoAcids.has_key((ord(toSave[0][i]))-ord('A')+1):
                    countAminoAcids[(ord(toSave[0][i]))-ord('A')+1]=(countAminoAcids[(ord(toSave[0][i]))-ord('A')+1])+1
                else:
                    countAminoAcids.update({(ord(toSave[0][i]))-ord('A')+1:1})
            count=0
            #energy=int(words[7].rstrip())
            if float(words[7].rstrip())>0:
                features140.write("1")
                result.write("1")
            else:
                features140.write("0")
                result.write("0")
            ascending=[]
            for i in countAminoAcids:
                ascending.append(i)
            ascending.sort()

            for i in range(0,len(countAminoAcids)):
                result.write(" "+str(ascending[i])+":"+str(countAminoAcids.get(ascending[i])))

            result.write(" 27:"+str(words[5])+" 28:"+str(words[6]+"\n"))


            wildType=ord(words[1])-ord('A')+1
            wildType=6*20+wildType
            mutation=ord(words[3])-ord('A')+1

            mutation=6*26+mutation
            for i in range(0,len(countAminoAcids)):
                features140.write(" "+str(ascending[i]+26*i)+":"+str(countAminoAcids.get(ascending[i])))
            if wildType < mutation:
                features140.write(" "+str(wildType)+":-1 "+str(mutation)+":1")
            else:
                features140.write(" "+str(mutation)+":1 "+str(wildType)+":-1")
            features140.write("\n")
        ###
        if count is 100:
            break
        ###
    print extra

    result.close()
    features140.close()



readFile(pathForBaseFile)