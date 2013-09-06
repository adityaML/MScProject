import math
import sys
import os
import re
from collections import defaultdict
from itertools import combinations

def checkForPPI(fid,templist,kp,dirname,flag):
	filenameA1=dirname+'/'+fid+'.a1'
	slist=['Gene','GeneFamily','mRNA','Operon','Promoter','Regulon','Site']
	fh=open(filenameA1,'r')
	tempDict={}
	for line in fh:
		line=line.strip()
		pat=re.compile("^T\d")
		res = pat.match(line)
		if res:
			tlist=line.split()
			tempDict[tlist[0]]=tlist[1]
	if flag==0:
		tl1={}
		for key1 in templist:
			ln=len(templist[key1])
			if ln > 1:
				if templist[key1][1] in kp:
					tl1[templist[key1][1]]=key1
			else:
				if templist[key1][0] in kp:
					tl1[templist[key1][0]]=key1
		ret=0
		for key2 in tl1:
			if tempDict[tl1[key2]] in slist:
				ret=1
	elif flag==1:
		tl1={}
		for key1 in templist:
			if key1 in kp:
				tl1[key1]=templist[key1][1]
		ret=0
		for key2 in tl1:
			if tempDict[tl1[key2]] in slist:
				ret=1
	return ret


def keywords(text1):
	templist=[]
	dirname=os.path.dirname(os.path.abspath(__file__))
	fn=dirname+'/keywords.txt'
	fh=open(fn,'r')
	for line in fh:
		line=line.strip()
		templist.append([line,0])
	fh.close
	for l in range(len(templist)):
		if text1.find(templist[l][0]) != -1:
			templist[l][1]=1
	return templist


def parseTxt(fid,list1,dict1,dict2,dirname):
	# Extract its corresponding names
	tempdict={}
	for key in dict2:
		if len(dict2[key]) > 1:
			if dict2[key][1] in list1:
				tempdict[dict2[key][1]]=dict2[key][0]
	# Extract text
	filename=dirname+'/'+fid+'.txt'
	fh = open(filename,'r')
	templist=[]
	for line in fh:
		templine=line.strip()
		templist.append(templine)
	fh.close()
	linetext=templist[0]
	# Lookup these names in the text file
	tags=[]
	for k in list1:
		temp=linetext.find(tempdict[k])
		tags.append([temp,tempdict[k]])
	# Direction feature
	if len(tags) < 2:
		return 0
	if tags[0][0] > tags[1][0]:
		direction = 0
	else:
		direction = 1
	# Distance feature , direction , 3,6,9 words features
	tags.sort()
	midtext=linetext[tags[0][0]:tags[1][0]]
	words_midtext=midtext.split()
	mid_dist=len(words_midtext)
	tm3 , tm6 , tm9=0, 0, 0
	if mid_dist > 3:
		tm3 = 1
	if mid_dist > 6:
		tm6 = 1
	if mid_dist > 9:
		tm9 = 1
	templ1=[direction,tm3,tm6,tm9]
	# templ2=bowFeatures(tags,linetext)
	# build patterns
	templ2=keywords(linetext)
	newkeylist=[]
	for some in templ2:
		if some[1] ==1:
			newkeylist.append(some[0])
	if len(newkeylist) > 0:
		pat=[0,0,0]
		for term in newkeylist:
			k=[tags[0][0],tags[1][0]]
			if linetext.find(term) <= min(k):
				pat[0]+=1
			if min(k) <= linetext.find(term) <= max(k):
				pat[1]+=1
			if linetext.find(term) >= max(k):
				pat[2]+=1
		val=pat.index(max(pat))
		pat=[0,0,0]
		pat[val]=1
		patt=','.join([str(vv) for vv in pat])
	else:
		patt='0,0,0'
	#keywords
	keywordsf= ','.join([str(t1[1]) for t1 in templ2])
	justkey=keywordsf
	justdirkey=str(direction)+','+keywordsf
	tempf=','.join([str(t) for t in templ1])
	tempret=tempf+','+keywordsf
	return tempret+','+patt


def parseA1Test(filename,dirname):
	fileNameA1=dirname+'/'+filename+'.a1'
	tempDict=defaultdict(lambda:'')
	parseDict={}
	fh=open(fileNameA1,'r+')
	for line in fh:
		line=line.strip()
		linelist=line.split()
		tempDict[linelist[0]]=line
	for val in tempDict:
		p4=re.compile("^N\d")
		res = p4.match(tempDict[val])
		if res:
			t1= tempDict[val].split()[2].split(':')[-1]
			t2= tempDict[val].split()[3].split(':')[-1]
			t3= tempDict[t1].split()[2]
			parseDict[t2]=[t3,t1]
	#strval= ''.join([str(i) for i in range(num)])
	#feature1str= ','.join([str(i) for i in featureEnt])
	#itr=combinations(strval,2)
	fh.close()
	return parseDict


def parseA1(fileName,dirname):
	#import pdb;pdb.set_trace()
	entity=['Gene','GeneFamily','mRNA','Operon','PolymeraseComplex','Promoter',
	'Protein','ProteinComplex','ProteinFamily','Regulon','Site']
	entityDict=defaultdict(lambda:0)
	fileNameA1=dirname+'/'+fileName+'.a1'
	tempDict=defaultdict(lambda:'')
	parseDict={}
	fh=open(fileNameA1,'r+')
	for line in fh:
		line=line.strip()
		linelist=line.split()
		tempDict[linelist[0]]=line
		p1=re.compile("^T\d")
		res1 = p1.match(line)
		if res1:
			entityDict[linelist[1]]+=0.01
	fh.close()
	for val in tempDict:
		p4=re.compile("^N\d")
		res = p4.match(tempDict[val])
		if res:
			t1= tempDict[val].split()[2].split(':')[-1]
			t2= tempDict[val].split()[3].split(':')[-1]
			t3= tempDict[t1].split()[2]
			parseDict[t2]=t3
	num= len(parseDict)
	featureEnt=[]
	for val in entity:
		if val in entityDict:
			featureEnt+=[entityDict[val]]
		else:
			featureEnt+=[0]
	strval= ''.join([str(i) for i in range(num)])
	feature1str= ','.join([str(i) for i in featureEnt])
	itr=combinations(strval,2)
	return [parseDict,feature1str]


def generateDictTrain(dirName):
	file_list=os.listdir(dirName)
	finalDict={}
	for l in file_list:
		lwhich=l.split('.')
		m1=l.find('.a1')
		if m1!=-1:
			fileName=l.split('.')[0]
			fileNameA1=dirName+'/'+fileName+'.a1'
			fileNameA2=dirName+'/'+fileName+'.a2'
			fileNameText=dirName+'/'+fileName+'.txt'
			fhtxt=open(fileNameText,'r+')
			temp=processA2File(fileNameA2,fileNameA1)
			finalDict[fileName]=temp
	return finalDict


def generateDictTest(dirName):
	file_list=os.listdir(dirName)
	finalDict={}
	for l in file_list:
		lwhich=l.split('.')
		m1=l.find('.a1')
		if m1!=-1:
			fileName=l.split('.')[0]
			fileNameA1=dirName+'/'+fileName+'.a1'
			fileNameText=dirName+'/'+fileName+'.txt'
			temp1=parseA1Test(fileName,dirName)
			temp2=parseA1(fileName,dirName)
			finalDict[fileName]=[temp1,temp2]
	return finalDict

# Generate feature vector from *.txt  *.a1 and *.a2 files
def generateFeatureTrain(dirname,flag):
	#import pdb;pdb.set_trace()
	trYlabels={'Regulation':1,'Activation':2,'Inhibition':3,'Requirement':4,
	'Binding':5,'Transcription':6,'Other':7}
	finalDict=generateDictTrain(dirname)
	if flag:
		fname=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/trainNLP.csv'
		fname1=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/dataGraphTrain.csv'
	else:
		fname=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/devNLP.csv'
		fname1=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/dataGraphDev.csv'
	trf = open(fname, 'w')
	trf1 = open(fname1, 'w')
	for key in finalDict:
		featureID=key
		tempDict=finalDict[key]
		for val in tempDict:
			templist= tempDict[val]
			# feature set 1 extracted from .a1 file by parseA1 function
			fset1=parseA1(featureID,dirname)
			for i in range(len(tempDict[val])/3):
				kp= [templist[3*i+0],templist[3*i+1]]
				# feature set 2 from parseTxt function
				fset2=parseTxt(featureID,kp,fset1[0],templist[3*i+2],dirname)
				nl=[eval(fset1[0][sval])/1000.0 for sval in kp]
				tempstr=','.join(map(str, nl))
				ret=checkForPPI(featureID,templist[3*i+2],kp,dirname,0)
				if fset2 != 0:
					if ret !=0:
						#finaltext= fset2+','+str(trYlabels[val])
						finaltext= tempstr+','+fset2+','+str(trYlabels[val])
						finaltext= tempstr+','+fset1[1]+','+fset2+','+str(trYlabels[val])
					else:
						#finaltext= fset2+','+str(trYlabels[val])
						#finaltext= tempstr+','+fset2+','+str(trYlabels[val])
						finaltext= tempstr+','+fset1[1]+','+fset2+','+str(trYlabels[val])
					print >>trf,finaltext
					str1=','.join([l for l in kp])
					print >>trf1,str1
	trf.close()
	trf1.close()

# Generate feature vector from *.a1 and .txt files
def generateFeatureTest(dirname):
	finalDict= generateDictTest(dirname)
	fname=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/testNLP.csv'
	fname1=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/dataGraph.csv'
	trf = open(fname, 'w')
	trf1 = open(fname1,'w')
	for key in finalDict:
		tempdict=finalDict[key][0]
		fileNameA1=dirname+'/'+key+'.a1'
		templist=[]
		for l in tempdict:
			p5=re.compile("[a-z]+[A-Z0-9]+")
			if p5.match(l):
				temptx=':'.join([str(s) for s in tempdict[l]])
				templist.append(temptx)
		num=len(templist)
		if num > 1:
			strval= ''.join([str(i) for i in range(num)])
			itr=combinations(strval,2)
			for tup in itr:
				nlist=[]
				for l in tup:
					nlist.append(l)
				#nl=[nlist]
				nl=[nlist,[nlist[1],nlist[0]]]
				for tup in nl:
					f1=templist[int(tup[0])].split(':')[-1]
					f2=templist[int(tup[1])].split(':')[-1]
					dict2= extractNames(f1,f2,fileNameA1)
					list1=[dict2[f1][1],dict2[f2][1]]
					dict1= finalDict[key][1][0]
					nl=[eval(dict1[k])/1000.0 for k in list1]
					t1=','.join(map(str, nl))
					t2 = finalDict[key][1][1]
					t3= parseTxt(key,list1,dict1,dict2,dirname)
					ret=checkForPPI(key,tempdict,list1,dirname,1)
					if t3 != 0:
						#finaltext=t3
						#finaltext=t1+','+t3
						finaltext=t1+','+t2+','+t3
						print >>trf,finaltext
						str1=','.join([l for l in list1])
						print >>trf1,str1
	trf.close()
	trf1.close()



# Extract entities and features from *.a1 files to be used for generating labels
def extractNames(f1,f2,fileNameA1):
	nameDict=defaultdict(list)
	fha1=open(fileNameA1,'r+')
	for line in fha1:
		line=line.strip()
		matchObj = re.match( r'^N\d', line)
		tempListA1=line.split('\t')
		if tempListA1[0]== f1:
			nameDict[f1].append(tempListA1[2])
		if tempListA1[0]== f2:
			nameDict[f2].append(tempListA1[2])
		if matchObj:
			tempRef=tempListA1[1].split(':')[1].split()[0]
			tempGeneName=tempListA1[1].split(':')[-1]
			nameDict[tempRef].append(tempGeneName)
	fha1.close()
	return nameDict


# Extract relations from *.a2 files.
def processA2File(fileNameA2,fileNameA1):
	tempDict=defaultdict(lambda:'')
	finalDict=defaultdict(lambda:[])
	fha2=open(fileNameA2,'r+')
	for line in fha2:
		line=line.strip()
		p1=re.compile("^R\d")
		p2=re.compile("^E\d")
		p3=re.compile("^T\d")
		p4=re.compile("^Interaction")
		newlist1=line.split('\t')
		tempDict[newlist1[0]]=newlist1[1:-1]+[newlist1[-1]]
	tempRDict,tempTDict,tempEDict,tempRDict1={},{},{},{}
	for key in tempDict:
		if p1.match(key):
			tempRDict[key]=tempDict[key]
			tempR=tempDict[key][0].split()
			if p4.match(tempR[0]):
				tempRAgent=tempR[1].split(':')[1]
				tempRTarget=tempR[2].split(':')[1]
				tempRval=[tempRAgent,tempRTarget]
				tempRDict1[key]=tempRval
		if p2.match(key):
			tempE=tempDict[key][0].split()
			tempEval=tempE[1].split(':')[1]
			tempEDict[key]=tempEval
		if p3.match(key):
			tempTDict[key]=tempDict[key]
	for r1 in tempRDict:
		[tempRel,tempTarget,tempAgent] = tempRDict[r1][0].split()
		AgentType=tempAgent.split(':')[1]
		TargetType=tempTarget.split(':')[1]
		tempAgent=AgentType
		tempTarget=TargetType
		if p4.match(tempRel):
			while True:
				if p2.match(AgentType):
					AgentType=tempEDict[AgentType]
				if p1.match(AgentType):
					AgentType=tempRDict1[AgentType][0]
				if p2.match(TargetType):
					TargetType=tempEDict[TargetType]
				if p1.match(TargetType):
					TargetType=tempRDict1[TargetType][1]
				if p3.match(AgentType) and p3.match(TargetType):
					break
			k=extractNames(AgentType,TargetType,fileNameA1)
			if len(k[AgentType]) > 1 and len(k[TargetType]) >1:
				finalDict[tempRel.split('.')[1]]+=[k[AgentType][1] , 
				k[TargetType][1],k]
	fha2.close()
	return finalDict#print >>trf, ','.join(map(str, features))

def main():
	# Example file
	#dirnametrain='/Users/adityajitta/Desktop/BioNLPS/sample'
	dirnametrain=os.path.dirname(os.path.abspath(__file__))+'/train'
	dirnamedev=os.path.dirname(os.path.abspath(__file__))+'/dev'
	dirnametest=os.path.dirname(os.path.abspath(__file__))+'/test'
	# Generate training feature vector
	generateFeatureTrain(dirnametrain,1)
	# Generate dev set vector
	generateFeatureTrain(dirnamedev,0)
	# Generate Features for testing
	generateFeatureTest(dirnametest)
	#print parseTxt()

if __name__ == '__main__':
	main()




