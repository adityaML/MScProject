import os,sys,re,math
from collections import defaultdict
from itertools import combinations

def genpredfile(infile1,infile2):
	fi1=open(infile1,'r')
	fi2=open(infile2,'w')
	for line in fi1:
		line=line.strip()
		print >>fi2,str(line.split(',')[-1])
	fi2.close()
	fi1.close()


def gengraphfile(outfile,infile,infile2):
	k={1:'[color="0.325 0.850 0.850"]',2:'[color="0.647 0.702 0.702"]',3:'[color="0.002 0.999 0.999"]',
	4:'[color="0.499 0.386 1.000"]',5:'[color="0.647 0.702 0.702"]',6:'[color="0.201 0.753 1.000"]'}
	inittext=['digraph prof {','size="6,4";','ratio = fill;','node [style=filled];']
	fo=open(outfile,'w')
	fi1=open(infile,'r')
	fi2=open(infile2,'r')
	tempDict={}
	for l1 in inittext:
		print >>fo,l1
	templist=[]
	for line1,line2 in zip(fi1,fi2):
		temp1=line1.strip()
		temp2=line2.strip()
		[a1,a2]=temp1.split(',')
		for t1 in [a1,a2]:
			tempDict[a1]='[color="0.650 0.200 1.000"]'
			tempDict[a2]='[color="0.650 0.200 1.000"]'
		a1='"'+a1+'"'
		a2='"'+a2+'"'
		temptext=a1+' '+'->'+' '+a2+' '+k[int(temp2)]+';'
		print >>fo,temptext
	for key in tempDict:
		tempstr= '"'+key+'"'+' '+tempDict[key]+';'
		print >>fo,tempstr
	print >>fo,'}'
	fo.close()
	fi1.close()
	fi2.close()

def gengraphsif(siffile,outfile):
	yl={'Regulation':1,'Activation':2,'Inhibition':3,'Requirement':4,
	'Binding':5,'Transcription':6}
	k={1:'[color="0.325 0.850 0.850"]',2:'[color="0.647 0.702 0.702"]',3:'[color="0.002 0.999 0.999"]',
	4:'[color="0.499 0.386 1.000"]',5:'[color="0.647 0.702 0.702"]',6:'[color="0.201 0.753 1.000"]'}
	inittext=['digraph prof {','size="6,4";','ratio = fill;','node [style=filled];']
	fi=open(siffile,'r')
	fo=open(outfile,'w')
	tempDict={}
	for l1 in inittext:
		print >>fo,l1
	for line in fi:
		line=line.strip()
		llist=line.split()
		[a1,a2]=[llist[0],llist[2]]
		temp2=llist[1].split('.')[1]
		for t1 in [a1,a2]:
			tempDict[a1]='[color="0.650 0.200 1.000"]'
			tempDict[a2]='[color="0.650 0.200 1.000"]'
		a1='"'+a1+'"'
		a2='"'+a2+'"'
		temptext=a1+' '+'->'+' '+a2+' '+k[yl[temp2]]+';'
		print >>fo,temptext
	for key in tempDict:
		tempstr= '"'+key+'"'+' '+tempDict[key]+';'
		print >>fo,tempstr
	print >>fo,'}'
	fo.close()
	fi.close()



if __name__ == '__main__':
	infile1=os.path.dirname(os.path.abspath(__file__))+'/trainNLP.csv'
	infile2=os.path.dirname(os.path.abspath(__file__))+'/pred.csv'
	infile=os.path.dirname(os.path.abspath(__file__))+'/dataGraphDev.csv'
	siffile=os.path.dirname(os.path.abspath(__file__))+'/output_thresh.sif'
	#genpredfile(infile1,infile2)
	outfile=os.path.dirname(os.path.abspath(__file__))+'/finalgraph2.gv'
	#gengraphfile(outfile,infile,infile2)
	gengraphsif(siffile,outfile)