#!/usr/bin/python

from numpy import loadtxt
from sklearn.svm import SVC
from sklearn import svm
from sklearn import multiclass
from sklearn.cross_validation import LeaveOneOut
import matplotlib.pylab as pl
import os,re,sys
import numpy as np
from gurobipy import *

def mipsMain(thresh):
	tex=loadtxt('/home/aditya/Project/csvFiles/trainNLP.csv',delimiter=',')
	tex_dev=loadtxt('/home/aditya/Project/csvFiles/devNLP.csv',delimiter=',')
	tex_test=loadtxt('/home/aditya/Project/csvFiles/testNLP.csv',delimiter=',')
	vals= tex.shape
	vals_dev= tex_dev.shape
	vals_test= tex_test.shape
	ty=tex[:,vals[1]-1]
	tx=tex[:,0:vals[1]-2]
	ty_dev=tex_dev[:,vals_dev[1]-1]
	tx_dev=tex_dev[:,0:vals_dev[1]-2]
	tx_test=tex_test[:,0:vals_test[1]-1]
	clf=SVC(C=5.0, cache_size=200, coef0=0.0, degree=3, gamma=5.0, kernel='linear',probability=True, scale_C=False, shrinking=True,tol=0.0001)
	clf.fit(tx,ty)
	#prev_pred=clf.predict(tx_test)
	#temp = clf.predict_proba(tx_test)
	prev_pred=clf.predict(tx_dev)
	temp = clf.predict_proba(tx_dev)
	(ninst,nclass)=temp.shape
	obj_coeff=temp.tolist()
	# Mixed Integer Programming
	m=Model("mips")
	# declare Variables
	vars={}
	for i in range(ninst):
		for j in range(nclass):
			vars[i,j] = m.addVar(obj=-1*obj_coeff[i][j], vtype=GRB.BINARY,name='x'+str(i)+'_'+str(j))
	m.update()
	# Add Constraints
	# everything is a Regulation:1,
	# Activation:2 and Inhibition:3 do not occur together
	# Requirement:4 participates in Activation and Inhibition
	# Binding:5 may or may not coexist Transcription:6
	for i in range(ninst):
		num_labels=sum(np.array([obj_coeff[i][jo] for jo in range(nclass)])> thresh)
		# constraint on number of labels predicted
		m.addConstr(quicksum(vars[i,j] for j in range(nclass)) == num_labels)
		# making sure Regulation always exist
		# m.addConstr(vars[i,0]==1)
		# making sure Activation and Inhibition do not occur together
		m.addConstr(vars[i,1]+vars[i,2] <= 1)
	m.update()
	# Optimize Model
	m._vars = vars
	m.optimize()
	final=[]
	for i in range(ninst):
		final.append([float(j+1) for j in range(nclass) if m._vars[i,j].x == 1.0])
	return [final,prev_pred]
	# Get Results

def genpredfile(temp,infile2):
	fh1=open(infile2,'w')
	for k in temp:
		if type(k) is list:
			for v in k:
				v=int(v)
				print >>fh1,str(v)
		else:
			k=int(k)
			print >>fh1,str(k)
	fh1.close()

if __name__ == '__main__':
	klab={1:'Regulation',2:'Activation',3:'Inhibition',4:'Requirement',
	5:'Binding',6:'Transcription',7:'Regulation'}
	thresh=0.3
	[k,pred]=mipsMain(thresh)
	outfile=os.path.dirname(os.path.abspath(__file__))+'/output_thresh.sif'
	outfiledev=os.path.dirname(os.path.abspath(__file__))+'/outputdev.sif'
	infile=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/dataGraphDev.csv'
	#infile=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/dataGraph.csv'
	infile1=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/devNLP.csv'
	infile2=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/pred.csv'
	genpredfile(pred,infile2)
	#genpredfile(k,infile2)
	fh=open(infile,'r')
	fo=open(outfile,'w')
	fo1=open(outfiledev,'w')
	for l,v in zip(fh,k):
		l=l.strip()
		[a1,a2]=l.split(',')
		for vs in v:
			temptext1=a1+'\t'+'Interaction.'+klab[int(vs)]+'\t'+a2
			print >>fo,temptext1
	fh.close()
	fo.close()
	fi3=open(infile2,'r')
	fi4=open(infile,'r')
	print k
	print pred
	for l1,l2 in zip(fi4,fi3):
		temp1=l1.strip()
		temp2=l2.strip()
		[a1,a2]=temp1.split(',')
		temptext1=a1+'\t'+'Interaction.'+klab[int(temp2)]+'\t'+a2
		print >>fo1,temptext1
	fo1.close()
	fi3.close()
	fi4.close()
	newtext='./GRN.py --a1-dir /home/aditya/Project/dev --a2-dir /home/aditya/Project/dev --pred-sif /home/aditya/Project/output_thresh.sif /home/aditya/Project/dev/PMID-*.txt'
	stream=os.popen(newtext)
	for lal in stream:
		lal=lal.strip()
		print lal
	












