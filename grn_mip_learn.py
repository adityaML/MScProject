#!/usr/bin/python

from numpy import loadtxt
from sklearn.svm import SVC
from sklearn.cross_validation import LeaveOneOut
import matplotlib.pylab as pl
import os,re,sys
import numpy as np
from gurobipy import *
def mipsMain(thresh):
    print thresh
    tex=loadtxt('/home/aditya/Project/csvFiles/trainNLP.csv',delimiter=',')
    tex_dev=loadtxt('/home/aditya/Project/csvFiles/devNLP.csv',delimiter=',')
    vals= tex.shape
    vals_dev= tex_dev.shape
    ty=tex[:,vals[1]-1]
    tx=tex[:,0:vals[1]-2]
    ty_dev=tex_dev[:,vals_dev[1]-1]
    tx_dev=tex_dev[:,0:vals_dev[1]-2]
    clf=SVC(C=4.0, cache_size=200, coef0=0.0, degree=3, gamma=9.0, kernel='linear',probability=True, scale_C=False, shrinking=True,tol=0.0001)
    clf.fit(tx,ty)
    #prev_pred=clf.predict(tx_dev[:,0:vals[1]-2])
    #temp = clf.predict_proba(tx_dev[:,0:vals[1]-2])
    prev_pred=clf.predict(tx[:,0:vals[1]-2])
    temp = clf.predict_proba(tx[:,0:vals[1]-2])
    (ninst,nclass)=temp.shape
    print temp.shape
    obj_coeff=temp.tolist()
    # Mixed Integer Programming
    m=Model()
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
        # constraint on number of labels predicted
        num_labels=sum(np.array([obj_coeff[i][jo] for jo in range(nclass)])> thresh)
        if num_labels > 0:
			# making sure Regulation always exist
            # m.addConstr(vars[i,0]==1)
            m.addConstr(quicksum(vars[i,j] for j in range(nclass)) == num_labels)
            # making sure Activation and Inhibition do not occur together
            m.addConstr(vars[i,1]+vars[i,2] <= 1)
            m.addConstr((vars[i,1]>0 and vars[i,0] >=0) or (vars[i,1] <= 0) )
            # adding a constraint restricting the co-existance of Inhibition and Transcription
            # m.addConstr(vars[i,2]+vars[i,5] <= 1)
            # m.addConstr(vars[i,1]+vars[i,3] == 2)
    m.update()
    # Optimize Model
    m._vars = vars
    m.optimize()
    final=[]
    for i in range(ninst):
        final.append([float(j+1) for j in range(nclass) if m._vars[i,j].x == 1.0])
    return [final,prev_pred]



def genpredfile(temp,infile2):
    fh1=open(infile2,'w')
    for k in temp:
        k=int(k)
        print >>fh1,str(k)
    fh1.close()


def mips_run(thresh):
    klab={1:'Regulation',2:'Activation',3:'Inhibition',4:'Requirement',5:'Binding',6:'Transcription',7:'Regulation'}
    [k,pred]=mipsMain(thresh)
    print k
    outfile=os.path.dirname(os.path.abspath(__file__))+'/output.sif'
    infile=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/dataGraphTrain.csv'
    infile1=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/devNLP.csv'
    infile2=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/pred.csv'
    genpredfile(pred,infile2)
    fh=open(infile,'r')
    fo=open(outfile,'w')
    for l,v in zip(fh,k):
        l=l.strip()
        [a1,a2]=l.split(',')
        for vs in v:
            temptext1=a1+'\t'+'Interaction.'+klab[int(vs)]+'\t'+a2
            print >>fo,temptext1
    fh.close()
    fo.close()

