from numpy import loadtxt
from sklearn.svm import SVC
from sklearn import multiclass
from sklearn.cross_validation import LeaveOneOut
import matplotlib.pylab as pl
import os,re,sys, gengraph 
import features
import numpy as np
from collections import defaultdict


def gensif(temp,f):
	infile2=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/pred.csv'
	fh1=open(infile2,'w')
	for k in temp:
		k=int(k)
		print >>fh1,str(k)
	fh1.close()
	gengraph.main(f)

features.main()
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
#clf=svm.SVC()
#clf.fit(tx,ty)/home/aditya/Project/dev
newtext_train='./GRN.py --a1-dir /home/aditya/Project/train --a2-dir /home/aditya/Project/train --pred-sif /home/aditya/Project/output_train.sif /home/aditya/Project/train/PMID-*.txt'
newtext_dev='./GRN.py --a1-dir /home/aditya/Project/dev --a2-dir /home/aditya/Project/dev --pred-sif /home/aditya/Project/output_dev.sif /home/aditya/Project/dev/PMID-*.txt'
print '....'
print 'Regulation:1, Activation:2, Inhibition:3, Requirement:4, Binding:5, Transcription:6'
print '\n'
print 'Training.......'

cvfile=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/final_cv_scores.txt'
cvfh=open(cvfile,'w')
loo=LeaveOneOut(len(ty))
tupvals=[]
tupsetting=[]
l1=[0.1+k*0.002 for k in range(2000)]
for myC in l1:
	for mygamma in [5]:
		final=[]
		estimator=SVC(C=myC, cache_size=200, coef0=0.0, degree=3, gamma=mygamma, kernel='linear',probability=True, scale_C=False, shrinking=True,tol=0.0001)
		tupsetting+=[(myC,mygamma)]
		clsfy=multiclass.OneVsRestClassifier(estimator)
		clsfy.fit(tx,ty)
		ty_dev_pred=clsfy.predict(tx_dev)
		temp=np.array([ty_dev_pred])
		final=np.append(final,temp)
		gensif(final,0)
		s=os.popen(newtext_dev)
		for k in s:
			print >>cvfh, k
			k= k.strip()
			ps=re.compile("^Slot")
			ress = ps.match(k)
			if ress:
				temp1=k.split()[3]
				tupvals.append(temp1) 
		s.close()
		tempval=map(str,[myC,mygamma])
		temptext=','.join([l for l in tempval])
		print str(temp1),tempval
		print >>cvfh,str(temp1),temptext
cvfh.close()



# Pick the optimal values from LOOC and test them on dev set
val=tupvals.index(min(tupvals))
myc,myGamma=tupsetting[val]
print myc,myGamma
estimator=SVC(C=myC, cache_size=200, coef0=0.0, gamma=mygamma, kernel='linear',probability=True, scale_C=False, shrinking=True,tol=0.0001)
clsfy=multiclass.OneVsRestClassifier(estimator)
clsfy.fit(tx,ty)
devpred=clsfy.predict(tx_dev)
gensif(devpred,0)
s=os.popen(newtext_dev)
for k in s:
	print k
print devpred
print ty_dev
pred_list={}
pred_dict=defaultdict(lambda:0)
true_dict=defaultdict(lambda:0)
for val1,val2 in zip(ty_dev,devpred):
	pred_dict[int(val2)]+=1
	true_dict[int(val1)]+=1
print pred_dict
print true_dict
print '1:Regulation, 2:Activation, 3:Inhibition, 4:Requirement, 5:Binding, 6:Transcription'
