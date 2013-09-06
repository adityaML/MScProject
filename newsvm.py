from numpy import loadtxt
#from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn import multiclass
from sklearn import preprocessing 
import matplotlib.pylab as pl
import os
import gengraph
import features
import re

def gensif(temp):
	infile2=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/pred.csv'
	fh1=open(infile2,'w')
	for k in temp:
		k=int(k)
		print >>fh1,str(k)
	fh1.close()
	gengraph.main()

lb=preprocessing.L
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
newtext='./GRN.py --a1-dir /home/aditya/Project/dev --a2-dir /home/aditya/Project/dev --pred-sif /home/aditya/Project/output.sif /home/aditya/Project/dev/PMID-*.txt'
print '....'
print '\n'
print 'Regulation:1, Activation:2, Inhibition:3, Requirement:4, Binding:5, Transcription:6'
print '\n'
print 'Training.......'
#clsfy=OneVsRestClassifier(SVC(kernel='linear'))
tempdict={}
for myC in [0.001,0.01,0.1,0.5,1.0,5.0,10.0]:
	for mygamma in [0.0,0.001,0.01,0.1,1.0,5.0,10.0]:
		for mydegree in [2,3,4,5]:
			estimator=SVC(C=myC, cache_size=200, coef0=0.0, degree=mydegree, gamma=mygamma, 
			kernel='rbf',probability=True, scale_C=False, shrinking=True,tol=0.0001)
			clsfy=multiclass.OneVsRestClassifier(estimator)
			clsfy.fit(tx,ty)
			estimator.fit(tx,ty)
			ty_dev_pred=clsfy.predict(tx_dev)
			sometup=(myC,mygamma,mydegree)
			gensif(ty_dev_pred)
			stream=os.popen(newtext)
			templis=[]
			for lal in stream:
				lal= lal.strip()
				pm=re.compile("^M")
				resm = pm.match(lal)
				ps=re.compile("^Slot")
				ress = ps.match(lal)
				if resm:
					templis.append(lal.split()[1])
				if ress:
					templis.append(lal.split()[3])
			stream.close()
			tempdict[sometup]=templis
			print sometup,templis
			# Output the SER score
estimator=SVC(C=1.0, cache_size=200, coef0=0.0, degree=3, gamma=1.0, kernel='linear',probability=True, scale_C=False, shrinking=True,tol=0.0001)
clsfy=multiclass.OneVsRestClassifier(estimator)
clsfy.fit(tx,ty)
estimator.fit(tx,ty)
ty_dev_pred=clsfy.predict(tx_dev)
print 1.0*(ty_dev == ty_dev_pred).sum()/ty_dev_pred.shape[0]
print '\n\n'
print 'Checking on dev set...'
print 'True Labels for dev set'
print ty_dev
print 'Predictions on dev set'
ty_dev_pred=clsfy.predict(tx_dev)
print ty_dev_pred
print 1.0*(ty_dev == ty_dev_pred).sum()/ty_dev_pred.shape[0]
print '\n\n'
print 'Predictions on test set'
print clsfy.predict(tx_test)
print estimator.predict_proba(tx_test)
temp=clsfy.predict(tx_test)
temp=ty_dev_pred

for l in tempdict:
	print l,tempdict[l]


def gensif():
	infile2=os.path.dirname(os.path.abspath(__file__))+'/csvFiles/pred.csv'
	fh1=open(infile2,'w')
	for k in temp:
		k=int(k)
		print >>fh1,str(k)
	fh1.close()
	gengraph.main()
#print clf.decision_function(tex[0,0:vals[1]-2])






