#!/usr/bin/python

import grn_mip_learn as gmipsl
import os
import re

newtext='./GRN.py --a1-dir /home/aditya/Project/dev --a2-dir /home/aditya/Project/dev --pred-sif /home/aditya/Project/output.sif /home/aditya/Project/dev/PMID-*.txt'

newtext1='./GRN.py --a1-dir /home/aditya/Project/train --a2-dir /home/aditya/Project/train --pred-sif /home/aditya/Project/output.sif /home/aditya/Project/train/PMID-*.txt'

fil1='/home/aditya/Project/csvFiles/Outcsv/mip_all_c.txt'
fh1=open(fil1,'w')
# It performs Cross Validation over these values of k or threshold
l1=[0.1+k*0.02 for k in range(20)]
for l in l1:
    print l
    gmipsl.mips_run(l)
    stream=os.popen(newtext1)
    templis=[]
    print >>fh1,l
    for lal in stream:
        lal=lal.strip()
        print >>fh1,lal 
    stream.close()
    #print templis
