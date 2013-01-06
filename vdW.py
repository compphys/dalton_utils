#!/usr/bin/env python

import sys,utils
import utils.gauss_legendre 

# polarizabilities alpha(iw) for A and B
A=[]
B=[]

fileA=sys.argv[1]
if len(sys.argv)==2:
    fileB=fileA
elif len(sys.argv)==3:
    fileB=sys.argv[2]

l=[]
for file in [fileA,fileB]:
    dalout = open(file,'r')
    s=dalout.readline()
    while s:
        p=s.find('Averaged value')
        if p != -1:
            l.append(float(s[p+30:p+50]))
        s=dalout.readline()
# we assume that the first frequency is 0.0 which is not 
# to be included in the quadratue integration
n=len(l)/2-1 # number quadrature points
A=l[1:n+1]; B=l[n+2:2*n+3]; 
print "alpha(0;0): A=", l[0], "  B=", l[n+1]

if n==8:
    utils.gauss_legendre.gauss_legendre_8pt(A,B)
elif n==10:
    utils.gauss_legendre.gauss_legendre_10pt(A,B)
elif n==12:
    utils.gauss_legendre.gauss_legendre_12pt(A,B)
else:
    print "vdW.py: incorrect number of alpha values"
    print "alpha=",A
    sys.exit(2)
