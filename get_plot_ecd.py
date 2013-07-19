#!/usr/bin/env python

import sys
from numpy import *
from pylab import *
from utils import utils
from utils.constants import a0,alpha
from utils.conversion import au2ev

filename=sys.argv[1]

try:
    with open(filename): pass
except IOError:
    print 'error: program argument must be a valid filename'
    print '       provided filename =',filename
    sys.exit(2)

grey=(0.5,0.5,0.5)
broadening=0.12399

fig1 = figure(1,figsize=(8,6))

def line_skip(fileptr,n):
    for i in range(n):
        fileptr.readline()

def mv_to_line(fileptr,string):
    found=False
    while not found:
        l=fileptr.readline()
        found=(l.find(string)>=0)

def get_rot_strength(filename):
    maxstates=100
    E=zeros(maxstates)
    R=zeros(maxstates)
    infile=open(filename,'r')
    mv_to_line(infile,'Oscillator and Scalar Rotational Strengths')
    line_skip(infile,8)
    i=0
    l=infile.readline()
    while l.find('--------------') == -1:
        E[i] = float(l.split()[2]) 
        R[i] = float(l.split()[7])
        i+=1
        l=infile.readline()
    return i,E[:i],R[:i]

nstates,E,R = get_rot_strength(filename)

# convert energy (eV) to wavelength (nm)
l = 2*pi/(alpha*E/au2ev)*a0*1e9
for i in range(nstates):
    plot([l[i],l[i]],[0,R[i]] ,'b-',linewidth=1.7)

x,y=utils.lorentzian(E,R,broadening,3.5,6.0,1000)
# convert energy (eV) to wavelength (nm)
lx = 2*pi/(alpha*x/au2ev)*a0*1e9
plot(lx,y,'k-', linewidth=1.0)

setp(gca(),xlim=(210,320))
xlabel('Photon energy (eV)')
ylabel(r'Rotational strength ($10^{-40}$ esu$^2$ cm$^2$)')
figtext(0.35,0.93,'file: ' + filename)
figtext(0.15,0.85,str(nstates)+' states')
#------------------------------------------------------------------------

savefig('ecd.pdf')
