#!/usr/bin/env python 

"""molecule [options] filename

purpose:
analyze Dalton molecule input file written in nonformatted
keyword-based format

options:
-h             prints this information
-e atom        define atom as ECP center
-o outfile     prints molecular coordinates in xyz-format to be
               read by e.g. molden (*.com type filenames will produce
               output file in gaussian format to be read by e.g. gview)
-v             verbose"""

import sys,math,getopt
import atomic_data, read_molinp
from conversion import au2ang

# process list of arguments
outfile = ''
ecp = False
ecp_atoms = []
verbose = False

try:
    opts, args = getopt.getopt(sys.argv[1:],'e:ho:v')
except getopt.error, msg:
    print msg
    print 'for help use -h'
    sys.exit(2)
else:
    if len(sys.argv) == 1:
        print __doc__
        sys.exit(0)
    else:
        filename = sys.argv[-1]        

for o, a in opts:
    if o == '-e':
        ecp = True
        ecp_atoms.append(a)
    elif o == '-h':
        print __doc__
        sys.exit(0)
    elif o == '-o':
        outfile = a
    elif o == '-v':
        verbose = True

# program start
def atomdist(a,b):
    dist=math.sqrt( (a[1]-b[1])**2 + \
                    (a[2]-b[2])**2 + \
                    (a[3]-b[3])**2 )
    return dist

def addsym_atoms(atomlist,symgenlist):
    for symgen in symgenlist:
        templist = atomlist[:]
        for tempatom in templist:
            symatom = [tempatom[0],\
                       tempatom[1]*symgen[0],\
                       tempatom[2]*symgen[1],\
                       tempatom[3]*symgen[2]]
            addatom = True
            for atom in atomlist:
                if atomdist(atom,symatom) < 0.1:
                    addatom = False
            if addatom:
                atomlist.append(symatom)

    if verbose:
        print 'After addition of symmetry dependent atoms:'
        for atom in atomlist:
            print atom

def nuclear_repulsion(atomlist):
    nucrep=0.0; total_charge=0.0
    for i in range(len(atomlist)):
        
        if atomlist[i][0] in ecp_atoms:
            qi = atomic_data.ecp_charge(atomlist[i][0])
        else:
            qi = atomic_data.charge(atomlist[i][0])

        total_charge += qi

        for j in range(i+1,len(atomlist)):

            if atomlist[j][0] in ecp_atoms:
                qj = atomic_data.ecp_charge(atomlist[j][0])
            else:
                qj = atomic_data.charge(atomlist[j][0])

            nucrep += qi*qj / atomdist(atomlist[i],atomlist[j])

    return nucrep, total_charge

def center_of_mass(atomlist):
    total_mass=0.0
    rcm=['CM', 0.0, 0.0, 0.0]
    for atom in atomlist:
        total_mass += atomic_data.mass(atom[0])
        rcm[1] += atomic_data.mass(atom[0])*atom[1]
        rcm[2] += atomic_data.mass(atom[0])*atom[2]
        rcm[3] += atomic_data.mass(atom[0])*atom[3]
    for i in range(1,4):
        rcm[i] = rcm[i]/total_mass
        if rcm[i]<1e-10:
            rcm[i]=0.0
    return rcm

def nuclear_dipole_moment(atomlist):
    rcc=['CC', 0.0, 0.0, 0.0]
    for atom in atomlist:
        rcc[1] += atomic_data.charge(atom[0])*atom[1]
        rcc[2] += atomic_data.charge(atom[0])*atom[2]
        rcc[3] += atomic_data.charge(atom[0])*atom[3]
    for i in range(1,4):
        if rcc[i]<1e-10:
            rcc[i]=0.0
    return rcc[1:]

def cavity_radius(rcm,atomlist):
    cavrad = 0.0
    for atom in atomlist:
        dist = atomdist(rcm,atom) + atomic_data.vdW(atom[0])/au2ang
        if dist > cavrad:
            cavrad = dist
            saveatom = atom
            saveatomdist=atomdist(rcm,atom)
    if verbose:
        print 'Outermost atom in cavity:'
        print saveatom
        print 'atom_to_mass-center dist=',saveatomdist,'a.u.'
        print 'cavity radius=',cavrad,'a.u.'
    return cavrad

atomlist, symgenlist = read_molinp.get_atoms(filename)
addsym_atoms(atomlist,symgenlist)

natoms=len(atomlist)
rcm=center_of_mass(atomlist)
nucrep,total_charge=nuclear_repulsion(atomlist)
nucdip=nuclear_dipole_moment(atomlist)
cavrad=cavity_radius(rcm,atomlist)

print '%25s %5d' % ('Number of atoms =', natoms)
print '%25s %7.1f' % ('Total nuclear charge =', total_charge)
print '%25s %14.8f' % ('Nuclear repulsion =', nucrep)
print '%25s %14.8f' % ('Cavity radius =', cavrad)
print '%25s' % ('Nuclear dipole moment ='), nucdip
print '%25s' % ('Center of mass ='), rcm[1:]

if outfile.find('.com')>0:
    out=open(outfile,'w')
    out.write('%NProcLinda=8'+'\n')
    out.write('%NProcShared=4'+'\n')
    out.write('#p opt b3lyp/cc-pVDZ Symm=Loose'+'\n'*2)
    out.write('Comment line' +'\n'*2)
    out.write('0 1\n')
    for atom in atomlist:
        l = '%2s%14.6f%12.6f%12.6f\n' % \
            (atom[0],atom[1]*au2ang,atom[2]*au2ang,atom[3]*au2ang)
        out.write(l)
    out.write('\n'+'#end of file\n')
    out.close
elif outfile:
    out=open(outfile,'w')
    out.write(str(natoms)+'\n'*2)
    for atom in atomlist:
        l = '%2s%14.6f%12.6f%12.6f\n' % \
            (atom[0],atom[1]*au2ang,atom[2]*au2ang,atom[3]*au2ang)
        out.write(l)
    out.close

# End of file
