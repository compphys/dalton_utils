#!/usr/bin/env python

def get_atoms(filename,verbose=False):
    
    import sys
    from utils import atomic_data
    from utils.conversion import au2ang

    try:
        dalinp=open(filename,'r')
    except:
        print 'file not found:', filename
        sys.exit(2)

    atomlist=[]
    symgenlist=[]

    fact = 1.0
    ntypes_count=0
    line=dalinp.readline()
    while line:
        if line.find('Atomtypes')>-1:
            pos=line.find('Atomtypes')+len('Atomtypes')+1
            ntypes=int(line[pos:pos+2])
        if line.find('Angstrom')>-1:
            fact=1.0/au2ang
        if line.find('Generators')>-1:
            if line.find(' X ')>0:
                symgenlist.append([-1, 1, 1])
            if line.find(' Y ')>0:
                symgenlist.append([1, -1, 1])
            if line.find(' Z ')>0:
                symgenlist.append([1, 1, -1])
            if line.find(' XY ')>0:
                symgenlist.append([-1, -1, 1])
            if line.find(' XZ ')>0:
                symgenlist.append([-1, 1, -1])
            if line.find(' YZ ')>0:
                symgenlist.append([1, -1, -1])
            if line.find('XYZ')>0:
                symgenlist.append([-1, -1, -1])
        if line.find('Atoms')>-1:
            ntypes_count+=1
            pos=line.find('Atoms')+len('Atoms')+1
            natoms=int(line[pos:].split()[0])
            pos=line.find('Charge')+len('Charge')+1
            charge=float(line[pos:pos+4])
            label=''
            for l in atomic_data.labels:
                if not label:
                    if atomic_data.charge(l)==charge:
                        label=l
            if not label:
                print 'atomic label not found for charge=',charge
                sys.exit()
            for i in range(0,natoms):
                line=dalinp.readline()
                field=line.split()
                atom=[label,\
                      float(field[1])*fact,\
                      float(field[2])*fact,\
                      float(field[3])*fact]
                atomlist.append(atom)
        line=dalinp.readline()

    if verbose:
        if fact == 1.0:
            print 'Input file coordinates were given in a.u.'
        else:
            print 'Input file coordinates were given in Angstrom'
        print 'Atoms specified in input file:'
        for atom in atomlist:
            print atom
        print 'Symmetry generators specified in input file:'
        for symgen in symgenlist:
            print symgen

    if ntypes!=ntypes_count:
        print 'Incorrect specification of number of atom types.'
        sys.exit(2)
        
    return atomlist,symgenlist

if __name__ == '__main__':
    print dir()
