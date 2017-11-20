import numpy
import string


def write(mol, fname,vel=None):
    '''
    Write lammpstrj to visualize GCMD runs, write lambda into velocities
    :Parameters:
        -fname  (str): name of the xyz file
        -mol    (obj): instance of a molclass
    '''
    natoms = mol.natoms 
    f = open(fname,"w")
    #### timestep header, not sure if necessary
    f.write('ITEM: TIMESTEP\n0.1\n')

    #### number of atoms with header
    f.write('ITEM: NUMBER OF ATOMS\n')
    f.write (str(natoms)+'\n')
    
    #### cell parameters with header
    # there is a problem here: either we wrap in box or we have to shift cell xyzminmax by the respective half 
    # as a test: go for the first one!
    mol.wrap_in_box()
    f.write('ITEM: BOX BOUNDS pp pp pp\n')
    f.write ('0.00 %f \n' % mol.cell[0,0])
    f.write ('0.00 %f \n' % mol.cell[1,1])
    f.write ('0.00 %f \n' % mol.cell[2,2])
    #### hope everything is orthorombic! otherwise tis aint gonna work! maybe irellevant anyway for  visualization ;)

    #### atom positions (and velocities) with header
    # i think vmd needs all, vx,vy,vz in order to read it properly!
    f.write('ITEM: ATOMS id type x y z vx vy vz\n')
    if vel is None: vel = numpy.zeros((natoms,3))
    for i in range(natoms):
        f.write("%i %2s %f %f %f %f %f %f \n" % (i,mol.elems[i], mol.xyz[i][0], mol.xyz[i][1], mol.xyz[i][2]))
        #f.write("%2s %12.6f %12.6f %12.6f\n" % (mol.elems[i], mol.xyz[i,0], mol.xyz[i,1], mol.xyz[i,2]))
    f.close()
    return


def write_raw(f,stepcount,natoms,cell,elems,xyz,lamb):
    '''
    Write lammpstrj to visualize GCMD runs, write lambda into velocities
    :Parameters:
        -fname  (str): name of the xyz file
        -mol    (obj): instance of a molclass
    '''
    #f = open(fname,"w")
    #### timestep header, not sure if necessary
    f.write('ITEM: TIMESTEP\n%12.1f\n' % float(stepcount))

    #### number of atoms with header
    f.write('ITEM: NUMBER OF ATOMS\n')
    f.write (str(natoms)+'\n')
    
    #### cell parameters with header
    # there is a problem here: either we wrap in box or we have to shift cell xyzminmax by the respective half 
    # as a test: go for the first one!
    f.write('ITEM: BOX BOUNDS pp pp pp\n')
    f.write ('0.00 %f \n' % cell[0,0])
    f.write ('0.00 %f \n' % cell[1,1])
    f.write ('0.00 %f \n' % cell[2,2])
    #### hope everything is orthorombic! otherwise tis aint gonna work! maybe irellevant anyway for  visualization ;)

    #### atom positions (and velocities) with header
    # i think vmd needs all, vx,vy,vz in order to read it properly!
    f.write('ITEM: ATOMS id type x y z vx vy vz\n')
    for i in range(natoms):
        f.write("%i %2s %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n" % (i,elems[i], xyz[i,0], xyz[i,1], xyz[i,2],lamb[i], 0.0,0.0))
        #f.write("%2s %12.6f %12.6f %12.6f\n" % (mol.elems[i], mol.xyz[i,0], mol.xyz[i,1], mol.xyz[i,2]))
    #f.close()
    return
