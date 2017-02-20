import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization! 
#fig = plt.figure()
import numpy


class plotter(object):
    
    def __init__(self,mol):
        self.mol = mol
        
    def add_cell(self,ax):
        mol = self.mol
        def axplt(ax,xx,yy):
            ax.plot([xx[0],yy[0]],[xx[1],yy[1]],[xx[2],yy[2]],color='black',linewidth=3)
            
        cell = mol.get_cell()
        zero = numpy.zeros(3)
        x,y,z=cell[:,0],cell[:,1],cell[:,2]
        axplt(ax,zero,x)
        axplt(ax,zero,y)
        axplt(ax,zero,z)
        #axplt(ax,,)
        axplt(ax,x,x+y)
        axplt(ax,y,x+y)
        axplt(ax,z+x,x+y+z)
        axplt(ax,x+y,x+y+z)
        axplt(ax,z+y,x+y+z)
        axplt(ax,z,z+x)
        axplt(ax,z,z+y)
        axplt(ax,x,x+z)
        axplt(ax,y,y+z)
        
        
        return ax
        
    def plot(self,scell=False,bonds=False,labels=False):
        mol = self.mol
        col = ['r','g','b','m','c','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k']+['k']*200
        fig = plt.figure(figsize=plt.figaspect(1.0)*1.5)
        #ax = fig.add_subplot(111, projection='3d')
        ax = Axes3D(fig)
        ax = self.add_cell(ax)
        atd = {}
        for i,aa in enumerate(list(set(mol.atypes))):
            atd.update({aa:col[i]})
        print atd
        if bonds:
            for i in range(mol.natoms):
                conn = mol.conn[i]
                for j in range(len(conn)):
                    if mol.pconn:
                        if numpy.sum(numpy.abs(mol.pconn[i][j])) == 0:
                            ax.plot([mol.xyz[i][0],mol.xyz[conn[j]][0]],[mol.xyz[i][1],mol.xyz[conn[j]][1]],[mol.xyz[i][2],mol.xyz[conn[j]][2]],color='black')
                        else:
                            xyznew = mol.get_image(mol.xyz[conn[j]],mol.pconn[i][j])
                            ax.scatter(xyznew[0],xyznew[1],xyznew[2],color='orange')
                            ax.plot([mol.xyz[i][0],xyznew[0]],[mol.xyz[i][1],xyznew[1]],[mol.xyz[i][2],xyznew[2]],color='green')
                    else:
                        ax.plot([mol.xyz[i][0],mol.xyz[conn[j]][0]],[mol.xyz[i][1],mol.xyz[conn[j]][1]],[mol.xyz[i][2],mol.xyz[conn[j]][2]],color=atd[mol.atypes[i]])

        if labels:
            for i in range(mol.natoms):
                label = str(i)+'-'+str(mol.atypes[i]) +'-'+str(len(mol.conn[i]))
                ax.text(mol.xyz[i][0], mol.xyz[i][1], mol.xyz[i][2]+0.005, label, color='k',fontsize=9)
        if scell:
            xyz3 = mol.make_333(out=True)
            xyz3 =  numpy.array(xyz3)
            ax.scatter(xyz3[:,0],xyz3[:,1],xyz3[:,2],color='r',alpha=0.5)
        xyz=numpy.array(mol.xyz)
        for i,xx in enumerate(xyz):

            ax.scatter(xx[0],xx[1],xx[2],color=atd[mol.atypes[i]])
        minbound = numpy.min([numpy.min(xyz[:,0]),numpy.min(xyz[:,1]),numpy.min(xyz[:,2])])
        maxbound = numpy.max([numpy.max(xyz[:,0]),numpy.max(xyz[:,1]),numpy.max(xyz[:,2])])
        ax.auto_scale_xyz([0.0, maxbound], [0.0, maxbound], [0.0, maxbound])
        #ax.scatter(xyz1[:,0],xyz1[:,1],xyz1[:,2],color='k')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()