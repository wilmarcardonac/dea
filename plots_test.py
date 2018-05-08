import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py


#######################################################################################
#######################################################################################
# Loading data for numerical and analytical solutions
#######################################################################################
#######################################################################################

a2,g1,g2,g3,g4,g5,g6 = np.loadtxt('./output/analytical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5,6])

am4,dmm4,dem4,vmm4,vdem4,phim4,psim4,dmprime = np.loadtxt('./output/numerical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5,6,7])

fig = py.figure()

labels_dmprime = (r'$ f \sigma_8 $')

py.plot(1./am4-1.,0.8*am4*dmprime/dmm4[-1])

py.xlabel(r'$a$',fontsize='large')
py.xlim(0.,2.)
py.legend(labels_dmprime,loc=0,ncol=2)

py.savefig('./output/fsigma8.pdf')
py.close()

fig = py.figure()

labels = (r'$\phi^{num}$',r'$\psi^{num}$',r'$\phi$',r'$\psi$')

py.loglog(am4,abs(phim4))
py.loglog(am4,abs(psim4))
py.loglog(a2,abs(g5),linestyle='dotted')
py.loglog(a2,abs(g6),linestyle='dashed')

#py.ylim(-3,2)
#py.xlim(9.e-5,2.e-4)
# Label for x-axis

py.xlabel(r'$a$',fontsize='large')

py.legend(labels,loc=0,ncol=2)

py.savefig('./output/potentials.pdf')
py.close()

fig = py.figure()

#labels = (r'$|\delta_{m}|$',r'$|V_{m}|$',r'$|\delta_{de}|$',r'$|V_{de}|$',
labels = (r'$|\delta_{de}^{num}|$',r'$|V_{de}^{num}|$',r'$|\delta_{m}^{num}|$',r'$|V_{m}^{num}|$',r'$|\delta^{sub}|$',r'$|V^{sub}|$')#,r'$v_{de}^{sup-hor}$',r'$v_{de}^{sub-sound}$')

#py.loglog(a,abs(dm),linestyle='dotted')
#py.loglog(a,abs(Vm),linestyle='dashed')
#py.loglog(a,abs(dde),linestyle='dotted')
#py.loglog(a,abs(Vde),linestyle='dashed')
py.loglog(am4,abs(dem4))
py.loglog(am4,abs(vdem4))
py.loglog(am4,abs(dmm4))
py.loglog(am4,abs(vmm4))
py.loglog(a2,abs(g1),linestyle='dotted',color='blue')
py.loglog(a2,abs(g2),linestyle='dotted',color='green')
py.loglog(a2,abs(g3),linestyle='dotted',color='red')
py.loglog(a2,abs(g4),linestyle='dotted',color='cyan')

#py.ylim(-3,2)
#py.xlim(9.e-5,2.e-4)
# Label for x-axis

py.xlabel(r'$a$',fontsize='large')

#######################################################################################
# Placing the legend 
#######################################################################################

py.legend(labels,loc=0,ncol=2)
#py.legend(labels,bbox_to_anchor=(-1.2,2.4,2.2,.102),loc=3,ncol=6,mode='expand',borderaxespad=0.)

#######################################################################################
# Saving the figure 
#######################################################################################

#py.savefig('example.pdf',bbox_inches='tight')

#######################################################################################
# Show the figure
#######################################################################################

#py.tight_layout()
py.savefig('./output/matter_perturbations.pdf')
py.close()


exit()

