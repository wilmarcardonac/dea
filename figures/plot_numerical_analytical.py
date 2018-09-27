import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py


#######################################################################################
#######################################################################################
# Loading data for numerical and analytical solutions
#######################################################################################
#######################################################################################

a,dm,d,Vm,V,phi,psi = np.loadtxt('../output/numerical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5,6])

x,de,Vde,dmm,Vmm,phii,psii,cs2 = np.loadtxt('../output/analytical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5,6,7])

# \delta_m

fig = py.figure()

labels_delta_m = (r'$\delta_m^{num}$',r'$\delta_m^{ana}$')

py.loglog(a,abs(dm),linestyle='-',c='b')
py.loglog(x,abs(dmm),linestyle='dashed',c='r')
py.xlabel(r'$a$',fontsize='large')

py.legend(labels_delta_m,loc=0)

py.savefig('comparison_delta_m.pdf')

py.close()

# V_m

fig = py.figure()

labels_V_m = (r'$V_m^{num}$',r'$V_m^{ana}$')

py.loglog(a,abs(Vm),linestyle='-',c='b')
py.loglog(x,abs(Vmm),linestyle='dashed',c='r')
py.xlabel(r'$a$',fontsize='large')

py.legend(labels_V_m,loc=0)

py.savefig('comparison_V_m.pdf')

py.close()

# \delta

fig = py.figure()

labels_delta = (r'$\delta^{num}$',r'$\delta^{ana}$')

py.loglog(a,abs(d),linestyle='-',c='b')
py.loglog(x,abs(de),linestyle='dashed',c='r')
py.xlabel(r'$a$',fontsize='large')

py.legend(labels_delta,loc=0)

py.savefig('comparison_delta.pdf')

py.close()

# V

fig = py.figure()

labels_V = (r'$V^{num}$',r'$V^{ana}$')

py.loglog(a,abs(V),linestyle='-',c='b')
py.loglog(x,abs(Vde),linestyle='dashed',c='r')
py.xlabel(r'$a$',fontsize='large')

py.legend(labels_V,loc=0)

py.savefig('comparison_V.pdf')

py.close()

# \phi

fig = py.figure()

labels_phi = (r'$\phi^{num}$',r'$\phi^{ana}$')

py.loglog(a,abs(phi),linestyle='-',c='b')
py.loglog(x,abs(phii),linestyle='dashed',c='r')
py.xlabel(r'$a$',fontsize='large')

py.legend(labels_phi,loc=0)

py.savefig('comparison_phi.pdf')

py.close()

# \psi

fig = py.figure()

labels_psi = (r'$\psi^{num}$',r'$\psi^{ana}$')

py.loglog(a,abs(psi),linestyle='-',c='b')
py.loglog(x,abs(psii),linestyle='dashed',c='r')
py.xlabel(r'$a$',fontsize='large')

py.legend(labels_psi,loc=0)

py.savefig('comparison_psi.pdf')

py.close()

exit()

