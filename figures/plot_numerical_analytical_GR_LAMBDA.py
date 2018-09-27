import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py


#######################################################################################
#######################################################################################
# Loading data for numerical and analytical solutions
#######################################################################################
#######################################################################################

a,dm,Vm,phi_plus,chi = np.loadtxt('../output/numerical_solution.txt',unpack=True,usecols=[0,1,2,3,4])

x,dmm,Vmm,phii_plus,chii = np.loadtxt('../output/analytical_solution.txt',unpack=True,usecols=[0,1,2,3,4])

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

# \phi_plus

fig = py.figure()

labels_phi_plus = (r'$\phi_{+}^{num}$',r'$\phi_{+}^{ana}$')

py.loglog(a,abs(phi_plus),linestyle='-',c='b')
py.loglog(x,abs(phii_plus),linestyle='dashed',c='r')
py.xlabel(r'$a$',fontsize='large')

py.legend(labels_phi_plus,loc=0)

py.savefig('comparison_phi_plus.pdf')

py.close()

# \csi

fig = py.figure()

labels_chi = (r'$\chi^{num}$',r'$\chi^{ana}$')

py.plot(a,chi,linestyle='-',c='b')
py.plot(x,chii,linestyle='dashed',c='r')
py.xlabel(r'$a$',fontsize='large')

py.legend(labels_chi,loc=0)

py.savefig('comparison_chi.pdf')

py.close()

exit()

