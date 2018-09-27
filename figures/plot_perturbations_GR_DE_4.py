import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

aa,adm,aVm,aphiplus,achi = np.loadtxt('../output/analytical_solution.txt',unpack=True,usecols=[0,1,2,3,4])

# \delta_m

fig = py.figure()

labels_delta_m = (r'$|\delta_m^{ana}(a)|$',r'$|\delta_m^{num}(a)|$')

py.loglog(aa,abs(adm),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_delta_m,loc=0)

py.savefig('./perturbations/delta_m.pdf')

py.close()

# V_m

fig = py.figure()

labels_V_m = (r'$|V_m^{ana}(a)|$',r'$|V_m^{num}(a)|$')

py.loglog(aa,abs(aVm),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_V_m,loc=0)

py.savefig('./perturbations/V_m.pdf')

py.close()

# phi_plus

fig = py.figure()

labels_phi_plus = (r'$|\phi_{+}^{ana}(a)|$',r'$|\phi_{+}^{num}(a)|$')

py.loglog(aa,abs(aphiplus),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_phi_plus,loc=0)

py.savefig('./perturbations/phi_plus.pdf')

py.close()

# chi

fig = py.figure()

labels_chi = (r'$|\chi^{ana}(a)|$',r'$|\chi^{num}(a)|$')

py.semilogx(aa,abs(achi),color='red')

py.xlabel(r'$a$',fontsize='large')

#py.ylim(-0.005,0.005)

py.legend(labels_chi,loc=0)

py.savefig('./perturbations/chi.pdf')

py.close()

exit()


