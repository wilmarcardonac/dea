import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

#a,d,V,dm,Vm,aphi,apsi,cs2_MD,cs2_MD_prime,angular_frequency = np.loadtxt('../output/analytical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9])
a,d,V,dm,Vm,aphi,apsi = np.loadtxt('../output/analytical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5,6])

# \delta

fig = py.figure()

labels_delta = (r'$|\delta^{ana}(a)|$',r'$|\delta_m^{num}(a)|$')

py.loglog(a,abs(d),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_delta,loc=0)

py.savefig('./perturbations/delta.pdf')

py.close()

# V

fig = py.figure()

labels_V = (r'$|V^{ana}(a)|$',r'$|V_m^{num}(a)|$')

py.loglog(a,abs(V),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_V,loc=0)

py.savefig('./perturbations/V.pdf')

py.close()

# \delta_m

fig = py.figure()

labels_delta_m = (r'$|\delta_m^{ana}(a)|$',r'$|\delta_m^{num}(a)|$')

py.loglog(a,abs(dm),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_delta_m,loc=0)

py.savefig('./perturbations/delta_m.pdf')

py.close()

# V_m

fig = py.figure()

labels_V_m = (r'$|V_m^{ana}(a)|$',r'$|V_m^{num}(a)|$')

py.loglog(a,abs(Vm),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_V_m,loc=0)

py.savefig('./perturbations/V_m.pdf')

py.close()

# phi

fig = py.figure()

labels_phi = (r'$|\phi^{ana}(a)|$',r'$|\phi_{+}^{num}(a)|$')

py.loglog(a,abs(aphi),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_phi,loc=0)

py.savefig('./perturbations/phi.pdf')

py.close()

# psi

fig = py.figure()

labels_psi = (r'$|\psi^{ana}(a)|$',r'$|\chi^{num}(a)|$')

py.loglog(a,abs(apsi),color='red')

py.xlabel(r'$a$',fontsize='large')

#py.ylim(-0.005,0.005)

py.legend(labels_psi,loc=0)

py.savefig('./perturbations/psi.pdf')

py.close()

exit()
# cs2_MD

fig = py.figure()

labels_cs2_MD = (r'$|c_s^2(a)|$',r'$|\chi^{num}(a)|$')

py.loglog(a,abs(cs2_MD),color='red')

py.xlabel(r'$a$',fontsize='large')

#py.ylim(-0.005,0.005)

py.legend(labels_cs2_MD,loc=0)

py.savefig('./perturbations/cs2_MD.pdf')

py.close()

# cs2_MD_prime

fig = py.figure()

labels_cs2_MD_prime = (r'$|c_s^{2,prime}(a)|$',r'$|\chi^{num}(a)|$')

py.loglog(a,abs(cs2_MD_prime),color='red')

py.xlabel(r'$a$',fontsize='large')

#py.ylim(-0.005,0.005)

py.legend(labels_cs2_MD_prime,loc=0)

py.savefig('./perturbations/cs2_MD_prime.pdf')

py.close()

# angular_frequency

fig = py.figure()

labels_angular_frequency = (r'$|\omega^{2}(a)|$',r'$|\chi^{num}(a)|$')

py.loglog(a,abs(angular_frequency),color='red')

py.xlabel(r'$a$',fontsize='large')

#py.ylim(-0.005,0.005)

py.legend(labels_angular_frequency,loc=0)

py.savefig('./perturbations/angular_frequency.pdf')

py.close()

exit()


