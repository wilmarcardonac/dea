import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

a,H,H_prime,H_double_prime,ricci_scalar,ricci_scalar_prime,fMG,fMG_R,fMG_R_prime,fMG_R_double_prime,DE_density,w_DE,fMG_RR,fMG_RR_prime,fMG_RR_double_prime,Omega_M,Omega_DE,Omega_DE_prime = np.loadtxt('../output/background_functions.txt',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])

# H

fig = py.figure()

labels_H = (r'$|H(a)|$')

py.loglog(a,abs(H),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_H,loc=0)

py.savefig('conformal_Hubble_parameter.pdf')

py.close()

# H prime

fig = py.figure()

labels_H_prime = (r'$|Hprime(a)|$')

py.loglog(a,abs(H_prime),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_H_prime,loc=0)

py.savefig('derivative_conformal_Hubble_parameter.pdf')

py.close()

# H double prime

fig = py.figure()

labels_H_double_prime = (r'$|Hdoubleprime(a)|$')

py.loglog(a,abs(H_double_prime),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_H_double_prime,loc=0)

py.savefig('second_derivative_conformal_Hubble_parameter.pdf')

py.close()

# Ricci scalar

fig = py.figure()

labels_ricci_scalar = (r'$|ricciscalar(a)|$')

py.loglog(a,abs(ricci_scalar),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_ricci_scalar,loc=0)

py.savefig('ricci_scalar.pdf')

py.close()

# Ricci scalar prime

fig = py.figure()

labels_ricci_scalar_prime = (r'$|ricciscalarprime(a)|$')

py.loglog(a,abs(ricci_scalar_prime),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_ricci_scalar_prime,loc=0)

py.savefig('ricci_scalar_prime.pdf')

py.close()

# f(R(a))

fig = py.figure()

labels_fMG = (r'$|f(R(a))|$')

py.loglog(a,abs(fMG),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_fMG,loc=0)

py.savefig('fMG.pdf')

py.close()

# f_R(R(a))

fig = py.figure()

labels_fMG_R = (r'$|f_R(R(a))|$')

py.loglog(a,abs(fMG_R),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_fMG_R,loc=0)

py.savefig('fMG_R.pdf')

py.close()

# f_R_prime(R(a))

fig = py.figure()

labels_fMG_R_prime = (r'$|f_R_prime(R(a))|$')

py.loglog(a,abs(fMG_R_prime),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_fMG_R_prime,loc=0)

py.savefig('fMG_R_prime.pdf')

py.close()

# f_R_double_prime(R(a))

fig = py.figure()

labels_fMG_R_double_prime = (r'$|f_R_double_prime(R(a))|$')

py.loglog(a,abs(fMG_R_double_prime),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_fMG_R_double_prime,loc=0)

py.savefig('fMG_R_double_prime.pdf')

py.close()

# DE_density(a) times kappa 

fig = py.figure()

labels_DE_density = (r'$|DE_density(a)|$')

py.loglog(a,abs(DE_density),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_DE_density,loc=0)

py.savefig('DE_density.pdf')

py.close()

# w_DE(a)  

fig = py.figure()

labels_w_DE = (r'$|w_DE(a)|$')

py.loglog(a,abs(w_DE),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_w_DE,loc=0)

py.savefig('w_DE.pdf')

py.close()

# f_RR(R(a))

fig = py.figure()

labels_fMG_RR = (r'$|f_RR(R(a))|$')

py.loglog(a,abs(fMG_RR),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_fMG_RR,loc=0)

py.savefig('fMG_RR.pdf')

py.close()

# f_RR_prime(R(a))

fig = py.figure()

labels_fMG_RR_prime = (r'$|f_RR_prime(R(a))|$')

py.loglog(a,abs(fMG_RR_prime),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_fMG_RR_prime,loc=0)

py.savefig('fMG_RR_prime.pdf')

py.close()

# f_RR_double_prime(R(a))

fig = py.figure()

labels_fMG_RR_double_prime = (r'$|f_RR_double_prime(R(a))|$')

py.loglog(a,abs(fMG_RR_double_prime),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_fMG_RR_double_prime,loc=0)

py.savefig('fMG_RR_double_prime.pdf')

py.close()

# Omega_M(a)

fig = py.figure()

labels_Omega_M = (r'$|\Omega_m(a)|$')

py.loglog(a,abs(Omega_M),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_Omega_M,loc=0)

py.savefig('Omega_M.pdf')

py.close()

# Omega_DE(a)

fig = py.figure()

labels_Omega_DE = (r'$|\Omega_DE(a)|$')

py.loglog(a,abs(Omega_DE),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_Omega_DE,loc=0)

py.savefig('Omega_DE.pdf')

py.close()

# Omega_DE_prime(a)

fig = py.figure()

labels_Omega_DE_prime = (r'$|\Omega_DE_prime(a)|$')

py.loglog(a,abs(Omega_DE_prime),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_Omega_DE_prime,loc=0)

py.savefig('Omega_DE_prime.pdf')

py.close()

exit()


