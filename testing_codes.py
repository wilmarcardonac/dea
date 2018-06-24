import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

a,w,wprime,H,cs2,ceff2,Omegam,OmegaDE,dprho,dm,Vm,dde,Vde,pi,GeffGN,Qeff,ca2,x = np.loadtxt('./output/functions.txt',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])

#am4,dmm4,dem4,vmm4,vdem4,phim4,psim4,dmprime = np.loadtxt('./output/numerical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5,6,7])

fig = py.figure()

labels_perturbations = (r'$\delta_m$',r'$|\delta|$')

# Perturbations

py.loglog(a,dm,color='blue')
py.loglog(a,abs(dde),color='black')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_perturbations,loc=0)

py.savefig('./output/testing_codes_matter_perturbations.pdf')

py.close()

fig = py.figure()

labels_terms_pressure_perturbation = ('pressure perturbation over density')

# Terms

py.loglog(a,abs(dprho),color='red') 

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_terms_pressure_perturbation,loc=0)

py.savefig('./output/testing_codes_pressure_perturbation.pdf')

py.close()

fig = py.figure()

labels_terms_velocity_perturbation = (r'$|V_m|$',r'$|V|$')

# Terms
py.loglog(a,abs(Vm),color='green')
py.loglog(a,abs(Vde),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_terms_velocity_perturbation,loc=0)

py.savefig('./output/testing_codes_velocity_perturbation.pdf')

py.close()

fig = py.figure()

labels_functions_fluid = (r'$w(a)$',r'$w^{prime}(a)$')

# w 

py.semilogx(a,w,color='red')
py.semilogx(a,wprime,color='blue')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_functions_fluid,loc=0)

py.savefig('./output/testing_codes_equation_of_state.pdf')

py.close()

fig = py.figure()

labels_H = (r'$H(a)/H_0$')

# H

py.loglog(a,abs(H),color='red')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_H,loc=0)

py.savefig('./output/testing_codes_H.pdf')

py.close()

fig = py.figure()

labels_velocities = (r'$c_{s}^2$',r'$c_{eff}^2$',r'$c_a^2$')

# Velocities

py.semilogx(a,cs2,color='red')
py.semilogx(a,ceff2,color='blue')
py.semilogx(a,ca2,color='black')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_velocities,loc=0)

py.savefig('./output/testing_codes_speeds.pdf')

py.close()


fig = py.figure()

labels_density_parameters = (r'$\Omega_m$',r'$ \Omega_{DE}$')

# Density parameters

py.loglog(a,abs(Omegam))

py.loglog(a,abs(OmegaDE))

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_density_parameters,loc=0)

py.savefig('./output/testing_codes_density_parameters.pdf')

py.close()

fig = py.figure()

labels_anisotropic_stress = (r'$\Pi(a)$')

# Anisotropic stress

py.loglog(a,abs(pi))

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_anisotropic_stress,loc=0)

py.savefig('./output/testing_codes_anisotropic_stress.pdf')

py.close()

fig = py.figure()

labels_Geff = (r'$G_{eff}/G_N$')

# Geff

py.loglog(a,abs(GeffGN))

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_Geff,loc=0)

py.savefig('./output/testing_codes_Geff.pdf')

py.close()

fig = py.figure()

labels_Qeff = (r'$Q_{eff}$')

# Qeff

py.loglog(a,abs(Qeff))

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_Qeff,loc=0)

py.savefig('./output/testing_codes_Qeff.pdf')

py.close()

fig = py.figure()

labels_X = (r'$X$')

# X

py.loglog(a,abs(x))

py.xlabel(r'$a$',fontsize='large')

#py.ylim(0.1,2)

py.legend(labels_X,loc=0)

py.savefig('./output/testing_codes_X.pdf')

py.close()

exit()

