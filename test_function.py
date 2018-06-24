import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

a,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26,f27,f28 = np.loadtxt('./output/functions.txt',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28])

a2,g1,g2 = np.loadtxt('./output/analytical_solution.txt',unpack=True,usecols=[0,1,2])

fig = py.figure()

labels_terms_pressure_perturbation = ('t1','t2','t3','t4','t1+t2','t1+t2+t3','t1+t2+t3+t4')

# Terms

py.loglog(a,abs(f19),color='red') # t1  
py.loglog(a,abs(f20),linestyle='dashed',color='red') #t2
py.loglog(a,abs(f26),color='red',linestyle='-.') # t3
py.loglog(a,abs(f27),color='red',linestyle='dotted') # t4

#py.loglog(a,abs(f21),color='blue',linestyle='dotted') #t3
#py.loglog(a,abs(f22),color='blue') #t4
#py.loglog(a,abs(f23),color='blue',linestyle='dashed') #t5
py.loglog(a,f19+f20,color='green') #t1+t2
#py.loglog(a,abs(f21+f22+f23),color='magenta') # t3+t4+t5
#py.loglog(a,abs(f24),color='black') # denominator
py.loglog(a,f19+f20+f26,color='black') #t1+t2+t3
py.loglog(a,f19+f20+f26+f27,color='blue') #t1+t2+t3+t4

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_terms_pressure_perturbation,loc=0)

py.savefig('./output/terms_comparison_pressure_perturbation.pdf')

py.close()

fig = py.figure()

labels_terms_velocity_perturbation = ('t5','t6','t7','t5+t6','t5+t6+t7')

# Terms

py.loglog(a,abs(f22),color='red') #t5
py.loglog(a,abs(f21),color='red',linestyle='dotted') #t6
py.loglog(a,abs(f23),color='red',linestyle='dashed') #t7
py.loglog(a,f22+f21,color='green') #t5+t6
py.loglog(a,abs(f22+f21+f23),color='magenta') # t5+t6+t7
#py.loglog(a,abs(f24),color='black') # denominator
#py.loglog(a,f19+f20+f26,color='black') #t1+t2+t6
#py.loglog(a,f19+f20+f26+f27,color='blue') #t1+t2+t6+t7

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_terms_velocity_perturbation,loc=0)

py.savefig('./output/terms_comparison_velocity_perturbation.pdf')

py.close()

fig = py.figure()

labels_functions_fluid = (r'$c_s^2$',r'$w(a)$',r'$\pi$',r'$w^{prime}(a)$')

# w cs2 pi

py.semilogx(a,f1,linestyle='dotted')
py.semilogx(a,f17,linestyle='dashed')
py.semilogx(a,f11)
py.semilogx(a,f28)

py.xlabel(r'$a$',fontsize='large')
py.ylim(-2.,1.)
py.legend(labels_functions_fluid,loc=0)

py.savefig('./output/fluid_constant_cs_w_pi.pdf')

py.close()

fig = py.figure()

labels_potentials = (r'$\phi$',r'$\psi$')

# Potentials

py.loglog(a,abs(f14),linestyle='dotted')
py.loglog(a,abs(f15),linestyle='dashed')
py.xlabel(r'$a$',fontsize='large')

py.legend(labels_potentials,loc=0)

py.savefig('./output/potentials_test.pdf')

py.close()

fig = py.figure()

labels_velocities = (r'$\frac{\delta P}{\delta \rho}$',r'$\frac{2 \pi}{3 \delta}$',r'$c_{eff}^2$',r'$c$',r'$c_{eff,full}^2$')

# Velocities

py.semilogx(a,f1,linestyle='dotted')
py.semilogx(a,f2,linestyle='dashed')
py.semilogx(a,f3,color='red')
py.semilogx(a,f4,color='magenta')
py.semilogx(a,f25,color='red',linestyle='dashed') 

py.xlabel(r'$a$',fontsize='large')
py.ylim(-1.,2.)
py.legend(labels_velocities,loc=0)

py.savefig('./output/velocities.pdf')

py.close()

fig = py.figure()

labels_perturbations = (r'$\delta_m$',r'$|V_m|$',r'$|\delta^{sub-sound}|$',r'$|V^{sub-sound}|$',r'$|\delta^{super-sound}|$',r'$|V^{super-sound}|$',r'$|\delta^{sub}|$',r'$|V^{sub}|$')

# Perturbations

py.loglog(a,f5)
py.loglog(a,abs(f6),color='green')
py.loglog(a,abs(f7))
py.loglog(a,abs(f8))
py.loglog(a,abs(f9),color='blue',linestyle='dashed')
py.loglog(a,abs(f10),color='green',linestyle='dashed')
py.loglog(a2,abs(g1),linestyle='dotted',color='blue')
py.loglog(a2,abs(g2),linestyle='dotted',color='green')

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_perturbations,loc=0)

py.savefig('./output/perturbations.pdf')

py.close()

fig = py.figure()

labels_perturbations = (r'$\pi(a)$',r'$ \delta P(a)$',r"\pi^{'}(a)")

# Perturbations

py.loglog(a,abs(f11))

py.loglog(a,abs(f12))

py.loglog(a,abs(f13))

py.xlabel(r'$a$',fontsize='large')

py.legend(labels_perturbations,loc=0)

py.savefig('./output/anisotropic_stress_pressure_perturbation.pdf')

py.close()

exit()

