import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py


#######################################################################################
#######################################################################################
# Loading data for numerical and analytical solutions
#######################################################################################
#######################################################################################

am4,dmm4,dem4,vmm4,vdem4,Z = np.loadtxt('../output/numerical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5])

#xm4,dea1m4,vdea1m4,dm,Vm = np.loadtxt('../output/analytical_solution.txt',unpack=True,usecols=[0,1,2,3,4])

# \Delta_m

fig = py.figure()
labels = (r'$\Delta_{m}^{num}$')#,r'$\delta_{m}^{ana}$')
py.loglog(am4,abs(dmm4),color='r')
#py.loglog(xm4,abs(dm),color='b',linestyle='dashed')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('comparison_Delta_m.pdf')
py.close()

# Delta_de

fig = py.figure()
labels = (r'$\Delta_{de}^{num}$')#,r'$V_{m}^{ana}$')
#py.loglog(am4,abs(dem4),color='r')
py.plot(am4,abs(dem4),color='r')
#py.loglog(xm4,abs(Vm),color='b',linestyle='dashed')
py.xscale('log')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('comparison_Delta_de.pdf')
py.close()

# Theta_m

fig = py.figure()
labels = (r'$\Theta_{m}^{num}$')#,r'$\delta^{ana}$')
py.loglog(am4,abs(vmm4),color='r')
#py.loglog(xm4,abs(dea1m4),color='b',linestyle='dashed')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('comparison_Theta_m.pdf')
py.close()

# Theta_de

fig = py.figure()
labels = (r'$\Theta_{de}^{num}$')#,r'$V^{ana}$')
#py.loglog(am4,abs(vdem4),color='r')
py.plot(am4,abs(vdem4),color='r')
#py.loglog(xm4,abs(vdea1m4),color='b',linestyle='dashed')
py.xscale('log')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('comparison_Theta_de.pdf')
py.close()

# Z

fig = py.figure()
labels = (r'$Z^{num}$')#,r'$V^{ana}$')
py.plot(am4,-Z/Z[-1],color='r')
#py.loglog(xm4,abs(vdea1m4),color='b',linestyle='dashed')
py.xscale('log')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('comparison_Z.pdf')
py.close()

exit()



