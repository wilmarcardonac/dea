import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py


#######################################################################################
#######################################################################################
# Loading data for numerical and analytical solutions
#######################################################################################
#######################################################################################

am4,dmm4,dem4,vmm4,vdem4 = np.loadtxt('../output/numerical_solution.txt',unpack=True,usecols=[0,1,2,3,4])

xm4,dea1m4,vdea1m4,dea2m4,vdea2m4,dm,Vm = np.loadtxt('../output/analytical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5,6])

# \delta_m

fig = py.figure()
labels = (r'$\delta_{m}^{num}$',r'$\delta_{m}^{ana}$')
py.loglog(am4,abs(dmm4),color='r')
py.loglog(xm4,abs(dm),color='b',linestyle='dashed')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('comparison_delta_m.pdf')
py.close()

# V_m

fig = py.figure()
labels = (r'$V_{m}^{num}$',r'$V_{m}^{ana}$')
py.loglog(am4,abs(vmm4),color='r')
py.loglog(xm4,abs(Vm),color='b',linestyle='dashed')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('comparison_V_m.pdf')
py.close()

# \delta

fig = py.figure()
labels = (r'$\delta^{num}$',r'$\delta^{sup}$',r'$\delta^{sub}$')
py.loglog(am4,abs(dem4),color='r')
py.loglog(xm4,abs(dea1m4),color='b',linestyle='dashed')
py.loglog(xm4,abs(dea2m4),color='g',linestyle='dotted')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('comparison_delta.pdf')
py.close()

# V

fig = py.figure()
labels = (r'$V^{num}$',r'$V^{sup}$',r'$V^{sub}$')
py.loglog(am4,abs(vdem4),color='r')
py.loglog(xm4,abs(vdea1m4),color='b',linestyle='dashed')
py.loglog(xm4,abs(vdea2m4),color='g',linestyle='dotted')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('comparison_V.pdf')
py.close()

exit()



