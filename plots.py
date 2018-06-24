import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py


#######################################################################################
#######################################################################################
# Loading data for numerical and analytical solutions
#######################################################################################
#######################################################################################

# Model 1

#am1,dmm1,dem1,vmm1,vdem1,phim1,psim1 = np.loadtxt('numericalmodel1.dat',unpack=True,usecols=[0,1,2,3,4,5,6])

#xm1,dea1m1,vdea1m1,dea2m1,vdea2m1 = np.loadtxt('analyticalmodel1.dat',unpack=True,usecols=[0,1,2,3,4])

# Model 2

#am2,dmm2,dem2,vmm2,vdem2,phim2,psim2 = np.loadtxt('numericalmodel2.dat',unpack=True,usecols=[0,1,2,3,4,5,6])

#xm2,dea1m2,vdea1m2,dea2m2,vdea2m2 = np.loadtxt('analyticalmodel2.dat',unpack=True,usecols=[0,1,2,3,4])

# Model 3

#am3,dmm3,dem3,vmm3,vdem3,phim3,psim3 = np.loadtxt('numericalmodel3.dat',unpack=True,usecols=[0,1,2,3,4,5,6])

#xm3,dea1m3,vdea1m3,dea2m3,vdea2m3 = np.loadtxt('analyticalmodel3.dat',unpack=True,usecols=[0,1,2,3,4])

# Model 4

am4,dmm4,dem4,vmm4,vdem4,phim4,psim4 = np.loadtxt('./output/numerical_solution.txt',unpack=True,usecols=[0,1,2,3,4,5,6])

xm4,dea1m4,vdea1m4,dea2m4,vdea2m4 = np.loadtxt('./output/analytical_solution.txt',unpack=True,usecols=[0,1,2,3,4])


#######################################################################################
#######################################################################################
# Plotting dark energy perturbations
#######################################################################################
#######################################################################################

fig = py.figure()
#fig.subplots_adjust(bottom=0.1, left=0.06, top=0.85, right=0.97,hspace=0.3)

#py.subplot(3,1,1)

#######################################################################################
# Model 1
#######################################################################################

#py.subplot(3,1,1)
#labels = (r'$\delta_{de}$',r'$\delta_{de}^{sup-hor}$',r'$\delta_{de}^{sub-sound}$',r'$v_{de}$',r'$v_{de}^{sup-hor}$',r'$v_{de}^{sub-sound}$')

#py.subplot(2,2,1)
#line1, = py.loglog(am1,abs(dem1),label=r'$\delta_{de}$')
#line2, = py.loglog(xm1,abs(dea1m1),label=r'$\delta_{de}^{sup-hor}$',linestyle='--')
#line3, = py.loglog(xm1,abs(dea2m1),label=r'$\delta_{de}^{sub-sound}$',linestyle='--')
#line4, = py.loglog(am1,abs(vdem1),label=r'$v_{de}$')
#line5, = py.loglog(xm1,abs(vdea1m1),label=r'$v_{de}^{sup-hor}$',linestyle='-.')
#line6, = py.loglog(xm1,abs(vdea2m1),label=r'$v_{de}^{sub-sound}$',linestyle='-.',color=(0,0,0))

# Scale factor at which the mode crosses the horizon

#py.vlines(7.85*10**(-5),10**(-6),10**3,color='k',linestyles='dotted')

# Scale factor at which the mode crosses the effective sound horizon 

#py.vlines(7.85*10**(-4),10**(-6),10**3,color='k',linestyles='dashed')

# Parameters

#py.text(10**(-3),10**(-4),r'$e_{\pi}=10^{-1},\, f_{\pi} = 0,$',fontsize='large')
#py.text(10**(-2),10**(-5),r'$g_{\pi} = 0$',fontsize='large')

# Label for x-axis

#py.xlabel(r'$a$',fontsize='large')

#py.legend(bbox_to_anchor=(0.38,0.35),loc=2,ncol=2)#,prop={'size':16})

#######################################################################################
# Model 2
#######################################################################################

#py.subplot(2,2,2)
#py.loglog(am2,abs(dem2))
#py.loglog(xm2,abs(dea1m2),linestyle='dashed')
#py.loglog(xm2,abs(dea2m2),linestyle='dashed')
#py.loglog(am2,abs(vdem2))
#py.loglog(xm2,abs(vdea1m2),linestyle='-.')
#py.loglog(xm2,abs(vdea2m2),linestyle='-.',color=(0,0,0))

# Scale factor at which the mode crosses the horizon

#py.vlines(7.85*10**(-5),10**(-3),10**2,color='k',linestyles='dotted')

# Scale factor at which the mode crosses the effective sound horizon 

#py.vlines(2.36*10**(-3),10**(-3),10**2,color='k',linestyles='dashed')

# Parameters

#py.text(4*10**(-3),2*10,r'$e_{\pi} = 0,\, f_{\pi} = 10^{-1},$',fontsize='large')
#py.text(2*10**(-2),6,r'$g_{\pi} = 10^{-2}$',fontsize='large')

# Label for x-axis

#py.xlabel(r'$a$',fontsize='large')

#######################################################################################
# Model 3
#######################################################################################

#py.subplot(2,2,3)
#py.loglog(am3,abs(dem3))
#py.loglog(xm3,abs(dea1m3),linestyle='dashed')
#py.loglog(xm3,abs(dea2m3),linestyle='dashed')
#py.loglog(am3,abs(vdem3))
#py.loglog(xm3,abs(vdea1m3),linestyle='-.')
#py.loglog(xm3,abs(vdea2m3),linestyle='-.',color=(0,0,0))

# Scale factor at which the mode crosses the horizon

#py.vlines(7.85*10**(-5),10**(-3),10,color='k',linestyles='dotted')

# Scale factor at which the mode crosses the effective sound horizon 

#py.vlines(2.36*10**(-3),10**(-3),10,color='k',linestyles='dashed')

# Title

#py.text(2*10**(-4),3,r'$e_{\pi} = 0$',fontsize='large')
#py.text(3*10**(-3),3,r'$f_{\pi} = 10^{-1},\, g_{\pi} = 1$',fontsize='large')

# Label for x-axis

#py.xlabel(r'$a$',fontsize='large')

#######################################################################################
# Model 4
#######################################################################################

#py.subplot(2,2,4)
py.loglog(am4,abs(dem4))
py.loglog(xm4,abs(dea1m4),linestyle='dashed')
py.loglog(xm4,abs(dea2m4),linestyle='dashed')
py.loglog(am4,abs(vdem4))
py.loglog(xm4,abs(vdea1m4),linestyle='-.')
py.loglog(xm4,abs(vdea2m4),linestyle='-.',color=(0,0,0))

# Scale factor at which the mode crosses the horizon

#py.vlines(7.85*10**(-5),10**(-5),10,color='k',linestyles='dotted')

# Scale factor at which the mode crosses the effective sound horizon 

#py.vlines(2.36*10**(-3),10**(-5),10,color='k',linestyles='dashed')

# Title 

#py.text(3*10**(-3),8*10**(-4),r'$e_{\pi}=0,\, f_{\pi} = 10^{-1},$',fontsize='large')
#py.text(10**(-2),10**(-4),r'$g_{\pi} = 10$',fontsize='large')

# Label for x-axis

py.xlabel(r'$a$',fontsize='large')

#######################################################################################
# Placing the legend 
#######################################################################################
labels = (r'$\delta_{de}$',r'$V_{de}$')
py.legend(labels,bbox_to_anchor=(-1.2,2.4,2.2,.102),loc=3,ncol=6,mode='expand',borderaxespad=0.)

#######################################################################################
# Saving the figure 
#######################################################################################

#py.savefig('example.pdf',bbox_inches='tight')

#######################################################################################
# Show the figure
#######################################################################################

#py.tight_layout()
py.savefig('./output/test.pdf')
py.close()
exit()

