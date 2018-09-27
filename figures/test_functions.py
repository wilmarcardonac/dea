import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

a,f1,f2,f3,f4,f5 = np.loadtxt('../output/functions.txt',unpack=True,usecols=[0,1,2,3,4,5])

fig = py.figure()

# anisotorpic stress

labels = ('MD','sub-MD','sub-MD-2')

py.semilogx(a,f1,linestyle='dotted',color='b')
py.semilogx(a,f2,linestyle='-.',color='g')
py.semilogx(a,f3,linestyle='dashed',color='r')
#py.semilogx(a,0.*f1,color='black')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.ylim(-1.e-8,1.e-8)
py.savefig('anisotropic_stress_comparison.pdf')

py.close()

fig = py.figure()

# sound speed squared

labels = (r'$c_s^{2,MD}$',r'$c_s^{2,sub-MD}$')

py.semilogx(a,f4,linestyle='dotted',color='b')
py.semilogx(a,f5,linestyle='-.',color='g')
py.xlabel(r'$a$',fontsize='large')
py.legend(labels,loc=0)
py.savefig('sound_speed_comparison.pdf')

py.close()

exit()

