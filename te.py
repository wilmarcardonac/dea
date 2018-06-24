import sh
import os

#######################################################################################
# Compiling the fortran code to compute dark energy and matter perturbations
#######################################################################################

sh.make('clean')

sh.make('def')

os.system('./depmd')





