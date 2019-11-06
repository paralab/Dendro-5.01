###############################################################################
# Apr.14.2019
# Script to generate rhs equations for quadratic gravity
# Here, we use diagonalize GH system
###############################################################################
import dendro
from sympy import *
from sympy.physics.vector.vector import Vector
from sympy.printing.dot import dotprint
###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')

# declare variables
RccSc   = dendro.scalar("RccSc", "[pp]")

Vaux  = dendro.vec3("Vaux", "[pp]")

gm  = dendro.sym_3x3("gm", "[pp]")
Rcct  = dendro.sym_3x3("Rcct", "[pp]")

# Lie derivative weight
weight = -Rational(2,3)
weight_Gt = Rational(2,3)

# specify the functions for computing first and second derivatives
d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction
kod = dendro.set_kreiss_oliger_dissipation('kograd')

d2 = dendro.d2

#f = Function('f')

# generate metric related quantities
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)
###################################################################
# evolution equations
###################################################################



###################################################################
# generate code
###################################################################

outs = [RccSc_rhs, gm_rhs, Rcct_rhs, Vaux_rhs, h3id_rhs]
vnames = ['RccSc_rhs', 'gm_rhs','Rcct_rhs','Vaux_rhs','h3id_rhs' ]
#dendro.generate_debug(outs, vnames)
#dendro.generate(outs, vnames, '[pp]')
#numVars=len(outs)
#for i in range(0,numVars):
#    dendro.generate_separate([outs[i]],[vnames[i]],'[pp]')

#dendro.generate_separate([Gt_rhs],['Gt_rhs'],'[pp]')
#dendro.generate([CalGt, Gt_rhs_s1, Gt_rhs_s2, Gt_rhs_s3, Gt_rhs_s4, Gt_rhs_s5, Gt_rhs_s6, Gt_rhs_s7, Gt_rhs,B_rhs],['CalGt', 'Gt_rhs_s1', 'Gt_rhs_s2', 'Gt_rhs_s3', 'Gt_rhs_s4', 'Gt_rhs_s5', 'Gt_rhs_s6', 'Gt_rhs_s7', 'Gt_rhs', 'B_rhs'],'[pp]')
