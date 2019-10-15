#ifndef _DERVS_H
#define _DERVS_H

#include <cmath>
#include "TreeNode.h"

#define IDX(i, j, k) ((i) + nx * ((j) + ny * (k)))

void deriv42_z(double *const Dzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);
void deriv42_y(double *const Dyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_x(double *const Dxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);

void deriv42_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);
void deriv42_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);

void ko_deriv42_z(double *const Du, const double *const u, const double dz, const unsigned *sz, unsigned bflag);
void ko_deriv42_y(double *const Du, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void ko_deriv42_x(double *const Du, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);

void deriv22_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);
void deriv22_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv22_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);


#endif
