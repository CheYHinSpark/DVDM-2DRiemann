#ifndef REDEQUIL_H
#define REDEQUIL_H

#include <cstring>

#include "settings.h"
#include "2DRPFunctions.h"

using namespace std;

void RedEquil_DVD(double rho, double u, double v, double E,
	double* alpha, double* rhoeq, double* ueq, double* sig2eq);

// ���L>0�����ص�rhoeq��g
double RedEquil_DVDDVM(double rho, double u, double E, double* rhoeq);

// ���L>0�����ص�rhoeq��g
double RedEquil_DVM(double rho, double u, double v, double E,
	double* alpha, double* rhoeq);

#endif