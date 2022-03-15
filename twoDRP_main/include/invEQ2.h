#ifndef INVEQ2_H
#define INVEQ2_H

extern const double M_PI_3; //, Rg, CFL;
extern const int N1, N2, N_grid;			// Spatial grid
extern const int Ndir, Nnod, Nmom;			// Numbers of directions (3), nodes (2) and moments (5)

/* Moment inversion algorithm for 2-node Gaussian-EQMOM */
// mpara里面的顺序为：w1, u1, w2, u2, sig^2
void invEQ2(double m0, double m1, double m2, double m3, double m4, double mpara[]);

#endif