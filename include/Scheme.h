#include "CArray.h"

void D_wind_up(int X, int Y, double dt);

void DUGKS(int X, int Y, double dt);

void DUGKS_gh(int X, int Y, double dt);

double DUGKS_9_interpol(double fmm, double fm0, double fmp,
    double f0m, double f00, double f0p,
    double fpm, double fp0, double fpp, double r, double s);

double Van_Leer_Limiter(double s_1, double s_2);