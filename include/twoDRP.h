#include <iostream>
#include <string>
#include <omp.h>    // ²¢ÐÐ

#define _USE_MATH_DEFINES
#include <math.h>

#include <time.h>

#include <stdio.h>

#include "settings.h"
#include "CArray.h"
#include "invEQ2.h"
#include "2DRPFunctions.h"
#include "RedEquil.h"
#include "Scheme.h"


void twoDRP_EQMOM(string init_setting, string compute_setting);
void twoDRP_DVDDVM(string init_setting, string compute_setting);
void twoDRP_DVM(string init_setting, string compute_setting);

void twoDRP_DVDDVM_gh(string init_setting, string compute_setting);
void twoDRP_DVM_gh(string init_setting, string compute_setting);
