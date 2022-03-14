#include <iostream>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

#include "../include/RedEquil.h"

using namespace std;
double RedEquil_DVDDVM(double rho, double u, double E, double* rhoeq)
{

    // ʹ��BFGS��ţ�ٷ�����Goldstein - Armijo ruleѡ�񲽳�
    // 2022/03/04������ʹ��GA���������ͣ�����������ѣ��������ã�
    // ʹ�ü򵥵Ĳ������򣺴�1��ʼ����Jû���½��������̲���
    // GA�Ĳ���
    double mGA[2] = { 0.25, 0.75 };
    double rGA[2] = { 0.5, 1.5 };   // ���̺ͷŴ󲽳���ϵ����ʵ���ϵڶ�������Ч���ǷŴ�rGA1* rGA2��

#pragma region Check realizability
/*
* ��ȥ����
*/
#pragma endregion

#pragma region ����׼��
    double* B_inv = new double[3 * 3]{};
    double* B_inv_old = new double[3 * 3]{};
    for (int i = 0;i < 3;++i)
    {
        B_inv[i * 3 + i] = 1;
    }
    double rt = 2 * E - u * u;
    double* a_old = new double[3]{ -0.5 * log(2 * M_PI * rt) - u * u / 2 / rt,u / rt,-1 / rt };
    double* a_new = new double[3]{};
    double* dlt_a = new double[3]{};    // ǰ������x�Ĳ�
    double* gradJ = new double[3]{};
    double* g_old = new double[3]{};
    double* dlt_g = new double[3]{};    // ǰ������gradJ�Ĳ�
    double* dir_t = new double[3]{};    // ��ǰ�����½�����
    double J, J_new{}, n2gradJ, expabc; // һЩ������Ҫ�õ���ʱ����
    double err = 1;
    int iter = 0;
    double lr = 1;
#pragma endregion

    // ��һ�μ����J��gradJ
    double d = pow(2 * M_PI / -a_old[2], L / 2.);
    for (int m = 0;m < M;++m)
    {
        expabc = d * exp(a_old[0] + (a_old[1] + 0.5 * a_old[2] * dvm[m]) * dvm[m]);
        gradJ[0] += expabc;
        gradJ[1] += expabc * dvm[m];
        gradJ[2] += expabc * 0.5 * dvm[m] * dvm[m];
    }
    J = gradJ[0] - (a_old[0] + a_old[1] * u + a_old[2] * E);
    // ���gradJ
    gradJ[2] = gradJ[2] - L * 0.5 * gradJ[0] / a_old[2] - E;
    gradJ[0] = gradJ[0] - 1;
    gradJ[1] = gradJ[1] - u;
    // ��ʼ����
    while (err > REDEQUIL_ERROR)
    {
        iter++;// ����

        // �����ⷽ��
        MtimesV(B_inv, gradJ, 3, dir_t);

        n2gradJ = pow(gradJ[0], 2) + pow(gradJ[1], 2) + pow(gradJ[2], 2);

#pragma region ѧϰ�ʵ���
        // ѧϰ�ʵ���
        lr = 1;
        bool needloop = true;
        while (needloop && err > REDEQUIL_ERROR)
        {
            iter++;
            needloop = false;   // �����Ѿ�����Ҫ��ѭ����
            for (int j = 0;j < 3;++j)
            {
                a_new[j] = a_old[j] - lr * dir_t[j];
            }

            if (L > 0 && a_new[2] >= 0)
            {
                needloop = true;
                lr = lr * rGA[0];
                continue;
            }

            // �����J_new
            J_new = -(a_new[0] + a_new[1] * u + a_new[2] * E);
            d = pow(2 * M_PI / -a_new[2], L / 2.);
            for (int m = 0;m < M;++m)
            {
                J_new += d * exp(a_new[0] + (a_new[1] + 0.5 * a_new[2] * dvm[m]) * dvm[m]);
            }

            if (0 > J - J_new)
            {
                needloop = true;
                lr = lr * rGA[0];
            }
          /*  if (lr * mGA[0] * n2gradJ > J - J_new)
            {
                needloop = true;
                lr = lr * rGA[0];
            }
            else if (lr * mGA[1] * n2gradJ < J - J_new)
            {
                needloop = true;
                lr = lr * rGA[1];
            }*/
            err = lr * sqrt(n2gradJ);
        }
#pragma endregion

        // ���һ�ε���
        // ����ɵ�ֵ
        swapPtr(&a_old, &a_new);
        swapPtr(&g_old, &gradJ);
        // ��ʼ���µ�gradJ��J
        memset(gradJ, 0, sizeof(double) * 3);
        d = pow(2 * M_PI / -a_old[2], L / 2.);
        for (int m = 0;m < M;++m)
        {
            expabc = d * exp(a_old[0] + (a_old[1] + 0.5 * a_old[2] * dvm[m]) * dvm[m]);
            gradJ[0] += expabc;
            gradJ[1] += expabc * dvm[m];
            gradJ[2] += expabc * 0.5 * dvm[m] * dvm[m];
        }
        J = gradJ[0] - (a_old[0] + a_old[1] * u + a_old[2] * E);
        // ���gradJ
        gradJ[2] = gradJ[2] - L * 0.5 * gradJ[0] / a_old[2] - E;
        gradJ[0] = gradJ[0] - 1;
        gradJ[1] = gradJ[1] - u;
        // BFGS���Ĳ��֣�����B_inverse
        for (int i = 0;i < 3;++i)
        {
            dlt_g[i] = gradJ[i] - g_old[i];
            dlt_a[i] = a_old[i] - a_new[i]; // ע���������ڽ�������������old-new
        }
        swapPtr(&B_inv, &B_inv_old);
        double x_dot_g_inv = 1 / VdotV(dlt_a, dlt_g, 3);

        for (int i = 0;i < 3;++i)
        {
            for (int j = 0;j < 3;++j)
            {
                B_inv[i * 3 + j] = dlt_a[i] * dlt_a[j] * x_dot_g_inv;
                for (int m = 0;m < 3;++m)
                {
                    double im = (double)(i == m) - dlt_a[i] * dlt_g[m] * x_dot_g_inv;
                    for (int n = 0;n < 3;++n)
                    {
                        B_inv[i * 3 + j] += im * B_inv_old[m * 3 + n]
                            * ((double)(j == n) - dlt_a[j] * dlt_g[n] * x_dot_g_inv);
                    }
                }
            }
        }
    }

    double sigeq2 = -1 / a_old[2];
    // ����rhoeq��ֵ
    d = pow(2 * M_PI / -a_old[2], L / 2.);
    double dmass = d * rho / (J_new + (a_old[0] + a_old[1] * u + a_old[2] * E));
    for (int m = 0;m < M;++m)
    {
        rhoeq[m] = exp(a_old[0] + (a_old[1] + 0.5 * a_old[2] * dvm[m]) * dvm[m]) * dmass;
    }

#pragma endregion

    delete[] B_inv;
    delete[] B_inv_old;
    delete[] a_old;
    delete[] a_new;
    delete[] dlt_a;
    delete[] gradJ;
    delete[] g_old;
    delete[] dlt_g;
    delete[] dir_t;

    //cout << iter << endl;


   /* u = rho * u;
    E = rho * E;
    for (int n = 0;n < M;++n)
    {
        rho -= rhoeq[n];
        u -= rhoeq[n] * dvm[n];
        E -= rhoeq[n] * (dvm[n] * dvm[n] + L * sigeq2) * 0.5;
    }
    cout << abs(rho) + abs(u) + abs(E) << endl;*/
    
    return sigeq2;
}

