#include <iostream>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

#include "../include/RedEquil.h"

using namespace std;
double RedEquil_DVM(double rho, double u, double v, double E, double* alpha, double* rhoeq)
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
    double* B_inv = new double[4 * 4]{};
    double* B_inv_old = new double[4 * 4]{};
    for (int i = 0;i < 4;++i)
    {
        B_inv[i * 4 + i] = 1;
    }
    double* a_old = new double[4]{};
    memcpy(a_old, alpha, sizeof(double) * 4);
    double* a_new = new double[4]{};
    double* dlt_a = new double[4]{};    // ǰ������x�Ĳ�
    double* gradJ = new double[4]{};
    double* g_old = new double[4]{};
    double* dlt_g = new double[4]{};    // ǰ������gradJ�Ĳ�
    double* dir_t = new double[4]{};    // ��ǰ�����½�����
    double J, J_new{}, n2gradJ, expabc; // һЩ������Ҫ�õ���ʱ����
    double err = 1;
    int iter = 0;
    double lr = 1;
#pragma endregion

    // ��һ�μ����J��gradJ
    double d = pow(2 * M_PI / -a_old[3], L / 2.);
    for (int k = M * N - 1;k >= 0;--k)
    {
        expabc = d * exp(a_old[0] + a_old[1] * uk[k] + a_old[2] * vk[k] + 0.5 * a_old[3] * U2k[k]);
        gradJ[0] += expabc;
        gradJ[1] += expabc * uk[k];
        gradJ[2] += expabc * vk[k];
        gradJ[3] += expabc * 0.5 * U2k[k];
    }
    J = gradJ[0] - (a_old[0] + a_old[1] * u + a_old[2] * v + a_old[3] * E);
    // ���gradJ
    gradJ[3] = gradJ[3] - L * 0.5 * gradJ[0] / a_old[3] - E;
    gradJ[0] = gradJ[0] - 1;
    gradJ[1] = gradJ[1] - u;
    gradJ[2] = gradJ[2] - v;
    // ��ʼ����
    while (err > REDEQUIL_ERROR)
    {
        iter++;// ����

        // �����ⷽ��
        MtimesV(B_inv, gradJ, 4, dir_t);

        n2gradJ = pow(dir_t[0], 2) + pow(dir_t[1], 2) + pow(dir_t[2], 2) + pow(dir_t[3], 2);

#pragma region ѧϰ�ʵ���
        // ѧϰ�ʵ���
        bool needloop = true;
        lr = 1;
        while (needloop && err > REDEQUIL_ERROR)
        {
            iter++;
            needloop = false;   // �����Ѿ�����Ҫ��ѭ����
            for (int j = 0;j < 4;++j)
            {
                a_new[j] = a_old[j] - lr * dir_t[j];
            }

            if (L > 0 && a_new[3] >= 0)
            {
                needloop = true;
                lr = lr * rGA[0];
                continue;
            }

            // �����J_new
            J_new = -(a_new[0] + a_new[1] * u + a_new[2] * v + a_new[3] * E);
            d = pow(2 * M_PI / -a_new[3], L / 2.);
            for (int k = M * N - 1;k >= 0;--k)
            {
                J_new += d * exp(a_new[0] + a_new[1] * uk[k] + a_new[2] * vk[k] + 0.5 * a_new[3] * U2k[k]);
            }
            if (0 > J - J_new)
            {
                needloop = true;
                lr = lr * rGA[0];
            }
           /* if (lr * mGA[0] * n2gradJ > J - J_new)
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
            //cout << lr << endl;
        }
#pragma endregion

        // ���һ�ε���
        // ����ɵ�ֵ
        swapPtr(&a_old, &a_new);
        swapPtr(&g_old, &gradJ);
        // ��ʼ���µ�gradJ��J
        memset(gradJ, 0, sizeof(double) * 4);
        d = pow(2 * M_PI / -a_old[3], L / 2.);
        for (int k = M * N - 1;k >= 0;--k)
        {
            expabc = d * exp(a_old[0] + a_old[1] * uk[k] + a_old[2] * vk[k] + 0.5 * a_old[3] * U2k[k]);
            gradJ[0] += expabc;
            gradJ[1] += expabc * uk[k];
            gradJ[2] += expabc * vk[k];
            gradJ[3] += expabc * 0.5 * U2k[k];
        }
        J = gradJ[0] - (a_old[0] + a_old[1] * u + a_old[2] * v + a_old[3] * E);
        // ���gradJ
        gradJ[3] = gradJ[3] - L * 0.5 * gradJ[0] / a_old[3] - E;
        gradJ[0] = gradJ[0] - 1;
        gradJ[1] = gradJ[1] - u;
        gradJ[2] = gradJ[2] - v;
        // BFGS���Ĳ��֣�����B_inverse
        for (int i = 0;i < 4;++i)
        {
            dlt_g[i] = gradJ[i] - g_old[i];
            dlt_a[i] = a_old[i] - a_new[i]; // ע���������ڽ�������������old-new
        }
        swapPtr(&B_inv, &B_inv_old);
        double x_dot_g_inv = 1 / VdotV(dlt_a, dlt_g, 4);
       
        for (int i = 0;i < 4;++i)
        {
            for (int j = 0;j < 4;++j)
            {
                B_inv[i * 4 + j] = dlt_a[i] * dlt_a[j] * x_dot_g_inv;
                for (int m = 0;m < 4;++m)
                {
                    double im = (double)(i == m) - dlt_a[i] * dlt_g[m] * x_dot_g_inv;
                    for (int n = 0;n < 4;++n)
                    {
                        B_inv[i * 4 + j] += im * B_inv_old[m * 4 + n]
                            * ((double)(j == n) - dlt_a[j] * dlt_g[n] * x_dot_g_inv);
                    }
                }
            }
        }
    }

    memcpy(alpha, a_old, sizeof(double) * 4);

    double sigeq2 = -1 / a_old[3];
    // ����rhoeq��ֵ
    d = pow(2 * M_PI / -a_old[3], L / 2.);
    double dmass = d* rho / (J_new + (a_old[0] + a_old[1] * u + a_old[2] * v + a_old[3] * E));
    for (int k = M * N - 1;k >= 0;--k)
    {
        rhoeq[k] = exp(a_old[0] + a_old[1] * uk[k] + a_old[2] * vk[k] + 0.5 * a_old[3] * U2k[k]) * dmass;
    }

    delete[] B_inv;
    delete[] B_inv_old;
    delete[] a_old;
    delete[] a_new;
    delete[] dlt_a;
    delete[] gradJ;
    delete[] g_old;
    delete[] dlt_g;
    delete[] dir_t;


  /*  u = rho * u;
    v = rho * v;
    E = rho * E;
    for (int k = 0;k < N*M;++k)
    {
        rho -= rhoeq[k];
        u -= rhoeq[k] * uk[k];
        v -= rhoeq[k] * vk[k];
        E -= rhoeq[k] * (U2k[k] + L * sigeq2) * 0.5;
    }
    cout << abs(rho) + abs(u) + abs(v) + abs(E) << endl;*/

    return sigeq2;
}
