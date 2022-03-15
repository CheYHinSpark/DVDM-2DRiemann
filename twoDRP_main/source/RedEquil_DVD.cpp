#include <iostream>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

#include "../include/RedEquil.h"

using namespace std;
void RedEquil_DVD(double rho, double u, double v, double E,
    double* alpha, double* rhoeq, double* ueq, double* sig2eq)
{
    /*
    ����gamma�ķ�ʽ��gamma��һ����N����������ʾ�����Ƕȡ�
    ʹ���Ż�alpha���㷨
    ʹ��BFGS��ţ�ٷ�����Goldstein - Armijo ruleѡ�񲽳�
    // 2022/03/04������ʹ��GA���������ͣ�����������ѣ��������ã�
    // ʹ�ü򵥵Ĳ������򣺴�1��ʼ����Jû���½��������̲���
    */

    // GA�Ĳ���
    double mGA[2] = { 0.25, 0.75 };
    double rGA[2] = { 0.5, 3.6 };   // ���̺ͷŴ󲽳���ϵ����ʵ���ϵڶ�������Ч���ǷŴ�rGA1* rGA2��

#pragma region Check realizability
   //������
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
    double J, J_new{}, n2gradJ, blk, expblk;
    double err = 1;
    int iter = 0;
    double lr = 1;
#pragma endregion

    // ��һ�μ����J��gradJ
    // ����ϵ��
    double d = pow(2. * M_PI / -a_old[3], (L + 1) / 2.) * exp(a_old[0]);
    for (int i = 0;i < N;++i)
    {
        blk = a_old[1] * x_n[i] + a_old[2] * y_n[i];
        expblk = d * exp(-blk * blk * 0.5 / a_old[3]);
        // ����gradJ��ע��J��gradJ��һ�����Ĺ�ϵ
        gradJ[0] += expblk;
        blk /= a_old[3];
        gradJ[1] += -blk * x_n[i] * expblk;
        gradJ[2] += -blk * y_n[i] * expblk;
        gradJ[3] += pow(blk, 2) * expblk;
    }
    // gradJ4����
    gradJ[3] = 0.5 * (gradJ[3] - (L + 1) * gradJ[0] / a_old[3]);
    // J
    J = gradJ[0] - (a_old[0] + a_old[1] * u + a_old[2] * v + a_old[3] * E);
    // ���gradJ
    gradJ[0] = gradJ[0] - 1.;
    gradJ[1] = gradJ[1] - u;
    gradJ[2] = gradJ[2] - v;
    gradJ[3] = gradJ[3] - E;
    // ��ʼ����
    while (err > REDEQUIL_ERROR)
    {
        iter++;// ����

        // �����ⷽ��
        MtimesV(B_inv, gradJ, 4, dir_t);

        n2gradJ = pow(gradJ[0], 2) + pow(gradJ[1], 2) + pow(gradJ[2], 2) + pow(gradJ[3], 2);

#pragma region ѧϰ�ʵ���
        // ѧϰ�ʵ���
        lr = 1;
        bool needloop = true;
        bool if_out = false;
        while ((needloop || if_out) && err > REDEQUIL_ERROR)
        {
            needloop = false;   // �����Ѿ�����Ҫ��ѭ����
            lr = lr * rGA[0];
            for (int j = 0;j < 4;++j)
            {
                a_new[j] = a_old[j] - lr * dir_t[j];
            }
            if (a_new[3] >= 0)
            {
                if_out = true;
                continue;
            }
            // �����J_new
            J_new = 0;
            for (int i = 0;i < N;++i)
            {
                J_new += exp(pow(a_new[1] * x_n[i] + a_new[2] * y_n[i], 2) * 0.5 / -a_new[3]);
            }
            J_new = pow(2 * M_PI / -a_new[3], (L + 1) / 2.) * exp(a_new[0]) * J_new
                - (a_new[0] + a_new[1] * u + a_new[2] * v + a_new[3] * E);
            if (0 > J - J_new)
            {
                needloop = true;
            }
            //if (lr * mGA[0] * n2gradJ > J - J_new)
            //{
            //    needloop = true;
            //}
            //else if (lr * mGA[1] * n2gradJ < J - J_new && !if_out)
            //{
            //    // ����Ҫ�ϴβ��ǳ��磬�ŷŴ�lr��Ҳ��Ϊ�˷�ֹ����ѭ��
            //    needloop = true;
            //    lr = lr * rGA[1];
            //}
            if_out = false; // �����һȦ��û�г���
            err = lr * sqrt(n2gradJ);
        }
#pragma endregion

        // ���һ�ε���
        // ����ɵ�ֵ
        swapPtr(&a_old, &a_new);
        swapPtr(&g_old, &gradJ);
        // �����J��gradJ
        memset(gradJ, 0, sizeof(double) * 4);
        // ����ϵ��
        d = pow(2. * M_PI / -a_old[3], (L + 1) / 2.) * exp(a_old[0]);
        for (int i = 0;i < N;++i)
        {
            blk = a_old[1] * x_n[i] + a_old[2] * y_n[i];
            expblk = d * exp(-blk * blk * 0.5 / a_old[3]);
            blk /= a_old[3];
            // ����gradJ��ע��J��gradJ��һ�����Ĺ�ϵ
            gradJ[0] += expblk;
            gradJ[1] += -blk * x_n[i] * expblk;
            gradJ[2] += -blk * y_n[i] * expblk;
            gradJ[3] += pow(blk, 2) * expblk;
        }
        // gradJ4����
        gradJ[3] = 0.5 * (gradJ[3] - (L + 1) * gradJ[0] / a_old[3]);
        // J
        J = gradJ[0] - (a_old[0] + a_old[1] * u + a_old[2] * v + a_old[3] * E);
        // ���gradJ
        gradJ[0] = gradJ[0] - 1.;
        gradJ[1] = gradJ[1] - u;
        gradJ[2] = gradJ[2] - v;
        gradJ[3] = gradJ[3] - E;
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

    // ����rhoeq�ĳ���ֵ
    *sig2eq = -1 / a_old[3];
    d = pow(2 * M_PI * *sig2eq, (L + 1) / 2.);
    double mass = rho / (J_new + (a_old[0] + a_old[1] * u + a_old[2] * v + a_old[3] * E));
    for (int i = 0;i < N;++i)
    {
        blk = x_n[i] * a_old[1] + y_n[i] * a_old[2];
        ueq[i] = blk * *sig2eq;
        rhoeq[i] = d * exp(a_old[0] + blk * blk * 0.5 * *sig2eq) * mass;
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


    //u = rho * u;
    //v = rho * v;
    //E = rho * E;
    //for (int n = 0;n < N;++n)
    //{
    //    rho -= rhoeq[n];
    //    u -= rhoeq[n] * ueq[n] * x_n[n];
    //    v -= rhoeq[n] * ueq[n] * y_n[n];
    //    E -= rhoeq[n] * (ueq[n] * ueq[n] + (L + 1) * *sig2eq) * 0.5;
    //}
    //cout << abs(rho) + abs(u) + abs(v) + abs(E) << endl;

    return;
}
