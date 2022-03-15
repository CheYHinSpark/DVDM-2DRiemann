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
    输入gamma的方式，gamma是一个长N的向量，表示各个角度。
    使用优化alpha的算法
    使用BFGS拟牛顿法，以Goldstein - Armijo rule选择步长
    // 2022/03/04：由于使用GA法则出现了停不下来的困难，最终弃用，
    // 使用简单的步长规则：从1开始，若J没有下降，则缩短步长
    */

    // GA的参数
    double mGA[2] = { 0.25, 0.75 };
    double rGA[2] = { 0.5, 3.6 };   // 缩短和放大步长的系数，实际上第二个数的效果是放大到rGA1* rGA2倍

#pragma region Check realizability
   //不做了
#pragma endregion


#pragma region 迭代准备
    double* B_inv = new double[4 * 4]{};
    double* B_inv_old = new double[4 * 4]{};
    for (int i = 0;i < 4;++i)
    {
        B_inv[i * 4 + i] = 1;
    }
    double* a_old = new double[4]{};
    memcpy(a_old, alpha, sizeof(double) * 4);
    double* a_new = new double[4]{};
    double* dlt_a = new double[4]{};    // 前后两次x的差
    double* gradJ = new double[4]{};
    double* g_old = new double[4]{};
    double* dlt_g = new double[4]{};    // 前后两次gradJ的差
    double* dir_t = new double[4]{};    // 当前步的下降方向
    double J, J_new{}, n2gradJ, blk, expblk;
    double err = 1;
    int iter = 0;
    double lr = 1;
#pragma endregion

    // 第一次计算出J和gradJ
    // 存下系数
    double d = pow(2. * M_PI / -a_old[3], (L + 1) / 2.) * exp(a_old[0]);
    for (int i = 0;i < N;++i)
    {
        blk = a_old[1] * x_n[i] + a_old[2] * y_n[i];
        expblk = d * exp(-blk * blk * 0.5 / a_old[3]);
        // 计算gradJ，注意J和gradJ第一个数的关系
        gradJ[0] += expblk;
        blk /= a_old[3];
        gradJ[1] += -blk * x_n[i] * expblk;
        gradJ[2] += -blk * y_n[i] * expblk;
        gradJ[3] += pow(blk, 2) * expblk;
    }
    // gradJ4补完
    gradJ[3] = 0.5 * (gradJ[3] - (L + 1) * gradJ[0] / a_old[3]);
    // J
    J = gradJ[0] - (a_old[0] + a_old[1] * u + a_old[2] * v + a_old[3] * E);
    // 完成gradJ
    gradJ[0] = gradJ[0] - 1.;
    gradJ[1] = gradJ[1] - u;
    gradJ[2] = gradJ[2] - v;
    gradJ[3] = gradJ[3] - E;
    // 开始迭代
    while (err > REDEQUIL_ERROR)
    {
        iter++;// 计数

        // 求出求解方向
        MtimesV(B_inv, gradJ, 4, dir_t);

        n2gradJ = pow(gradJ[0], 2) + pow(gradJ[1], 2) + pow(gradJ[2], 2) + pow(gradJ[3], 2);

#pragma region 学习率调整
        // 学习率调整
        lr = 1;
        bool needloop = true;
        bool if_out = false;
        while ((needloop || if_out) && err > REDEQUIL_ERROR)
        {
            needloop = false;   // 假设已经不需要再循环了
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
            // 计算出J_new
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
            //    // 必须要上次不是出界，才放大lr，也是为了防止无限循环
            //    needloop = true;
            //    lr = lr * rGA[1];
            //}
            if_out = false; // 完成了一圈，没有出界
            err = lr * sqrt(n2gradJ);
        }
#pragma endregion

        // 完成一次迭代
        // 储存旧的值
        swapPtr(&a_old, &a_new);
        swapPtr(&g_old, &gradJ);
        // 计算出J和gradJ
        memset(gradJ, 0, sizeof(double) * 4);
        // 存下系数
        d = pow(2. * M_PI / -a_old[3], (L + 1) / 2.) * exp(a_old[0]);
        for (int i = 0;i < N;++i)
        {
            blk = a_old[1] * x_n[i] + a_old[2] * y_n[i];
            expblk = d * exp(-blk * blk * 0.5 / a_old[3]);
            blk /= a_old[3];
            // 计算gradJ，注意J和gradJ第一个数的关系
            gradJ[0] += expblk;
            gradJ[1] += -blk * x_n[i] * expblk;
            gradJ[2] += -blk * y_n[i] * expblk;
            gradJ[3] += pow(blk, 2) * expblk;
        }
        // gradJ4补完
        gradJ[3] = 0.5 * (gradJ[3] - (L + 1) * gradJ[0] / a_old[3]);
        // J
        J = gradJ[0] - (a_old[0] + a_old[1] * u + a_old[2] * v + a_old[3] * E);
        // 完成gradJ
        gradJ[0] = gradJ[0] - 1.;
        gradJ[1] = gradJ[1] - u;
        gradJ[2] = gradJ[2] - v;
        gradJ[3] = gradJ[3] - E;
        // BFGS核心部分，更新B_inverse
        for (int i = 0;i < 4;++i)
        {
            dlt_g[i] = gradJ[i] - g_old[i];
            dlt_a[i] = a_old[i] - a_new[i]; // 注意这里由于交换过，所以是old-new
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

    // 计算rhoeq的初步值
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
