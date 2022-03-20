#include "../include/Scheme.h"
#include "../include/settings.h"
#include "../include/2DRPFunctions.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <iostream>
#include "../include/RedEquil.h"


using namespace std;

void D_wind_up(int X, int Y, double dt)
{
#pragma omp parallel for
    for (int i = 0;i < X;++i)
    {
        for (int j = 0;j < Y;++j)
        {
            int i_p1 = min(i + 1, X - 1);
            int i_m1 = max(i - 1, 0);
            int j_p1 = min(j + 1, Y - 1);
            int j_m1 = max(j - 1, 0);

            // 计算本地温度
            double rt_local = macro(i, j, 3) - 0.5 * (pow(macro(i, j, 1), 2) + pow(macro(i, j, 2), 2));
            // 此时本轮计算宏观量已经没有用了
            macro(i, j, 0) = 0;
            macro(i, j, 1) = 0;
            macro(i, j, 2) = 0;
            macro(i, j, 3) = 0;
            for (int k = M * N - 1;k >= 0;--k)
            {
                double dl_term = 0;


                // 得到下一时刻MF
                if (flag_collide == 1 && kappa == 0)
                {
                    // 迎风格式
                    if (uk[k] > 0) { dl_term += (dt / dx) * uk[k] * (Meq(i, j, 0, k) - Meq(i_m1, j, 0, k)); }
                    else { dl_term += (dt / dx) * uk[k] * (Meq(i_p1, j, 0, k) - Meq(i, j, 0, k)); }
                    if (vk[k] > 0) { dl_term += (dt / dy) * vk[k] * (Meq(i, j, 0, k) - Meq(i, j_m1, 0, k)); }
                    else { dl_term += (dt / dy) * vk[k] * (Meq(i, j_p1, 0, k) - Meq(i, j, 0, k)); }

                    // 用平衡态计算
                    MF(i, j, 0, k) = Meq(i, j, 0, k) - dl_term;
                }
                else
                {
                    // 迎风格式
                    if (uk[k] > 0) { dl_term += (dt / dx) * uk[k] * (MF_c(i, j, 0, k) - MF_c(i_m1, j, 0, k)); }
                    else { dl_term += (dt / dx) * uk[k] * (MF_c(i_p1, j, 0, k) - MF_c(i, j, 0, k)); }
                    if (vk[k] > 0) { dl_term += (dt / dy) * vk[k] * (MF_c(i, j, 0, k) - MF_c(i, j_m1, 0, k)); }
                    else { dl_term += (dt / dy) * vk[k] * (MF_c(i, j_p1, 0, k) - MF_c(i, j, 0, k)); }

                    if (flag_collide == 1)
                    {
                        // 补上碰撞项
                        // 计算出本地的tau
                        // TODO
                        double tau = kappa / sqrt(rt_local);

                        // 使用指数型半隐式格式，注意0到4阶都要加
                        MF(i, j, 0, k) = exp(-dt / tau) * (
                            MF_c(i, j, 0, k) - dl_term
                            ) + (1 - exp(-dt / tau)) * Meq(i, j, 0, k);
                    }
                    else
                    {
                        // 无碰撞
                        MF(i, j, 0, k) = MF_c(i, j, 0, k) - dl_term;
                    }
                }

#pragma region 计算下一时刻的宏观量
                // 计算下一时刻的宏观量
                macro(i, j, 0) += MF(i, j, 0, k);                           // rho
                macro(i, j, 1) += MF(i, j, 0, k) * uk[k];  // rho u
                macro(i, j, 2) += MF(i, j, 0, k) * vk[k];  // rho v
                macro(i, j, 3) += MF(i, j, 0, k) * 0.5 * U2k[k];   // rho E
            }
            macro(i, j, 1) /= macro(i, j, 0);           // u
            macro(i, j, 2) /= macro(i, j, 0);           // v
            macro(i, j, 3) /= macro(i, j, 0);           // E
#pragma endregion

        }
    }
}


void DUGKS(int X, int Y, double dt)
{
    memset(&macro(0, 0, 0), 0, sizeof(double) * 4 * X * Y);
#pragma omp parallel for
    for (int i = 0;i < X;++i)
    {
        for (int j = 0;j < Y;++j)
        {
            int i_p1 = min(i + 1, X - 1);
            int i_m1 = max(i - 1, 0);
            int j_p1 = min(j + 1, Y - 1);
            int j_m1 = max(j - 1, 0);
            int i_p2 = min(i + 2, X - 1);
            int i_m2 = max(i - 2, 0);
            int j_p2 = min(j + 2, Y - 1);
            int j_m2 = max(j - 2, 0);

#pragma region 无碰撞情形
            if (flag_collide == 0)
            {
                // 创建四条边上的\Delta t / 2处的值
                double f_ip_j_h, f_im_j_h, f_i_jp_h, f_i_jm_h;
                double xr, ys;
                for (int k = M * N - 1;k >= 0;--k)
                {
#pragma region 九点插值
                    /*
                    // 使用9点插值方法
                    xr = uk[k] * 0.5 * (dt / dx);
                    ys = vk[k] * 0.5 * (dt / dy);
                    if (xr > 0)
                    {
                        // 应该使用靠左的点
                        f_im_j_h = DUGKS_9_interpol(
                            MF_c(i_m2, j_m1, 0, k), MF_c(i_m2, j, 0, k), MF_c(i_m2, j_p1, 0, k),
                            MF_c(i_m1, j_m1, 0, k), MF_c(i_m1, j, 0, k), MF_c(i_m1, j_p1, 0, k),
                            MF_c(i, j_m1, 0, k), MF_c(i, j, 0, k), MF_c(i, j_p1, 0, k),
                            0.5 - xr, -ys);
                        f_ip_j_h = DUGKS_9_interpol(
                            MF_c(i_m1, j_m1, 0, k), MF_c(i_m1, j, 0, k), MF_c(i_m1, j_p1, 0, k),
                            MF_c(i, j_m1, 0, k), MF_c(i, j, 0, k), MF_c(i, j_p1, 0, k),
                            MF_c(i_p1, j_m1, 0, k), MF_c(i_p1, j, 0, k), MF_c(i_p1, j_p1, 0, k),
                            0.5 - xr, -ys);
                    }
                    else
                    {
                        f_im_j_h = DUGKS_9_interpol(
                            MF_c(i_m1, j_m1, 0, k), MF_c(i_m1, j, 0, k), MF_c(i_m1, j_p1, 0, k),
                            MF_c(i, j_m1, 0, k), MF_c(i, j, 0, k), MF_c(i, j_p1, 0, k),
                            MF_c(i_p1, j_m1, 0, k), MF_c(i_p1, j, 0, k), MF_c(i_p1, j_p1, 0, k),
                            -0.5 - xr, -ys);
                        f_ip_j_h = DUGKS_9_interpol(
                            MF_c(i, j_m1, 0, k), MF_c(i, j, 0, k), MF_c(i, j_p1, 0, k),
                            MF_c(i_p1, j_m1, 0, k), MF_c(i_p1, j, 0, k), MF_c(i_p1, j_p1, 0, k),
                            MF_c(i_p2, j_m1, 0, k), MF_c(i_p2, j, 0, k), MF_c(i_p2, j_p1, 0, k),
                            -0.5 - xr, -ys);
                    }
                    if (ys > 0)
                    {
                        f_i_jm_h = DUGKS_9_interpol(
                            MF_c(i_m1, j_m2, 0, k), MF_c(i_m1, j_m1, 0, k), MF_c(i_m1, j, 0, k),
                            MF_c(i, j_m2, 0, k), MF_c(i, j_m1, 0, k), MF_c(i, j, 0, k),
                            MF_c(i_p1, j_m2, 0, k), MF_c(i_p1, j_m1, 0, k), MF_c(i_p1, j, 0, k),
                            -xr, 0.5 - ys);
                        f_i_jp_h = DUGKS_9_interpol(
                            MF_c(i_m1, j_m1, 0, k), MF_c(i_m1, j, 0, k), MF_c(i_m1, j_p1, 0, k),
                            MF_c(i, j_m1, 0, k), MF_c(i, j, 0, k), MF_c(i, j_p1, 0, k),
                            MF_c(i_p1, j_m1, 0, k), MF_c(i_p1, j, 0, k), MF_c(i_p1, j_p1, 0, k),
                            -xr, 0.5 - ys);
                    }
                    else
                    {
                        f_i_jm_h = DUGKS_9_interpol(
                            MF_c(i_m1, j_m1, 0, k), MF_c(i_m1, j, 0, k), MF_c(i_m1, j_p1, 0, k),
                            MF_c(i, j_m1, 0, k), MF_c(i, j, 0, k), MF_c(i, j_p1, 0, k),
                            MF_c(i_p1, j_m1, 0, k), MF_c(i_p1, j, 0, k), MF_c(i_p1, j_p1, 0, k),
                            -xr, -0.5 - ys);
                        f_i_jp_h = DUGKS_9_interpol(
                            MF_c(i_m1, j, 0, k), MF_c(i_m1, j_p1, 0, k), MF_c(i_m1, j_p2, 0, k),
                            MF_c(i, j, 0, k), MF_c(i, j_p1, 0, k), MF_c(i, j_p2, 0, k),
                            MF_c(i_p1, j, 0, k), MF_c(i_p1, j_p1, 0, k), MF_c(i_p1, j_p2, 0, k),
                            -xr, -0.5 - ys);
                    }
                    */
#pragma endregion

#pragma region VanLeer限制器5点格式
                    // 使用Van Leer限制器的格式
                    xr = uk[k] * 0.5 * (dt / dx);
                    ys = vk[k] * 0.5 * (dt / dy);
                    double sx1, sx2, sy1, sy2;  // 分别表示左右x格点差，下上y格点差
                    // 省略除以dx，dy的步骤
                    if (xr > 0)
                    {
                        // 应该使用靠左的点
                        sx1 = MF_c(i_m1, j, 0, k) - MF_c(i_m2, j, 0, k);
                        sx2 = MF_c(i, j, 0, k) - MF_c(i_m1, j, 0, k);
                        sy1 = MF_c(i_m1, j, 0, k) - MF_c(i_m1, j_m1, 0, k);
                        sy2 = MF_c(i_m1, j_p1, 0, k) - MF_c(i_m1, j, 0, k);

                        f_im_j_h = MF_c(i_m1, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (0.5 - xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0. - ys);

                        sx1 = sx2;
                        sx2 = MF_c(i_p1, j, 0, k) - MF_c(i, j, 0, k);
                        sy1 = MF_c(i, j, 0, k) - MF_c(i, j_m1, 0, k);
                        sy2 = MF_c(i, j_p1, 0, k) - MF_c(i, j, 0, k);

                        f_ip_j_h = MF_c(i, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (0.5 - xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0. - ys);
                    }
                    else
                    {
                        // 应该使用靠右的点
                        sx1 = MF_c(i, j, 0, k) - MF_c(i_m1, j, 0, k);
                        sx2 = MF_c(i_p1, j, 0, k) - MF_c(i, j, 0, k);
                        sy1 = MF_c(i, j, 0, k) - MF_c(i, j_m1, 0, k);
                        sy2 = MF_c(i, j_p1, 0, k) - MF_c(i, j, 0, k);

                        f_im_j_h = MF_c(i, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-0.5 - xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0. - ys);

                        sx1 = sx2;
                        sx2 = MF_c(i_p2, j, 0, k) - MF_c(i_p1, j, 0, k);
                        sy1 = MF_c(i_p1, j, 0, k) - MF_c(i_p1, j_m1, 0, k);
                        sy2 = MF_c(i_p1, j_p1, 0, k) - MF_c(i_p1, j, 0, k);

                        f_ip_j_h = MF_c(i_p1, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-0.5 - xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0. - ys);
                    }
                    if (ys > 0)
                    {
                        // 应该使用靠下的点
                        sx1 = MF_c(i, j_m1, 0, k) - MF_c(i_m1, j_m1, 0, k);
                        sx2 = MF_c(i_p1, j_m1, 0, k) - MF_c(i, j_m1, 0, k);
                        sy1 = MF_c(i, j_m1, 0, k) - MF_c(i, j_m2, 0, k);
                        sy2 = MF_c(i, j, 0, k) - MF_c(i, j_m1, 0, k);

                        f_i_jm_h = MF_c(i, j_m1, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0.5 - ys);

                        sx1 = MF_c(i, j, 0, k) - MF_c(i_m1, j, 0, k);
                        sx2 = MF_c(i_p1, j, 0, k) - MF_c(i, j, 0, k);
                        sy1 = sy2;
                        sy2 = MF_c(i, j_p1, 0, k) - MF_c(i, j, 0, k);

                        f_i_jp_h = MF_c(i, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0.5 - ys);
                    }
                    else
                    {
                        // 应该使用靠上的点
                        sx1 = MF_c(i, j, 0, k) - MF_c(i_m1, j, 0, k);
                        sx2 = MF_c(i_p1, j, 0, k) - MF_c(i, j, 0, k);
                        sy1 = MF_c(i, j, 0, k) - MF_c(i, j_m1, 0, k);
                        sy2 = MF_c(i, j_p1, 0, k) - MF_c(i, j, 0, k);

                        f_i_jm_h = MF_c(i, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-xr)
                            + Van_Leer_Limiter(sy1, sy2) * (-0.5 - ys);

                        sx1 = MF_c(i, j_p1, 0, k) - MF_c(i_m1, j_p1, 0, k);
                        sx2 = MF_c(i_p1, j_p1, 0, k) - MF_c(i, j_p1, 0, k);
                        sy1 = sy2;
                        sy2 = MF_c(i, j_p2, 0, k) - MF_c(i, j_p1, 0, k);

                        f_i_jp_h = MF_c(i, j_p1, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-xr)
                            + Van_Leer_Limiter(sy1, sy2) * (-0.5 - ys);
                    }
#pragma endregion

                    MF(i, j, 0, k) = MF_c(i, j, 0, k)
                        - (dt / dx) * uk[k] * (f_ip_j_h - f_im_j_h)
                        - (dt / dy) * vk[k] * (f_i_jp_h - f_i_jm_h);

                    // 计算下一时刻的宏观量
                    macro(i, j, 0) += MF(i, j, 0, k);               // rho
                    macro(i, j, 1) += MF(i, j, 0, k) * uk[k];       // rho u
                    macro(i, j, 2) += MF(i, j, 0, k) * vk[k];       // rho v
                    macro(i, j, 3) += MF(i, j, 0, k) * 0.5 * U2k[k];// rho E
                }
            }
#pragma endregion
            else if (kappa == 0)
            {
                // 创建四条边上的\Delta t / 2处的值
                double f_ip_j_h, f_im_j_h, f_i_jp_h, f_i_jm_h;
                double xr, ys;
                for (int k = M * N - 1;k >= 0;--k)
                {
#pragma region 九点插值
                    /*
                    // 使用9点插值方法
                    xr = uk[k] * 0.5 * (dt / dx);
                    ys = vk[k] * 0.5 * (dt / dy);
                    if (xr > 0)
                    {
                        // 应该使用靠左的点
                        f_im_j_h = DUGKS_9_interpol(
                            Meq(i_m2, j_m1, 0, k), Meq(i_m2, j, 0, k), Meq(i_m2, j_p1, 0, k),
                            Meq(i_m1, j_m1, 0, k), Meq(i_m1, j, 0, k), Meq(i_m1, j_p1, 0, k),
                            Meq(i, j_m1, 0, k), Meq(i, j, 0, k), Meq(i, j_p1, 0, k),
                            0.5 - xr, -ys);
                        f_ip_j_h = DUGKS_9_interpol(
                            Meq(i_m1, j_m1, 0, k), Meq(i_m1, j, 0, k), Meq(i_m1, j_p1, 0, k),
                            Meq(i, j_m1, 0, k), Meq(i, j, 0, k), Meq(i, j_p1, 0, k),
                            Meq(i_p1, j_m1, 0, k), Meq(i_p1, j, 0, k), Meq(i_p1, j_p1, 0, k),
                            0.5 - xr, -ys);
                    }
                    else
                    {
                        f_im_j_h = DUGKS_9_interpol(
                            Meq(i_m1, j_m1, 0, k), Meq(i_m1, j, 0, k), Meq(i_m1, j_p1, 0, k),
                            Meq(i, j_m1, 0, k), Meq(i, j, 0, k), Meq(i, j_p1, 0, k),
                            Meq(i_p1, j_m1, 0, k), Meq(i_p1, j, 0, k), Meq(i_p1, j_p1, 0, k),
                            -0.5 - xr, -ys);
                        f_ip_j_h = DUGKS_9_interpol(
                            Meq(i, j_m1, 0, k), Meq(i, j, 0, k), Meq(i, j_p1, 0, k),
                            Meq(i_p1, j_m1, 0, k), Meq(i_p1, j, 0, k), Meq(i_p1, j_p1, 0, k),
                            Meq(i_p2, j_m1, 0, k), Meq(i_p2, j, 0, k), Meq(i_p2, j_p1, 0, k),
                            -0.5 - xr, -ys);
                    }
                    if (ys > 0)
                    {
                        f_i_jm_h = DUGKS_9_interpol(
                            Meq(i_m1, j_m2, 0, k), Meq(i_m1, j_m1, 0, k), Meq(i_m1, j, 0, k),
                            Meq(i, j_m2, 0, k), Meq(i, j_m1, 0, k), Meq(i, j, 0, k),
                            Meq(i_p1, j_m2, 0, k), Meq(i_p1, j_m1, 0, k), Meq(i_p1, j, 0, k),
                            -xr, 0.5 - ys);
                        f_i_jp_h = DUGKS_9_interpol(
                            Meq(i_m1, j_m1, 0, k), Meq(i_m1, j, 0, k), Meq(i_m1, j_p1, 0, k),
                            Meq(i, j_m1, 0, k), Meq(i, j, 0, k), Meq(i, j_p1, 0, k),
                            Meq(i_p1, j_m1, 0, k), Meq(i_p1, j, 0, k), Meq(i_p1, j_p1, 0, k),
                            -xr, 0.5 - ys);
                    }
                    else
                    {
                        f_i_jm_h = DUGKS_9_interpol(
                            Meq(i_m1, j_m1, 0, k), Meq(i_m1, j, 0, k), Meq(i_m1, j_p1, 0, k),
                            Meq(i, j_m1, 0, k), Meq(i, j, 0, k), Meq(i, j_p1, 0, k),
                            Meq(i_p1, j_m1, 0, k), Meq(i_p1, j, 0, k), Meq(i_p1, j_p1, 0, k),
                            -xr, -0.5 - ys);
                        f_i_jp_h = DUGKS_9_interpol(
                            Meq(i_m1, j, 0, k), Meq(i_m1, j_p1, 0, k), Meq(i_m1, j_p2, 0, k),
                            Meq(i, j, 0, k), Meq(i, j_p1, 0, k), Meq(i, j_p2, 0, k),
                            Meq(i_p1, j, 0, k), Meq(i_p1, j_p1, 0, k), Meq(i_p1, j_p2, 0, k),
                            -xr, -0.5 - ys);
                    }
                    */
#pragma endregion


#pragma region VanLeer限制器5点格式
                    // 使用Van Leer限制器的格式
                    xr = uk[k] * 0.5 * (dt / dx);
                    ys = vk[k] * 0.5 * (dt / dy);
                    double sx1, sx2, sy1, sy2;  // 分别表示左右x格点差，下上y格点差
                    // 省略除以dx，dy的步骤
                    if (xr > 0)
                    {
                        // 应该使用靠左的点
                        sx1 = Meq(i_m1, j, 0, k) - Meq(i_m2, j, 0, k);
                        sx2 = Meq(i, j, 0, k) - Meq(i_m1, j, 0, k);
                        sy1 = Meq(i_m1, j, 0, k) - Meq(i_m1, j_m1, 0, k);
                        sy2 = Meq(i_m1, j_p1, 0, k) - Meq(i_m1, j, 0, k);

                        f_im_j_h = Meq(i_m1, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (0.5 - xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0. - ys);

                        sx1 = sx2;
                        sx2 = Meq(i_p1, j, 0, k) - Meq(i, j, 0, k);
                        sy1 = Meq(i, j, 0, k) - Meq(i, j_m1, 0, k);
                        sy2 = Meq(i, j_p1, 0, k) - Meq(i, j, 0, k);

                        f_ip_j_h = Meq(i, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (0.5 - xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0. - ys);
                    }
                    else
                    {
                        // 应该使用靠右的点
                        sx1 = Meq(i, j, 0, k) - Meq(i_m1, j, 0, k);
                        sx2 = Meq(i_p1, j, 0, k) - Meq(i, j, 0, k);
                        sy1 = Meq(i, j, 0, k) - Meq(i, j_m1, 0, k);
                        sy2 = Meq(i, j_p1, 0, k) - Meq(i, j, 0, k);

                        f_im_j_h = Meq(i, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-0.5 - xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0. - ys);

                        sx1 = sx2;
                        sx2 = Meq(i_p2, j, 0, k) - Meq(i_p1, j, 0, k);
                        sy1 = Meq(i_p1, j, 0, k) - Meq(i_p1, j_m1, 0, k);
                        sy2 = Meq(i_p1, j_p1, 0, k) - Meq(i_p1, j, 0, k);

                        f_ip_j_h = Meq(i_p1, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-0.5 - xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0. - ys);
                    }
                    if (ys > 0)
                    {
                        // 应该使用靠下的点
                        sx1 = Meq(i, j_m1, 0, k) - Meq(i_m1, j_m1, 0, k);
                        sx2 = Meq(i_p1, j_m1, 0, k) - Meq(i, j_m1, 0, k);
                        sy1 = Meq(i, j_m1, 0, k) - Meq(i, j_m2, 0, k);
                        sy2 = Meq(i, j, 0, k) - Meq(i, j_m1, 0, k);

                        f_i_jm_h = Meq(i, j_m1, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0.5 - ys);

                        sx1 = Meq(i, j, 0, k) - Meq(i_m1, j, 0, k);
                        sx2 = Meq(i_p1, j, 0, k) - Meq(i, j, 0, k);
                        sy1 = sy2;
                        sy2 = Meq(i, j_p1, 0, k) - Meq(i, j, 0, k);

                        f_i_jp_h = Meq(i, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-xr)
                            + Van_Leer_Limiter(sy1, sy2) * (0.5 - ys);
                    }
                    else
                    {
                        // 应该使用靠上的点
                        sx1 = Meq(i, j, 0, k) - Meq(i_m1, j, 0, k);
                        sx2 = Meq(i_p1, j, 0, k) - Meq(i, j, 0, k);
                        sy1 = Meq(i, j, 0, k) - Meq(i, j_m1, 0, k);
                        sy2 = Meq(i, j_p1, 0, k) - Meq(i, j, 0, k);

                        f_i_jm_h = Meq(i, j, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-xr)
                            + Van_Leer_Limiter(sy1, sy2) * (-0.5 - ys);

                        sx1 = Meq(i, j_p1, 0, k) - Meq(i_m1, j_p1, 0, k);
                        sx2 = Meq(i_p1, j_p1, 0, k) - Meq(i, j_p1, 0, k);
                        sy1 = sy2;
                        sy2 = Meq(i, j_p2, 0, k) - Meq(i, j_p1, 0, k);

                        f_i_jp_h = Meq(i, j_p1, 0, k)
                            + Van_Leer_Limiter(sx1, sx2) * (-xr)
                            + Van_Leer_Limiter(sy1, sy2) * (-0.5 - ys);
                    }
#pragma endregion

                    MF(i, j, 0, k) = Meq(i, j, 0, k)
                        - (dt / dx) * uk[k] * (f_ip_j_h - f_im_j_h)
                        - (dt / dy) * vk[k] * (f_i_jp_h - f_i_jm_h);

                    // 计算下一时刻的宏观量
                    macro(i, j, 0) += MF(i, j, 0, k);               // rho
                    macro(i, j, 1) += MF(i, j, 0, k) * uk[k];       // rho u
                    macro(i, j, 2) += MF(i, j, 0, k) * vk[k];       // rho v
                    macro(i, j, 3) += MF(i, j, 0, k) * 0.5 * U2k[k];// rho E
                }
            }
            macro(i, j, 1) /= macro(i, j, 0);           // u
            macro(i, j, 2) /= macro(i, j, 0);           // v
            macro(i, j, 3) /= macro(i, j, 0);           // E
        }
    }
}


void DUGKS_gh(int X, int Y, double dt, double mode)
{
    memset(&macro(0, 0, 0), 0, sizeof(double) * 4 * X * Y);
    if (flag_collide == 0 || kappa == 0)
    {
#pragma region 无碰撞情形或者连续情形
        CArray* fs{};
        if (flag_collide == 0)
        {
            fs = &MF_c;
        }
        else if (kappa == 0)
        {
            fs = &Meq;
        }
        CArray mf = *fs;   // 注意这个mf是无法改变的，只能读取，原理不懂

#pragma omp parallel for
        for (int i = 0;i < X;++i)
        {
            for (int j = 0;j < Y;++j)
            {
                int i_p1 = min(i + 1, X - 1);
                int i_m1 = max(i - 1, 0);
                int j_p1 = min(j + 1, Y - 1);
                int j_m1 = max(j - 1, 0);
                int i_p2 = min(i + 2, X - 1);
                int i_m2 = max(i - 2, 0);
                int j_p2 = min(j + 2, Y - 1);
                int j_m2 = max(j - 2, 0);

                // 创建四条边上的\Delta t / 2处的值
                double f_ip_j_h, f_im_j_h, f_i_jp_h, f_i_jm_h;
                double xr, ys;
                int ik = 0;
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    ik = k % (M * N); // 得到小ik
#pragma region VanLeer限制器5点格式
                    // 使用Van Leer限制器的格式
                    xr = uk[ik] * 0.5 * (dt / dx);
                    ys = vk[ik] * 0.5 * (dt / dy);
                    // 省略除以dx，dy的步骤
                    if (xr > 0)
                    {
                        // 应该使用靠左的点
                        f_im_j_h = DUGKS_5_vnLr(mf(i_m2, j, 0, 0, k), mf(i_m1, j, 0, 0, k),
                            mf(i, j, 0, 0, k), mf(i_m1, j_m1, 0, 0, k), mf(i_m1, j_p1, 0, 0, k),
                            0.5 - xr, -ys);
                        f_ip_j_h = DUGKS_5_vnLr(mf(i_m1, j, 0, 0, k), mf(i, j, 0, 0, k),
                            mf(i_p1, j, 0, 0, k), mf(i, j_m1, 0, 0, k), mf(i, j_p1, 0, 0, k),
                            0.5 - xr, -ys);
                    }
                    else
                    {
                        // 应该使用靠右的点
                        f_im_j_h = DUGKS_5_vnLr(mf(i_m1, j, 0, 0, k), mf(i, j, 0, 0, k),
                            mf(i_p1, j, 0, 0, k), mf(i, j_m1, 0, 0, k), mf(i, j_p1, 0, 0, k),
                            -0.5 - xr, -ys);
                        f_ip_j_h = DUGKS_5_vnLr(mf(i, j, 0, 0, k), mf(i_p1, j, 0, 0, k),
                            mf(i_p2, j, 0, 0, k), mf(i_p1, j_m1, 0, 0, k), mf(i_p1, j_p1, 0, 0, k),
                            -0.5 - xr, -ys);
                    }
                    if (ys > 0)
                    {
                        // 应该使用靠下的点
                        f_i_jm_h = DUGKS_5_vnLr(mf(i_m1, j_m1, 0, 0, k), mf(i, j_m1, 0, 0, k),
                            mf(i_p1, j_m1, 0, 0, k), mf(i, j_m2, 0, 0, k), mf(i, j, 0, 0, k),
                            -xr, 0.5 - ys);
                        f_i_jp_h = DUGKS_5_vnLr(mf(i_m1, j, 0, 0, k), mf(i, j, 0, 0, k),
                            mf(i_p1, j, 0, 0, k), mf(i, j_m1, 0, 0, k), mf(i, j_p1, 0, 0, k),
                            -xr, 0.5 - ys);
                    }
                    else
                    {
                        // 应该使用靠上的点
                        f_i_jm_h = DUGKS_5_vnLr(mf(i_m1, j, 0, 0, k), mf(i, j, 0, 0, k),
                            mf(i_p1, j, 0, 0, k), mf(i, j_m1, 0, 0, k), mf(i, j_p1, 0, 0, k),
                            -xr, -0.5 - ys);
                        f_i_jp_h = DUGKS_5_vnLr(mf(i_m1, j_p1, 0, 0, k), mf(i, j_p1, 0, 0, k),
                            mf(i_p1, j_p1, 0, 0, k), mf(i, j, 0, 0, k), mf(i, j_p2, 0, 0, k),
                            -xr, -0.5 - ys);
                    }
#pragma endregion
                    MF(i, j, 0, 0, k) = mf(i, j, 0, 0, k)
                        - (dt / dx) * uk[ik] * (f_ip_j_h - f_im_j_h)
                        - (dt / dy) * vk[ik] * (f_i_jp_h - f_i_jm_h);
                }
                for (int k = M * N - 1;k >= 0;--k)
                {
                    // 计算下一时刻的宏观量
                    macro(i, j, 0) += MF(i, j, 0, 0, k);        // rho
                    macro(i, j, 1) += MF(i, j, 0, 0, k) * uk[k];// rho u
                    macro(i, j, 2) += MF(i, j, 0, 0, k) * vk[k];// rho v
                    macro(i, j, 3) += (MF(i, j, 0, 0, k) * U2k[k]
                        + MF(i, j, 1, 0, k)) * 0.5;             // rho E
                }

                macro(i, j, 1) /= macro(i, j, 0);           // u
                macro(i, j, 2) /= macro(i, j, 0);           // v
                macro(i, j, 3) /= macro(i, j, 0);           // E
            }
        }
#pragma endregion
    }
    else
    {
#pragma region 一般情形
        CArray* fs = &MF_c;
        CArray mf = *fs;   // 注意这个mf是无法改变的，只能读取，原理不懂
        int finish_count;
        int grid_num;
#pragma region 第一轮
        // 第一轮，算出f(i,j+-1/2,n+1/2)，共有X*(Y+1)个，存到MF_m中
        finish_count = 0;
        grid_num = X * (Y + 1);
#pragma omp parallel for
        for (int i = 0;i < X;++i)
        {
            for (int j = 0;j < Y + 1;++j)
            {
                // 算出\bar f(i,j+-1/2,n+1/2)，共有X*(Y+1)个，存到MF_m中
                int i_p = min(i + 1, X - 1);
                int i_m = max(i - 1, 0);
                int j_p;
                int j_c;
                int j_m;

                double xr, ys;
                int ik = 0;
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    ik = k % (M * N); // 得到小ik
#pragma region VanLeer限制器5点格式
                    // 使用Van Leer限制器的格式
                    xr = uk[ik] * 0.5 * (dt / dx);
                    ys = vk[ik] * 0.5 * (dt / dy);
                    // 省略除以dx，dy的步骤
                    if (ys > 0)
                    {
                        // 应该使用靠下的点
                        j_p = min(j, Y - 1);
                        j_c = max(j - 1, 0);
                        j_m = max(j - 2, 0);
                        xr = -xr;
                        ys = 0.5 - ys;
                    }
                    else
                    {
                        // 应该使用靠上的点
                        j_p = min(j + 1, Y - 1);
                        j_c = min(j, Y - 1);
                        j_m = max(j - 1, 0);
                        xr = -xr;
                        ys = -0.5 - ys;
                    }
                    MF_m(i, j, 0, 0, k) = DUGKS_5_vnLr(mf(i_m, j_c, 0, 0, k),
                        mf(i, j_c, 0, 0, k), mf(i_p, j_c, 0, 0, k), mf(i, j_m, 0, 0, k),
                        mf(i, j_p, 0, 0, k), xr, ys);
#pragma endregion
                }
                // 计算1/2处宏观量
                double rho{}, u{}, v{}, E{}, rt{};
                for (int k = M * N - 1;k >= 0;--k)
                {
                    rho += MF_m(i, j, 0, 0, k);         // rho
                    u += MF_m(i, j, 0, 0, k) * uk[k];   // rho u
                    v += MF_m(i, j, 0, 0, k) * vk[k];   // rho v
                    E += (MF_m(i, j, 0, 0, k) * U2k[k]
                        + MF_m(i, j, 1, 0, k)) * 0.5;   // rho E
                }
                u /= rho;   // u
                v /= rho;   // v
                E /= rho;   // E
                rt = (2 * E - u * u - v * v) / (D + L);
                // 准备转化为f(i,j+-1/2,n+1/2)
                // 计算平衡态
                double* req_m = new double[2 * M * N]{};
                double a[4] = { -1.5 * log(2 * M_PI * rt) - (u * u + v * v) / 2 / rt, u / rt, v / rt, -1 / rt };
                if (mode == 3)
                {
                    double sig2 = RedEquil_DVM(rho, u, v, E, a, req_m);
                    for (int k = M * N - 1;k >= 0;--k)
                    {
                        req_m[M * N + k] = L * sig2 * req_m[k];
                    }
                }
                else if (mode == 4)
                {
                    double* req = new double[N] {};
                    double* ueq = new double[N] {};
                    double sig2eq = 0;
                    RedEquil_DVD(rho, u, v, E, a, req, ueq, &sig2eq);
                    for (int n = 0;n < N;++n)
                    {
                        double sig2 = RedEquil_DVDDVM(req[n], ueq[n], 0.5 * (ueq[n] * ueq[n] + (L + 1) * sig2eq),
                            &req_m[n * M]);
                        for (int m = 0;m < M;++m)
                        {
                            req_m[(n + N) * M + m] = L * sig2 * req_m[n * M + m];
                        }
                    }
                    delete[] req;
                    delete[] ueq;
                }
                // 转化为f(i,j+-1/2,n+1/2)，Guo 2013文章，20式
                double e1 = (2 * kappa) / (2 * kappa + dt * 0.5);
                double e2 = (dt * 0.5) / (2 * kappa + dt * 0.5);
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    MF_m(i, j, 0, 0, k) = e1 * MF_m(i, j, 0, 0, k) + e2 * req_m[k];
                }
                delete[] req_m;

                finish_count++;
                printf("计算中间平衡态X，已完成 %d / %d    \r", finish_count, grid_num);
            }
        }
        // 把这些东西全部加入MF
#pragma omp parallel for
        for (int i = 0;i < X;++i)
        {
            for (int j = 0;j < Y;++j)
            {
                for (int k = 2 * M * N;k >= 0;--k)
                {
                    MF(i, j, 0, 0, k) -= (dt / dy) * vk[k % (M * N)] 
                        * (MF_m(i, j + 1, 0, 0, k) - MF_m(i, j, 0, 0, k));
                }
            }
        }
#pragma endregion

#pragma region 第二轮
        // 第二轮，算出f(i+-1/2,j,n+1/2)，共有(X+1)*Y个，存到MF_m中
        finish_count = 0;
        grid_num = X * (Y + 1);
#pragma omp parallel for
        for (int i = 0;i < X + 1;++i)
        {
            for (int j = 0;j < Y;++j)
            {
                // 算出\bar f(i+-1/2,j,n+1/2)，共有(X+1)*Y个，存到MF_m中
                int j_p = min(j + 1, Y - 1);
                int j_m = max(j - 1, 0);
                int i_p;
                int i_c;
                int i_m;

                double xr, ys;
                int ik = 0;
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    ik = k % (M * N); // 得到小ik
#pragma region VanLeer限制器5点格式
                    // 使用Van Leer限制器的格式
                    xr = uk[ik] * 0.5 * (dt / dx);
                    ys = vk[ik] * 0.5 * (dt / dy);
                    // 省略除以dx，dy的步骤
                    if (xr > 0)
                    {
                        // 应该使用靠左的点
                        i_p = min(i, X - 1);
                        i_c = max(i - 1, 0);
                        i_m = max(i - 2, 0);
                        xr = 0.5 - xr;
                        ys = -ys;
                    }
                    else
                    {
                        // 应该使用靠右的点、
                        i_p = min(i+1, X - 1);
                        i_c = min(i, X - 1);
                        i_m = max(i - 1, 0);
                        xr = -0.5 - xr;
                        ys = -ys;
                    }
                    MF_m(i, j, 0, 0, k) = DUGKS_5_vnLr(mf(i_m, j, 0, 0, k),
                        mf(i_c, j, 0, 0, k), mf(i_p, j, 0, 0, k), mf(i_c, j_m, 0, 0, k),
                        mf(i_c, j_p, 0, 0, k), xr, ys);
#pragma endregion
                }
                // 计算1/2处宏观量
                double rho{}, u{}, v{}, E{}, rt{};
                for (int k = M * N - 1;k >= 0;--k)
                {
                    rho += MF_m(i, j, 0, 0, k);         // rho
                    u += MF_m(i, j, 0, 0, k) * uk[k];   // rho u
                    v += MF_m(i, j, 0, 0, k) * vk[k];   // rho v
                    E += (MF_m(i, j, 0, 0, k) * U2k[k]
                        + MF_m(i, j, 1, 0, k)) * 0.5;   // rho E
                }
                u /= rho;   // u
                v /= rho;   // v
                E /= rho;   // E
                rt = (2 * E - u * u - v * v) / (D + L);
                // 准备转化为f(i,j+-1/2,n+1/2)
                // 计算平衡态
                double* req_m = new double[2 * M * N]{};
                double a[4] = { -1.5 * log(2 * M_PI * rt) - (u * u + v * v) / 2 / rt, u / rt, v / rt, -1 / rt };
                if (mode == 3)
                {
                    double sig2 = RedEquil_DVM(rho, u, v, E, a, req_m);
                    for (int k = M * N - 1;k >= 0;--k)
                    {
                        req_m[M * N + k] = L * sig2 * req_m[k];
                    }
                }
                else if (mode == 4)
                {
                    double* req = new double[N] {};
                    double* ueq = new double[N] {};
                    double sig2eq = 0;
                    RedEquil_DVD(rho, u, v, E, a, req, ueq, &sig2eq);
                    for (int n = 0;n < N;++n)
                    {
                        double sig2 = RedEquil_DVDDVM(req[n], ueq[n], 0.5 * (ueq[n] * ueq[n] + (L + 1) * sig2eq),
                            &req_m[n * M]);
                        for (int m = 0;m < M;++m)
                        {
                            req_m[(n + N) * M + m] = L * sig2 * req_m[n * M + m];
                        }
                    }
                    delete[] req;
                    delete[] ueq;
                }
                // 转化为f(i+-1/2,j,n+1/2)，Guo 2013文章，20式
                double e1 = (2 * kappa) / (2 * kappa + dt * 0.5);
                double e2 = (dt * 0.5) / (2 * kappa + dt * 0.5);
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    MF_m(i, j, 0, 0, k) = e1 * MF_m(i, j, 0, 0, k) + e2 * req_m[k];
                }
                delete[] req_m;

                finish_count++;
                printf("计算中间平衡态Y，已完成 %d / %d    \r", finish_count, grid_num);
            }
        }
        // 把这些东西全部加入MF
#pragma omp parallel for
        for (int i = 0;i < X;++i)
        {
            for (int j = 0;j < Y;++j)
            {
                for (int k = 2 * M * N;k >= 0;--k)
                {
                    MF(i, j, 0, 0, k) -= (dt / dx) * uk[k % (M * N)]
                        * (MF_m(i + 1, j, 0, 0, k) - MF_m(i, j, 0, 0, k));
                }

                for (int k = M * N - 1;k >= 0;--k)
                {
                    // 计算下一时刻的宏观量
                    macro(i, j, 0) += MF(i, j, 0, 0, k);        // rho
                    macro(i, j, 1) += MF(i, j, 0, 0, k) * uk[k];// rho u
                    macro(i, j, 2) += MF(i, j, 0, 0, k) * vk[k];// rho v
                    macro(i, j, 3) += (MF(i, j, 0, 0, k) * U2k[k]
                        + MF(i, j, 1, 0, k)) * 0.5;             // rho E
                }

                macro(i, j, 1) /= macro(i, j, 0);           // u
                macro(i, j, 2) /= macro(i, j, 0);           // v
                macro(i, j, 3) /= macro(i, j, 0);           // E
            }
        }
#pragma endregion

#pragma endregion
    }
}

double DUGKS_5_vnLr(double fl, double fc, double fr, double fd, double fu,
    double r, double s)
{
    return fc + Van_Leer_Limiter(fc - fl, fr - fc) * r
        + Van_Leer_Limiter(fc - fd, fu - fc) * s;
}

double DUGKS_9_interpol(double fmm, double fm0, double fmp,
    double f0m, double f00, double f0p,
    double fpm, double fp0, double fpp, double r, double s)
{
    return 0.5 * r * (r - 1) * (0.5 * s * (s - 1) * fmm + (1 - s * s) * fm0 + 0.5 * s * (s + 1) * fmp)
        + (1 - r * r) * (0.5 * s * (s - 1) * f0m + (1 - s * s) * f00 + 0.5 * s * (s + 1) * f0p)
        + 0.5 * r * (r + 1) * (0.5 * s * (s - 1) * fpm + (1 - s * s) * fp0 + 0.5 * s * (s + 1) * fpp);
}

double Van_Leer_Limiter(double s_1, double s_2)
{
    if (sign(s_1) + sign(s_2) == 0) { return 0; }
    return 2 * s_1 * s_2 / (s_1 + s_2);
}