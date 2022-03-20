#include "../include/settings.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <iostream>

using namespace std;
/* 定义设置变量 */

#pragma region 各种设置

// 维度设置
extern int D = 3;
extern int L = 0;

// DVD相关
extern int N = 2;
extern double* x_n = new double[N];  // 各个方向x
extern double* y_n = new double[N];  // 各个方向y
extern double* z_n = new double[N];  // 各个方向z
// DVM相关
extern double dvm_bg = -6.25;   // dvm
extern double dvm_st = 0.5;
extern double* dvm = new double[2] {-1, 1};
extern int M = 2;
extern double* uk = new double[M * N];
extern double* vk = new double[M * N];
extern double* wk = new double[M * N];
extern double* U2k = new double[M * N];
// BGK相关
extern double flag_collide = 1;
extern double kappa = 0.0;
// 网格差分相关
extern double dx = 0.005;
extern double dy = 0.005;
extern double CFLMax = 0.5;

//2D黎曼问题初始设置
/*  |-------|(2*x_num-1,2*y_num-1)
    |0,1|1,1|
    |-------|
    |0,0|1,0|
    |-------|(2*x_num-1,0)  */
extern double r_init[2][2] = { {1.0,1.0}, {1.0,1.0} };
extern double u_init[2][2] = { {}, {} };
extern double v_init[2][2] = { {}, {} };
extern double p_init[2][2] = { {1.0,1.0}, {1.0,1.0} };
extern double E_init[2][2] = { {},{} };

// 计算总时间等
extern double t_current = 0.;
extern int step_count = -1;    // 从-1开始

// 计算精度
extern double REDEQUIL_ERROR = 5e-8;

// 差分格式
extern std::string SCHEME = "DUGKS";
// 离散点选取方法
extern std::string dvm_strategy = "common-difference";


// 保存用的数组
extern int save_max = 26;            // 最大储存次数
extern double save_interval = 0.01;  // 每隔0.01秒储存一次
extern int Save_Count = 1;

// 计算用大数组
extern CArray MF(1, 1);       //前一时刻数据(x，y，方向，dvm)
extern CArray MF_c(1, 1);     //复制数据(x，y，方向，dvm)
extern CArray MF_m(1, 1);     //DUGKS专用的中间态数据(x，y，方向，dvm)
extern CArray Meq(1, 1);      //平衡态
extern CArray macro(1, 1);    //宏观量(x，y，（rho,u,v,e）)

void SetDVM()
{
    if (dvm_strategy == "common-difference")
    {
        M = 0;
        double bg = -abs(dvm_bg);
        for (double db = bg;db <= abs(bg) + 1e-10;db += dvm_st) { ++M; }
        delete[] dvm;
        dvm = new double[M];
        for (int m = 0;m < M;++m)
        {
            dvm[m] = bg;
            bg += dvm_st;
        }
    }
    else if (dvm_strategy == "sqrt-central")
    {
        int md = floor(abs(dvm_bg));
        M = md * 2 + 1;
        delete[] dvm;
        dvm = new double[M] {};
        for (int i = 1;i <= md;++i)
        {
            dvm[md - i] = -sqrt((double)i) * dvm_st;
            dvm[md + i] = sqrt((double)i) * dvm_st;
        }
    }
    else if (dvm_strategy == "sqrt-symmetric")
    {
        int md = floor(abs(dvm_bg));
        M = md * 2;
        delete[] dvm;
        dvm = new double[M] {};
        for (int i = 1;i <= md;++i)
        {
            dvm[md - i] = -(sqrt((double)i) - 0.5) * dvm_st;
            dvm[md + i - 1] = (sqrt((double)i) - 0.5) * dvm_st;
        }
    }
    else if (dvm_strategy == "3sqrt-central")
    {
        int md = floor(abs(dvm_bg));
        M = md * 2 + 1;
        delete[] dvm;
        dvm = new double[M] {};
        for (int i = 1;i <= md;++i)
        {
            dvm[md - i] = -pow((double)i, 1. / 3.) * dvm_st;
            dvm[md + i] = pow((double)i, 1. / 3.) * dvm_st;
        }
    }
    else if (dvm_strategy == "3sqrt-symmetric")
    {
        int md = floor(abs(dvm_bg));
        M = md * 2;
        delete[] dvm;
        dvm = new double[M] {};
        for (int i = 1;i <= md;++i)
        {
            dvm[md - i] = -(pow((double)i, 1. / 3.) - 0.5) * dvm_st;
            dvm[md + i - 1] = (pow((double)i, 1. / 3.) - 0.5) * dvm_st;
        }
    }
}

// 发生在设置Direction之后
void SetUVW()
{
    delete[] uk;
    delete[] vk;
    delete[] wk;
    delete[] U2k;
    uk = new double[M * N];
    vk = new double[M * N];
    wk = new double[M * N];
    U2k = new double[M * N];
    for (int n = 0;n < N;++n)
    {
        for (int m = 0;m < M;++m)
        {
            uk[n * M + m] = x_n[n] * dvm[m];
            vk[n * M + m] = y_n[n] * dvm[m];
            wk[n * M + m] = z_n[n] * dvm[m];
            U2k[n * M + m] = dvm[m] * dvm[m];
        }
    }
}

//发生在UVW之前
void SetDirections()
{
    if (L > 0)
    {
        delete[] x_n;
        delete[] y_n;
        delete[] z_n;
        x_n = new double[N] {};
        y_n = new double[N] {};
        z_n = new double[N] {};
        for (int n = 0;n < N;++n)
        {
            x_n[n] = cos((n + 0.5) / N * M_PI);
            y_n[n] = sin((n + 0.5) / N * M_PI);
            z_n[n] = 0;
        }
        return;
    }
    int Nlen = N;
    int qN = N * N * N - (N - 1) * (N - 1) * (N - 1);   //  总方向数目的1/4
    N = 4 * qN;
    delete[] x_n;
    delete[] y_n;
    delete[] z_n;
    x_n = new double[N] {};
    y_n = new double[N] {};
    z_n = new double[N] {};
    double* qx_n = new double[qN] {};
    double* qy_n = new double[qN] {};
    double* qz_n = new double[qN] {};
    int k = 0;
    for (int x = 0;x < Nlen;++x)
    {
        for (int y = 0;y < Nlen;++y)
        {
            qx_n[k] = x + 0.5;
            qy_n[k] = y + 0.5;
            qz_n[k] = Nlen - 0.5;
            ++k;
        }
    }
    for (int z = 0;z < Nlen - 1;++z)
    {
        for (int x = 0;x < Nlen;++x)
        {
            qx_n[k] = x + 0.5;
            qy_n[k] = Nlen - 0.5;
            qz_n[k] = z + 0.5;
            ++k;
        }
        for (int y = 0;y < Nlen-1;++y)
        {
            qx_n[k] = Nlen - 0.5;
            qy_n[k] = y + 0.5;
            qz_n[k] = z + 0.5;
            ++k;
        }
    }
    for (int k1 = 0;k1 < qN;++k1)
    {
        Set3DNormalVec(&x_n[k1], &y_n[k1], &z_n[k1], qx_n[k1], qy_n[k1], qz_n[k1]);
    }
    for (int k2 = 0;k2 < qN;++k2)
    {
        Set3DNormalVec(&x_n[k2 + qN], &y_n[k2 + qN], &z_n[k2 + qN], qx_n[k2], -qy_n[k2], qz_n[k2]);
    }
    for (int k3 = 0;k3 < qN;++k3)
    {
        Set3DNormalVec(&x_n[k3 + 2 * qN], &y_n[k3 + 2 * qN], &z_n[k3 + 2 * qN], -qx_n[k3], qy_n[k3], qz_n[k3]);
    }
    for (int k4 = 0;k4 < qN;++k4)
    {
        Set3DNormalVec(&x_n[k4 + 3 * qN], &y_n[k4 + 3 * qN], &z_n[k4 + 3 * qN], -qx_n[k4], -qy_n[k4], qz_n[k4]);
    }
    delete[] qx_n;
    delete[] qy_n;
    delete[] qz_n;
}


void Set3DNormalVec(double* v1, double* v2, double* v3, double x, double y, double z)
{
    double md = sqrt(x * x + y * y + z * z);
    x /= md;
    y /= md;
    z /= md;
    *v1 = x;
    *v2 = y;
    *v3 = z;
}

#pragma endregion