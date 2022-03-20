#include "../include/Scheme.h"
#include "../include/settings.h"
#include "../include/2DRPFunctions.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <iostream>
#include "../include/RedEquil.h"


using namespace std;


void DUGKS_gh(int X, int Y, double dt, double mode)
{
    memset(&macro(0, 0, 0), 0, sizeof(double) * 4 * X * Y);
    if (flag_collide == 0 || kappa == 0)
    {
#pragma region ����ײ���λ�����������
        CArray* fs{};
        if (flag_collide == 0)
        {
            fs = &MF_c;
        }
        else if (kappa == 0)
        {
            fs = &Meq;
        }
        CArray mf = *fs;   // ע�����mf���޷��ı�ģ�ֻ�ܶ�ȡ��ԭ������

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

                // �����������ϵ�\Delta t / 2����ֵ
                double f_ip_j_h, f_im_j_h, f_i_jp_h, f_i_jm_h;
                double xr, ys;
                int ik = 0;
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    ik = k % (M * N); // �õ�Сik
#pragma region VanLeer������5���ʽ
                    // ʹ��Van Leer�������ĸ�ʽ
                    xr = uk[ik] * 0.5 * (dt / dx);
                    ys = vk[ik] * 0.5 * (dt / dy);
                    // ʡ�Գ���dx��dy�Ĳ���
                    if (xr > 0)
                    {
                        // Ӧ��ʹ�ÿ���ĵ�
                        f_im_j_h = Interpol_5_vnLr(mf(i_m2, j, 0, 0, k), mf(i_m1, j, 0, 0, k),
                            mf(i, j, 0, 0, k), mf(i_m1, j_m1, 0, 0, k), mf(i_m1, j_p1, 0, 0, k),
                            0.5 - xr, -ys);
                        f_ip_j_h = Interpol_5_vnLr(mf(i_m1, j, 0, 0, k), mf(i, j, 0, 0, k),
                            mf(i_p1, j, 0, 0, k), mf(i, j_m1, 0, 0, k), mf(i, j_p1, 0, 0, k),
                            0.5 - xr, -ys);
                    }
                    else
                    {
                        // Ӧ��ʹ�ÿ��ҵĵ�
                        f_im_j_h = Interpol_5_vnLr(mf(i_m1, j, 0, 0, k), mf(i, j, 0, 0, k),
                            mf(i_p1, j, 0, 0, k), mf(i, j_m1, 0, 0, k), mf(i, j_p1, 0, 0, k),
                            -0.5 - xr, -ys);
                        f_ip_j_h = Interpol_5_vnLr(mf(i, j, 0, 0, k), mf(i_p1, j, 0, 0, k),
                            mf(i_p2, j, 0, 0, k), mf(i_p1, j_m1, 0, 0, k), mf(i_p1, j_p1, 0, 0, k),
                            -0.5 - xr, -ys);
                    }
                    if (ys > 0)
                    {
                        // Ӧ��ʹ�ÿ��µĵ�
                        f_i_jm_h = Interpol_5_vnLr(mf(i_m1, j_m1, 0, 0, k), mf(i, j_m1, 0, 0, k),
                            mf(i_p1, j_m1, 0, 0, k), mf(i, j_m2, 0, 0, k), mf(i, j, 0, 0, k),
                            -xr, 0.5 - ys);
                        f_i_jp_h = Interpol_5_vnLr(mf(i_m1, j, 0, 0, k), mf(i, j, 0, 0, k),
                            mf(i_p1, j, 0, 0, k), mf(i, j_m1, 0, 0, k), mf(i, j_p1, 0, 0, k),
                            -xr, 0.5 - ys);
                    }
                    else
                    {
                        // Ӧ��ʹ�ÿ��ϵĵ�
                        f_i_jm_h = Interpol_5_vnLr(mf(i_m1, j, 0, 0, k), mf(i, j, 0, 0, k),
                            mf(i_p1, j, 0, 0, k), mf(i, j_m1, 0, 0, k), mf(i, j_p1, 0, 0, k),
                            -xr, -0.5 - ys);
                        f_i_jp_h = Interpol_5_vnLr(mf(i_m1, j_p1, 0, 0, k), mf(i, j_p1, 0, 0, k),
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
                    // ������һʱ�̵ĺ����
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
#pragma region һ������
        CArray* fs = &MF_c;
        CArray mf = *fs;   // ע�����mf���޷��ı�ģ�ֻ�ܶ�ȡ��ԭ������
        int finish_count;
        int grid_num;
#pragma region ��һ��
        // ��һ�֣����f(i,j+-1/2,n+1/2)������X*(Y+1)�����浽MF_m��
        finish_count = 0;
        grid_num = X * (Y + 1);
#pragma omp parallel for
        for (int i = 0;i < X;++i)
        {
            for (int j = 0;j < Y + 1;++j)
            {
                // ���\bar f(i,j+-1/2,n+1/2)������X*(Y+1)�����浽MF_m��
                int i_p = min(i + 1, X - 1);
                int i_m = max(i - 1, 0);
                int j_p;
                int j_c;
                int j_m;

                double xr, ys;
                int ik = 0;
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    ik = k % (M * N); // �õ�Сik
#pragma region VanLeer������5���ʽ
                    // ʹ��Van Leer�������ĸ�ʽ
                    xr = uk[ik] * 0.5 * (dt / dx);
                    ys = vk[ik] * 0.5 * (dt / dy);
                    // ʡ�Գ���dx��dy�Ĳ���
                    if (ys > 0)
                    {
                        // Ӧ��ʹ�ÿ��µĵ�
                        j_p = min(j, Y - 1);
                        j_c = max(j - 1, 0);
                        j_m = max(j - 2, 0);
                        xr = -xr;
                        ys = 0.5 - ys;
                    }
                    else
                    {
                        // Ӧ��ʹ�ÿ��ϵĵ�
                        j_p = min(j + 1, Y - 1);
                        j_c = min(j, Y - 1);
                        j_m = max(j - 1, 0);
                        xr = -xr;
                        ys = -0.5 - ys;
                    }
                    MF_m(i, j, 0, 0, k) = Interpol_5_vnLr(mf(i_m, j_c, 0, 0, k),
                        mf(i, j_c, 0, 0, k), mf(i_p, j_c, 0, 0, k), mf(i, j_m, 0, 0, k),
                        mf(i, j_p, 0, 0, k), xr, ys);
#pragma endregion
                }
                // ����1/2�������
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
                // ׼��ת��Ϊf(i,j+-1/2,n+1/2)
                // ����ƽ��̬
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
                // ת��Ϊf(i,j+-1/2,n+1/2)��Guo 2013���£�20ʽ
                double e1 = (2 * kappa) / (2 * kappa + dt * 0.5);
                double e2 = (dt * 0.5) / (2 * kappa + dt * 0.5);
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    MF_m(i, j, 0, 0, k) = e1 * MF_m(i, j, 0, 0, k) + e2 * req_m[k];
                }
                delete[] req_m;

                finish_count++;
                printf("�����м�ƽ��̬X������� %d / %d    \r", finish_count, grid_num);
            }
        }
        // ����Щ����ȫ������MF
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

#pragma region �ڶ���
        // �ڶ��֣����f(i+-1/2,j,n+1/2)������(X+1)*Y�����浽MF_m��
        finish_count = 0;
        grid_num = X * (Y + 1);
#pragma omp parallel for
        for (int i = 0;i < X + 1;++i)
        {
            for (int j = 0;j < Y;++j)
            {
                // ���\bar f(i+-1/2,j,n+1/2)������(X+1)*Y�����浽MF_m��
                int j_p = min(j + 1, Y - 1);
                int j_m = max(j - 1, 0);
                int i_p;
                int i_c;
                int i_m;

                double xr, ys;
                int ik = 0;
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    ik = k % (M * N); // �õ�Сik
#pragma region VanLeer������5���ʽ
                    // ʹ��Van Leer�������ĸ�ʽ
                    xr = uk[ik] * 0.5 * (dt / dx);
                    ys = vk[ik] * 0.5 * (dt / dy);
                    // ʡ�Գ���dx��dy�Ĳ���
                    if (xr > 0)
                    {
                        // Ӧ��ʹ�ÿ���ĵ�
                        i_p = min(i, X - 1);
                        i_c = max(i - 1, 0);
                        i_m = max(i - 2, 0);
                        xr = 0.5 - xr;
                        ys = -ys;
                    }
                    else
                    {
                        // Ӧ��ʹ�ÿ��ҵĵ㡢
                        i_p = min(i + 1, X - 1);
                        i_c = min(i, X - 1);
                        i_m = max(i - 1, 0);
                        xr = -0.5 - xr;
                        ys = -ys;
                    }
                    MF_m(i, j, 0, 0, k) = Interpol_5_vnLr(mf(i_m, j, 0, 0, k),
                        mf(i_c, j, 0, 0, k), mf(i_p, j, 0, 0, k), mf(i_c, j_m, 0, 0, k),
                        mf(i_c, j_p, 0, 0, k), xr, ys);
#pragma endregion
                }
                // ����1/2�������
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
                // ׼��ת��Ϊf(i,j+-1/2,n+1/2)
                // ����ƽ��̬
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
                // ת��Ϊf(i+-1/2,j,n+1/2)��Guo 2013���£�20ʽ
                double e1 = (2 * kappa) / (2 * kappa + dt * 0.5);
                double e2 = (dt * 0.5) / (2 * kappa + dt * 0.5);
                for (int k = 2 * M * N - 1;k >= 0;--k)
                {
                    MF_m(i, j, 0, 0, k) = e1 * MF_m(i, j, 0, 0, k) + e2 * req_m[k];
                }
                delete[] req_m;

                finish_count++;
                printf("�����м�ƽ��̬Y������� %d / %d    \r", finish_count, grid_num);
            }
        }
        // ����Щ����ȫ������MF
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
                    // ������һʱ�̵ĺ����
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


void LBGK_gh(int X, int Y, double dt)
{
#pragma omp parallel for
    for (int i = 0;i < X;++i)
    {
        for (int j = 0;j < Y;++j)
        {
            int i_p = min(i + 1, X - 1);
            int i_m = max(i - 1, 0);
            int j_p = min(j + 1, Y - 1);
            int j_m = max(j - 1, 0);

            double xr, ys;
            int ik = 0;
            for (int k = 2 * M * N - 1;k >= 0;--k)
            {
                ik = k % (M * N); // �õ�Сik
#pragma region VanLeer������5���ʽ
                // ʹ��Van Leer�������ĸ�ʽ
                xr = uk[ik] * (dt / dx);
                ys = vk[ik] * (dt / dy);
                // ʡ�Գ���dx��dy�Ĳ���
                MF(i, j, 0, 0, k) = Interpol_5_vnLr(MF_c(i_m, j, 0, 0, k), MF_c(i, j, 0, 0, k),
                    MF_c(i_p, j, 0, 0, k), MF_c(i, j_m, 0, 0, k), MF_c(i, j_p, 0, 0, k),
                    -xr,-ys);
#pragma endregion
            }
            for (int k = M * N - 1;k >= 0;--k)
            {
                // ������һʱ�̵ĺ����
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

double Interpol_5_vnLr(double fl, double fc, double fr, double fd, double fu,
    double r, double s)
{
    return fc + Van_Leer_Limiter(fc - fl, fr - fc) * r
        + Van_Leer_Limiter(fc - fd, fu - fc) * s;
}

double Van_Leer_Limiter(double s_1, double s_2)
{
    if (sign(s_1) + sign(s_2) == 0) { return 0; }
    return 2 * s_1 * s_2 / (s_1 + s_2);
}