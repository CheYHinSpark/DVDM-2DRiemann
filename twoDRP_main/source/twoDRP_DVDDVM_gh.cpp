#include "../include/twoDRP.h"

void twoDRP_DVDDVM_gh(string init_setting, string compute_setting, int mode)
{
#pragma region 初始化计算参数，创建需要的数组
    L = 1;
    D = 2;
    // 设置各个方向，要避免平行于壁面的方向
    SetDirections();

    // 完成设置dvm
    SetDVM();

    // 以当前时间确定一个文件夹名
    string foldername;
    if (mode == 3)
    {
        foldername = "DVM" + to_string(N) + "x" + to_string(M)
            + init_setting + compute_setting + CurrentLocalTime() + SCHEME;
    }
    else if (mode == 4)
    {
        foldername = "DVDDVM" + to_string(N) + "x" + to_string(M)
            + init_setting + compute_setting + CurrentLocalTime() + SCHEME;
    }
    string command;
    command = "mkdir " + foldername;
    system(command.c_str());

    // 差分网格尺寸
    int x_num = int(0.5 / dx);    //  x轴格点数的一半
    int y_num = int(0.5 / dy);    //  y轴格点数的一半
    int grid_num = 4 * x_num * y_num;

    // 创建各种数组
    SetUVW();

    MF.reConstruct(2 * x_num, 2 * y_num, 2, N, M);      //前一时刻数据(x，y，gh, 方向，dvm)
    MF_c.reConstruct(2 * x_num, 2 * y_num, 2, N, M);    //复制数据(x，y，gh, 方向，dvm)
    Meq.reConstruct(2 * x_num, 2 * y_num, 2, N, M);     //平衡态
    CArray alpha(2 * x_num, 2 * y_num, 4);              //平衡态参数
    macro.reConstruct(2 * x_num, 2 * y_num, 4);         //宏观量
    if (flag_collide == 1 && kappa != 0)
    {
        MF_m.reConstruct(2 * x_num + 1, 2 * y_num + 1, 2, N, M);
    }

    // 保存用的
    double* save_t = new double[save_max] {};
    // 注意这里就是要把x和y反过来
    CArray save_r(2 * y_num, 2 * x_num);      // rho
    CArray save_u(2 * y_num, 2 * x_num);      // u
    CArray save_v(2 * y_num, 2 * x_num);      // v
    CArray save_e(2 * y_num, 2 * x_num);      // e
#pragma endregion

    if (mode == 3)
    {
        cout << "计算方法DVM_gh，";
    }
    else if (mode == 4)
    {
        cout << "计算方法DVDDVM_gh，";
    }
    cout << "方向数N=" << N << "，离散速度个数" << M << "，碰撞参数flag_collide=" << flag_collide <<
        "，kappa=" << kappa << "，网格数" << 2 * x_num << "*" << 2 * y_num << "，CFLMax=" << CFLMax << endl;
    cout << "将计算至t=" << (save_max - 1) * save_interval << endl;
    cout << "即将开始计算" << endl;

#pragma region 处理初始值
    // 初始能量
    for (int i = 0;i < 2;++i)
    {
        for (int j = 0;j < 2;++j)
        {
            E_init[i][j] = 0.5 * (pow(u_init[i][j], 2) + pow(v_init[i][j], 2) + 3 * p_init[i][j] / r_init[i][j]);
        }
    }
    double* MFi = new double[2 * M * N]{};
    double* req_i = new double[N] {};   // mode3不用
    double* ueq_i = new double[N] {};   // mode3不用
    double sig2eq_i = 0;
    // I,J大区坐标
    for (int I = 0;I < 2;++I)
    {
        for (int J = 0;J < 2;++J)
        {
            double u = u_init[I][J];
            double v = v_init[I][J];
            double rt = p_init[I][J] / r_init[I][J];
            double a[4] = { -1.5 * log(2 * M_PI * rt) - (u * u + v * v) / 2 / rt, u / rt, v / rt, -1 / rt };  // 分别为a, b1, b2, c

            if (mode == 3)
            {
                double sig2 = RedEquil_DVM(r_init[I][J], u_init[I][J], v_init[I][J], E_init[I][J],
                    a, MFi);
                for (int k = M * N - 1;k >= 0;--k)
                {
                    MFi[M * N + k] = MFi[k] * L * sig2;
                }
            }
            else if (mode == 4)
            {
                RedEquil_DVD(r_init[I][J], u_init[I][J], v_init[I][J], E_init[I][J],
                    a, req_i, ueq_i, &sig2eq_i);
                for (int n = 0;n < N;++n)
                {
                    double sig2 = RedEquil_DVDDVM(req_i[n], ueq_i[n], 0.5 * (ueq_i[n] * ueq_i[n] + (L + 1) * sig2eq_i),
                        &MFi[n * M]);
                    for (int m = 0;m < M;++m)
                    {
                        MFi[(N + n) * M + m] = MFi[n * M + m] * L * sig2;
                    }
                }
            }
            for (int i = 0;i < x_num;++i)
            {
                for (int j = 0;j < y_num;++j)
                {
                    memcpy(&alpha(i + I * x_num, j + J * y_num, 0), a, sizeof(double) * 4);
                    memcpy(&MF(i + I * x_num, j + J * y_num, 0, 0, 0), MFi, sizeof(double) * 2 * M * N);
                    // 初始宏观量
                    macro(i + I * x_num, j + J * y_num, 0) = r_init[I][J];
                    macro(i + I * x_num, j + J * y_num, 1) = u_init[I][J];
                    macro(i + I * x_num, j + J * y_num, 2) = v_init[I][J];
                    macro(i + I * x_num, j + J * y_num, 3) = E_init[I][J];

                    save_r(j + J * y_num, i + I * x_num) = r_init[I][J];
                    save_u(j + J * y_num, i + I * x_num) = u_init[I][J];
                    save_v(j + J * y_num, i + I * x_num) = v_init[I][J];
                    save_e(j + J * y_num, i + I * x_num) = E_init[I][J];
                }
            }
        }
    }
    delete[] MFi;
    delete[] req_i;
    delete[] ueq_i;

    // 第一次储存
    CSVWriteMatrix(foldername + "\\" + "r0.csv", &save_r(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "u0.csv", &save_u(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "v0.csv", &save_v(0, 0), 2 * y_num, 2 * x_num);
    CSVWriteMatrix(foldername + "\\" + "e0.csv", &save_e(0, 0), 2 * y_num, 2 * x_num);

#pragma endregion

    cout << "初始值处理完毕，主循环开始" << endl;

#pragma region 主循环
    double clock_1, clock_2, clock_3, clock_4, clock_5, clock_begin, clock_end;
    clock_begin = clock();
    while (Save_Count < save_max)
    {
        step_count++;
        clock_1 = clock();

#pragma region 第1步，确定时间步长
        // Guo: 可以无视\tau
        double dt = CFLMax * min(dx, dy) / abs(dvm_bg);
        t_current += dt;
#pragma endregion

        clock_2 = clock();

#pragma region 第2步，根据宏观量计算出各方向的平衡态
        // 如果是无碰撞情形可以跳过这一步
        if (flag_collide == 1)
        {
            int finish_count = 0;
            if (mode == 3)
            {
#pragma omp parallel for
                for (int i = grid_num - 1; i >= 0;i--)
                {
                    double sig2 = RedEquil_DVM(macro(0, i, 0), macro(0, i, 1), macro(0, i, 2), macro(0, i, 3),
                        &alpha(0, i, 0), &Meq(0, i, 0, 0, 0));
                    for (int j = M * N - 1;j >= 0;--j)
                    {
                        Meq(0, i, 1, 0, j) = L * sig2 * Meq(0, i, 0, 0, j);
                    }
                    finish_count++;

                    printf("计算平衡态，已完成 %d / %d\r", finish_count, grid_num);
                }
            }
            else if (mode == 4)
            {
#pragma omp parallel for
                for (int i = grid_num - 1; i >= 0;i--)
                {
                    double* req = new double[N] {};
                    double* ueq = new double[N] {};
                    double sig2eq = 0;
                    RedEquil_DVD(macro(0, i, 0), macro(0, i, 1), macro(0, i, 2), macro(0, i, 3),
                        &alpha(0, i, 0), req, ueq, &sig2eq);

                    for (int n = 0;n < N;++n)
                    {
                        double sig2 = RedEquil_DVDDVM(req[n], ueq[n], 0.5 * (ueq[n] * ueq[n] + (L + 1) * sig2eq),
                            &Meq(0, i, 0, n, 0));
                        for (int m = 0;m < M;++m)
                        {
                            Meq(0, i, 1, n, m) = L * sig2 * Meq(0, i, 0, n, m);
                        }
                    }
                    delete[] req;
                    delete[] ueq;
                    finish_count++;

                    printf("计算平衡态，已完成 %d / %d    \r", finish_count, grid_num);
                }
            }
        }
#pragma endregion

        clock_3 = clock();

#pragma region 第3步，复制

        if (flag_collide == 0 || kappa == 0)
        {
            swapPtr(&MF.ptr, &MF_c.ptr);
        }
        else
        {
            // 对于一般的情况MF存的其实是\tilde f，为了节约内存
            // 这里要将MF变成当前时刻的\tilde f+，
            // 而MF_c变成当前时刻的\bar f+，
            // Guo 2013文章21,22式
            double e1 = (2 * kappa - dt * 0.5) / (2 * kappa + dt);
            double e2 = (1.5 * dt) / (2 * kappa + dt);
#pragma omp parallel for
            for (int i = 8 * x_num * y_num * M * N - 1;i >= 0;--i)
            {
                MF_c(0, 0, 0, 0, i) = e1 * MF(0, 0, 0, 0, i) + e2 * Meq(0, 0, 0, 0, i);
                MF(0, 0, 0, 0, i) = 4.0 / 3.0 * MF_c(0, 0, 0, 0, i) - 1.0 / 3.0 * MF(0, 0, 0, 0, i);
            }
        }

#pragma endregion

        clock_4 = clock();

#pragma region 第4步，求解差分方程，更新双节点EQMOM参数，处理宏观量等

        // DUGKS格式
        DUGKS_gh(2 * x_num, 2 * y_num, dt, mode);

#pragma endregion

        clock_5 = clock();

        // 一次循环结束，输出当前时间
        cout << "第" << step_count + 1 << "步，" << "t=" << t_current
            << "，计算耗时" << (clock_5 - clock_1) / CLOCKS_PER_SEC << "秒"
            /* << "，算步长" << (clock_2 - clock_1) / CLOCKS_PER_SEC << "秒"
             << "，算平衡态" << (clock_3 - clock_2) / CLOCKS_PER_SEC << "秒"
             << "，算Flux" << (clock_4 - clock_3) / CLOCKS_PER_SEC << "秒"
             << "，求差分方程等" << (clock_5 - clock_4) / CLOCKS_PER_SEC << "秒"*/
            << endl;

#pragma region 看看是不是需要储存
        if (t_current >= Save_Count * save_interval)
        {
            save_t[Save_Count] = t_current;
#pragma omp parallel for
            for (int i = 0;i < 2 * x_num;++i)
            {
                for (int j = 0;j < 2 * y_num;++j)
                {
                    save_r(j, i) = macro(i, j, 0);
                    save_u(j, i) = macro(i, j, 1);
                    save_v(j, i) = macro(i, j, 2);
                    save_e(j, i) = macro(i, j, 3);
                }
            }

            CSVWriteMatrix(foldername + "\\" + "r" + to_string(Save_Count) + ".csv", &save_r(0, 0), 2 * y_num, 2 * x_num);
            CSVWriteMatrix(foldername + "\\" + "u" + to_string(Save_Count) + ".csv", &save_u(0, 0), 2 * y_num, 2 * x_num);
            CSVWriteMatrix(foldername + "\\" + "v" + to_string(Save_Count) + ".csv", &save_v(0, 0), 2 * y_num, 2 * x_num);
            CSVWriteMatrix(foldername + "\\" + "e" + to_string(Save_Count) + ".csv", &save_e(0, 0), 2 * y_num, 2 * x_num);

            clock_end = clock();
            cout << "第" << Save_Count << "次储存，当前t=" << t_current
                << "，计算已用时" << (clock_end - clock_begin) / CLOCKS_PER_SEC << "秒"
                << endl;
            Save_Count++;
        }
#pragma endregion
    }
    clock_end = clock();
#pragma endregion

    cout << "主循环结束，总用时" << (clock_end - clock_begin) / CLOCKS_PER_SEC << "秒" << endl;

#pragma region 储存数据

    CSVWriteArray(foldername + "\\" + "t.csv", &save_t[1], save_max - 1);
    CSVWriteArray(foldername + "\\" + "dvm.csv", &dvm[0], M);
    // 最后储存一些基本信息
    ofstream outFile;
    outFile.open(foldername + "\\" + "setting.csv", ios::out); // 打开模式可省略
    outFile << "N," << N << endl;
    outFile << "dvm_strategy," << dvm_strategy << endl;
    outFile << "dvm_bg," << dvm_bg << endl;
    outFile << "dvm_st," << dvm_st << endl;
    outFile << "dx," << dx << endl;
    outFile << "dy," << dy << endl;
    outFile << "flag_collide," << flag_collide << endl;
    outFile << "kappa," << kappa << endl;
    outFile.close();

#pragma endregion

    cout << "数据储存完毕" << endl;
}
